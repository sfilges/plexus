# ================================================================================
# Genome resource management — presets, registry, and init orchestration
#
# Manages reference genome resources: FASTA downloads, BLAST indexes, and
# the genome registry (~/.plexus/data/registry.json).
#
# Author: Stefan Filges (stefan.filges@pm.me)
# Copyright (c) 2026 Stefan Filges
# ================================================================================

from __future__ import annotations

import gzip
import hashlib
import json
import os
import subprocess
import urllib.request
from pathlib import Path

from loguru import logger

from plexus.utils.utils import (
    detect_chrom_naming_mismatch,
    get_fasta_contigs,
    get_vcf_contigs,
)

# ---------------------------------------------------------------------------
# Cache directory
# ---------------------------------------------------------------------------

DEFAULT_DATA_DIR = Path.home() / ".plexus" / "data"
ENV_DATA_DIR = "PLEXUS_DATA_DIR"


def get_cache_dir() -> Path:
    """Return the data cache directory, respecting ``$PLEXUS_DATA_DIR``."""
    env = os.environ.get(ENV_DATA_DIR)
    if env:
        return Path(env)
    return DEFAULT_DATA_DIR


# ---------------------------------------------------------------------------
# Genome presets
# ---------------------------------------------------------------------------

GENOME_PRESETS: dict[str, dict] = {
    "hg38": {
        "description": "Human GRCh38/hg38",
        # FASTA
        "fasta_url": (
            "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
        ),
        "fasta_filename": "hg38.fa",
        "fasta_size_note": "~900 MB compressed / ~3.1 GB uncompressed",
        # SNP VCF — gnomAD AF-only from GATK best-practices bucket
        "snp_vcf_url": (
            "https://storage.googleapis.com/gatk-best-practices/somatic-hg38/"
            "af-only-gnomad.hg38.vcf.gz"
        ),
        "snp_tbi_url": (
            "https://storage.googleapis.com/gatk-best-practices/somatic-hg38/"
            "af-only-gnomad.hg38.vcf.gz.tbi"
        ),
        "snp_vcf_filename": "af-only-gnomad.hg38.vcf.gz",
        "snp_tbi_filename": "af-only-gnomad.hg38.vcf.gz.tbi",
    },
}

_REGISTRY_FILENAME = "registry.json"

# ---------------------------------------------------------------------------
# Registry helpers
# ---------------------------------------------------------------------------


def _registry_path() -> Path:
    return get_cache_dir() / _REGISTRY_FILENAME


def _load_registry() -> dict:
    path = _registry_path()
    if not path.is_file():
        return {}
    with path.open() as f:
        return json.load(f)


def _save_registry(registry: dict) -> None:
    path = _registry_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        json.dump(registry, f, indent=2)


def _register_fasta(genome: str, fasta_path: Path) -> None:
    registry = _load_registry()
    registry[genome] = {"fasta": str(fasta_path)}
    _save_registry(registry)
    logger.info(f"Registered {genome} FASTA: {fasta_path}")


# ---------------------------------------------------------------------------
# SHA-256 checksums
# ---------------------------------------------------------------------------

_CHUNK_SIZE = 1024 * 1024  # 1 MB


def _compute_sha256(path: Path) -> str:
    """Compute the SHA-256 hex digest of a file, reading in 1 MB chunks."""
    h = hashlib.sha256()
    with path.open("rb") as f:
        while chunk := f.read(_CHUNK_SIZE):
            h.update(chunk)
    return h.hexdigest()


def _parse_checksums_file(path: Path) -> dict[str, str]:
    """Parse a ``sha256sum``-format checksums file.

    Expected format (one entry per line)::

        <hex_digest>  <filename>

    Returns a dict mapping filename (basename only) to hex digest.
    Raises ``ValueError`` on malformed lines.
    """
    checksums: dict[str, str] = {}
    with path.open() as f:
        for lineno, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            # sha256sum uses two spaces; some tools use one space + asterisk
            parts = line.replace("*", " ").split(None, 1)
            if len(parts) != 2:
                raise ValueError(
                    f"Malformed checksums file {path}, line {lineno}: {line!r}"
                )
            digest, filename = parts
            checksums[Path(filename).name] = digest
    return checksums


def _register_genome_resources(
    genome: str,
    fasta_path: Path,
    snp_vcf_path: Path | None = None,
    expected_checksums: dict[str, str] | None = None,
) -> None:
    """Register genome resources with SHA-256 checksums.

    If *expected_checksums* is provided (``{filename: hex_digest}``), the
    on-disk file is hashed and compared.  A mismatch raises ``ValueError``
    **before** anything is written to the registry.
    """
    entry: dict[str, str] = {"fasta": str(fasta_path)}

    # FASTA checksum
    logger.info(f"Computing SHA-256 for {fasta_path.name} ...")
    actual_fasta = _compute_sha256(fasta_path)
    if expected_checksums and fasta_path.name in expected_checksums:
        expected = expected_checksums[fasta_path.name]
        if actual_fasta != expected:
            raise ValueError(
                f"FASTA checksum mismatch for {fasta_path.name}: "
                f"expected {expected[:16]}…, got {actual_fasta[:16]}…"
            )
        logger.info(f"FASTA checksum verified: {actual_fasta[:16]}…")
    entry["fasta_sha256"] = actual_fasta

    # SNP VCF checksum
    if snp_vcf_path is not None:
        logger.info(f"Computing SHA-256 for {snp_vcf_path.name} ...")
        actual_vcf = _compute_sha256(snp_vcf_path)
        if expected_checksums and snp_vcf_path.name in expected_checksums:
            expected = expected_checksums[snp_vcf_path.name]
            if actual_vcf != expected:
                raise ValueError(
                    f"SNP VCF checksum mismatch for {snp_vcf_path.name}: "
                    f"expected {expected[:16]}…, got {actual_vcf[:16]}…"
                )
            logger.info(f"SNP VCF checksum verified: {actual_vcf[:16]}…")
        entry["snp_vcf"] = str(snp_vcf_path)
        entry["snp_vcf_sha256"] = actual_vcf

    registry = _load_registry()
    registry[genome] = entry
    _save_registry(registry)
    logger.info(f"Registered {genome} resources in registry.")


def verify_resource_checksums(genome: str) -> dict[str, bool | None]:
    """Verify stored checksums for *genome* against on-disk files.

    Returns ``{resource: True | False | None}`` where ``None`` means
    no checksum is stored for that resource.
    """
    registry = _load_registry()
    entry = registry.get(genome, {})
    results: dict[str, bool | None] = {}

    for resource, hash_key in [
        ("fasta", "fasta_sha256"),
        ("snp_vcf", "snp_vcf_sha256"),
    ]:
        path_str = entry.get(resource)
        stored_hash = entry.get(hash_key)
        if not path_str or not stored_hash:
            results[resource] = None
            continue
        path = Path(path_str)
        if not path.is_file():
            results[resource] = False
            continue
        actual = _compute_sha256(path)
        results[resource] = actual == stored_hash
        if not results[resource]:
            logger.warning(
                f"Checksum mismatch for {resource}: "
                f"expected {stored_hash[:16]}…, got {actual[:16]}…"
            )
    return results


# ---------------------------------------------------------------------------
# Global configuration (operational mode)
# ---------------------------------------------------------------------------

_GLOBAL_CONFIG_FILENAME = "config.json"


def _global_config_path() -> Path:
    # ~/.plexus/config.json (one level above the data dir)
    return get_cache_dir().parent / _GLOBAL_CONFIG_FILENAME


def _load_global_config() -> dict:
    path = _global_config_path()
    if not path.is_file():
        return {"mode": "research"}
    with path.open() as f:
        return json.load(f)


def _save_global_config(cfg: dict) -> None:
    path = _global_config_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        json.dump(cfg, f, indent=2)


def get_operational_mode() -> str:
    """Return the operational mode: ``'research'`` or ``'compliance'``."""
    return _load_global_config().get("mode", "research")


def set_operational_mode(mode: str) -> None:
    """Set the operational mode (``'research'`` or ``'compliance'``)."""
    if mode not in ("research", "compliance"):
        raise ValueError(f"Invalid mode '{mode}'. Must be 'research' or 'compliance'.")
    cfg = _load_global_config()
    cfg["mode"] = mode
    _save_global_config(cfg)
    logger.info(f"Operational mode set to '{mode}'.")


# ---------------------------------------------------------------------------
# Public query functions
# ---------------------------------------------------------------------------


def list_genomes() -> list[str]:
    """Return names of all supported genomes."""
    return list(GENOME_PRESETS.keys())


def get_registered_fasta(genome: str) -> Path | None:
    """Return the registered FASTA path for *genome*, or None if absent/missing."""
    registry = _load_registry()
    entry = registry.get(genome)
    if not entry:
        return None
    fasta_path = Path(entry["fasta"])
    if not fasta_path.is_file():
        logger.warning(f"Registered FASTA for '{genome}' not found at {fasta_path}")
        return None
    return fasta_path


def is_blast_db_ready(fasta_path: Path) -> bool:
    """Return True if a BLAST database exists adjacent to *fasta_path*."""
    db = str(fasta_path).rsplit(".", 1)[0]
    return all(Path(f"{db}{suf}").is_file() for suf in (".nhr", ".nin", ".nsq"))


def is_fai_ready(fasta_path: Path) -> bool:
    """Return True if a samtools/pysam .fai index exists for *fasta_path*."""
    return Path(str(fasta_path) + ".fai").is_file()


def genome_status(genome: str) -> dict:
    """Return readiness flags and checksums for all resources of *genome*.

    Returns
    -------
    dict with bool flags (fasta, fai, blast_db, snp_vcf) and optional
    checksum strings (fasta_sha256, snp_vcf_sha256).
    """
    from plexus.snpcheck.resources import is_resource_available

    registry = _load_registry()
    entry = registry.get(genome, {})
    fasta = get_registered_fasta(genome)
    return {
        "fasta": fasta is not None,
        "fai": fasta is not None and is_fai_ready(fasta),
        "blast_db": fasta is not None and is_blast_db_ready(fasta),
        "snp_vcf": is_resource_available(),
        "fasta_sha256": entry.get("fasta_sha256"),
        "snp_vcf_sha256": entry.get("snp_vcf_sha256"),
    }


# ---------------------------------------------------------------------------
# Download helpers
# ---------------------------------------------------------------------------


def _download_with_progress(url: str, dest: Path) -> None:
    """Download *url* to *dest* with a rich progress bar (atomic .part pattern)."""
    from rich.progress import (
        BarColumn,
        DownloadColumn,
        Progress,
        TextColumn,
        TimeRemainingColumn,
        TransferSpeedColumn,
    )

    part = dest.with_suffix(dest.suffix + ".part")

    try:
        response = urllib.request.urlopen(url)  # noqa: S310
        total = int(response.headers.get("Content-Length", 0))

        with (
            Progress(
                TextColumn("[bold blue]{task.fields[filename]}"),
                BarColumn(),
                DownloadColumn(),
                TransferSpeedColumn(),
                TimeRemainingColumn(),
            ) as progress,
            part.open("wb") as f,
        ):
            task = progress.add_task(
                "download",
                total=total or None,
                filename=dest.name,
            )
            while True:
                chunk = response.read(1024 * 256)  # 256 KB chunks
                if not chunk:
                    break
                f.write(chunk)
                progress.advance(task, len(chunk))

        part.rename(dest)
    except Exception:
        part.unlink(missing_ok=True)
        raise


def _download_fasta(genome: str) -> Path:
    """Download and decompress the FASTA for *genome* into the genome cache dir.

    Returns
    -------
    Path to the uncompressed .fa file.
    """
    from rich.progress import (
        BarColumn,
        Progress,
        TextColumn,
        TimeRemainingColumn,
        TransferSpeedColumn,
    )

    preset = GENOME_PRESETS[genome]
    genome_dir = get_cache_dir() / "genomes" / genome
    genome_dir.mkdir(parents=True, exist_ok=True)

    gz_dest = genome_dir / (preset["fasta_filename"] + ".gz")
    fa_dest = genome_dir / preset["fasta_filename"]

    # Download compressed FASTA
    logger.info(
        f"Downloading {genome} FASTA ({preset['fasta_size_note']}) from UCSC ..."
    )
    _download_with_progress(preset["fasta_url"], gz_dest)

    # Decompress with progress tracking on compressed bytes consumed
    gz_size = gz_dest.stat().st_size
    logger.info(f"Decompressing {gz_dest.name} ...")

    try:
        with (
            Progress(
                TextColumn("[bold blue]{task.fields[filename]}"),
                BarColumn(),
                TransferSpeedColumn(),
                TimeRemainingColumn(),
            ) as progress,
            gzip.open(gz_dest, "rb") as gz_f,
            fa_dest.open("wb") as out_f,
        ):
            task = progress.add_task(
                "decompress",
                total=gz_size,
                filename=gz_dest.name,
            )
            last_compressed_pos = 0
            while True:
                chunk = gz_f.read(1024 * 256)
                if not chunk:
                    break
                out_f.write(chunk)
                # Track progress via position in the compressed file
                current_pos = gz_f.fileobj.tell()
                progress.advance(task, current_pos - last_compressed_pos)
                last_compressed_pos = current_pos
    except Exception:
        fa_dest.unlink(missing_ok=True)
        raise

    gz_dest.unlink()
    return fa_dest


# ---------------------------------------------------------------------------
# Index builders
# ---------------------------------------------------------------------------


def build_blast_index(fasta_path: Path) -> None:
    """Build a BLAST nucleotide database adjacent to *fasta_path*."""
    from plexus.blast.blast_runner import _check_blast_tools

    _check_blast_tools()
    db_path = str(fasta_path).rsplit(".", 1)[0]
    logger.info(f"Building BLAST index: {db_path}")
    subprocess.run(
        [
            "makeblastdb",
            "-in",
            str(fasta_path),
            "-dbtype",
            "nucl",
            "-parse_seqids",
            "-out",
            db_path,
        ],
        check=True,
    )


# ---------------------------------------------------------------------------
# Top-level init orchestrator
# ---------------------------------------------------------------------------


def init_genome(
    genome: str,
    fasta: Path | None = None,
    snp_vcf: Path | None = None,
    force: bool = False,
    skip_blast: bool = False,
    skip_snp: bool = False,
    download: bool = False,
    mode: str | None = None,
    checksums: Path | None = None,
) -> None:
    """Register and index reference resources for *genome*.

    Parameters
    ----------
    genome     : Genome name, must be a key in GENOME_PRESETS.
    fasta      : Path to a local FASTA. Required unless *download* is True.
    snp_vcf    : Path to a local gnomAD VCF. Required unless *download* or
                 *skip_snp* is True.
    force      : Rebuild/re-download even if resources already exist.
    skip_blast : Skip BLAST index creation.
    skip_snp   : Skip SNP VCF step entirely.
    download   : Download FASTA and/or VCF from preset URLs when local
                 paths are not provided.
    mode       : If given, set the operational mode ('research' or 'compliance').
    checksums  : Path to a sha256sum-format checksums file. When provided,
                 on-disk files are verified against expected hashes before
                 registration. A mismatch raises ``ValueError``.
    """
    from plexus.snpcheck.resources import download_gnomad_vcf

    if genome not in GENOME_PRESETS:
        supported = ", ".join(GENOME_PRESETS)
        raise ValueError(f"Unknown genome '{genome}'. Supported genomes: {supported}")

    # ── Operational mode ─────────────────────────────────────────────────────
    if mode is not None:
        set_operational_mode(mode)

    # ── Parse checksums file (if provided) ───────────────────────────────────
    expected_checksums: dict[str, str] | None = None
    if checksums is not None:
        if not checksums.is_file():
            raise FileNotFoundError(f"Checksums file not found: {checksums}")
        expected_checksums = _parse_checksums_file(checksums)

    # ── Step 1: SNP VCF ──────────────────────────────────────────────────────
    resolved_vcf: Path | None = None
    if not skip_snp:
        if snp_vcf is not None:
            if not Path(snp_vcf).is_file():
                raise FileNotFoundError(f"SNP VCF not found: {snp_vcf}")
            logger.info(f"Using user-supplied SNP VCF: {snp_vcf}")
            resolved_vcf = Path(snp_vcf)
        elif download:
            logger.info("Downloading gnomAD AF-only VCF ...")
            resolved_vcf = download_gnomad_vcf(force=force)
        else:
            raise ValueError(
                "No SNP VCF provided. Pass --snp-vcf /path/to/snps.vcf.gz "
                "or use --download to fetch the gnomAD VCF automatically."
            )

    # ── Step 2: FASTA ────────────────────────────────────────────────────────
    if fasta is not None:
        if not fasta.is_file():
            raise FileNotFoundError(f"FASTA file not found: {fasta}")
        fasta_path = fasta
    elif download:
        fasta_path = _download_fasta(genome)
    else:
        raise ValueError(
            "No FASTA provided. Pass --fasta /path/to/genome.fa "
            "or use --download to fetch the genome automatically."
        )

    # ── Step 3: FAI index ────────────────────────────────────────────────────
    if not is_fai_ready(fasta_path) or force:
        import pysam

        logger.info(f"Building .fai index for {fasta_path.name} ...")
        pysam.faidx(str(fasta_path))
    else:
        logger.info(".fai index already present, skipping.")

    # ── Step 4: BLAST index ──────────────────────────────────────────────────
    if not skip_blast:
        if not is_blast_db_ready(fasta_path) or force:
            build_blast_index(fasta_path)
        else:
            logger.info("BLAST index already present, skipping.")

    # ── Step 5: Chromosome naming validation ─────────────────────────────────
    if resolved_vcf is not None:
        mismatch = None
        try:
            logger.info("Verifying chromosome naming consistency...")
            fasta_contigs = get_fasta_contigs(fasta_path)
            vcf_contigs = get_vcf_contigs(resolved_vcf)
            mismatch = detect_chrom_naming_mismatch(fasta_contigs, vcf_contigs)
        except Exception as e:
            logger.warning(f"Could not verify chromosome naming consistency: {e}")

        if mismatch:
            # In compliance/clinical mode, this is a fatal error.
            # In research mode, it is a loud warning.
            logger.warning(f"Resource naming conflict detected: {mismatch}")
            if get_operational_mode() == "compliance":
                logger.error(
                    "Clinical/Compliance mode active: initialization aborted "
                    "due to naming mismatch."
                )
                raise ValueError(mismatch)

    # ── Step 6: Register with checksums ──────────────────────────────────────
    _register_genome_resources(
        genome,
        fasta_path,
        snp_vcf_path=resolved_vcf,
        expected_checksums=expected_checksums,
    )
