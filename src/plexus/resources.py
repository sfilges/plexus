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
import json
import subprocess
import urllib.request
from pathlib import Path

from loguru import logger

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
    from plexus.snpcheck.resources import get_cache_dir

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


def genome_status(genome: str) -> dict[str, bool]:
    """Return readiness flags for all resources of *genome*.

    Returns
    -------
    dict with bool flags: fasta, fai, blast_db, snp_vcf
    """
    from plexus.snpcheck.resources import is_resource_available

    fasta = get_registered_fasta(genome)
    return {
        "fasta": fasta is not None,
        "fai": fasta is not None and is_fai_ready(fasta),
        "blast_db": fasta is not None and is_blast_db_ready(fasta),
        "snp_vcf": is_resource_available(),
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

    from plexus.snpcheck.resources import get_cache_dir

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
) -> None:
    """Download and index all resources for *genome*.

    Parameters
    ----------
    genome     : Genome name, must be a key in GENOME_PRESETS.
    fasta      : Path to a local FASTA. If None, downloads from preset URL.
    snp_vcf    : Path to a local gnomAD VCF. If None, downloads from preset URL.
    force      : Rebuild/re-download even if resources already exist.
    skip_blast : Skip BLAST index creation.
    skip_snp   : Skip SNP VCF step entirely.
    """
    from plexus.snpcheck.resources import download_gnomad_vcf

    if genome not in GENOME_PRESETS:
        supported = ", ".join(GENOME_PRESETS)
        raise ValueError(f"Unknown genome '{genome}'. Supported genomes: {supported}")

    # ── Step 1: SNP VCF ──────────────────────────────────────────────────────
    if not skip_snp:
        if snp_vcf is not None:
            logger.info(f"Using user-supplied SNP VCF: {snp_vcf}")
        else:
            logger.info("Downloading gnomAD AF-only VCF ...")
            download_gnomad_vcf(force=force)

    # ── Step 2: FASTA ────────────────────────────────────────────────────────
    if fasta is not None:
        if not fasta.is_file():
            raise FileNotFoundError(f"FASTA file not found: {fasta}")
        fasta_path = fasta
    else:
        fasta_path = _download_fasta(genome)

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

    # ── Step 5: Register in genome registry ──────────────────────────────────
    _register_fasta(genome, fasta_path)
