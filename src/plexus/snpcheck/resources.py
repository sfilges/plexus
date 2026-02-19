# ================================================================================
# SNP resource management — download, caching, and verification
#
# Manages the gnomAD AF-only VCF used for SNP checking. The file is downloaded
# explicitly via `plexus download-resources` and cached locally for reuse.
#
# Author: Stefan Filges (stefan.filges@pm.me)
# Copyright (c) 2025-2026 Stefan Filges
# ================================================================================

from __future__ import annotations

import os
import urllib.request
from pathlib import Path

from loguru import logger

# gnomAD AF-only VCF (hg38) from the GATK best-practices bucket
GNOMAD_VCF_URL = (
    "https://storage.googleapis.com/gatk-best-practices/somatic-hg38/"
    "af-only-gnomad.hg38.vcf.gz"
)
GNOMAD_TBI_URL = (
    "https://storage.googleapis.com/gatk-best-practices/somatic-hg38/"
    "af-only-gnomad.hg38.vcf.gz.tbi"
)

GNOMAD_VCF_FILENAME = "af-only-gnomad.hg38.vcf.gz"
GNOMAD_TBI_FILENAME = "af-only-gnomad.hg38.vcf.gz.tbi"

DEFAULT_CACHE_DIR = Path.home() / ".plexus" / "data"

ENV_DATA_DIR = "PLEXUS_DATA_DIR"
ENV_SNP_VCF = "PLEXUS_SNP_VCF"


def get_cache_dir() -> Path:
    """Return the cache directory, respecting ``$PLEXUS_DATA_DIR``."""
    env = os.environ.get(ENV_DATA_DIR)
    if env:
        return Path(env)
    return DEFAULT_CACHE_DIR


def get_cached_vcf_path() -> Path:
    """Expected path for the cached gnomAD VCF."""
    return get_cache_dir() / GNOMAD_VCF_FILENAME


def get_cached_tbi_path() -> Path:
    """Expected path for the cached gnomAD TBI index."""
    return get_cache_dir() / GNOMAD_TBI_FILENAME


def is_resource_available() -> bool:
    """Check whether both the cached VCF and its TBI index exist."""
    return get_cached_vcf_path().is_file() and get_cached_tbi_path().is_file()


def _download_with_progress(url: str, dest: Path) -> None:
    """Download *url* to *dest* with a rich progress bar.

    Uses an atomic ``.part`` file pattern — the final file only appears
    once the download completes successfully.
    """
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
        # Clean up partial file on failure
        part.unlink(missing_ok=True)
        raise


def download_gnomad_vcf(force: bool = False) -> Path:
    """Download the gnomAD AF-only VCF and TBI to the cache directory.

    Parameters
    ----------
    force : bool
        Re-download even if files already exist.

    Returns
    -------
    Path
        Path to the downloaded VCF.
    """
    vcf_path = get_cached_vcf_path()
    tbi_path = get_cached_tbi_path()

    if not force and is_resource_available():
        logger.info(f"gnomAD VCF already cached: {vcf_path}")
        return vcf_path

    cache_dir = get_cache_dir()
    cache_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"Downloading gnomAD VCF to {vcf_path}")
    _download_with_progress(GNOMAD_VCF_URL, vcf_path)

    logger.info(f"Downloading gnomAD TBI to {tbi_path}")
    _download_with_progress(GNOMAD_TBI_URL, tbi_path)

    logger.info("Download complete.")
    return vcf_path


def resource_status_message() -> str:
    """Return a human-readable status string for the SNP resource."""
    vcf_path = get_cached_vcf_path()
    if is_resource_available():
        size_mb = vcf_path.stat().st_size / (1024 * 1024)
        return f"gnomAD VCF: cached ({size_mb:.0f} MB) at {vcf_path}"
    return "gnomAD VCF: not downloaded (run `plexus init`)"
