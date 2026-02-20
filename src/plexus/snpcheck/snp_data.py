# ================================================================================
# SNP data acquisition and preparation
#
# Provides automatic SNP VCF resolution for primer checking. When the user
# doesn't provide a VCF, uses the cached gnomAD AF-only VCF and intersects
# it with panel target regions via bcftools to produce a small
# region-specific VCF.
#
# Author: Stefan Filges (stefan.filges@pm.me)
# Copyright (c) 2026 Stefan Filges
# ================================================================================

from __future__ import annotations

import os
import shutil
import tempfile
from pathlib import Path

from loguru import logger

from plexus.designer.multiplexpanel import MultiplexPanel
from plexus.snpcheck.resources import (
    ENV_SNP_VCF,
    get_cached_vcf_path,
    is_resource_available,
)
from plexus.utils.env import check_executable
from plexus.utils.utils import (
    detect_chrom_naming_mismatch,
    get_vcf_contigs,
    run_command,
)


def _check_bcftools() -> str:
    """Verify bcftools is available and return the path."""

    bcftools = shutil.which("bcftools")
    if bcftools is None or not check_executable("bcftools"):
        raise RuntimeError(
            "bcftools is required for SNP checking but was not found on PATH. "
            "Install with: conda install -c bioconda bcftools"
        )
    return bcftools


def _write_regions_bed(
    panel: MultiplexPanel,
    tmpdir: Path,
    padding: int = 200,
    valid_contigs: set[str] | None = None,
) -> Path:
    """Generate a BED file of panel regions for bcftools intersection.

    Parameters
    ----------
    panel : MultiplexPanel
        Panel with junctions that have design_start coordinates.
    tmpdir : Path
        Directory to write the BED file.
    padding : int
        Extra bases around each junction region.
    valid_contigs : set[str] | None
        Set of contig names that are valid in the target VCF. If provided,
        regions on other contigs will be skipped to avoid bcftools errors.

    Returns
    -------
    Path
        Path to the BED file.
    """
    bed_path = tmpdir / "panel_regions.bed"
    count = 0
    with bed_path.open("w") as f:
        for junction in panel.junctions:
            if not hasattr(junction, "design_start") or junction.design_start is None:
                continue
            chrom = junction.chrom

            if valid_contigs is not None and chrom not in valid_contigs:
                logger.debug(f"Skipping region on contig '{chrom}' (not in SNP VCF)")
                continue

            start = max(0, junction.design_start - padding)
            if hasattr(junction, "design_end") and junction.design_end is not None:
                end = junction.design_end + padding
            else:
                end = max(junction.start, junction.end) + padding

            f.write(f"{chrom}\t{start}\t{end}\n")
            count += 1

    if count == 0:
        logger.warning(
            "No regions from the panel were found in the SNP VCF contig list."
        )

    return bed_path


def intersect_vcf_with_regions(
    vcf_path: Path,
    panel: MultiplexPanel,
    output_dir: Path,
    padding: int = 200,
) -> Path:
    """Intersect a VCF with panel target regions.

    Uses ``bcftools view -R`` to extract only variants overlapping the
    panel's junction regions, producing a small output VCF.

    Parameters
    ----------
    vcf_path : Path
        Path to the source VCF (tabix-indexed).
    panel : MultiplexPanel
        Panel with junction coordinates.
    output_dir : Path
        Directory to write the region-specific VCF.
    padding : int
        Extra bases around each junction.

    Returns
    -------
    Path
        Path to the region-specific VCF.
    """
    bcftools = _check_bcftools()
    output_vcf = output_dir / "snpcheck_regions.vcf.gz"

    if output_vcf.exists():
        logger.info(f"Using existing region-intersected VCF: {output_vcf}")
        return output_vcf

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)

        # Pre-filter contigs to avoid bcftools crashing on missing contigs
        logger.debug(f"Getting contig list from {vcf_path}...")
        try:
            vcf_contigs = get_vcf_contigs(vcf_path)

            # DOC-02: detect chr-prefix mismatch before silently producing empty results
            panel_chroms = {
                junction.chrom
                for junction in panel.junctions
                if getattr(junction, "chrom", None) is not None
            }
            mismatch_msg = detect_chrom_naming_mismatch(panel_chroms, vcf_contigs)
            if mismatch_msg:
                logger.warning(mismatch_msg)

            bed_path = _write_regions_bed(
                panel, tmpdir_path, padding, valid_contigs=vcf_contigs
            )
        except Exception as e:
            logger.warning(
                f"Failed to get contig list from VCF: {e}. Proceeding without filtering."
            )
            bed_path = _write_regions_bed(panel, tmpdir_path, padding)

        logger.info("Extracting panel regions from SNP VCF...")

        run_command(
            [
                bcftools,
                "view",
                "-R",
                str(bed_path),
                str(vcf_path),
                "-Oz",
                "-o",
                str(output_vcf),
            ],
            check=True,
            retries=3,
        )

        # Index the output
        run_command(
            [bcftools, "index", "-t", str(output_vcf)],
            check=True,
        )

    logger.info(f"Created region-specific VCF: {output_vcf}")
    return output_vcf


def get_snp_vcf(
    panel: MultiplexPanel,
    output_dir: Path,
    user_vcf: str | Path | None = None,
    padding: int = 200,
) -> Path:
    """Get a SNP VCF ready for primer checking.

    Resolution order:
      1. User-provided VCF (``--snp-vcf``) — used directly
      2. ``$PLEXUS_SNP_VCF`` environment variable — used directly
      3. Cached gnomAD VCF — intersected with panel regions
      4. Error with actionable instructions

    Parameters
    ----------
    panel : MultiplexPanel
        Panel with junction coordinates.
    output_dir : Path
        Output directory for intermediate files.
    user_vcf : str | Path | None
        User-provided VCF path. If given, used directly.
    padding : int
        Extra bases around each junction region.

    Returns
    -------
    Path
        Path to a tabix-indexed VCF ready for SNP checking.
    """
    # 1. User-provided VCF (--snp-vcf flag)
    if user_vcf is not None:
        user_vcf = Path(user_vcf)
        if not user_vcf.exists():
            raise FileNotFoundError(f"SNP VCF not found: {user_vcf}")
        logger.info(f"Using user-provided SNP VCF: {user_vcf}")
        return user_vcf

    # 2. Environment variable
    env_vcf = os.environ.get(ENV_SNP_VCF)
    if env_vcf:
        env_path = Path(env_vcf)
        if not env_path.exists():
            raise FileNotFoundError(
                f"SNP VCF from ${ENV_SNP_VCF} not found: {env_path}"
            )
        logger.info(f"Using SNP VCF from ${ENV_SNP_VCF}: {env_path}")
        return env_path

    # 3. Cached gnomAD VCF → intersect with panel regions
    if is_resource_available():
        cached_vcf = get_cached_vcf_path()
        logger.info(f"Using cached gnomAD VCF: {cached_vcf}")
        return intersect_vcf_with_regions(
            cached_vcf, panel, output_dir, padding=padding
        )

    # 4. Nothing available — actionable error
    raise FileNotFoundError(
        "No SNP VCF available. Options:\n"
        "  1. Run `plexus init` to download the gnomAD VCF\n"
        "  2. Provide a VCF via --snp-vcf /path/to/snps.vcf.gz\n"
        "  3. Set $PLEXUS_SNP_VCF=/path/to/snps.vcf.gz\n"
        "  4. Skip SNP checking with --skip-snpcheck"
    )
