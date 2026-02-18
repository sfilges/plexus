# ================================================================================
# SNP overlap checking using tabix-indexed VCF/BCF files
#
# Queries a local VCF/BCF file (e.g. dbSNP) to find variants overlapping
# primer binding sites. Primers with SNP overlaps receive a penalty that
# feeds into the multiplex optimizer via pair_penalty.
#
# Author: Stefan Filges (stefan@simsendiagnostics.com)
# Copyright (c) 2025-2026 Simsen Diagnostics AB
# ================================================================================

from __future__ import annotations

from loguru import logger

from plexus.designer.multiplexpanel import MultiplexPanel


def _primer_genomic_coords(primer, junction):
    """Return (chrom, genomic_start, genomic_end) for a primer.

    Primer.start is 0-based within the design region.
    junction.design_start is 1-based genomic coordinate of the design region.
    """
    genomic_start = junction.design_start + primer.start
    genomic_end = genomic_start + primer.length
    return junction.chrom, genomic_start, genomic_end


def _get_allele_frequency(record):
    """Extract the maximum alternate allele frequency from a VCF record.

    Handles common INFO field names: AF, CAF, FREQ.
    Returns float or None.
    """
    if "AF" in record.info:
        af_values = record.info["AF"]
        if isinstance(af_values, list | tuple):
            return max(float(v) for v in af_values)
        return float(af_values)

    if "CAF" in record.info:
        caf = record.info["CAF"]
        if len(caf) > 1:
            alt_freqs = [float(x) for x in caf[1:] if x != "."]
            return max(alt_freqs) if alt_freqs else None

    return None


def _count_snps_in_region(vcf, chrom, start, end, af_threshold):
    """Count variants with AF >= af_threshold in a genomic region.

    Parameters
    ----------
    vcf : pysam.VariantFile
    chrom : str
    start : int
        1-based inclusive start.
    end : int
        1-based exclusive end.
    af_threshold : float

    Returns
    -------
    tuple[int, list[int]]
        (count, list of 1-based positions)
    """
    count = 0
    positions = []
    try:
        # pysam.VariantFile.fetch uses 0-based half-open coordinates
        for record in vcf.fetch(chrom, start - 1, end - 1):
            af = _get_allele_frequency(record)
            if af is not None and af >= af_threshold:
                count += 1
                positions.append(record.pos)
    except ValueError:
        # Contig not in VCF (e.g. alt contigs)
        pass
    return count, positions


def run_snp_check(
    panel: MultiplexPanel,
    vcf_path: str,
    af_threshold: float = 0.01,
    snp_penalty_weight: float = 5.0,
) -> None:
    """Check all primer pairs in the panel for SNP overlaps using a local VCF.

    Modifies Primer and PrimerPair objects in-place:
      - primer.snp_count
      - pair.snp_count, pair.snp_penalty
      - pair.pair_penalty (SNP penalty added)

    Parameters
    ----------
    panel : MultiplexPanel
    vcf_path : str
        Path to tabix-indexed VCF or BCF file.
    af_threshold : float
        Minimum allele frequency (default 0.01 = 1%).
    snp_penalty_weight : float
        Penalty per SNP (default 5.0).
    """
    import pysam

    logger.info(f"Running SNP check with VCF: {vcf_path}")
    logger.info(f"AF threshold: {af_threshold}, penalty weight: {snp_penalty_weight}")

    with pysam.VariantFile(vcf_path) as vcf:
        total_snps = 0
        primers_with_snps = 0

        for junction in panel.junctions:
            if not junction.primer_pairs:
                continue
            if not hasattr(junction, "design_start") or junction.design_start is None:
                logger.warning(
                    f"No design_start for junction {junction.name}, skipping SNP check"
                )
                continue

            for pair in junction.primer_pairs:
                pair_snp_count = 0

                for primer in (pair.forward, pair.reverse):
                    chrom, gstart, gend = _primer_genomic_coords(primer, junction)
                    count, positions = _count_snps_in_region(
                        vcf, chrom, gstart, gend, af_threshold
                    )

                    primer.snp_count = count
                    pair_snp_count += count

                    if count > 0:
                        primers_with_snps += 1
                        total_snps += count
                        logger.debug(
                            f"Primer {primer.name}: {count} SNPs at positions {positions}"
                        )

                pair.snp_count = pair_snp_count
                pair.snp_penalty = pair_snp_count * snp_penalty_weight

                if pair.snp_penalty > 0:
                    if pair.pair_penalty is not None:
                        pair.pair_penalty += pair.snp_penalty
                    else:
                        pair.pair_penalty = pair.snp_penalty

    logger.info(
        f"SNP check complete: {total_snps} SNPs found across {primers_with_snps} primers"
    )


def filter_snp_pairs(panel: MultiplexPanel) -> int:
    """Remove primer pairs that overlap any SNP (snp_count > 0).

    For each junction, pairs with ``snp_count == 0`` are kept.  If *all*
    pairs for a junction have SNPs, the single pair with the lowest
    ``snp_count`` is retained so the junction is not lost from the panel.

    Parameters
    ----------
    panel : MultiplexPanel
        Panel whose junctions have already been through ``run_snp_check``.

    Returns
    -------
    int
        Total number of primer pairs removed across all junctions.
    """
    total_removed = 0

    for junction in panel.junctions:
        if not junction.primer_pairs:
            continue

        clean = [p for p in junction.primer_pairs if p.snp_count == 0]
        dirty = [p for p in junction.primer_pairs if p.snp_count > 0]

        if not dirty:
            # Nothing to remove
            continue

        if clean:
            removed = len(dirty)
            junction.primer_pairs = clean
            logger.info(
                f"Junction {junction.name}: removed {removed} primer pair(s) "
                f"overlapping SNPs, {len(clean)} clean pair(s) remain"
            )
        else:
            # All pairs have SNPs — keep the least affected one
            best = min(junction.primer_pairs, key=lambda p: p.snp_count)
            removed = len(junction.primer_pairs) - 1
            junction.primer_pairs = [best]
            logger.warning(
                f"Junction {junction.name}: all pairs overlap SNPs; "
                f"keeping pair {best.pair_id} with lowest snp_count={best.snp_count}"
            )

        total_removed += removed

    return total_removed
