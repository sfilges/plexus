# ================================================================================
# SNP overlap checking using tabix-indexed VCF/BCF files
#
# Queries a local VCF/BCF file (e.g. dbSNP) to find variants overlapping
# primer binding sites. Primers with SNP overlaps receive a penalty that
# feeds into the multiplex optimizer via pair_penalty.
#
# Author: Stefan Filges (stefan.filges@pm.me)
# Copyright (c) 2025-2026 Stefan Filges
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
    """Find variants with AF >= af_threshold in a genomic region.

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
    tuple[int, list[tuple[int, float]]]
        (count, list of (1-based position, allele frequency) tuples)
    """
    count = 0
    snps = []
    try:
        # pysam.VariantFile.fetch uses 0-based half-open coordinates
        for record in vcf.fetch(chrom, start - 1, end - 1):
            af = _get_allele_frequency(record)
            if af is not None and af >= af_threshold:
                count += 1
                snps.append((record.pos, af))
    except ValueError:
        # Contig not in VCF (e.g. alt contigs)
        pass
    return count, snps


def _calc_weighted_snp_penalty(
    snps: list[tuple[int, float]],
    primer_genomic_start: int,
    primer_length: int,
    orientation: str,
    base_weight: float,
    three_prime_window: int,
    three_prime_multiplier: float,
) -> float:
    """Weight SNP penalties by proximity to the primer's 3' end.

    Coordinate logic (verified against design.py k-mer generation):
    - Forward primer: 5'->3' runs left-to-right on forward strand
      3' end genomic pos = genomic_start + length - 1
      dist_from_3prime = (genomic_start + length - 1) - snp_pos
    - Reverse primer: 5'->3' runs right-to-left on forward strand
      3' end genomic pos = genomic_start (the leftmost coordinate)
      dist_from_3prime = snp_pos - genomic_start
    """
    penalty = 0.0
    for snp_pos, _af in snps:
        if orientation == "forward":
            dist_from_3prime = (primer_genomic_start + primer_length - 1) - snp_pos
        else:
            dist_from_3prime = snp_pos - primer_genomic_start

        # Clamp to valid range (BLAST coords may be slightly off)
        dist_from_3prime = max(0, dist_from_3prime)

        multiplier = (
            three_prime_multiplier if dist_from_3prime < three_prime_window else 1.0
        )
        penalty += base_weight * multiplier
    return penalty


def run_snp_check(
    panel: MultiplexPanel,
    vcf_path: str,
    af_threshold: float = 0.01,
    *,
    snp_penalty_weight: float,
    snp_3prime_window: int = 5,
    snp_3prime_multiplier: float = 3.0,
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
        Base penalty per SNP.
    snp_3prime_window : int
        Number of bases from 3' end considered high-impact (default 5).
    snp_3prime_multiplier : float
        Penalty multiplier for SNPs within the 3' window (default 3.0).
    """
    import pysam

    logger.info(f"Running SNP check with VCF: {vcf_path}")
    logger.info(f"AF threshold: {af_threshold}, penalty weight: {snp_penalty_weight}")
    logger.info(
        f"3' window: {snp_3prime_window} bp, 3' multiplier: {snp_3prime_multiplier}x"
    )

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
                pair_snp_penalty = 0.0

                for primer in (pair.forward, pair.reverse):
                    chrom, gstart, gend = _primer_genomic_coords(primer, junction)
                    count, snps = _count_snps_in_region(
                        vcf, chrom, gstart, gend, af_threshold
                    )

                    primer.snp_count = count
                    pair_snp_count += count

                    if count > 0:
                        primers_with_snps += 1
                        total_snps += count
                        orientation = (
                            "forward" if primer.direction == "forward" else "reverse"
                        )
                        penalty = _calc_weighted_snp_penalty(
                            snps=snps,
                            primer_genomic_start=gstart,
                            primer_length=primer.length,
                            orientation=orientation,
                            base_weight=snp_penalty_weight,
                            three_prime_window=snp_3prime_window,
                            three_prime_multiplier=snp_3prime_multiplier,
                        )
                        pair_snp_penalty += penalty
                        positions = [pos for pos, _ in snps]
                        logger.debug(
                            f"Primer {primer.name}: {count} SNPs at positions {positions}"
                        )

                pair.snp_count = pair_snp_count
                pair.snp_penalty = pair_snp_penalty

                if pair.snp_penalty > 0:
                    if pair.pair_penalty is not None:
                        pair.pair_penalty += pair.snp_penalty
                    else:
                        pair.pair_penalty = pair.snp_penalty

    logger.info(
        f"SNP check complete: {total_snps} SNPs found across {primers_with_snps} primers"
    )


def filter_snp_pairs(panel: MultiplexPanel) -> tuple[int, list[str]]:
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
    tuple[int, list[str]]
        Total number of primer pairs removed across all junctions, and a list
        of junction names where the fallback (all pairs had SNPs) was triggered.
    """
    total_removed = 0
    fallback_junctions: list[str] = []

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
            fallback_junctions.append(junction.name)
            logger.warning(
                f"Junction {junction.name}: all pairs overlap SNPs; "
                f"keeping pair {best.pair_id} with lowest snp_count={best.snp_count}"
            )

        total_removed += removed

    return total_removed, fallback_junctions
