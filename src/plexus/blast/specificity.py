import os
import resource

import pandas as pd
from loguru import logger

from plexus.blast.annotator import BlastResultsAnnotator
from plexus.blast.blast_runner import BlastRunner
from plexus.blast.offtarget_finder import AmpliconFinder
from plexus.designer.multiplexpanel import MultiplexPanel


def _log_df_memory(name: str, df: pd.DataFrame) -> None:
    """Log row count and deep memory usage of a DataFrame."""
    mb = df.memory_usage(deep=True).sum() / 1_048_576
    logger.info(f"{name}: {len(df):,} rows, {mb:.1f} MB")


def _log_process_memory() -> None:
    """Log peak RSS of the current process."""
    ru = resource.getrusage(resource.RUSAGE_SELF)
    # macOS reports ru_maxrss in bytes; Linux in kilobytes
    peak_mb = (
        ru.ru_maxrss / 1_048_576
        if os.uname().sysname == "Darwin"
        else ru.ru_maxrss / 1024
    )
    logger.info(f"Process peak RSS: {peak_mb:.1f} MB")


def run_specificity_check(
    panel: MultiplexPanel,
    work_dir: str,
    genome_fasta: str,
    num_threads: int = 1,
    *,
    length_threshold: int = 15,
    evalue_threshold: float = 10.0,
    max_mismatches: int = 3,
    three_prime_tolerance: int = 3,
    max_amplicon_size: int = 2000,
    ontarget_tolerance: int = 5,
    blast_evalue: float = 30000.0,
    blast_word_size: int = 7,
    blast_reward: int = 1,
    blast_penalty: int = -1,
    blast_max_hsps: int = 100,
    blast_dust: str = "yes",
    max_bound_per_primer: int | None = 10000,
):
    """
    Run BLAST on all candidate primers in the panel to check for specificity
    and identify off-target amplicons.

    This assumes that the panel has already aggregated primers and assigned unique IDs.

    Args:
        panel: The MultiplexPanel object containing junctions and primer designs.
        work_dir: Directory to store temporary BLAST files.
        genome_fasta: Path to the reference genome FASTA file.
        num_threads: Number of BLAST threads.
        length_threshold: Minimum 3'-anchored alignment length (bp) to predict binding.
        evalue_threshold: E-value cutoff for predicted binding.
        max_mismatches: Maximum mismatches in a 3'-anchored alignment.
        three_prime_tolerance: Max unaligned bases at the 3' end for a hit to be
            considered 3'-anchored.  Compensates for BLAST clipping terminal mismatches.
        max_amplicon_size: Maximum amplicon size (bp) to consider.
        ontarget_tolerance: Coordinate tolerance (bp) for on-target classification.
        blast_evalue: E-value passed to blastn via -evalue to control search sensitivity.
        blast_word_size: Word size passed to blastn via -word_size.
        blast_reward: Match reward passed to blastn via -reward.
        blast_penalty: Mismatch penalty passed to blastn via -penalty.
    """
    logger.info("Starting specificity check (BLAST)...")
    os.makedirs(work_dir, exist_ok=True)

    # 1. Ensure primers are aggregated and mapped
    if not hasattr(panel, "unique_primer_map") or not panel.unique_primer_map:
        panel.aggregate_primers()

    if not panel.unique_primer_map:
        logger.warning("No primers found to blast.")
        return

    # 2. Write FASTA using the panel's method
    input_fasta = os.path.join(work_dir, "all_primers.fasta")
    panel.save_candidate_primers_to_fasta(input_fasta)

    # 3. Run BLAST
    blast_table = os.path.join(work_dir, "blast_table.txt")

    runner = BlastRunner(input_fasta, genome_fasta)
    runner.create_database()
    runner.run(
        output_table=blast_table,
        num_threads=num_threads,
        evalue=blast_evalue,
        word_size=blast_word_size,
        reward=blast_reward,
        penalty=blast_penalty,
        max_hsps=blast_max_hsps,
        dust=blast_dust,
    )
    blast_df = runner.get_dataframe()
    _log_df_memory("blast_df", blast_df)
    _log_process_memory()

    if blast_df.empty:
        logger.warning("BLAST returned no hits.")
        return

    # 4. Annotate Results
    target_map = getattr(panel, "primer_target_map", {})
    annotator = BlastResultsAnnotator(blast_df, target_map=target_map)
    annotator.build_annotation_dict(
        length_threshold=length_threshold,
        evalue_threshold=evalue_threshold,
        max_mismatches=max_mismatches,
        three_prime_tolerance=three_prime_tolerance,
    )
    annotator.add_annotations()

    # 5. Find Off-Target Amplicons
    bound_df = annotator.get_predicted_bound()
    _log_df_memory("bound_df", bound_df)

    if max_bound_per_primer is not None:
        pre_cap = len(bound_df)
        bound_df = (
            bound_df.sort_values("evalue")
            .groupby("qseqid", observed=True)
            .head(max_bound_per_primer)
            .reset_index(drop=True)
        )
        if len(bound_df) < pre_cap:
            n_capped = (
                bound_df.groupby("qseqid", observed=True)
                .size()
                .eq(max_bound_per_primer)
                .sum()
            )
            logger.info(
                f"Per-primer cap ({max_bound_per_primer}): "
                f"{pre_cap:,} → {len(bound_df):,} bound hits "
                f"({n_capped} primer(s) capped)"
            )
            _log_df_memory("bound_df (after cap)", bound_df)

    # Free the full BLAST results — only bound_df is needed from here on
    del blast_df, runner, annotator
    _log_process_memory()

    # 6. Build reverse lookup: (f_id, r_id) -> [(junction, pair), ...]
    seq_to_id = panel.unique_primer_map
    id_to_pairs: dict[tuple[str, str], list[tuple]] = {}
    for junction in panel.junctions:
        for pair in junction.primer_pairs:
            f_id = seq_to_id.get(pair.forward.seq)
            r_id = seq_to_id.get(pair.reverse.seq)
            if not f_id or not r_id:
                continue
            pair.specificity_checked = True
            pair.off_target_products = []
            pair.on_target_detected = False
            id_to_pairs.setdefault((f_id, r_id), []).append((junction, pair))
            id_to_pairs.setdefault((r_id, f_id), []).append((junction, pair))

    # 7. Process amplicons per chromosome to bound peak memory
    finder = AmpliconFinder(bound_df, target_map=target_map)
    total_amplicons = 0
    for _chrom, chrom_amp_df in finder.find_amplicons_by_chrom(
        max_size_bp=max_amplicon_size
    ):
        total_amplicons += len(chrom_amp_df)
        for (f_id, r_id), group_df in chrom_amp_df.groupby(
            ["F_primer", "R_primer"], observed=True
        ):
            for junction, pair in id_to_pairs.get((f_id, r_id), []):
                on_mask = _is_on_target_vec(
                    group_df, junction, pair, tolerance=ontarget_tolerance
                )
                if on_mask.any():
                    pair.on_target_detected = True
                off_df = group_df[~on_mask]
                if not off_df.empty:
                    pair.off_target_products.extend(off_df.to_dict("records"))

    if total_amplicons == 0:
        logger.info("No amplicons found (on- or off-target).")
    else:
        logger.info(f"Processed {total_amplicons:,} amplicons across all chromosomes")

    # Post-loop: log per-pair debug info and per-junction warnings
    for junction in panel.junctions:
        n_checked = sum(1 for p in junction.primer_pairs if p.specificity_checked)
        n_missing = 0
        for pair in junction.primer_pairs:
            if not pair.specificity_checked:
                continue
            if not pair.on_target_detected:
                n_missing += 1
                logger.debug(
                    f"Pair {pair.pair_id}: on-target amplicon not detected by BLAST."
                )
            if pair.off_target_products:
                logger.debug(
                    f"Pair {pair.pair_id} has {len(pair.off_target_products)} "
                    "off-target products."
                )
        if n_missing > 0:
            logger.warning(
                f"Junction {junction.name}: {n_missing}/{n_checked} pairs "
                "have no on-target amplicon detected by BLAST. "
                "Check BLAST sensitivity or junction coordinates."
            )

    _log_process_memory()
    logger.info("Specificity check complete.")


def filter_offtarget_pairs(panel: MultiplexPanel) -> tuple[int, list[str]]:
    """Remove primer pairs that have off-target amplicons.

    For each junction, pairs with no off-target products are kept.  If *all*
    pairs for a junction have off-targets, the single pair with the fewest
    off-target products is retained so the junction is not lost from the panel.

    Parameters
    ----------
    panel : MultiplexPanel
        Panel whose junctions have already been through ``run_specificity_check``.

    Returns
    -------
    tuple[int, list[str]]
        Total number of primer pairs removed across all junctions, and a list
        of junction names where the fallback (all pairs had off-targets) was
        triggered.
    """
    total_removed = 0
    fallback_junctions: list[str] = []

    for junction in panel.junctions:
        if not junction.primer_pairs:
            continue

        clean = [p for p in junction.primer_pairs if len(p.off_target_products) == 0]
        dirty = [p for p in junction.primer_pairs if len(p.off_target_products) > 0]

        if not dirty:
            continue

        if clean:
            removed = len(dirty)
            junction.primer_pairs = clean
            logger.info(
                f"Junction {junction.name}: removed {removed} primer pair(s) "
                f"with off-target products, {len(clean)} clean pair(s) remain"
            )
        else:
            # All pairs have off-targets — keep all with the fewest
            min_ot = min(len(p.off_target_products) for p in junction.primer_pairs)
            best_pairs = [
                p for p in junction.primer_pairs if len(p.off_target_products) == min_ot
            ]
            removed = len(junction.primer_pairs) - len(best_pairs)
            junction.primer_pairs = best_pairs
            fallback_junctions.append(junction.name)
            logger.warning(
                f"Junction {junction.name}: all pairs have off-target products; "
                f"keeping {len(best_pairs)} pair(s) with fewest "
                f"off-targets={min_ot}"
            )

        total_removed += removed

    return total_removed, fallback_junctions


def _is_on_target_vec(df, junction, pair, tolerance: int = 5):
    """Vectorized on-target classification for a DataFrame of amplicons."""
    design_start = getattr(junction, "design_start", None) or 0
    expected_fwd = design_start + pair.forward.start
    expected_rev = design_start + pair.reverse.start + pair.reverse.length - 1
    return (
        (df["chrom"] == junction.chrom)
        & ((df["F_start"] - expected_fwd).abs() <= tolerance)
        & ((df["R_start"] - expected_rev).abs() <= tolerance)
    )


def _is_on_target(prod: dict, junction, pair, tolerance: int = 5) -> bool:
    """Check if a BLAST amplicon overlaps the intended target region.

    Compares the BLAST hit genomic coordinates against the expected
    primer binding positions derived from the junction's design region.

    Coordinate conventions
    ----------------------
    - ``design_start`` is the 1-based genomic coordinate of the first base
      of the design region sequence.
    - ``pair.forward.start`` and ``pair.reverse.start`` are 0-based indices
      into the design region string; the forward primer's leftmost position
      is ``design_start + forward.start`` (1-based).
    - BLAST ``F_start`` (plus-strand ``sstart``) = leftmost (5') genomic
      position of the forward primer — matches ``design_start + forward.start``.
    - BLAST ``R_start`` (minus-strand ``sstart``) = **rightmost** (5') genomic
      position of the reverse primer — matches
      ``design_start + reverse.start + reverse.length - 1``.
    """
    if prod["chrom"] != junction.chrom:
        return False

    design_start = getattr(junction, "design_start", None) or 0

    # Forward primer: BLAST sstart (plus strand) = leftmost = 5' end
    expected_fwd_start = design_start + pair.forward.start

    # Reverse primer: BLAST sstart (minus strand) = rightmost = 5' end of
    # the reverse primer.  pair.reverse.start is the leftmost (3') position,
    # so the rightmost = start + length - 1.
    expected_rev_start = design_start + pair.reverse.start + pair.reverse.length - 1

    # Allow small tolerance for BLAST coordinate alignment differences.
    # Default 5 bp absorbs minor alignment shifts while still distinguishing
    # on-target hits from nearby off-target loci.
    fwd_match = abs(prod["F_start"] - expected_fwd_start) <= tolerance
    rev_match = abs(prod["R_start"] - expected_rev_start) <= tolerance

    return fwd_match and rev_match
