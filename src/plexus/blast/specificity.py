import os

from loguru import logger

from plexus.blast.annotator import BlastResultsAnnotator
from plexus.blast.blast_runner import BlastRunner
from plexus.blast.offtarget_finder import AmpliconFinder
from plexus.designer.multiplexpanel import MultiplexPanel


def run_specificity_check(panel: MultiplexPanel, work_dir: str, genome_fasta: str):
    """
    Run BLAST on all candidate primers in the panel to check for specificity
    and identify off-target amplicons.

    This assumes that the panel has already aggregated primers and assigned unique IDs.

    Args:
        panel: The MultiplexPanel object containing junctions and primer designs.
        work_dir: Directory to store temporary BLAST files.
        genome_fasta: Path to the reference genome FASTA file.
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
    blast_archive = os.path.join(work_dir, "blast_archive")
    blast_table = os.path.join(work_dir, "blast_table.txt")

    runner = BlastRunner(input_fasta, genome_fasta)
    runner.create_database()
    runner.run(output_archive=blast_archive, word_size=11)
    runner.reformat_output_as_table(blast_table)

    blast_df = runner.get_dataframe()

    if blast_df.empty:
        logger.warning("BLAST returned no hits.")
        return

    # 4. Annotate Results
    target_map = getattr(panel, "primer_target_map", {})
    annotator = BlastResultsAnnotator(blast_df, target_map=target_map)
    annotator.build_annotation_dict(length_threshold=15, evalue_threshold=10)
    annotator.add_annotations()

    # 5. Find Off-Target Amplicons
    bound_df = annotator.get_predicted_bound()
    finder = AmpliconFinder(bound_df, target_map=target_map)
    finder.find_amplicons(max_size_bp=2000)

    all_amplicons_df = finder.amplicon_df

    if all_amplicons_df is None or all_amplicons_df.empty:
        logger.info("No amplicons found (on- or off-target).")
        return

    # 6. Map results back to PrimerPairs
    amplicon_map = {}
    for _, row in all_amplicons_df.iterrows():
        f_id = row["F_primer"]  # SEQ_X
        r_id = row["R_primer"]  # SEQ_Y

        if (f_id, r_id) not in amplicon_map:
            amplicon_map[(f_id, r_id)] = []
        amplicon_map[(f_id, r_id)].append(row.to_dict())

    # Use the panel's unique map to look up IDs
    seq_to_id = panel.unique_primer_map

    for junction in panel.junctions:
        n_checked = 0
        n_missing_on_target = 0

        for pair in junction.primer_pairs:
            f_seq = pair.forward.seq
            r_seq = pair.reverse.seq

            f_id = seq_to_id.get(f_seq)
            r_id = seq_to_id.get(r_seq)

            if not f_id or not r_id:
                continue

            pair.specificity_checked = True
            n_checked += 1

            potential_products = amplicon_map.get((f_id, r_id), [])

            off_targets = []
            on_targets = []
            for prod in potential_products:
                if _is_on_target(prod, junction, pair):
                    on_targets.append(prod)
                else:
                    off_targets.append(prod)

            pair.off_target_products = off_targets
            pair.on_target_detected = len(on_targets) > 0

            if not pair.on_target_detected:
                n_missing_on_target += 1
                logger.debug(
                    f"Pair {pair.pair_id}: on-target amplicon not detected by BLAST."
                )

            if off_targets:
                logger.debug(
                    f"Pair {pair.pair_id} has {len(off_targets)} off-target products."
                )

        if n_missing_on_target > 0:
            logger.warning(
                f"Junction {junction.name}: {n_missing_on_target}/{n_checked} pairs "
                "have no on-target amplicon detected by BLAST. "
                "Check BLAST sensitivity or junction coordinates."
            )

    logger.info("Specificity check complete.")


def _is_on_target(prod: dict, junction, pair) -> bool:
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

    # Allow small tolerance for BLAST coordinate alignment differences
    tolerance = 5  # bp
    fwd_match = abs(prod["F_start"] - expected_fwd_start) <= tolerance
    rev_match = abs(prod["R_start"] - expected_rev_start) <= tolerance

    return fwd_match and rev_match
