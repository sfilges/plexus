# ================================================================================
# Primer design module — "plexus" design algorithm
# ================================================================================

from __future__ import annotations

import warnings

from loguru import logger

from plexus.designer.multiplexpanel import (
    MultiplexPanel,
    PrimerDesigns,
)
from plexus.designer.thal import (
    calculate_single_primer_thermodynamics,
)
from plexus.utils.utils import generate_kmers, reverse_complement

# ================================================================================
# Wrapper function for different primer design algorithms.
# ================================================================================


def design_primers(
    panel: MultiplexPanel,
    method: str = "plexus",
    on_junction_done: callable | None = None,
) -> MultiplexPanel:
    """
    Wrapper function to call the design algorithm.

    Args:
        panel: An instantiated MultiplexPanel object created with panel_factory.
        method: Design algorithm to use; defaults to "plexus".
        on_junction_done: Optional callback invoked after each junction is processed.

    Returns:
        A MultiplexPanel object with primer designs.
    """
    if method == "plexus":
        return design_multiplex_primers(panel, on_junction_done=on_junction_done)
    raise ValueError(f"Unknown design method: {method}")


# ================================================================================
# Custom primer design function based on primer3
# ================================================================================


def design_multiplex_primers(
    panel: MultiplexPanel,
    on_junction_done: callable | None = None,
) -> MultiplexPanel:
    """
    A function that picks individual primers left and right of the provided junctions.

    Args:
        panel: A MultiplexPanel object with loaded junctions
        on_junction_done: Optional callback invoked after each junction is processed.

    Returns:
        A MultiplexPanel object containing the left and right primer designs for each junction.
    """

    # =============================================
    # Common parameters for all junctions
    # =============================================

    # Get config sections for easier access
    singleplex = panel.config.singleplex_design_parameters
    pair_params = panel.config.primer_pair_parameters

    # Single primer parameters
    min_primer_length = singleplex.primer_min_length
    max_primer_length = singleplex.primer_max_length
    max_poly_X = singleplex.primer_max_poly_x
    max_poly_gc = singleplex.primer_max_poly_gc
    max_N = singleplex.primer_max_n
    min_gc = singleplex.primer_min_gc
    max_gc = singleplex.primer_max_gc
    gc_clamp = singleplex.primer_gc_clamp
    min_region_length = 2 * max_primer_length

    # Primer pair parameters
    min_amplicon_length = (
        pair_params.PRIMER_PRODUCT_MIN_INSERT_SIZE + 2 * min_primer_length
    )
    max_amplicon_length = pair_params.PRIMER_PRODUCT_MAX_SIZE
    max_primer_tm_difference = pair_params.PRIMER_PAIR_MAX_DIFF_TM
    pair_product_opt_size = pair_params.PRIMER_PRODUCT_OPT_SIZE

    # Weights for primer pair penalty calculations
    wt_pr_penalty = pair_params.PRIMER_PAIR_WT_PR_PENALTY
    wt_product_size_gt = pair_params.PRIMER_PAIR_WT_PRODUCT_SIZE_GT
    wt_product_size_lt = pair_params.PRIMER_PAIR_WT_PRODUCT_SIZE_LT
    wt_diff_tm = pair_params.PRIMER_PAIR_WT_DIFF_TM

    # Tail sequences (prepended to primers for dimer scoring and output)
    forward_tail = singleplex.forward_tail
    reverse_tail = singleplex.reverse_tail

    # =============================================
    # Start designing
    # =============================================

    # Pick primers to left and right of each junction
    for junction in panel.junctions:
        logger.info(
            f"#========================// {junction.name} //============================#"
        )

        try:
            # Get the design regions left and right of the junction, respectively.
            left_region = junction.design_region[0 : junction.jmin_coordinate]

            # Should the right region be reverse complemented here?
            right_region = reverse_complement(
                junction.design_region[
                    junction.jmax_coordinate : junction.junction_length
                ]
            )

            # This generates candidate primers, returns a list of Primer class objects,
            # which have been filtered based on length and basic sequence properties.
            # The offset is necessary as the start positions of the generated kmer
            # (on the design_region) are starting at position jmax_coordinate for the
            # reverse primer and at 0 for the forward primer.
            for orientation, region, offset in [
                ("forward", left_region, 0),
                ("reverse", right_region, junction.jmax_coordinate),
            ]:
                if len(region) < min_region_length:
                    error_msg = f"{orientation} primer region for junction {junction.name} is less than {min_region_length} (2x max primer length)."
                    logger.error(error_msg)
                    raise ValueError(error_msg)

                kmers = generate_kmers(
                    target_name=junction.name,
                    target_sequence=region,
                    orientation=orientation,
                    position_offset=offset,
                    k_min=min_primer_length,
                    k_max=max_primer_length,
                    max_poly_X=max_poly_X,
                    max_N=max_N,
                    min_gc=min_gc,
                    max_gc=max_gc,
                    gc_clamp=gc_clamp,
                    max_poly_gc=max_poly_gc,
                )
                if orientation == "forward":
                    left_kmers = kmers
                else:
                    right_kmers = kmers

            # Warn or raise error if too few kmers found.
            if len(left_kmers) == 0:
                logger.error("No left kmers found.")
                raise ValueError("No left kmers found.")
            if len(left_kmers) < 100:
                msg = "Fewer than 100 left kmers found."
                logger.warning(msg)
                warnings.warn(msg, stacklevel=2)

            if len(right_kmers) == 0:
                logger.error("No right kmers found.")
                raise ValueError("No right kmers found.")
            if len(right_kmers) < 100:
                msg = "Fewer than 100 right kmers found."
                logger.warning(msg)
                warnings.warn(msg, stacklevel=2)

            # Calculate thermodynamic properties of candidate primers and remove low quality primers based on config.
            left_primers, left_eval_string = calculate_single_primer_thermodynamics(
                left_kmers, panel.config, orientation="left"
            )
            right_primers, right_eval_string = calculate_single_primer_thermodynamics(
                right_kmers, panel.config, orientation="right"
            )

            primer_table = [left_primers, right_primers]

            eval_string = (
                f"LEFT PRIMERS: {left_eval_string}, RIGHT PRIMERS: {right_eval_string}"
            )

            # Save designs in junction object
            junction.primer_designs = PrimerDesigns(
                name=f"{junction.name}_designs",
                target=junction.name,
                design_region=junction.design_region,
                eval_string=eval_string,
                primer_table=primer_table,
            )

            logger.info(
                f"Primers retained post theromdynamic alignment: Left: {len(left_primers)}, right {len(right_primers)}"
            )

            # Find suitable primer pairs from initial designs and calculate primer pair penalties
            junction.primer_pairs = junction.find_primer_pairs(
                min_amplicon_length=min_amplicon_length,
                max_amplicon_length=max_amplicon_length,
                max_primer_tm_difference=max_primer_tm_difference,
                product_opt_size=pair_product_opt_size,
                wt_pr_penalty=wt_pr_penalty,
                wt_product_size_gt=wt_product_size_gt,
                wt_product_size_lt=wt_product_size_lt,
                wt_diff_tm=wt_diff_tm,
                forward_tail=forward_tail,
                reverse_tail=reverse_tail,
            )

            # Cap candidate pairs per junction to bound the search space
            max_pairs = pair_params.max_pairs_per_junction
            if max_pairs is not None and len(junction.primer_pairs) > max_pairs:
                original_count = len(junction.primer_pairs)
                junction.primer_pairs.sort(key=lambda p: p.pair_penalty)
                junction.primer_pairs = junction.primer_pairs[:max_pairs]
                logger.info(
                    f"Capped {junction.name} to {max_pairs} pairs "
                    f"(from {original_count})"
                )
        except Exception as e:
            logger.warning(
                f"Primer design failed for junction {junction.name}: {e}. Skipping."
            )
            junction.primer_pairs = []
            junction._design_error = str(e)
        finally:
            if on_junction_done:
                on_junction_done()

    # Separate failed junctions (no primer pairs)
    failed = [jn for jn in panel.junctions if not jn.primer_pairs]
    if failed:
        logger.warning(
            f"{len(failed)} junction(s) failed primer design and will be excluded."
        )
        for fj in failed:
            logger.warning(
                f"  - {fj.name}: {getattr(fj, '_design_error', 'no valid primer pairs')}"
            )
    panel.failed_junctions = failed
    panel.junctions = [jn for jn in panel.junctions if jn.primer_pairs]

    logger.info(
        f"Finished designing primers for {len(panel.junctions)} junctions in panel {panel.panel_name}."
    )
    return panel
