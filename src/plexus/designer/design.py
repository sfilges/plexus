# ================================================================================
# Primer design module — "plexus" design algorithm
# ================================================================================

from __future__ import annotations

import warnings

from loguru import logger

from plexus.config import DesignerConfig
from plexus.designer.multiplexpanel import (
    Junction,
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
# Helper: run thermodynamic filtering + pair finding for a single attempt
# ================================================================================


def _attempt_pair_design(
    junction: Junction,
    left_kmers: list,
    right_kmers: list,
    config: DesignerConfig,
) -> tuple[list, list, list, str]:
    """Run thal filtering and pair-finding with the given config.

    Returns (left_primers, right_primers, primer_pairs, eval_string).
    """
    singleplex = config.singleplex_design_parameters
    pair_params = config.primer_pair_parameters

    left_primers, left_eval = calculate_single_primer_thermodynamics(
        left_kmers, config, orientation="left"
    )
    right_primers, right_eval = calculate_single_primer_thermodynamics(
        right_kmers, config, orientation="right"
    )

    eval_string = f"LEFT PRIMERS: {left_eval}, RIGHT PRIMERS: {right_eval}"

    # Save designs so find_primer_pairs can access design_region
    junction.primer_designs = PrimerDesigns(
        name=f"{junction.name}_designs",
        target=junction.name,
        design_region=junction.design_region,
        eval_string=eval_string,
        primer_table=[left_primers, right_primers],
    )

    logger.info(
        f"Primers retained post thermodynamic filtering: "
        f"Left: {len(left_primers)}, Right: {len(right_primers)}"
    )

    min_primer_length = singleplex.primer_min_length
    min_amplicon_length = (
        pair_params.PRIMER_PRODUCT_MIN_INSERT_SIZE + 2 * min_primer_length
    )

    pairs = junction.find_primer_pairs(
        min_amplicon_length=min_amplicon_length,
        max_amplicon_length=pair_params.PRIMER_PRODUCT_MAX_SIZE,
        max_primer_tm_difference=pair_params.PRIMER_PAIR_MAX_DIFF_TM,
        product_opt_size=pair_params.PRIMER_PRODUCT_OPT_SIZE,
        wt_pr_penalty=pair_params.PRIMER_PAIR_WT_PR_PENALTY,
        wt_product_size_gt=pair_params.PRIMER_PAIR_WT_PRODUCT_SIZE_GT,
        wt_product_size_lt=pair_params.PRIMER_PAIR_WT_PRODUCT_SIZE_LT,
        wt_diff_tm=pair_params.PRIMER_PAIR_WT_DIFF_TM,
        forward_tail=singleplex.forward_tail,
        reverse_tail=singleplex.reverse_tail,
    )

    return left_primers, right_primers, pairs, eval_string


def _compute_min_amplicon(left_primers: list, right_primers: list) -> int | None:
    """Compute the smallest possible amplicon from surviving primers."""
    if not left_primers or not right_primers:
        return None
    max_left_end = max(p.start + p.length for p in left_primers)
    min_right_start = min(p.start for p in right_primers)
    min_right_length = min(
        p.length for p in right_primers if p.start == min_right_start
    )
    return (
        min_right_start
        + min_right_length
        - (
            max_left_end
            - max(p.length for p in left_primers if p.start + p.length == max_left_end)
        )
    )


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

    config = panel.config
    singleplex = config.singleplex_design_parameters
    pair_params = config.primer_pair_parameters

    # K-mer generation parameters (unchanged across rescue tiers)
    min_primer_length = singleplex.primer_min_length
    max_primer_length = singleplex.primer_max_length
    max_poly_X = singleplex.primer_max_poly_x
    max_poly_gc = singleplex.primer_max_poly_gc
    max_N = singleplex.primer_max_n
    min_gc = singleplex.primer_min_gc
    max_gc = singleplex.primer_max_gc
    gc_clamp = singleplex.primer_gc_clamp
    min_region_length = 2 * max_primer_length

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

            right_region = reverse_complement(
                junction.design_region[
                    junction.jmax_coordinate : junction.junction_length
                ]
            )

            # Generate candidate k-mers (once — filters don't change across rescue tiers)
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

            # --- Attempt default design ---
            left_primers, right_primers, pairs, eval_string = _attempt_pair_design(
                junction,
                left_kmers,
                right_kmers,
                config,
            )
            rescue_tier = 0

            # --- Rescue loop for failed junctions ---
            if not pairs and config.enable_rescue:
                for tier_idx, tier in enumerate(config.rescue_tiers, start=1):
                    logger.info(
                        f"Rescue tier {tier_idx} for {junction.name}: {tier.description}"
                    )
                    rescue_config = config.apply_rescue_tier(tier_idx - 1)
                    left_primers, right_primers, pairs, eval_string = (
                        _attempt_pair_design(
                            junction,
                            left_kmers,
                            right_kmers,
                            rescue_config,
                        )
                    )
                    if pairs:
                        rescue_tier = tier_idx
                        logger.info(
                            f"Rescue tier {tier_idx} succeeded for {junction.name}: "
                            f"{len(pairs)} pairs found"
                        )
                        break
                    logger.info(
                        f"Rescue tier {tier_idx} for {junction.name}: " f"still 0 pairs"
                    )

            junction.primer_pairs = pairs
            junction._rescue_tier = rescue_tier

            # Tag all pairs with their rescue tier
            for pair in junction.primer_pairs:
                pair.rescue_tier = rescue_tier

            # Loose pre-filter cap to limit memory / BLAST workload
            max_pre = pair_params.max_pairs_pre_filter
            if max_pre is not None and len(junction.primer_pairs) > max_pre:
                original_count = len(junction.primer_pairs)
                junction.primer_pairs.sort(key=lambda p: p.pair_penalty)
                junction.primer_pairs = junction.primer_pairs[:max_pre]
                logger.info(
                    f"Pre-filter cap: {junction.name} to {max_pre} pairs "
                    f"(from {original_count})"
                )

            # Diagnostics for still-failed junctions
            if not junction.primer_pairs:
                min_amp = _compute_min_amplicon(left_primers, right_primers)
                max_allowed = pair_params.PRIMER_PRODUCT_MAX_SIZE
                if config.enable_rescue and config.rescue_tiers:
                    last_tier = config.rescue_tiers[-1]
                    if last_tier.PRIMER_PRODUCT_MAX_SIZE is not None:
                        max_allowed = last_tier.PRIMER_PRODUCT_MAX_SIZE

                if min_amp is not None:
                    junction._design_error = (
                        f"no valid primer pairs (closest achievable amplicon: "
                        f"{min_amp}bp, max allowed: {max_allowed}bp)"
                    )
                else:
                    junction._design_error = "no valid primer pairs"

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

    # Log rescued junctions
    rescued = [
        jn
        for jn in panel.junctions
        if jn.primer_pairs and getattr(jn, "_rescue_tier", 0) > 0
    ]
    if rescued:
        logger.info(f"{len(rescued)} junction(s) rescued with relaxed parameters.")
        for rj in rescued:
            tier_idx = rj._rescue_tier
            tier = config.rescue_tiers[tier_idx - 1]
            logger.info(f"  - {rj.name}: {tier.description}")

    panel.failed_junctions = failed
    panel.junctions = [jn for jn in panel.junctions if jn.primer_pairs]

    logger.info(
        f"Finished designing primers for {len(panel.junctions)} junctions in panel {panel.panel_name}."
    )
    return panel
