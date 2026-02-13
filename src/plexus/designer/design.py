# ================================================================================
# Two functions for primer design:
#       - simsen_design_primers: Custom function
#       - primer3py_design_primers: Uses primer3-py bindings for primer3
#       - primer3_design_primers: Runs local primer3 (Not currently implemented)
# ================================================================================

import json
import os
import warnings

import primer3
from loguru import logger

from plexus.designer.multiplexpanel import (
    MultiplexPanel,
    PrimerDesigns,
)
from plexus.designer.thal import (
    calculate_single_primer_thermodynamics,
    calculate_single_primer_thermodynamics_parallel,
)
from plexus.utils.utils import generate_kmers, reverse_complement

# ================================================================================
# Wrapper function for different primer design algorithms.
# ================================================================================


def design_primers(
    panel: MultiplexPanel, method: str = "simsen", parallel: bool = False
) -> MultiplexPanel:
    """
    Wrapper function to call the various design algorithms.

    Args:
        panel: An instantiated MultiplexPanel object created with panel_factory.
        method: Design algorithm to use; defaults to "simsen".
            - "simsen": Custom thermodynamic-aware design (recommended, production-ready).
            - "primer3py": EXPERIMENTAL — Uses primer3-py bindings. Does not yet
              produce PrimerPair objects, so downstream multiplex optimization will
              receive no candidates. Under active development.
            - "primer3": NOT IMPLEMENTED — Placeholder for local primer3 execution.
        parallel: Boolean. If true, uses parallelized functions. Default is False.

    Returns:
        A MultiplexPanel object with primer designs.
    """
    if method == "simsen":
        return design_multiplex_primers(panel, parallel=parallel)
    if method == "primer3py":
        logger.warning(
            "The 'primer3py' design method is experimental and does not yet "
            "produce PrimerPair objects. Downstream multiplex optimization will "
            "receive no candidates. Use 'simsen' for a complete pipeline run."
        )
        return primer3py_design_primers(panel)
    if method == "primer3":
        raise NotImplementedError(
            "The 'primer3' design method is not yet implemented. "
            "Use 'simsen' (recommended) or 'primer3py' (experimental)."
        )
    raise ValueError(f"Unknown design method: {method}")


# ================================================================================
# Custom primer design function based on primer3
# ================================================================================


def design_multiplex_primers(
    panel: MultiplexPanel, parallel: bool = False
) -> MultiplexPanel:
    """
    A function that picks individual primers left and right of the provided junctions.

    Args:
        panel: A MultiplexPanel object with loaded junctions

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
    max_N = singleplex.primer_max_n
    min_gc = singleplex.primer_min_gc
    max_gc = singleplex.primer_max_gc
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

    # =============================================
    # Start designing
    # =============================================

    j = 0
    # Pick primers to left and right of each junction
    for junction in panel.junctions:
        logger.info(
            f"#========================// {junction.name} //============================#"
        )
        j += 1

        # Get the design regions left and right of the junction, respectively.
        left_region = junction.design_region[1 : junction.jmin_coordinate]

        # Should the right region be reverse complemented here?
        right_region = reverse_complement(
            junction.design_region[junction.jmax_coordinate : junction.junction_length]
        )

        # This generates candidate primers, returns a list of Primer class objects,
        # which have been filtered based on length and basic sequence properties.
        # The offset is necessary as the start positions of the generated kmer
        # (on the design_region) are starting at position jmax_coordinate for the
        # reverse primer and at 1 for the forward primer.
        for orientation, region, offset in [
            ("forward", left_region, 1),
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
            )
            if orientation == "forward":
                left_kmers = kmers
            else:
                right_kmers = kmers

        # Warn or raise error if too few kmers found.
        if len(left_kmers) < 100:
            msg = "Fewer than 100 left kmers found."
            logger.warning(msg)
            warnings.warn(msg, stacklevel=2)
        elif len(left_kmers) == 0:
            logger.error("No left kmers found.")
            raise ValueError("No left kmers found.")
        if len(right_kmers) < 100:
            logger.warning("Fewer than 100 right kmers found.")
            warnings.warn("Fewer than 100 right kmers found.", stacklevel=2)
        elif len(right_kmers) == 0:
            logger.error("No right kmers found.")
            raise ValueError("No right kmers found.")

        # Calculate thermodynamic properties of candidate primers and remove low quality primers based on config.
        if parallel:
            # Run both left and right primers at the same time
            (
                left_primers,
                right_primers,
            ) = calculate_single_primer_thermodynamics_parallel(
                left_kmers=left_kmers,
                right_kmers=right_kmers,
                config=panel.config,
            )
        else:
            left_primers, left_eval_string = calculate_single_primer_thermodynamics(
                left_kmers, panel.config, orientation="left"
            )
            right_primers, right_eval_string = calculate_single_primer_thermodynamics(
                right_kmers, panel.config, orientation="right"
            )

        # TODO: Actually generate table of primers to output
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
        )

    logger.info(
        f"Finished designing primers for {j} junctions in panel {panel.panel_name}."
    )
    return panel


# ================================================================================
# EXPERIMENTAL: primer3-py bindings design method
#
# This method uses the primer3-py Python bindings to design primer pairs.
# It currently stores raw primer3 output on each junction but does NOT create
# PrimerPair objects, so downstream steps (multiplex optimization, BLAST) will
# have no candidates to work with. Under active development.
# ================================================================================


def primer3py_design_primers(
    panel: MultiplexPanel, thal: int = 1, save_designs: bool = False, outdir: str = None
) -> MultiplexPanel:
    """
    EXPERIMENTAL: Run primer3 design on the junctions and retain the top primers.

    .. warning::
        This method does not yet produce ``PrimerPair`` objects. The raw primer3
        output is stored as ``junction.primer3_designs`` but the ``junction.primer_pairs``
        list remains empty. Downstream multiplex optimization will therefore receive
        no candidates.

    Returns a MultiplexPanel object.
    """
    # Get config sections for easier access
    singleplex = panel.config.singleplex_design_parameters
    pair_params = panel.config.primer_pair_parameters
    pcr = panel.config.pcr_conditions

    num_expected = singleplex.PRIMER_NUM_RETURN

    for junction in panel.junctions:
        min_product_length = (
            2 * singleplex.primer_min_length
            + junction.jmax_coordinate
            - junction.jmin_coordinate
        )
        max_product_length = pair_params.PRIMER_PRODUCT_MAX_SIZE

        # Set global design arguments for primer3. For details check the [manual](https://primer3.org/manual.html)
        global_args = {
            "PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT": thal,
            "PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT": thal,
            "PRIMER_TM_FORMULA": 1,
            "PRIMER_SALT_CORRECTIONS": 1,
            "PRIMER_NUM_RETURN": num_expected,
            "PRIMER_OPT_SIZE": singleplex.PRIMER_OPT_SIZE,
            "PRIMER_MIN_SIZE": singleplex.primer_min_length,
            "PRIMER_MAX_SIZE": singleplex.primer_max_length,
            "PRIMER_PRODUCT_OPT_SIZE": pair_params.PRIMER_PRODUCT_OPT_SIZE,
            "PRIMER_OPT_TM": singleplex.PRIMER_OPT_TM,
            "PRIMER_MIN_TM": singleplex.PRIMER_MIN_TM,
            "PRIMER_MAX_TM": singleplex.PRIMER_MAX_TM,
            "PRIMER_PAIR_MAX_DIFF_TM": pair_params.PRIMER_PAIR_MAX_DIFF_TM,
            "PRIMER_SALT_MONOVALENT": pcr.mv_concentration,
            "PRIMER_SALT_DIVALENT": pcr.dv_concentration,
            "PRIMER_DNA_CONC": pcr.primer_concentration,
            "PRIMER_DNTP_CONC": pcr.dntp_concentration,
            "PRIMER_DMSO_CONC": pcr.dmso_concentration,
            "PRIMER_FORMAMIDE_CONC": pcr.formamide_concentration,
            "PRIMER_OPT_GC_PERCENT": singleplex.PRIMER_OPT_GC_PERCENT,
            "PRIMER_MIN_GC": singleplex.primer_min_gc,
            "PRIMER_MAX_GC": singleplex.primer_max_gc,
            "PRIMER_GC_CLAMP": singleplex.primer_gc_clamp,
            "PRIMER_MAX_END_GC": 4,
            "PRIMER_MAX_END_STABILITY": singleplex.PRIMER_MAX_END_STABILITY,
            "PRIMER_MAX_POLY_X": singleplex.primer_max_poly_x,
            "PRIMER_MAX_NS_ACCEPTED": singleplex.primer_max_n,
            "PRIMER_MAX_SELF_ANY": 8.0,
            "PRIMER_MAX_SELF_END": 3.0,
            "PRIMER_PAIR_MAX_COMPL_ANY": 8.0,
            "PRIMER_PAIR_MAX_COMPL_END": 8.0,
            "PRIMER_MAX_SELF_ANY_TH": singleplex.PRIMER_MAX_SELF_ANY_TH,
            "PRIMER_MAX_SELF_END_TH": singleplex.PRIMER_MAX_SELF_END_TH,
            "PRIMER_PAIR_MAX_COMPL_ANY_TH": 47.0,  # Default primer3 value
            "PRIMER_PAIR_MAX_COMPL_END_TH": 47.0,  # Default primer3 value
            "PRIMER_MAX_HAIRPIN_TH": singleplex.PRIMER_MAX_HAIRPIN_TH,
            "PRIMER_MAX_TEMPLATE_MISPRIMING_TH": singleplex.PRIMER_MAX_TEMPLATE_MISPRIMING_TH,
            "PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH": 70.0,  # Default primer3 value
            "PRIMER_WT_SIZE_LT": singleplex.PRIMER_WT_SIZE_LT,
            "PRIMER_WT_SIZE_GT": singleplex.PRIMER_WT_SIZE_GT,
            "PRIMER_WT_GC_PERCENT_LT": singleplex.PRIMER_WT_GC_PERCENT_LT,
            "PRIMER_WT_GC_PERCENT_GT": singleplex.PRIMER_WT_GC_PERCENT_GT,
            "PRIMER_WT_SELF_ANY": 5.0,
            "PRIMER_WT_SELF_ANY_TH": singleplex.PRIMER_WT_SELF_ANY_TH,
            "PRIMER_WT_HAIRPIN_TH": singleplex.PRIMER_WT_HAIRPIN_TH,
            "PRIMER_WT_TEMPLATE_MISPRIMING_TH": 1.0,
            "PRIMER_PAIR_WT_COMPL_ANY": 5.0,
            "PRIMER_PAIR_WT_COMPL_ANY_TH": 5.0,
            "PRIMER_PAIR_WT_COMPL_END_TH": 5.0,
            "PRIMER_PAIR_WT_TEMPLATE_MISPRIMING_TH": 1.0,
            "PRIMER_PAIR_WT_PRODUCT_SIZE_LT": pair_params.PRIMER_PAIR_WT_PRODUCT_SIZE_LT,
            "PRIMER_PAIR_WT_PRODUCT_SIZE_GT": pair_params.PRIMER_PAIR_WT_PRODUCT_SIZE_GT,
            "PRIMER_PRODUCT_SIZE_RANGE": [[min_product_length, max_product_length]],
        }

        # This tag allows detailed specification of possible locations of left and right primers in primer pairs.
        # The associated value must be a semicolon-separated list of
        # <left_start>,<left_length>,<right_start>,<right_length>

        ok_regions_list = [
            1,
            junction.jmin_coordinate,
            junction.jmax_coordinate,
            len(junction.design_region) - junction.jmax_coordinate,
        ]

        # Set the sequence arguments for primer3
        seq_args = {
            "SEQUENCE_ID": junction.name,
            "SEQUENCE_TEMPLATE": junction.design_region,
            "SEQUENCE_PRIMER_PAIR_OK_REGION_LIST": [ok_regions_list],
        }

        logger.info(
            f"Running primer3 design engine for target junction: {junction.name}"
        )

        # Store primer3 designs in junction object
        junction.primer3_designs = primer3.bindings.design_primers(
            seq_args=seq_args, global_args=global_args
        )

        # Log complete primer3 output as a single row.
        # logger.info(junction.primer3_designs)

        # Formated primer3 log:
        logger.info(f"Design region size (nt): {len(junction.design_region)}")
        logger.info(
            f"Okay regions (<left_start>,<left_length>,<right_start>,<right_length>): {ok_regions_list}"
        )
        logger.info(
            f"LEFT EXPLAIN ({junction.name}): {junction.primer3_designs['PRIMER_LEFT_EXPLAIN']}"
        )
        logger.info(
            f"RIGHT EXPLAIN ({junction.name}): {junction.primer3_designs['PRIMER_RIGHT_EXPLAIN']}"
        )
        logger.info(
            f"PAIR EXPLAIN ({junction.name}): {junction.primer3_designs['PRIMER_PAIR_EXPLAIN']}"
        )
        logger.info(
            f"PRIMER_LEFT_RETURNED ({junction.name}): {junction.primer3_designs['PRIMER_LEFT_NUM_RETURNED']}"
        )
        logger.info(
            f"PRIMER_RIGHT_RETURNED ({junction.name}): {junction.primer3_designs['PRIMER_RIGHT_NUM_RETURNED']}"
        )
        logger.info(
            f"PRIMER_PAIR_RETURNED ({junction.name}): {junction.primer3_designs['PRIMER_PAIR_NUM_RETURNED']}"
        )

        if junction.primer3_designs["PRIMER_PAIR_NUM_RETURNED"] == 0:
            logger.warning(
                f"No suitable primer pairs found for junction: {junction.name}"
            )
            tm_range = (
                f"Tm range: {singleplex.PRIMER_MIN_TM} - {singleplex.PRIMER_MAX_TM}"
            )
            logger.info(f"{tm_range}. Consider increasing the Tm range.")
        elif junction.primer3_designs["PRIMER_PAIR_NUM_RETURNED"] < num_expected:
            logger.warning(f"Fewer primer pairs found than desired: {num_expected}")

        # Save design to file
        if save_designs:
            outfile = f"{junction.name}_primer3_out.json"
            # Determine the output directory
            if outdir:
                output_dir = outdir
            else:
                output_dir = os.path.expanduser("~")  # User's home directory

            # Ensure the output directory exists
            os.makedirs(output_dir, exist_ok=True)

            # Construct the full output path
            outfile_path = os.path.join(output_dir, outfile)

            # Save the file
            with open(outfile_path, "w") as f:
                json.dump(junction.primer3_designs, f, indent=4)

            logger.info(f"Saving primer3 designs to file: {outfile_path}")

    return panel


# ================================================================================
# NOT IMPLEMENTED: Local primer3 execution
#
# Placeholder for running a local primer3 binary. This requires primer3 to be
# installed on the system PATH and is not yet implemented.
# ================================================================================


def primer3_design_primers(panel: MultiplexPanel) -> MultiplexPanel:
    """
    NOT IMPLEMENTED: Run primer3 locally.

    Raises ``NotImplementedError``. Use ``design_method="simsen"`` instead.
    """
    raise NotImplementedError(
        "Local primer3 execution is not yet implemented. "
        "Use 'simsen' (recommended) or 'primer3py' (experimental)."
    )
