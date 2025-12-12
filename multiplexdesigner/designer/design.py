import json
import primer3
import warnings
from multiplexdesigner.designer.thal import calculate_single_primer_thermodynamics, calculate_single_primer_thermodynamics_parallel
from multiplexdesigner.designer.multiplexpanel import PrimerDesigns
from multiplexdesigner.utils.utils import generate_kmers, filter_kmers, reverse_complement

# ================================================================================
# Two functions for primer design:
#       - simsen_design_primers: Custom function
#       - primer3py_design_primers: Uses primer3-py bindings for primer3
#       - primer3_design_primers: Runs local primer3
# ================================================================================

# TODO: All designer options (custom, primer3py, primer3) should return data in the same format for downstram processing.

def design_primers(multiplexpanel, parallel: bool = False):
    """
    A function that picks individual primers left and right of the provided junctions.

    Args:
        multiplexpanel - A MultiPlex Panel object with loaded junctions
    
    Returns:
        A MultiplexPanel object containing the left and right primer designs for each junction.
    """
    
    # Pretty CLI interface
    #display_welcome()

    # Get logger and panel object
    logger = multiplexpanel[0]
    panel = multiplexpanel[1]

    # Common parameters for all junctions
    min_primer_length = panel.config['singleplex_design_parameters']['primer_min_length']
    max_primer_length = panel.config['singleplex_design_parameters']['primer_max_length']
    max_poly_X = panel.config['singleplex_design_parameters']['primer_max_poly_x']
    max_N = panel.config['singleplex_design_parameters']['primer_max_n']

    # Pick primers to left and right of each junction
    for junction in panel.junctions:
        logger.info("#=============================================================#")
        logger.info(junction.name)

        # Get the design regions left and right of the junction, respectively.
        left_region = junction.design_region[1:junction.jmin_coordinate]

        # TODO: Should the right region be reverse complemented here?
        right_region = reverse_complement(junction.design_region[junction.jmax_coordinate:junction.junction_length])

        # Check if design regions allow enough space to evaluate primers, expect around 50-60 bases for typical values of max_primer_length
        # TODO this should probably be checked upstream when design regions are calculated.
        if len(left_region) < 2 * max_primer_length:
            logger.error("Left primer region less than 2x max_primer_length")
            raise ValueError("Left primer region less than 2x max_primer_length")
        if len(right_region) < 2 * max_primer_length:
            logger.error("Right primer region less than 2x max_primer_length")
            raise ValueError("Right primer region less than 2x max_primer_length")

        # Generate k-mer solutions for left and right target regions (all possible sequences within min/max primer length)
        # TODO: these should generate Primer objects instead
        left_kmers = generate_kmers(left_region, min_primer_length, max_primer_length)
        right_kmers = generate_kmers(right_region, min_primer_length, max_primer_length)
        logger.info(f'kmers to evaluate: Left: {len(left_kmers)}; Right: {len(right_kmers)}')

        # Filter kmers based on basic sequence properties
        left_kmers = filter_kmers(left_kmers, max_poly_X, max_N)
        right_kmers = filter_kmers(right_kmers, max_poly_X, max_N)
        logger.info(f'kmers after basic filtering: Left: {len(left_kmers)}; Right: {len(right_kmers)}')

        # Warn or raise error if too few kmers found.
        # TODO adjust the values and improve handling.
        if len(left_kmers) < 100:
            msg = "Fewer than 100 left kmers found."
            logger.warning(msg)
            warnings.warn(msg)
        elif len(left_kmers) == 0:
            logger.error("No left kmers found.")
            raise ValueError("No left kmers found.")
        if len(right_kmers) < 100:
            logger.warning("Fewer than 100 right kmers found.") 
            warnings.warn("Fewer than 100 right kmers found.") 
        elif len(right_kmers) == 0:
            logger.error("No right kmers found.")
            raise ValueError("No right kmers found.")

        if parallel:
            # Run both left and right primers at the same time using 
            left_primers, right_primers = calculate_single_primer_thermodynamics_parallel(left_kmers = left_kmers, right_kmers = right_kmers, config = panel.config, logger = logger)
        else:
            left_primers, left_eval_string = calculate_single_primer_thermodynamics(left_kmers, panel.config, logger, orientation = "left")
            right_primers, right_eval_string = calculate_single_primer_thermodynamics(right_kmers, panel.config, logger, orientation = "right")

        # TODO: Actually generate table or primers to output
        primer_table = [left_primers, right_primers]

        eval_string = f"LEFT PRIMERS: {left_eval_string}, RIGHT PRIMERS: {right_eval_string}"

        junction.primer_designs = PrimerDesigns(
            name = f"{junction.name}_designs",
            target = junction.name,
            design_region = junction.design_region,
            eval_string = eval_string, 
            primer_table = primer_table
        )

        logger.info(f'primers retained post theromdynamic alignment: Left: {len(left_primers)}, right {len(right_primers)}')

    return(panel)


def primer3py_design_primers(panel_logger, thal = 1, save_designs = False):
    """
    Create a panel object and calculate junctions. Then run primer3 design on the junctions
    and retain the top primers for each juction.

    Returns a MultiplexPanel object.
    """
    panel = panel_logger[1]
    logger = panel_logger[0]

    num_expected = panel.config['singleplex_design_parameters']['PRIMER_NUM_RETURN']

    for junction in panel.junctions:

        min_product_length = 2 * panel.config['singleplex_design_parameters']['primer_min_length'] + junction.jmax_coordinate - junction.jmin_coordinate
        max_product_length = panel.config['singleplex_design_parameters']['max_amplicon_length']

        # Set global design arguments for primer3. For details check the [manual](https://primer3.org/manual.html)
        global_args={
            'PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT': thal,     # If set to 1, use thermodynamic values (parameters ending in "_TH")
            'PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT': thal,
            'PRIMER_TM_FORMULA': 1,
            'PRIMER_SALT_CORRECTIONS': 1,
            'PRIMER_NUM_RETURN': num_expected,
            'PRIMER_OPT_SIZE': panel.config['singleplex_design_parameters']['PRIMER_OPT_SIZE'],
            'PRIMER_MIN_SIZE': panel.config['singleplex_design_parameters']['primer_min_length'],
            'PRIMER_MAX_SIZE': panel.config['singleplex_design_parameters']['primer_max_length'],
            'PRIMER_PRODUCT_OPT_SIZE': panel.config['singleplex_design_parameters']['PRIMER_PRODUCT_OPT_SIZE'],
            'PRIMER_OPT_TM': panel.config['singleplex_design_parameters']['PRIMER_OPT_TM'],
            'PRIMER_MIN_TM': panel.config['singleplex_design_parameters']['PRIMER_MIN_TM'],
            'PRIMER_MAX_TM': panel.config['singleplex_design_parameters']['PRIMER_MAX_TM'],
            'PRIMER_PAIR_MAX_DIFF_TM': panel.config['singleplex_design_parameters']['primer_pair_max_tm_difference'],
            'PRIMER_SALT_MONOVALENT': panel.config['pcr_conditions']['mv_concentration'],           # The millimolar (mM) concentration of monovalent salt cations (usually KCl) in the PCR.
            'PRIMER_SALT_DIVALENT': panel.config['pcr_conditions']['dv_concentration'],             # The millimolar concentration of divalent salt cations (usually MgCl^(2+)) in the PCR. 
            'PRIMER_DNA_CONC': panel.config['pcr_conditions']['primer_concentration'],
            'PRIMER_DNTP_CONC': panel.config['pcr_conditions']['dntp_concentration'],               # Millimolar concentration of the sum of all deoxyribonucleotide triphosphates (e.g. 4x 0.2=0.8)
            'PRIMER_DMSO_CONC': panel.config['pcr_conditions']['dmso_concentration'],
            'PRIMER_FORMAMIDE_CONC': panel.config['pcr_conditions']['formamide_concentration'],
            'PRIMER_OPT_GC_PERCENT': panel.config['singleplex_design_parameters']['PRIMER_OPT_GC_PERCENT'],     # Optimum GC percent. This parameter influences primer selection only if PRIMER_WT_GC_PERCENT_GT or PRIMER_WT_GC_PERCENT_LT are non-0.
            'PRIMER_MIN_GC': panel.config['singleplex_design_parameters']['primer_min_gc'],
            'PRIMER_MAX_GC': panel.config['singleplex_design_parameters']['primer_max_gc'],
            'PRIMER_GC_CLAMP': 0,
            'PRIMER_MAX_END_GC': 4,
            'PRIMER_MAX_END_STABILITY': panel.config['singleplex_design_parameters']['PRIMER_MAX_END_STABILITY'],
            'PRIMER_MAX_POLY_X': panel.config['singleplex_design_parameters']['primer_max_poly_x'],
            'PRIMER_MAX_NS_ACCEPTED': 0,
            'PRIMER_MAX_SELF_ANY': 8.0,
            'PRIMER_MAX_SELF_END': 3.0,
            'PRIMER_PAIR_MAX_COMPL_ANY': 8.0,
            'PRIMER_PAIR_MAX_COMPL_END': 8.0,
            'PRIMER_MAX_SELF_ANY_TH': panel.config['singleplex_design_parameters']['PRIMER_MAX_SELF_ANY_TH'],
            'PRIMER_MAX_SELF_END_TH': panel.config['singleplex_design_parameters']['PRIMER_MAX_SELF_END_TH'],
            'PRIMER_PAIR_MAX_COMPL_ANY_TH': panel.config['singleplex_design_parameters']['PRIMER_PAIR_MAX_COMPL_ANY_TH'],
            'PRIMER_PAIR_MAX_COMPL_END_TH': panel.config['singleplex_design_parameters']['PRIMER_PAIR_MAX_COMPL_END_TH'],
            'PRIMER_MAX_HAIRPIN_TH': panel.config['singleplex_design_parameters']['PRIMER_MAX_HAIRPIN_TH'],
            'PRIMER_MAX_TEMPLATE_MISPRIMING_TH': panel.config['singleplex_design_parameters']['PRIMER_MAX_TEMPLATE_MISPRIMING_TH'],
            'PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH': panel.config['singleplex_design_parameters']['PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH'],
            'PRIMER_WT_SIZE_LT': panel.config['singleplex_design_parameters']['primer_length_penalty'],       # Penalty weight for primers shorter than PRIMER_OPT_SIZE.
            'PRIMER_WT_SIZE_GT': panel.config['singleplex_design_parameters']['primer_length_penalty'],       # Penalty weight for primers longer than PRIMER_OPT_SIZE.
            'PRIMER_WT_GC_PERCENT_LT': panel.config['singleplex_design_parameters']['PRIMER_WT_GC_PERCENT_LT'],                                          # Penaly weight for lower than opt GC content
            'PRIMER_WT_GC_PERCENT_GT': panel.config['singleplex_design_parameters']['PRIMER_WT_GC_PERCENT_GT'],                                          # Penaly weight for greater than opt GC content
            'PRIMER_WT_SELF_ANY': 5.0,
            'PRIMER_WT_SELF_ANY_TH': 5.0,
            'PRIMER_WT_HAIRPIN_TH': 1.0,
            'PRIMER_WT_TEMPLATE_MISPRIMING_TH': 1.0,
            'PRIMER_PAIR_WT_COMPL_ANY': 5.0,
            'PRIMER_PAIR_WT_COMPL_ANY_TH': 5.0,
            'PRIMER_PAIR_WT_COMPL_END_TH': 5.0,
            'PRIMER_PAIR_WT_TEMPLATE_MISPRIMING_TH': 1.0,
            'PRIMER_PAIR_WT_PRODUCT_SIZE_LT': panel.config['singleplex_design_parameters']['PRIMER_PAIR_WT_PRODUCT_SIZE_LT'],
            'PRIMER_PAIR_WT_PRODUCT_SIZE_GT': panel.config['singleplex_design_parameters']['PRIMER_PAIR_WT_PRODUCT_SIZE_GT'],
            'PRIMER_PRODUCT_SIZE_RANGE': [[min_product_length, max_product_length]]
        }

        # This tag allows detailed specification of possible locations of left and right primers in primer pairs.
        # The associated value must be a semicolon-separated list of
        # <left_start>,<left_length>,<right_start>,<right_length>

        ok_regions_list = [1,junction.jmin_coordinate,junction.jmax_coordinate,len(junction.design_region)-junction.jmax_coordinate]
    
        # Set the sequence arguments for primer3
        seq_args={
            'SEQUENCE_ID': junction.name,
            'SEQUENCE_TEMPLATE': junction.design_region,
            'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST': [ok_regions_list]
        }

        logger.info(f"Running primer3 design engine for target junction: {junction.name}")

        # Store primer3 designs in junction object
        junction.primer3_designs = primer3.bindings.design_primers(
            seq_args=seq_args,
            global_args=global_args
        )

        # Log complete primer3 output as a single row.
        #logger.info(junction.primer3_designs)
        
        # Formated primer3 log:
        logger.info(f"Design region size (nt): {len(junction.design_region)}")
        logger.info(f"Okay regions (<left_start>,<left_length>,<right_start>,<right_length>): {ok_regions_list}")
        logger.info(f"LEFT EXPLAIN ({junction.name}): {junction.primer3_designs["PRIMER_LEFT_EXPLAIN"]}")
        logger.info(f"RIGHT EXPLAIN ({junction.name}): {junction.primer3_designs["PRIMER_RIGHT_EXPLAIN"]}")
        logger.info(f"PAIR EXPLAIN ({junction.name}): {junction.primer3_designs["PRIMER_PAIR_EXPLAIN"]}")
        logger.info(f"PRIMER_LEFT_RETURNED ({junction.name}): {junction.primer3_designs["PRIMER_LEFT_NUM_RETURNED"]}")
        logger.info(f"PRIMER_RIGHT_RETURNED ({junction.name}): {junction.primer3_designs["PRIMER_RIGHT_NUM_RETURNED"]}")
        logger.info(f"PRIMER_PAIR_RETURNED ({junction.name}): {junction.primer3_designs["PRIMER_PAIR_NUM_RETURNED"]}")

        if junction.primer3_designs["PRIMER_PAIR_NUM_RETURNED"] == 0:
            logger.warning(f"No suitable primer pairs found for junction: {junction.name}")
            tm_range = f"Tm range: {panel.config['singleplex_design_parameters']['PRIMER_MIN_TM']} - {panel.config['singleplex_design_parameters']['PRIMER_MAX_TM']}"
            logger.info(f"{tm_range}. Consider increasing the Tm range.")
        elif junction.primer3_designs["PRIMER_PAIR_NUM_RETURNED"] < num_expected:
            logger.warning(f"Fewer primer pairs found than desired: {num_expected}")

        # Save design to file
        if save_designs:
            outfile = f"{junction.name}_primer3_out.json"
            with open(outfile, 'w') as f:
                json.dump(junction.primer3_designs, f, indent = 4)
            logger.info(f"Saving primer3 designs to file: {outfile}")

    return(panel)