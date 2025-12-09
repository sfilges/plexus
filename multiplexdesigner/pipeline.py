# Import modules
from multiplexdesigner import calculate_primer_design_regions, extract_genomic_region, load_parameters, read_junction_coordinates 
from multiplexdesigner import primer3_design_primers

#------// Workflow //------
# 1) Prepare inputs
# 2) Design single-plex solutions for each target
# 3) Check the specificity of each target
# 4) Pick the optimal multiplex solution
# 5) Output multiplex

# Example usage:
# candidate_sequence = Seq("TTGTGAGTTTTTGAAATCTCTGTGA")
# hairpin_result = calculate_hairpin(candidate_sequence, pcr_conditions)
# tm = calculate_melting_temperature(candidate_sequence, pcr_conditions)
# homodimer_result = calculate_homodimer(candidate_sequence, pcr_conditions)
# heterodimer_result = calculate_heterodimer(
#     candidate_sequence, 
#     candidate_sequence.reverse_complement(), 
#     pcr_conditions
# )

# Load configuration files
design_parameters = load_parameters("conf/design.json")
pcr_conditions = load_parameters("conf/pcr.config")
penalties = load_parameters("conf/penalties.config")

#-----// (1) Prepare inputs //---------

# Target genome
panel_name  = "test_panel"
genome      = "hg38"
sample_type = "cfDNA"
protocol    = "simsen"

# Get the mutation list (junctions) around which to design primer pairs
junctions = read_junction_coordinates('./data/design_regions.csv')

# Import genome reference
# TODO If genome is not available, automatically download
if genome == "hg38":
    fasta_file  = "data/genome.fa"
else:
    print("Please select a valid genome.")

# Define "tail-sequences" or adapters to be appended to each FP or RP, respectively
# TODO Add default support for specific protocols.
if protocol == "simsen":
    fp_tail = "GGACACTCTTTCCCTACACGACGCTCTTCCGATCTAAAAAAAAAAAAAAAAAAAATGGGAAAGAGTGTCC"
    rp_tail = "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"
else:
    fp_tail = design_parameters['5_primer_tail']
    rp_tail = design_parameters['3_primer_tail']

# Design primers for each junction in the junction list.
for junction in junctions:
    all_designs = []

    # Define and merge junctions
    # Each junction has a le
    # ft (j_min) and right (j_max) position

    j_min = junction['jmin']
    j_max = junction['jmax']

    left_primer_region, right_primer_region = calculate_primer_design_regions(j_min, j_max, design_parameters)

    # Use pysam to get design region from reference fasta
    # https://pysam.readthedocs.io/en/latest/api.html
    sequence = extract_genomic_region(
        fasta_file = fasta_file,
        chromosome = junction['chrom'],  # or just "1" depending on your FASTA format
        start = left_primer_region[0],   # 0-based start position
        end = right_primer_region[1],    # 0-based end position (exclusive)
    )

    seq_args={
        'SEQUENCE_ID': junction['name'],
        'SEQUENCE_TEMPLATE': sequence,
        'SEQUENCE_INCLUDED_REGION': junction['junction']
        }

    global_args={
        'PRIMER_OPT_SIZE': 20,
        'PRIMER_PICK_INTERNAL_OLIGO': 1,
        'PRIMER_INTERNAL_MAX_SELF_END': 8,
        'PRIMER_MIN_SIZE': 18,
        'PRIMER_MAX_SIZE': 30,
        'PRIMER_OPT_TM': 60.0,
        'PRIMER_MIN_TM': 57.0,
        'PRIMER_MAX_TM': 63.0,
        'PRIMER_MIN_GC': 20.0,
        'PRIMER_MAX_GC': 80.0,
        'PRIMER_MAX_POLY_X': 100,
        'PRIMER_INTERNAL_MAX_POLY_X': 100,
        'PRIMER_SALT_MONOVALENT': 50.0,
        'PRIMER_DNA_CONC': 50.0,
        'PRIMER_MAX_NS_ACCEPTED': 0,
        'PRIMER_MAX_SELF_ANY': 12,
        'PRIMER_MAX_SELF_END': 3.0,
        'PRIMER_PAIR_MAX_COMPL_ANY': 12,
        'PRIMER_PAIR_MAX_COMPL_END': 8,
        'PRIMER_MAX_HAIRPIN_TH': 47.0,
        'PRIMER_PICK_ANYWAY': 1,
        'PRIMER_PRODUCT_SIZE_RANGE': [[40,100]]
    }

    # Run primer3 to get primer solutions for each junction
    # TODO how to handle no primers found?
    junction_primers = primer3_design_primers(seq_args=seq_args, global_args=global_args)
    # Collect all primers for a given target
    all_designs += junction_primers