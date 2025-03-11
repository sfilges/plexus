import os
import primer3
#import pysam
import itertools
from Bio import SeqIO
from Bio.Seq import Seq
from src.hello import say_hello_to
import pandas as pd
import json

say_hello_to("Stefan")


# Load configuration files
design_parameters = load_parameters("conf/design.json")
pcr_conditions = load_parameters("conf/pcr.config")
penalties = load_parameters("conf/penalties.config")

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

#------// Workflow //------
# 1) Prepare inputs
# 2) Design single-plex solutions for each target
# 3) Check the specificity of each target
# 4) Pick the optimal multiplex solution
# 5) Output multiplex

#-----// (1) Prepare inputs //---------
# Get the mutation list (junctions) around which to design primer pairs

# Target genome
panel_name  = "test_panel"
genome      = "hg38"
sample_type = "cfDNA"
protocol    = "simsen"

# Define "tail-sequences" or adapters to be appended to each FP or RP, respectively

fp_tail = ""
rp_tail = ""

# Define and merge junctions
# Each junction has a le
# ft (j_min) and right (j_max) position
j_min = 123456
j_max = 123456

left_primer_region, right_primer_region = calculate_primer_design_regions(j_min, j_max, design_parameters)

# Import genome reference
# import data/genome.fa

# Use pysam to get design region from reference fasta
# https://pysam.readthedocs.io/en/latest/api.html
