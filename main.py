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

def calculate_primer_design_regions(j_min, j_max, design_parameters):
    """
    Calculate the forward and reverse primer design regions based on junction positions
    and design parameters.
    
    Args:
        j_min (int): Left position of the junction
        j_max (int): Right position of the junction
        design_parameters (dict): Dictionary containing design parameters including:
            - min_primer_length (int): Minimum primer length
            - max_primer_length (int): Maximum primer length
            - min_amplicon_gap (int): Minimum amplicon gap
            - max_amplicon_gap (int): Maximum amplicon gap
            - padding_bases (int): Number of padding bases
            
    Returns:
        tuple: Two lists containing [start, end] positions for left and right primer regions
    """
    
    # Extract parameters
    max_primer_length = design_parameters['max_primer_length']
    max_amplicon_gap = design_parameters['max_amplicon_gap']
    padding_bases = design_parameters['padding_bases']
    
    # Calculate forward primer region
    fp_start = j_max - max_primer_length - max_amplicon_gap + padding_bases
    fp_end = j_min + padding_bases
    
    # Calculate reverse primer region
    rp_start = j_max + padding_bases
    rp_end = j_min + max_primer_length + max_amplicon_gap - padding_bases
    
    # Return both regions as a tuple
    return [fp_start, fp_end], [rp_start, rp_end]

def depth_first_search_with_pruning(node, depth, alpha, beta, maximizing_player):
    # Base case: reached maximum depth or leaf node
    if depth == 0 or is_leaf_node(node):
        return evaluate_node(node)
    
    if maximizing_player:
        max_eval = float('-inf')
        for child in get_children(node):
            eval = depth_first_search_with_pruning(child, depth - 1, alpha, beta, False)
            max_eval = max(max_eval, eval)
            alpha = max(alpha, eval)
            
            # Pruning step: beta cut-off
            if beta <= alpha:
                break  # Cut off this branch
        return max_eval
    else:
        min_eval = float('inf')
        for child in get_children(node):
            eval = depth_first_search_with_pruning(child, depth - 1, alpha, beta, True)
            min_eval = min(min_eval, eval)
            beta = min(beta, eval)
            
            # Pruning step: alpha cut-off
            if beta <= alpha:
                break  # Cut off this branch
        return min_eval
    
# TODO how to handle overlapping junctions? There can
# be different cases:
# (1) Junctions very close together -> merge immediately
# (2) Moderately far apart junctions -> allow up to 14 nt 5' overlap
# (3) Identify junctions which cannot be merged -> report and/or ignore one of them

def read_junction_coordinates(csv_path):
    """
    Read junction coordinates from CSV file and convert to list of pairs.
    
    Args:
        csv_path (str): Path to the CSV file containing junction coordinates
        
    Returns:
        list: List of [min, max] coordinate pairs
    """
    df = pd.read_csv(csv_path)
    
    junctions = []
    for _, row in df.iterrows():
        chrom = row['Chrom']
        five_prime = int(row['Five_Prime_Coordinate'])
        three_prime = int(row['Three_Prime_Coordinate'])
        junctions.append([chrom, min(five_prime, three_prime), max(five_prime, three_prime)])
    
    return junctions

def sort_junctions(junctions):
    """
    Sort junctions first by chromosome and then by their start position.
    
    Args:
        junctions (list): List of [chrom, min, max] coordinate pairs
        
    Returns:
        list: Sorted list of junction pairs
    """
    return sorted(junctions, key=lambda x: (x[0], x[1]))

def merge_adjacent_junctions(sorted_junctions, max_gap):
    """
    Merge junctions that are within the specified maximum gap.
    
    Args:
        sorted_junctions (list): Sorted list of [chrom, min, max] coordinate pairs
        max_gap (int): Maximum allowed gap between junctions to merge
        
    Returns:
        list: List of merged junction coordinates
    """
    if not sorted_junctions:
        return []
        
    merged_junctions = []
    current = sorted_junctions[0]
    
    for next_junction in sorted_junctions[1:]:
        # Only merge if on same chromosome and within max_gap
        if (current[0] == next_junction[0] and  # same chromosome
            # check if 5' position of current is within max_gap of 3' position of next junction
            current[2] + max_gap >= next_junction[1]):  # within gap distance
            # Merge the junctions
            current = [
                current[0],  # keep chromosome
                min(current[1], next_junction[1]),  # min position
                max(current[2], next_junction[2])   # max position
            ]
        else:
            merged_junctions.append(current)
            current = next_junction
    
    merged_junctions.append(current)
    return merged_junctions

def process_and_merge_junctions(csv_path, design_parameters):
    """
    Import junctions from CSV file, process them, and merge where appropriate.
    
    Args:
        csv_path (str): Path to the CSV file containing junction coordinates
        design_parameters (dict): Dictionary containing design parameters including max_amplicon_gap
        
    Returns:
        list: List of merged junction coordinates as [min, max] pairs
    """
    # Read and convert coordinates
    junctions = read_junction_coordinates(csv_path)
    
    # Sort junctions
    sorted_junctions = sort_junctions(junctions)
    
    # Merge adjacent junctions
    merged_junctions = merge_adjacent_junctions(
        sorted_junctions, 
        design_parameters['max_amplicon_gap']
    )
    
    return merged_junctions

def load_parameters(config_path):
    """
    Load primer design parameters from a JSON configuration file.
    
    Args:
        config_path (str): Path to the JSON configuration file
        
    Returns:
        dict: Dictionary containing design parameters
    """
    try:
        with open(config_path, 'r') as f:
            return json.load(f)
    except FileNotFoundError:
        raise FileNotFoundError(f"Design parameters file not found at: {config_path}")
    except json.JSONDecodeError:
        raise ValueError(f"Invalid JSON format in design parameters file: {config_path}")

# loaded configuration files
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

# Define and merge junctions
# Each junction has a le
# ft (j_min) and right (j_max) position
j_min = 123456
j_max = 123456

left_primer_region, right_primer_region = calculate_primer_design_regions(j_min, j_max, design_parameters)


# Use pysam to get region from fasta
# https://pysam.readthedocs.io/en/latest/api.html


# Define "tail-sequences" or adapters to be appended to each FP or RP, respectively

fp_adapter = "forward adapter"
rp_adapter = "reverse adapter"

"""
>hg38_dna range=chr3:179195014-179195414 5'pad=200 3'pad=200 strand=+ repeatMasking=none
AAAAAAAAAAAAAAAAACACAGGGTAGAGACAGATAACAAGGGATTCCTG
ACACCTACGCTGTAGATAGCTATAACATTTCAATAGGAATCTTGGGAATC
TGGGGAACTGGAAGCAGCTGGGAAAGCACCGGCATCTGAGCCACTCAGTA
CTCTTCCAGGCTTGTGTGTGAGCCTTGCCACCTTCAGGTATTAGCACTTG
AAATCTAACTTCTTTATGAAGCTCCTTATTTACTTGCCTTCTCGGTGAAA
AAAAAAAACAACAACAAAAACCGTTTCCTTCCCCATCTGGTGCCAGCACT
GTATAATTCCTAATAAGCTTGAAAAGATAATGTTGGCAATACTTTACAAG
TTTATTGCTAATGTTTAACATTTATTAGGTGCTTGCTGTGTGCTGATCTC
T
"""


#-----// (2) Designer //---------

# DESIGNER is the core engine for creating all the primer and probe candidates and computing all the
# thermodynamic aspects of design such as target unimolecular folding, primer folding, bimolecular
# hybridization, and solving the multi-state coupled equilibria for the amount bound for the desired bimolecular
# duplex. The amount bound is directly related to the amount of signal generated in a diagnostic assay.
# DESIGNER also analyzes each primer and probe design for a series of heuristic properties such as
# sequence complexity, polyG test, oligo length penalty, amplicon length penalty, etc. Each of the scoring effects
# are multiplied by weighting factors and combined into an overall score for each primer/ probe set.

# In December 2018 functionality was added to PanelPlex is support for dbSNP (for human genome
# version GRCh38 and mouse genome version mm10). PanelPlex automatically detects all the positions with
# alternative alleles with CAF > 0.01 within each of the design regions and then automatically designs the
# primers to AVOID those SNP sites.

# Hybridization conditions


# 1. Generate potential solutions TODO Generate k-mers as initial solutions
# 2. Score each solutions TODO define scoring algorithm, including SNP penalty
# 3. Select top candidates for each target (init_solutions)
# 4. BLAST against amplicons in panel and discard matches TODO implement BLAST search
# 5. Keep top N solutions (top_solutions_to_keep)
# 6. Perform multiplex picking in N^X space for N solutions and X targets TODO implement multiplex picker algorithm


#-----// Generate k-mer solutions for targets

region = "CACAGGGTAGAGACAGATAACAAGGGATTCCTGACACCAAAAAAAAAAATACGCTGTAGATAGCTATAACATTTCAATAGGAATCTTGGGAATC"

k_min = min_primer_length
k_max = max_primer_length

def generate_kmers(k_min, k_max, sequence):
    kmers = []
    for k in range (k_min, k_max):
        # max position to search
        max_pos = len(sequence) + 1 - k
        
        for x in range(max_pos):
            kmer = sequence[x:x+k]
            kmers.append(kmer)

    return(kmers)


obj = generate_kmers(k_min, k_max, region)
print(obj)
print(len(obj))
