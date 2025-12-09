import json
import pysam
import pandas as pd

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

def merge_adjacent_junctions(sorted_junctions, max_amplicon_gap):
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
            current[2] + max_amplicon_gap >= next_junction[1]):  # within gap distance
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

def extract_genomic_region(fasta_file, chromosome, start, end, output_file=None):
    """
    Extract a genomic region from a FASTA file using pysam
    
    Args:
        fasta_file (str): Path to the FASTA file (e.g., 'genome.fa')
        chromosome (str): Chromosome/contig name (e.g., 'chr1', '1')
        start (int): Start position (0-based, inclusive)
        end (int): End position (0-based, exclusive)
        output_file (str, optional): Path to save the extracted sequence
    
    Returns:
        str: The extracted DNA sequence
    """
    
    # Open the FASTA file
    fasta = pysam.FastaFile(fasta_file)
    
    try:
        # Extract the sequence
        sequence = fasta.fetch(chromosome, start, end)
        
        # Optionally save to file
        if output_file:
            with open(output_file, 'w') as f:
                f.write(f">{chromosome}:{start}-{end}\n")
                f.write(sequence + "\n")
        
        return sequence
        
    except KeyError:
        print(f"Chromosome '{chromosome}' not found in FASTA file")
        return None
    except Exception as e:
        print(f"Error extracting region: {e}")
        return None
    finally:
        fasta.close()
