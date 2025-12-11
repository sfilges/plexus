import pandas as pd
import warnings
import logging
from datetime import datetime

def setup_logger():
    # Create a timestamp for the log file
    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    log_filename = f"multiplexdesigner_{timestamp}.log"

    # Configure the logger
    logging.basicConfig(
        filename=log_filename,
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s"
    )

    return logging.getLogger()


def gc_content(sequence: str) -> float:
    """
    Calculate the GC content of a DNA sequence.
    
    Parameters:
        sequence (str): DNA sequence consisting of A, T, G, C.
    
    Returns:
        float: GC content as a percentage.
    """
    if not sequence:
        return 0.0  # Handle empty string safely
    
    sequence = sequence.upper()  # Ensure case-insensitivity
    gc_count = sequence.count("G") + sequence.count("C")
    return (gc_count / len(sequence)) * 100

def reverse_complement(dna):
    """
    Returns the reverse complement of a DNA sequence.
    
    Args:
        dna (str): DNA sequence string (A, T, G, C)
    
    Returns:
        str: Reverse complement of the input DNA sequence
    """
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    
    # Handle lowercase and raise error for invalid bases
    dna = dna.upper()
    
    try:
        reverse_comp = ''.join(complement[base] for base in reversed(dna))
        return reverse_comp
    except KeyError as e:
        raise ValueError(f"Invalid DNA base: {e.args[0]}")


# Convert different sections to DataFrames
def convert_primer_data_to_tables(data):
    tables = {}
    
    # 1. Primer Pairs Summary
    if 'PRIMER_PAIR' in data:
        primer_pairs_df = pd.DataFrame(data['PRIMER_PAIR'])
        primer_pairs_df.index.name = 'Pair_ID'
        primer_pairs_df = primer_pairs_df.round(3)  # Round decimals
        tables['primer_pairs'] = primer_pairs_df
    
    # 2. Left Primers
    if 'PRIMER_LEFT' in data:
        left_primers_df = pd.DataFrame(data['PRIMER_LEFT'])
        left_primers_df.index.name = 'Left_Primer_ID'
        left_primers_df = left_primers_df.round(3)
        tables['left_primers'] = left_primers_df
    
    # 3. Right Primers  
    if 'PRIMER_RIGHT' in data:
        right_primers_df = pd.DataFrame(data['PRIMER_RIGHT'])
        right_primers_df.index.name = 'Right_Primer_ID'
        right_primers_df = right_primers_df.round(3)
        tables['right_primers'] = right_primers_df
    
    # 4. Internal Primers
    if 'PRIMER_INTERNAL' in data:
        internal_primers_df = pd.DataFrame(data['PRIMER_INTERNAL'])
        internal_primers_df.index.name = 'Internal_Primer_ID'
        internal_primers_df = internal_primers_df.round(3)
        tables['internal_primers'] = internal_primers_df
    
    # 5. Summary Statistics
    summary_data = {
        'Metric': ['Left Primers Returned', 'Right Primers Returned', 'Internal Primers Returned', 'Primer Pairs Returned'],
        'Count': [
            data.get('PRIMER_LEFT_NUM_RETURNED', 0),
            data.get('PRIMER_RIGHT_NUM_RETURNED', 0), 
            data.get('PRIMER_INTERNAL_NUM_RETURNED', 0),
            data.get('PRIMER_PAIR_NUM_RETURNED', 0)
        ]
    }
    tables['summary'] = pd.DataFrame(summary_data)
    
    return tables

def create_primer_dataframe(primer_data):
    """
    Convert primer design results to organized pandas DataFrames
    """
    
    # Method 1: Create a comprehensive primer pairs DataFrame
    primer_pairs = []
    
    num_pairs = primer_data.get('PRIMER_PAIR_NUM_RETURNED', 0)
    
    for i in range(num_pairs):
        pair_info = {
            'pair_id': i,
            'pair_penalty': primer_data.get(f'PRIMER_PAIR_{i}_PENALTY'),
            'product_size': primer_data.get(f'PRIMER_PAIR_{i}_PRODUCT_SIZE'),
            'product_tm': primer_data.get(f'PRIMER_PAIR_{i}_PRODUCT_TM'),
            'compl_any_th': primer_data.get(f'PRIMER_PAIR_{i}_COMPL_ANY_TH'),
            'compl_end_th': primer_data.get(f'PRIMER_PAIR_{i}_COMPL_END_TH'),
            'template_misprime_th': primer_data.get(f'PRIMER_PAIR_{i}_TEMPLATE_MISPRIMING_TH'),
            
            # Left primer info
            'left_sequence': primer_data.get(f'PRIMER_LEFT_{i}_SEQUENCE'),
            'left_coords': primer_data.get(f'PRIMER_LEFT_{i}'),
            'left_tm': primer_data.get(f'PRIMER_LEFT_{i}_TM'),
            'left_gc_percent': primer_data.get(f'PRIMER_LEFT_{i}_GC_PERCENT'),
            'left_penalty': primer_data.get(f'PRIMER_LEFT_{i}_PENALTY'),
            'left_self_any_th': primer_data.get(f'PRIMER_LEFT_{i}_SELF_ANY_TH'),
            'left_self_end_th': primer_data.get(f'PRIMER_LEFT_{i}_SELF_END_TH'),
            'left_hairpin_th': primer_data.get(f'PRIMER_LEFT_{i}_HAIRPIN_TH'),
            'left_end_stability': primer_data.get(f'PRIMER_LEFT_{i}_END_STABILITY'),
            ''
            
            # Right primer info
            'right_sequence': primer_data.get(f'PRIMER_RIGHT_{i}_SEQUENCE'),
            'right_coords': primer_data.get(f'PRIMER_RIGHT_{i}'),
            'right_tm': primer_data.get(f'PRIMER_RIGHT_{i}_TM'),
            'right_gc_percent': primer_data.get(f'PRIMER_RIGHT_{i}_GC_PERCENT'),
            'right_penalty': primer_data.get(f'PRIMER_RIGHT_{i}_PENALTY'),
            'right_self_any_th': primer_data.get(f'PRIMER_RIGHT_{i}_SELF_ANY_TH'),
            'right_self_end_th': primer_data.get(f'PRIMER_RIGHT_{i}_SELF_END_TH'),
            'right_hairpin_th': primer_data.get(f'PRIMER_RIGHT_{i}_HAIRPIN_TH'),
            'right_end_stability': primer_data.get(f'PRIMER_RIGHT_{i}_END_STABILITY'),
        }
        primer_pairs.append(pair_info)
    
    df_pairs = pd.DataFrame(primer_pairs)

    return(df_pairs)


def generate_kmers(sequence, k_min: int = 18, k_max: int = 25):
    """
    Generate k-mers as candidate primers.

    Args:
        sequence - Region from which to extract kmers
        k_min - Minimum length of kmers
        k_max - Maximum length of kmers
    """
    if k_min >= k_max:
        raise ValueError(f'Min kmer length ({k_min}) must be smaller than max kmer length ({k_max})')

    if k_min < 10 or k_min > 20:
        warnings.warn(f'Provided value for kmin, {k_min}, is outside the expected range: 10 - 20')

    if k_max < 20 or k_max > 30:
        warnings.warn(f'Provided value for kmax, {k_max}, is outside the expected range: 20 - 30')

    kmers = []
    for k in range(k_min, k_max):
        # max position to search
        max_pos = len(sequence) + 1 - k
        
        for x in range(max_pos):
            kmer = sequence[x:x+k]
            kmers.append(kmer)
    
    return(kmers)

def filter_kmers(kmers, max_poly_X = 4, max_N = 0):
    """
    Filter kmers (putative primers).
    """
    filtered_kmers = []
    for kmer in kmers:
        # Check if kmers contains more than n N-bases.
        if check_N_in_kmers(kmer, max_N):
            continue

        # Exclude kmers with too high or low GC-content.
        kmer_gc = gc_content(kmer)
        if kmer_gc < 30 or kmer_gc > 70:
            continue
        
        # Exclude kmers with too many repeated bases.
        if find_max_poly_X(kmer, max_poly_X):
            continue

        filtered_kmers.append(kmer)

    return(filtered_kmers)

def check_N_in_kmers(kmer, n=0):
    """
    Check if a k-mer contains more than n 'N' bases.
    
    Args:
        kmer (str): The k-mer sequence to check
        n (int): The maximum allowed number of N bases (default is 0)
    
    Returns:
        bool: True if the k-mer contains more than n 'N' bases, False otherwise
    """
    return kmer.upper().count('N') > n

def find_max_poly_X(kmer, n):
    """
    Check if a DNA sequence has more than n consecutive identical bases.
    
    Args:
        dna_sequence (str): The DNA sequence to check
        n (int): The threshold for consecutive bases (default is 4)
    
    Returns:
        bool: True if there are more than n consecutive identical bases, False otherwise
    """
    if len(kmer) <= n:
        return False
    
    current_base = kmer[0]
    count = 1
    
    for i in range(1, len(kmer)):
        if kmer[i] == current_base:
            count += 1
            if count > n:
                return True
        else:
            current_base = kmer[i]
            count = 1
    
    return False