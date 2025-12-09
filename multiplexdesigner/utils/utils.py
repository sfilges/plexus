import pandas as pd
import json
from datetime import datetime
import uuid
from typing import Dict, List, Tuple, Optional
from Bio import SeqIO
from Bio.Seq import Seq
import warnings

class Primer3:
    """Class to represent primer3 output"""
    
    def __init__(self, name: str, junction_name: str, chrom: str, five_prime: int, three_prime: int):
        self.panel_name = name
        self.junction_name = junction_name
        self.chrom = chrom
        self.five_prime = five_prime
        self.three_prime = three_prime
        self.design_region = None
        self.design_start = None
        self.design_end = None
        self.primer_pairs_table = None
        self.left_primer_table = None
        self.right_primer_table = None

class Junction:
    """Class to represent a genomic junction/mutation position"""
    
    def __init__(self, name: str, chrom: str, five_prime: int, three_prime: int):
        self.name = name
        self.chrom = chrom
        self.five_prime = five_prime
        self.three_prime = three_prime
        self.design_region = None
        self.design_start = None
        self.design_end = None
        self.junction_length = None
        self.jmin_coordinate = None
        self.jmax_coordinate = None
        self.primer_designs = None
        self.left_primer_table = None
        self.right_primer_table = None
    
    def __repr__(self):
        return f"Junction({self.name}, {self.chrom}:{self.five_prime}-{self.three_prime})"

class MultiplexPanel:
    """Main class for managing multiplex PCR panel design"""
    
    def __init__(self, panel_name: str, genome: str = "hg38"):
        self.panel_name = panel_name
        self.genome = genome
        self.date = datetime.now().isoformat()
        self.panel_uuid = str(uuid.uuid4())
        self.junctions = []
        self.config = None
        self.junction_df = None
        
    def load_config(self, config_path: str = None, config_dict: dict = None):
        """Load design and PCR configuration parameters"""
        if config_dict:
            config = config_dict
        elif config_path:
            with open(config_path, 'r') as f:
                config = json.load(f)
        else:
            # Default configuration
            config = {
                "singleplex_design_parameters": {
                    "primer_length_range": [18, 30],
                    "amplicon_gap_range": [20, 60],
                    "max_amplicon_length": 100,
                    "junction_padding_bases": 3,
                    "initial_solutions": 100,
                    "top_solutions_to_keep": 4,
                    "minimum_plexity": 5,
                    "maximum_plexity": 20,
                    "target_plexity": 20,
                    "force_plexity": False,
                    "variant_threshold": 0.01,
                    "5_primer_tail": "",
                    "3_prime_tail": "",
                    "snp_penalty": 1.0,
                    "primer_length_penalty": 1.0,
                    "primer_complexity_penalty": 1.0,
                    "polyA_penalty": 5,
                    "polyT_penalty": 5,
                    "polyC_penalty": 10,
                    "polyG_penalty": 10,
                    "amplicon_length_penalty": 1.0
                },
                "pcr": {
                    "annealing_temperature": 60,
                    "primer_concentration": 50,
                    "dntp_concentration": 0.6,
                    "dna_concentration": 50,
                    "mv_concentration": 50,
                    "dv_concentration": 1.5,
                    "dmso_concentration": 0.0,
                    "dmso_fact": 0.6,
                    "formamide_concentration": 0.8
                }
            }
        
        self.config = config
        
    def import_junctions_csv(self, file_path: str):
        """Import junctions from CSV file using pandas"""
        try:
            df = pd.read_csv(file_path)
            required_cols = ['Name', 'Chrom', 'Five_Prime_Coordinate', 'Three_Prime_Coordinate']
            
            if not all(col in df.columns for col in required_cols):
                raise ValueError(f"CSV must contain columns: {required_cols}")
            
            self.junction_df = df.copy()
            self._create_junction_objects_from_df()
            print(f"Successfully imported {len(self.junctions)} junctions from {file_path}")
            
        except Exception as e:
            print(f"Error importing CSV file: {e}")
            raise
    
    def _create_junction_objects_from_df(self):
        """Create Junction objects from DataFrame"""
        self.junctions = []
        for _, row in self.junction_df.iterrows():
            junction = Junction(
                name=row['Name'],
                chrom=row['Chrom'],
                five_prime=int(row['Five_Prime_Coordinate']),
                three_prime=int(row['Three_Prime_Coordinate'])
            )
            self.junctions.append(junction)
    
    def merge_close_junctions(self):
        """Merge junctions that are close together based on max_amplicon_gap from config"""
        if not self.junctions or not self.config:
            print("No junctions to merge or no design config loaded")
            return
        
        # Get max_amplicon_length from config (use as merge distance)
        max_amplicon_gap = self.config.get('max_amplicon_length', 100)
        
        print(f"Merging junctions within {max_amplicon_gap} bp on same chromosome...")
        
        if self.junction_df is not None:
            # Work with DataFrame for easier manipulation
            merged_df = self._merge_junctions_df(self.junction_df.copy(), max_amplicon_gap)
            self.junction_df = merged_df
            self._create_junction_objects_from_df()
        else:
            # Fallback to junction objects
            self._merge_junctions_objects(max_amplicon_gap)
        
        print(f"After merging: {len(self.junctions)} junctions remain")
    
    def _merge_junctions_df(self, df: pd.DataFrame, max_gap: int) -> pd.DataFrame:
        """Merge junctions in DataFrame based on distance threshold"""
        df = df.sort_values(['Chrom', 'Five_Prime_Coordinate']).reset_index(drop=True)
        merged_rows = []
        i = 0
        
        while i < len(df):
            current_row = df.iloc[i].copy()
            merge_group = [current_row]
            j = i + 1
            
            # Look for adjacent junctions to merge
            while j < len(df):
                next_row = df.iloc[j]
                
                # Check if same chromosome
                if next_row['Chrom'] != current_row['Chrom']:
                    break
                
                # Calculate distance between junction regions
                current_end = max(current_row['Five_Prime_Coordinate'], current_row['Three_Prime_Coordinate'])
                next_start = min(next_row['Five_Prime_Coordinate'], next_row['Three_Prime_Coordinate'])
                distance = next_start - current_end
                
                # If within merge distance, add to group
                if distance <= max_gap:
                    merge_group.append(next_row)
                    # Update current_row to encompass the merged region
                    current_row['Three_Prime_Coordinate'] = max(
                        current_row['Three_Prime_Coordinate'],
                        next_row['Five_Prime_Coordinate'],
                        next_row['Three_Prime_Coordinate']
                    )
                    current_row['Five_Prime_Coordinate'] = min(
                        current_row['Five_Prime_Coordinate'],
                        next_row['Five_Prime_Coordinate'],
                        next_row['Three_Prime_Coordinate']
                    )
                    j += 1
                else:
                    break
            
            # Create merged junction
            if len(merge_group) > 1:
                # Merge names
                names = [row['Name'] for row in merge_group]
                merged_name = "_".join(names)
                
                # Get coordinate bounds
                all_coords = []
                for row in merge_group:
                    all_coords.extend([row['Five_Prime_Coordinate'], row['Three_Prime_Coordinate']])
                
                merged_row = current_row.copy()
                merged_row['Name'] = merged_name
                merged_row['Five_Prime_Coordinate'] = min(all_coords)
                merged_row['Three_Prime_Coordinate'] = max(all_coords)
                
                print(f"Merged {len(merge_group)} junctions: {merged_name}")
            else:
                merged_row = current_row
            
            merged_rows.append(merged_row)
            i = j if j > i + 1 else i + 1
        
        return pd.DataFrame(merged_rows).reset_index(drop=True)
    
    def _merge_junctions_objects(self, max_gap: int):
        """Fallback method to merge junction objects directly"""
        # Sort junctions by chromosome and position
        sorted_junctions = sorted(self.junctions, 
                                key=lambda x: (x.chrom, min(x.five_prime, x.three_prime)))
        
        merged_junctions = []
        i = 0
        
        while i < len(sorted_junctions):
            current_junction = sorted_junctions[i]
            merge_group = [current_junction]
            j = i + 1
            
            # Look for adjacent junctions to merge
            while j < len(sorted_junctions):
                next_junction = sorted_junctions[j]
                
                # Check if same chromosome
                if next_junction.chrom != current_junction.chrom:
                    break
                
                # Calculate distance
                current_end = max(current_junction.five_prime, current_junction.three_prime)
                next_start = min(next_junction.five_prime, next_junction.three_prime)
                distance = next_start - current_end
                
                if distance <= max_gap:
                    merge_group.append(next_junction)
                    j += 1
                else:
                    break
            
            # Create merged junction
            merged_junction = self._merge_junction_group(merge_group)
            merged_junctions.append(merged_junction)
            
            i = j if j > i + 1 else i + 1
        
        self.junctions = merged_junctions
    
    def _merge_junction_group(self, junction_group: List[Junction]) -> Junction:
        """Merge a group of junctions into a single junction"""
        if len(junction_group) == 1:
            return junction_group[0]
        
        # Create merged name
        names = [j.name for j in junction_group]
        merged_name = "_".join(names)
        
        # Get chromosome (should be same for all)
        chrom = junction_group[0].chrom
        
        # Get min and max coordinates
        all_coords = []
        for j in junction_group:
            all_coords.extend([j.five_prime, j.three_prime])
        
        merged_five_prime = min(all_coords)
        merged_three_prime = max(all_coords)
        
        return Junction(merged_name, chrom, merged_five_prime, merged_three_prime)
    
    def extract_design_regions_from_fasta(self, fasta_file: str, padding: int = 200):
        """Extract genomic sequences for design regions with padding from FASTA file"""
        print(f"Extracting design regions from {fasta_file} with {padding}bp padding...")
        
        regions_extracted = 0
        for junction in self.junctions:
            try:
                # Calculate design region coordinates
                junction_start = min(junction.five_prime, junction.three_prime)
                junction_end = max(junction.five_prime, junction.three_prime)
                
                design_start = junction_start - padding
                design_end = junction_end + padding
                
                # Find matching chromosome in genome
                design_sequence = self.seqio_extract_genomic_sequence(fasta_file, junction.chrom, design_start, design_end)
                
                # Store in junction object
                junction.design_region = design_sequence.upper()
                junction.design_start = design_start
                junction.design_end = design_end
                regions_extracted += 1
                
            except Exception as e:
                print(f"Error extracting region for {junction.name}: {e}")
        
        print(f"Successfully extracted {regions_extracted} design regions")

    def seqio_extract_genomic_sequence(self, fasta_file, sequence_id, start, end, strand='+'):
        """Memory-efficient extraction for large files."""
        with open(fasta_file, 'r') as handle:
            for record in SeqIO.parse(handle, "fasta"):
                if record.id == sequence_id:
                    seq = record.seq[start-1:end]
                    if strand == '-':
                        seq = seq.reverse_complement()
                    return str(seq)
        raise ValueError(f"Sequence ID '{sequence_id}' not found")
    
    def calculate_junction_coordinates_in_design_region(self):
        """Calculate junction coordinates within design regions with proper logic"""
        if not self.config:
            print("Warning: No design config loaded, using default padding of 3bp")
            padding = 3
        else:
            padding = self.config.get('junction_padding_bases', 3)
        
        coordinates_calculated = 0
        
        for junction in self.junctions:
            if junction.design_region is None or junction.design_start is None:
                print(f"Warning: No design region found for {junction.name}")
                continue
            
            try:
                # Length of the design region sequence
                junction.junction_length = len(junction.design_region)
                
                # Calculate junction coordinates relative to design region start (0-based)
                # junction positions are 1-based genomic coordinates
                # design_start is 0-based genomic coordinate
                
                # Convert junction coordinates to 0-based relative to design region
                junction_five_rel = junction.five_prime - junction.design_start - 1  # Convert to 0-based
                junction_three_rel = junction.three_prime - junction.design_start - 1  # Convert to 0-based
                
                # Add padding around the junction region
                jmin_coordinate = min(junction_five_rel, junction_three_rel) - padding
                jmax_coordinate = max(junction_five_rel, junction_three_rel) + padding
                
                # Ensure coordinates are within design region bounds
                jmin_coordinate = max(0, jmin_coordinate)
                jmax_coordinate = min(junction.junction_length - 1, jmax_coordinate)
                
                junction.jmin_coordinate = jmin_coordinate
                junction.jmax_coordinate = jmax_coordinate
                coordinates_calculated += 1
                
                # Verify the calculation makes sense
                if jmin_coordinate >= jmax_coordinate:
                    print(f"Warning: Invalid junction coordinates for {junction.name}")
                
            except Exception as e:
                print(f"Error calculating coordinates for {junction.name}: {e}")
        
        print(f"Calculated junction coordinates for {coordinates_calculated} junctions")
    
    def verify_junction_coordinates(self):
        """Verify that junction coordinate calculations make sense"""
        print("\nVerifying junction coordinates...")
        
        for junction in self.junctions:
            if all(x is not None for x in [junction.jmin_coordinate, junction.jmax_coordinate, 
                                         junction.design_region, junction.design_start]):
                
                print(f"\nJunction: {junction.name}")
                print(f"  Genomic coordinates: {junction.chrom}:{junction.five_prime}-{junction.three_prime}")
                print(f"  Design region: {junction.chrom}:{junction.design_start}-{junction.design_end}")
                print(f"  Design region length: {junction.junction_length}")
                print(f"  Junction in design region: {junction.jmin_coordinate}-{junction.jmax_coordinate}")
                
                # Extract the junction sequence for verification
                if 0 <= junction.jmin_coordinate < junction.jmax_coordinate <= len(junction.design_region):
                    junction_seq = junction.design_region[junction.jmin_coordinate:junction.jmax_coordinate+1]
                    print(f"  Junction sequence: {junction_seq[:50]}...")
                else:
                    print(f"  Warning: Junction coordinates out of bounds!")
        
        print("\nCoordinate verification complete.")
    
    def create_junction_table(self) -> pd.DataFrame:
        """Create pandas DataFrame with all junction information"""
        data = []
        for junction in self.junctions:
            row = {
                'Name': junction.name,
                'Chrom': junction.chrom,
                'Five_Prime_Coordinate': junction.five_prime,
                'Three_Prime_Coordinate': junction.three_prime,
                'Design_Region': junction.design_region if junction.design_region else '',
                'Design_Start': junction.design_start if junction.design_start else '',
                'Design_End': junction.design_end if junction.design_end else '',
                'Junction_Length': junction.junction_length if junction.junction_length else '',
                'Jmin_Coordinate': junction.jmin_coordinate if junction.jmin_coordinate else '',
                'Jmax_Coordinate': junction.jmax_coordinate if junction.jmax_coordinate else ''
            }
            data.append(row)
        
        self.junction_df = pd.DataFrame(data)
        return self.junction_df


# Helper function for complete workflow
def process_multiplex_panel(csv_file: str, fasta_file: str, config_file: str = None, 
                           panel_name: str = "multiplex_panel", output_dir: str = ".",
                           padding: int = 200):
    """Complete workflow for processing a multiplex panel"""
    
    print(f"Processing multiplex panel: {panel_name}")
    
    # Initialize panel
    panel = MultiplexPanel(panel_name, "hg38")
    
    # Load configuration
    if config_file:
        panel.load_config(config_file)
    else:
        panel.load_config()  # Use defaults
    
    # Step 1: Import junctions
    print("\nStep 1: Importing junctions...")
    panel.import_junctions_csv(csv_file)
    
    # Step 2: Merge close junctions
    print("\nStep 2: Merging close junctions...")
    panel.merge_close_junctions()
    
    # Step 3: Extract design regions
    print("\nStep 3: Extracting design regions...")
    panel.extract_design_regions_from_fasta(fasta_file, padding=padding)
    
    # Step 4: Calculate junction coordinates
    print("\nStep 4: Calculating junction coordinates...")
    panel.calculate_junction_coordinates_in_design_region()
    
    # Step 5: Verify coordinates
    print("\nStep 5: Verifying coordinates...")
    panel.verify_junction_coordinates()
    
    # Step 6: Export results
    print("\nStep 6: Exporting results...")
    import os
    junction_table_file = os.path.join(output_dir, f"{panel_name}_junctions.csv")
    panel_data_file = os.path.join(output_dir, f"{panel_name}_complete.json")
    
    panel.export_junction_table(junction_table_file)
    panel.save_panel(panel_data_file)
    
    print(f"\nPanel processing complete!")
    print(f"Junction table saved to: {junction_table_file}")
    print(f"Complete panel data saved to: {panel_data_file}")
    
    return panel


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