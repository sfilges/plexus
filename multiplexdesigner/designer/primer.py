import json
import uuid
import pandas as pd
from typing import List
from Bio import SeqIO
from dataclasses import dataclass, field
from datetime import datetime

# TODO: import primer pairs from primer3 output. Wouldn't it be better to design forward and reverse primers
# separately instead of in pairs, and only consider pairing while also checking for primer dimers?

def load_primer_pairs_from_primer3_output(primer3_output, add_target=None):
    """
    Given an output file from primer3, return a list of
    PrimerPair objects

    params
        primer3_output_path: str
            Path to an output file produced by primer3. This will
            contain information a series of primer pairs, in a format
            <key>=<value>.
        add_target: Target [optional]
            A Target object, containing information about on which target
            primer3 run.

    """

    primer_pairs = []



    return primer_pairs


@dataclass
class Primer:
    """
    Define a single primer.
    """
    seq: str
    direction: str
    start: int
    length: int
    tm: float
    gc: float
    name: str=""

    def name_primer(self):
        """
        Name the primer
        """
        pass
    
    def add_tails(self):
        """
        Add tail sequences to primer.
        """
        pass

@dataclass
class PrimerPair:
    """
    Define a pair of primers.
    """
    forward: Primer
    reverse: Primer
    product_length: int
    pair_id: str = field(default="", repr=False)


@dataclass
class Primer3:
    """Class to represent primer3 output"""
    panel_name: str
    junction_name: str
    chrom: str
    five_prime: int
    three_prime: int
    design_region: str = None
    design_start: int = None
    design_end: int = None
    primer_pairs_table: object = None
    left_primer_table: object = None
    right_primer_table: object = None

@dataclass
class Junction:
    """Class to represent a genomic junction/mutation position"""
    name: str
    chrom: str
    five_prime: int
    three_prime: int
    design_region: str = None
    design_start: int = None
    design_end: int = None
    junction_length: int = None
    jmin_coordinate: int = None
    jmax_coordinate: int = None
    primer3_designs: object = None
    left_primer_table: object = None
    right_primer_table: object = None

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
        
    def load_config(self, logger, config_path: str = None, config_dict: dict = None):
        """Load design and PCR configuration parameters"""
        if config_dict:
            config = config_dict
        elif config_path:
            with open(config_path, 'r') as f:
                config = json.load(f)
            logger.info("Config imported successfully.")
        else:
            logger.error("MultiplexPanel.load_config - No suitable config file or dict provided.")
            raise ValueError("MultiplexPanel.load_config - No suitable config file or dict provided.")

        self.config = config
        
    def import_junctions_csv(self, file_path: str, logger):
        """Import junctions from CSV file using pandas"""
        try:
            df = pd.read_csv(file_path)
            required_cols = ['Name', 'Chrom', 'Five_Prime_Coordinate', 'Three_Prime_Coordinate']
            
            if not all(col in df.columns for col in required_cols):
                raise ValueError(f"CSV must contain columns: {required_cols}")
            
            self.junction_df = df.copy()
            self._create_junction_objects_from_df()
            print(f"Successfully imported {len(self.junctions)} junctions from {file_path}")
            logger.info(f"Successfully imported {len(self.junctions)} junctions from {file_path}")
            
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
    
    def merge_close_junctions(self, logger):
        """Merge junctions that are close together based on max_amplicon_gap from config"""
        if not self.junctions or not self.config:
            print("No junctions to merge or no design config loaded")
            logger.warning("No junctions to merge or no design config loaded")
            return
        
        # Get max_amplicon_length from config (use as merge distance)
        max_amplicon_gap = self.config.get('max_amplicon_length', 100)
        
        print(f"Merging junctions within {max_amplicon_gap} bp on same chromosome...")
        logger.info(f"Merging junctions within {max_amplicon_gap} bp on same chromosome...")

        if self.junction_df is not None:
            # Work with DataFrame for easier manipulation
            merged_df = self._merge_junctions_df(self.junction_df.copy(), max_amplicon_gap)
            self.junction_df = merged_df
            self._create_junction_objects_from_df()
        else:
            # Fallback to junction objects
            self._merge_junctions_objects(max_amplicon_gap)
        
        print(f"After merging: {len(self.junctions)} junctions remain")
        logger.info(f"After merging: {len(self.junctions)} junctions remain")
    
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
    
    def extract_design_regions_from_fasta(self, fasta_file: str, logger, padding: int = 200):
        """
        Extract genomic sequences for design regions with padding from FASTA file.

        Args:
            self:
            fasta_file: str
            logger: logging.getLogger() object
            padding: int = 200
        
        Return:
            Junction object
        """

        logger.info(f"Extracting design regions from {fasta_file} with {padding}bp padding...")
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
        logger.info(f"Successfully extracted {regions_extracted} design regions")

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
