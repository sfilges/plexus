# ================================================================================
# Multiplexpanel classes and associated functions
#
# Author: Stefan Filges (stefan@simsendiagnostics.com)
# Copyright (c) 2025 Simsen Diagnostics AB
# ================================================================================

import json
import os
import uuid
from dataclasses import dataclass
from datetime import datetime

import pandas as pd
from Bio import SeqIO
from loguru import logger

from multiplexdesigner.designer.primer import PrimerPair
from multiplexdesigner.utils.root_dir import ROOT_DIR

# ================================================================================
# A class to hold primer designs and design metrics.
# ================================================================================


@dataclass
class PrimerDesigns:
    """
    A class to hold primer designs and metrics.
    """

    name: str
    target: str
    design_region: str = None
    eval_string: str = None
    primer_table: object = None


# ================================================================================
# Define a single multiplex PCR target/junction
# ================================================================================


@dataclass
class Junction:
    """
    Class to represent a genomic junction/mutation position.

    Args:
        - name: Name of of the junction provided by the user
        - chrom: chromosome of the junction. Must be in the same format as the reference genome
        - start: 5' coordinate of the junction
        - end: 3' coordinate of the junction
        - design region: Computed template sequence for primer design
        - design_start
        - design_end
        - junction_length
        - jmin_coordinate
        - jmax_coordinate
        - primer_designs: PrimerDesigns object
    """

    name: str
    chrom: str
    start: int
    end: int
    design_region: str = None
    design_start: int = None
    design_end: int = None
    junction_length: int = None
    jmin_coordinate: int = None
    jmax_coordinate: int = None
    primer_designs: object = None
    primer_pairs: list = None

    def __repr__(self):
        return f"Junction({self.name}, {self.chrom}:{self.start}-{self.end})"

    def get_primer_designs(self):
        return self.primer_designs.primer_table

    # TODO: Write method to print/save self.primer_designs.primer_table
    def save_primer_designs(self):
        pass

    def find_primer_pairs(
        self,
        min_amplicon_length: int = 40,
        max_amplicon_length: int = 100,
        max_primer_tm_difference: int = 3,
        product_opt_size: int = 0,
        wt_pr_penalty: float = 1.0,
        wt_product_size_gt: float = 1.0,
        wt_product_size_lt: float = 1.0,
        wt_diff_tm: float = 1.0,
    ) -> list[tuple[object, object]]:
        """
        Optimized version using sorted primers for better performance.

        Args:
            self: A junction object
            min_amplicon_length: Minimal length of the entire amplicon, incl. primers.
            max_amplicon_length: Maximal length of the entire amplicon, incl. primers.
            max_primer_tm_difference: Maximum difference in melting temperature between primers.
            product_opt_size: Optimal PCR product size

        Returns:
            A list of tuples with (forward, reverse) primer objects.
        """
        valid_pairs = []

        logger.info(f"finding suitable primer pairs for junction: {self.name}")

        left_primers = self.primer_designs.primer_table[0]
        right_primers = self.primer_designs.primer_table[1]

        # Sort primers by start position
        left_sorted = sorted(left_primers, key=lambda p: p.start)
        right_sorted = sorted(right_primers, key=lambda p: p.start)

        # Init counters
        n_too_short = 0
        n_too_long = 0
        primer_tm_diff_too_large = 0
        primer_dimer_too_stable = 0
        primer_three_prime_too_stable = 0

        for left_primer in left_sorted:
            for right_primer in right_sorted:
                # Skip right primers that start before left primer ends
                if right_primer.start < left_primer.start + left_primer.length:
                    continue

                primer_pair_tm_dif = abs(right_primer.tm - left_primer.tm)
                if primer_pair_tm_dif > max_primer_tm_difference:
                    primer_tm_diff_too_large += 1
                    continue

                # Amplicon length includes primer sequences
                amplicon_length = (
                    right_primer.start + right_primer.length
                ) - left_primer.start

                # If amplicon is too short, continue to next right primer
                if amplicon_length < min_amplicon_length:
                    n_too_short += 1
                    continue

                # If amplicon is too long, no need to check further right primers
                if amplicon_length > max_amplicon_length:
                    n_too_long += 1
                    break

                # Get amplicon sequence and calculate product Tm
                design_region = self.design_region
                amplicon_start = left_primer.start
                amplicon_end = right_primer.start + right_primer.length
                amplicon_sequence = design_region[amplicon_start:amplicon_end]

                insert_start = left_primer.start + left_primer.length
                insert_end = right_primer.start
                insert_size = len(design_region[insert_start:insert_end])

                pair = PrimerPair(
                    forward=left_primer,
                    reverse=right_primer,
                    insert_size=insert_size,
                    amplicon_sequence=amplicon_sequence,
                    amplicon_length=amplicon_length,
                    pair_id=f"{left_primer.name}_{right_primer.name}",
                )

                # TODO: Calculate the thermodyanamics of primer-to-primer binding in each pair. This can also be
                # done at the primer pair selection step? Should these be calculated here?

                # logger.info(f"Left seq: {left_primer.seq}; Right seq: {right_primer.seq}")
                # logger.info(f"Left penalty: {left_primer.penalty}; Right penalty: {right_primer.penalty}")

                # Calculate primer pair penalty (similar to primer3)
                pair.pair_penalty = pair.calculate_primer_pair_penalty_th(
                    primer_left_penalty=left_primer.penalty,
                    primer_right_penalty=right_primer.penalty,
                    primer_left_tm=left_primer.tm,
                    primer_right_tm=right_primer.tm,
                    product_size=amplicon_length,
                    product_opt_size=product_opt_size,
                    wt_pr_penalty=wt_pr_penalty,
                    wt_product_size_gt=wt_product_size_gt,
                    wt_product_size_lt=wt_product_size_lt,
                    wt_diff_tm=wt_diff_tm,
                )

                valid_pairs.append(pair)

        logger.info(
            f"Considered {len(left_sorted)} left and {len(right_sorted)} right primers. Found: {len(valid_pairs)} valid primer pairs."
        )
        logger.info(
            f"Amplicon too short: {n_too_short}; Amplicon too long: {n_too_long}; Primer Tm difference larger than {max_primer_tm_difference}: {primer_tm_diff_too_large}"
        )
        logger.info(
            f"Primer dimers too stable: {primer_dimer_too_stable}; 3'-ends too stable: {primer_three_prime_too_stable}"
        )

        return valid_pairs


# ================================================================================
# Main class for panel design
# ================================================================================


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

    def load_config(
        self,
        preset: str = "default",
        config_path: str = None,
        config_dict: dict = None,
    ):
        """
        Load design and PCR configuration parameters.

        Args:
            self: A MultiplexPanel object
            preset: Which standard config file to use? Either "default" or "lenient".
            config_path: User-provided config from a file. If none is provided, the config is chosen based on the preset value.
            config_dict: User-provided config from a dictionary. If none is provided, the config is chosen based on the preset value.

        Returns:
            A MultiplexPanel object with loaded configuration.
        """
        if config_dict:
            config = config_dict
            logger.info("Config loaded from dictionary.")
        elif config_path:
            with open(config_path) as f:
                config = json.load(f)
            logger.info(f"Config loaded from: {config_path}")
        else:
            logger.info(
                f"No config file provided, using preset configuration: {preset}"
            )
            if preset == "default":
                config_path = f"{ROOT_DIR}/config/designer_default_config.json"
            elif preset == "lenient":
                config_path = f"{ROOT_DIR}/config/designer_lenient_config.json"
            else:
                config_path = f"{ROOT_DIR}/config/designer_default_config.json"

                warn_msg = f"Preset value `{preset}` must be either 'default' or 'lenient'. Using default instead."
                logger.warning(warn_msg)

            with open(config_path) as f:
                config = json.load(f)
            logger.info(f"Config loaded from: {config_path}")

        self.config = config

    def import_junctions_csv(self, file_path: str):
        """Import junctions from CSV file using pandas"""
        try:
            df = pd.read_csv(file_path)
            required_cols = [
                "Name",
                "Chrom",
                "Five_Prime_Coordinate",
                "Three_Prime_Coordinate",
            ]

            if not all(col in df.columns for col in required_cols):
                raise ValueError(f"CSV must contain columns: {required_cols}")

            self.junction_df = df.copy()
            self._create_junction_objects_from_df()
            logger.info(
                f"Successfully imported {len(self.junctions)} junctions from {file_path}"
            )

        except Exception as e:
            logger.error(f"Error importing CSV file: {e}")
            raise

    def _create_junction_objects_from_df(self):
        """Create Junction objects from DataFrame"""
        self.junctions = []
        for _, row in self.junction_df.iterrows():
            junction = Junction(
                name=row["Name"],
                chrom=row["Chrom"],
                start=int(row["Five_Prime_Coordinate"]),
                end=int(row["Three_Prime_Coordinate"]),
            )
            self.junctions.append(junction)

    def merge_close_junctions(self):
        """Merge junctions that are close together based on max_amplicon_gap from config"""
        if not self.junctions or not self.config:
            logger.warning("No junctions to merge or no design config loaded")
            return

        # Get max_amplicon_length from config (use as merge distance)
        max_amplicon_gap = self.config.get("max_amplicon_length", 100)

        logger.info(
            f"Merging junctions within {max_amplicon_gap} bp on same chromosome..."
        )

        if self.junction_df is not None:
            # Work with DataFrame for easier manipulation
            merged_df = self._merge_junctions_df(
                self.junction_df.copy(), max_amplicon_gap
            )
            self.junction_df = merged_df
            self._create_junction_objects_from_df()
        else:
            # Fallback to junction objects
            self._merge_junctions_objects(max_amplicon_gap)

        logger.info(f"After merging: {len(self.junctions)} junctions remain")

    def _merge_junctions_df(self, df: pd.DataFrame, max_gap: int) -> pd.DataFrame:
        """Merge junctions in DataFrame based on distance threshold"""
        df = df.sort_values(["Chrom", "Five_Prime_Coordinate"]).reset_index(drop=True)
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
                if next_row["Chrom"] != current_row["Chrom"]:
                    break

                # Calculate distance between junction regions
                current_end = max(
                    current_row["Five_Prime_Coordinate"],
                    current_row["Three_Prime_Coordinate"],
                )
                next_start = min(
                    next_row["Five_Prime_Coordinate"],
                    next_row["Three_Prime_Coordinate"],
                )
                distance = next_start - current_end

                # If within merge distance, add to group
                if distance <= max_gap:
                    merge_group.append(next_row)
                    # Update current_row to encompass the merged region
                    current_row["Three_Prime_Coordinate"] = max(
                        current_row["Three_Prime_Coordinate"],
                        next_row["Five_Prime_Coordinate"],
                        next_row["Three_Prime_Coordinate"],
                    )
                    current_row["Five_Prime_Coordinate"] = min(
                        current_row["Five_Prime_Coordinate"],
                        next_row["Five_Prime_Coordinate"],
                        next_row["Three_Prime_Coordinate"],
                    )
                    j += 1
                else:
                    break

            # Create merged junction
            if len(merge_group) > 1:
                # Merge names
                names = [row["Name"] for row in merge_group]
                merged_name = "_".join(names)

                # Get coordinate bounds
                all_coords = []
                for row in merge_group:
                    all_coords.extend(
                        [row["Five_Prime_Coordinate"], row["Three_Prime_Coordinate"]]
                    )

                merged_row = current_row.copy()
                merged_row["Name"] = merged_name
                merged_row["Five_Prime_Coordinate"] = min(all_coords)
                merged_row["Three_Prime_Coordinate"] = max(all_coords)

                logger.debug(f"Merged {len(merge_group)} junctions: {merged_name}")
            else:
                merged_row = current_row

            merged_rows.append(merged_row)
            i = j if j > i + 1 else i + 1

        return pd.DataFrame(merged_rows).reset_index(drop=True)

    def _merge_junctions_objects(self, max_gap: int):
        """Fallback method to merge junction objects directly"""
        # Sort junctions by chromosome and position
        sorted_junctions = sorted(
            self.junctions, key=lambda x: (x.chrom, min(x.start, x.end))
        )

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
                current_end = max(current_junction.start, current_junction.end)
                next_start = min(next_junction.start, next_junction.end)
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

    def _merge_junction_group(self, junction_group: list[Junction]) -> Junction:
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
            all_coords.extend([j.start, j.end])

        merged_five_prime = min(all_coords)
        merged_three_prime = max(all_coords)

        return Junction(merged_name, chrom, merged_five_prime, merged_three_prime)

    def extract_design_regions_from_fasta(self, fasta_file: str, padding: int = 200):
        """
        Extract genomic sequences for design regions with padding from FASTA file.

        Args:
            fasta_file: Path to the FASTA file
            padding: Number of bases to pad on each side of the junction

        Return:
            None (modifies junction objects in place)
        """
        logger.info(
            f"Extracting design regions from {fasta_file} with {padding}bp padding..."
        )

        regions_extracted = 0
        for junction in self.junctions:
            try:
                # Calculate design region coordinates
                junction_start = min(junction.start, junction.end)
                junction_end = max(junction.start, junction.end)

                design_start = junction_start - padding
                design_end = junction_end + padding

                # Find matching chromosome in genome
                design_sequence = self.seqio_extract_genomic_sequence(
                    fasta_file, junction.chrom, design_start, design_end
                )

                # Store in junction object
                junction.design_region = design_sequence.upper()
                junction.design_start = design_start
                junction.design_end = design_end
                regions_extracted += 1

            except Exception as e:
                logger.error(f"Error extracting region for {junction.name}: {e}")

        logger.info(f"Successfully extracted {regions_extracted} design regions")

    def seqio_extract_genomic_sequence(
        self, fasta_file, sequence_id, start, end, strand="+"
    ):
        """Memory-efficient extraction for large files."""
        with open(fasta_file) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                if record.id == sequence_id:
                    seq = record.seq[start - 1 : end]
                    if strand == "-":
                        seq = seq.reverse_complement()
                    return str(seq)
        raise ValueError(f"Sequence ID '{sequence_id}' not found")

    def calculate_junction_coordinates_in_design_region(self):
        """Calculate junction coordinates within design regions with proper logic"""
        if not self.config:
            logger.warning("No design config loaded, using default padding of 3bp")
            padding = 3
        else:
            padding = self.config.get("junction_padding_bases", 3)

        coordinates_calculated = 0

        for junction in self.junctions:
            if junction.design_region is None or junction.design_start is None:
                logger.warning(f"No design region found for {junction.name}")
                continue

            try:
                # Length of the design region sequence
                junction.junction_length = len(junction.design_region)

                # Calculate junction coordinates relative to design region start (0-based)
                # junction positions are 1-based genomic coordinates
                # design_start is 0-based genomic coordinate

                # Convert junction coordinates to 0-based relative to design region
                junction_five_rel = (
                    junction.start - junction.design_start - 1
                )  # Convert to 0-based
                junction_three_rel = (
                    junction.end - junction.design_start - 1
                )  # Convert to 0-based

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
                    logger.warning(f"Invalid junction coordinates for {junction.name}")

            except Exception as e:
                logger.error(f"Error calculating coordinates for {junction.name}: {e}")

        logger.info(
            f"Calculated junction coordinates for {coordinates_calculated} junctions"
        )

    def verify_junction_coordinates(self):
        """Verify that junction coordinate calculations make sense"""
        logger.info("Verifying junction coordinates...")

        for junction in self.junctions:
            if all(
                x is not None
                for x in [
                    junction.jmin_coordinate,
                    junction.jmax_coordinate,
                    junction.design_region,
                    junction.design_start,
                ]
            ):
                logger.debug(f"Junction: {junction.name}")
                logger.debug(
                    f"Genomic coordinates: {junction.chrom}:{junction.start}-{junction.end}"
                )
                logger.debug(
                    f"Design region: {junction.chrom}:{junction.design_start}-{junction.design_end}"
                )
                logger.debug(f"Design region length: {junction.junction_length}")
                logger.debug(
                    f"Junction in design region: {junction.jmin_coordinate}-{junction.jmax_coordinate}"
                )

                # Extract the junction sequence for verification
                if (
                    0
                    <= junction.jmin_coordinate
                    < junction.jmax_coordinate
                    <= len(junction.design_region)
                ):
                    junction_seq = junction.design_region[
                        junction.jmin_coordinate : junction.jmax_coordinate + 1
                    ]
                    logger.debug(f"Junction sequence: {junction_seq[:50]}...")
                else:
                    logger.warning(
                        f"Junction coordinates out of bounds for {junction.name}!"
                    )

        logger.info("Coordinate verification complete.")

    def create_junction_table(self) -> pd.DataFrame:
        """
        Create pandas DataFrame with all junction information.

        Args:
            self: A MultiplexPanel object containing junctions.
        """
        data = []
        for junction in self.junctions:
            row = {
                "Name": junction.name,
                "Chrom": junction.chrom,
                "Five_Prime_Coordinate": junction.start,
                "Three_Prime_Coordinate": junction.end,
                "Design_Region": junction.design_region
                if junction.design_region
                else "",
                "Design_Start": junction.design_start if junction.design_start else "",
                "Design_End": junction.design_end if junction.design_end else "",
                "Junction_Length": junction.junction_length
                if junction.junction_length
                else "",
                "Jmin_Coordinate": junction.jmin_coordinate
                if junction.jmin_coordinate
                else "",
                "Jmax_Coordinate": junction.jmax_coordinate
                if junction.jmax_coordinate
                else "",
            }
            data.append(row)

        self.junction_df = pd.DataFrame(data)
        return self.junction_df

    def write_junctions_to_csv(self, file_path: str):
        """
        Write the junctions from the MultiplexPanel object to a CSV file.

        Args:
            file_path (str): Path to the output CSV file.

        Raises:
            ValueError: If no junctions are available or if the file path is invalid.
            IOError: If there is an issue writing to the file.
        """
        # Check if there are any junctions to write
        if not hasattr(self, "junctions") or not self.junctions:
            error_msg = "No junctions available to write to CSV."
            logger.error(error_msg)
            raise ValueError(error_msg)

        # Create a DataFrame with all junction information
        try:
            junction_df = self.create_junction_table()
        except Exception as e:
            error_msg = f"Failed to create junction table: {e}"
            logger.error(error_msg)
            raise RuntimeError(error_msg) from e

        # Check if the DataFrame is empty
        if junction_df.empty:
            error_msg = "No junction data available to write to CSV."
            logger.error(error_msg)
            raise ValueError(error_msg)

        # Check if the file path is valid
        if not file_path:
            error_msg = "File path cannot be empty."
            logger.error(error_msg)
            raise ValueError(error_msg)

        # Ensure the directory exists
        dir_path = os.path.dirname(file_path)
        if dir_path and not os.path.exists(dir_path):
            error_msg = f"Directory does not exist: {dir_path}"
            logger.error(error_msg)
            raise OSError(error_msg)

        # Write the DataFrame to a CSV file
        try:
            junction_df.to_csv(file_path, index=False)
            logger.info(f"Junctions successfully written to {file_path}")
        except Exception as e:
            error_msg = f"Failed to write junctions to {file_path}: {e}"
            logger.error(error_msg)
            raise OSError(error_msg) from e

    def aggregate_primers(self) -> dict:
        """
        Collect all unique primers from all primer pairs in the panel.

        Returns:
            dict: A dictionary mapping primer sequences to a list of primer names (or IDs).
                  Also creates a simplified mapping of sequence -> synthetic ID (SEQ_0, SEQ_1, ...)
                  stored in self.unique_primer_map for downstream use.
        """
        seq_to_names = {}
        count_pairs = 0

        for junction in self.junctions:
            if not junction.primer_pairs:
                continue
            for pair in junction.primer_pairs:
                count_pairs += 1
                # Forward
                if pair.forward.seq not in seq_to_names:
                    seq_to_names[pair.forward.seq] = []
                seq_to_names[pair.forward.seq].append(pair.forward.name)

                # Reverse
                if pair.reverse.seq not in seq_to_names:
                    seq_to_names[pair.reverse.seq] = []
                seq_to_names[pair.reverse.seq].append(pair.reverse.name)

        logger.info(
            f"Aggregated {len(seq_to_names)} unique primer sequences from {count_pairs} candidate pairs."
        )

        # Create a stable synthetic ID map for these unique sequences
        self.unique_primer_map = {}
        for i, seq in enumerate(sorted(seq_to_names.keys())):
            self.unique_primer_map[seq] = f"SEQ_{i}"

        return seq_to_names

    def save_candidate_primers_to_fasta(self, file_path: str):
        """
        Save all unique candidate primers to a FASTA file.
        Uses synthetic IDs (SEQ_0, SEQ_1, ...) to avoid duplicates.

        Args:
            file_path: Path to the output FASTA file.
        """
        if not hasattr(self, "unique_primer_map") or not self.unique_primer_map:
            self.aggregate_primers()

        try:
            with open(file_path, "w") as f:
                for seq, uid in self.unique_primer_map.items():
                    f.write(f">{uid}\n{seq}\n")
            logger.info(f"Saved unique primers to {file_path}")
        except Exception as e:
            logger.error(f"Failed to save primers to FASTA: {e}")
            raise

    def save_candidate_pairs_to_csv(self, file_path: str):
        """
        Save all candidate primer pairs to a CSV file.

        Args:
            file_path: Path to the output CSV file.
        """
        data = []
        for junction in self.junctions:
            if not junction.primer_pairs:
                continue
            for pair in junction.primer_pairs:
                row = {
                    "Junction": junction.name,
                    "Pair_ID": pair.pair_id,
                    "Forward_Name": pair.forward.name,
                    "Reverse_Name": pair.reverse.name,
                    "Forward_Seq": pair.forward.seq,
                    "Reverse_Seq": pair.reverse.seq,
                    "Amplicon_Length": pair.amplicon_length,
                    "Penalty": pair.pair_penalty,
                    "Insert_Size": pair.insert_size,
                }
                data.append(row)

        if not data:
            logger.warning("No candidate pairs to save.")
            return

        try:
            df = pd.DataFrame(data)
            df.to_csv(file_path, index=False)
            logger.info(f"Saved {len(df)} candidate pairs to {file_path}")
        except Exception as e:
            logger.error(f"Failed to save pairs to CSV: {e}")
            raise


def panel_factory(
    name: str,
    genome: str,
    design_input_file: str,
    fasta_file: str,
    preset: str = "default",
    config_file: str = None,
    padding: int = 200,
) -> "MultiplexPanel":
    """
    Create and configure a MultiplexPanel ready for primer design.

    Args:
        name: Required. Name of the panel.
        genome: Required. Name of the reference genome to pull.
        design_input_file: Required. Samplesheet containing junction around which to design primers.
        fasta_file: Required. Reference fasta file.
        preset: Optional. Which preset configuration to use ("default" or "lenient").
        config_file: Optional. Configuration file to use, otherwise uses the default config.
        padding: Optional. How many bases around the junction should be extracted for the design region? Default is 200.

    Returns:
        A configured MultiplexPanel object.
    """
    logger.info(f"Creating multiplex panel: {name}")

    # Create and configure the multiplex panel object
    panel = MultiplexPanel(name, genome)

    # Load config file. This handles defaults internally if none is provided by the user.
    # User provided configs overwrite all defaults.
    panel.load_config(preset=preset, config_path=config_file)

    # Retrieve design regions and calculate junctions
    panel.import_junctions_csv(file_path=design_input_file)
    panel.merge_close_junctions()
    panel.extract_design_regions_from_fasta(fasta_file, padding=padding)
    panel.calculate_junction_coordinates_in_design_region()

    return panel
