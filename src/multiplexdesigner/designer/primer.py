# ================================================================================
# Primer classes and associated functions
#
# Author: Stefan Filges (stefan@simsendiagnostics.com)
# Copyright (c) 2025 Simsen Diagnostics AB
# ================================================================================

from dataclasses import dataclass, field

from multiplexdesigner.designer.thal import seqtm
from multiplexdesigner.utils.utils import create_primer_dataframe, gc_content


@dataclass
class Primer:
    """
    Define a single primer.
    """

    name: str
    seq: str
    direction: str
    start: int
    length: int
    bound: float = None
    tm: float = None
    tm_primer3: float = None
    gc: float = None
    penalty: float = None
    self_any_th: float = None
    self_end_th: float = None
    hairpin_th: float = None
    end_stability: float = None
    engine: str = None

    def add_tail(self, tail_seq: str, tail_direction: str = "five_prime"):
        """
        Add tail sequence to primer.
        """
        if tail_direction == "five_prime":
            self.seq = f"{tail_seq}{self.seq}"
        elif tail_direction == "three_prime":
            self.seq = f"{self.seq}{tail_seq}"


@dataclass
class PrimerPair:
    """
    Define a pair of primers.

    Args:
        forward: Forward Primer object
        reverse: Reverse Primer object
        insert_size: Length of the inter-primer region
        amplicon_sequence: Sequence of the amplicon
        amplicon_length: Length of the total amplicon, incl. primer sequences
        pair_penalty: Total primer pair penalty
        pair_id: Unique ID of the primer pair
    """

    forward: object
    reverse: object
    insert_size: int
    amplicon_sequence: str
    amplicon_length: int
    pair_penalty: float = None
    pair_id: str = field(default="", repr=False)
    specificity_checked: bool = False
    off_target_products: list = field(default_factory=list)
    dimer_score: float = None
    selected: bool = False

    @staticmethod
    def calculate_primer_pair_penalty_th(
        primer_left_penalty,
        primer_right_penalty,
        primer_left_tm,
        primer_right_tm,
        product_size,
        product_opt_size: int = 0,
        wt_pr_penalty=1.0,
        wt_product_size_gt=0.0,
        wt_product_size_lt=0.0,
        wt_diff_tm=0.0,
    ):
        """
        Calculate primer pair penalty using thermodynamic approach. Based
        on primer3 docs.

        Parameters:
        -----------
        primer_left_penalty : float
            Penalty for left primer
        primer_right_penalty : float
            Penalty for right primer
        primer_left_tm : float
            Melting temperature of left primer
        primer_right_tm : float
            Melting temperature of right primer
        product_tm : float
            Product melting temperature
        product_size : int
            Product size in base pairs
        compl_any_th : float
            Thermodynamic complementarity (any position)
        compl_end_th : float
            Thermodynamic complementarity (end position)
        template_mispriming_th : float
            Thermodynamic template mispriming value
        product_opt_tm : float
            Optimal product melting temperature
        product_opt_size : int
            Optimal product size
        wt_* : float
            Various weight parameters (default values shown)
        use_internal_oligo : bool
            Whether internal oligo is used

        Returns:
        --------
        float : Total penalty for the primer pair
        """

        penalty = 0.0

        # Add single primer penalties
        penalty += wt_pr_penalty * (primer_left_penalty + primer_right_penalty)

        # Product size penalties. This parameter influences primer pair selection
        # only if wt_product_size_gt or wt_product_size_lt are non-0. Standard
        # primer3 does not use this, i.e. product_opt_size = 0. But some ctDNA
        # it should be beneficial to select short amplicons and short amplicons
        # should be favoured.
        if product_opt_size > 0:
            if product_size > product_opt_size:
                penalty += wt_product_size_gt * (product_size - product_opt_size)
            elif product_size < product_opt_size:
                penalty += wt_product_size_lt * (product_opt_size - product_size)

        # Difference in Tm between primers
        penalty += wt_diff_tm * abs(primer_left_tm - primer_right_tm)

        return penalty


# TODO: import primer pairs from primer3 output. Wouldn't it be better to design
# forward and reverse primers separately instead of in pairs, and only consider pairing
# while also checking for primer dimers?


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
    directions = ["LEFT", "RIGHT"]

    for row in primer3_output:
        for d in directions:
            pair = {}
            pair[d] = Primer(
                seq=row["LEFT"],
                direction="F" if d == "LEFT" else "R",
                start=1,
                length=2,
                tm=1,
                gc=1,
            )
        primer_pair = PrimerPair(F=pair["LEFT"], R=pair["RIGHT"])
        primer_pairs.append(primer_pair)

    return primer_pairs


def get_primer_dict(junction):
    """
    Extract unique primer pairs from a junction.

    Parameters:
    df (pandas.DataFrame): DataFrame containing primer pair data generated with primer3.

    Returns:
    tuple: A tuple containing:
        - list of unique PrimerPair objects
        - dictionary mapping primer sequences and directions to Primer objects
    """

    df = create_primer_dataframe(junction.primer3_designs)

    primer_dict = {}  # Maps (sequence, direction) to Primer object

    for index, row in df.iterrows():
        # FORWARD
        left_seq = row["left_sequence"]
        left_direction = "forward"
        left_key = (left_seq, left_direction)
        if left_key not in primer_dict:
            left_tm_bound = seqtm(left_seq)
            left_primer = Primer(
                name=f"{junction.name}_{index}_forward",
                seq=left_seq,
                direction=left_direction,
                start=row["left_coords"][0],
                length=row["left_coords"][1],
                tm_primer3=round(row["left_tm"], 2),
                tm=round(left_tm_bound.Tm, 2),
                bound=round(left_tm_bound.bound, 2),
                gc=round(gc_content(left_seq), 2),
                penalty=row["left_penalty"],
            )
            primer_dict[left_key] = left_primer

        # REVERSE
        right_seq = row["right_sequence"]
        right_direction = "reverse"
        right_key = (right_seq, right_direction)
        if right_key not in primer_dict:
            right_tm_bound = seqtm(right_seq)
            right_primer = Primer(
                name=f"{junction.name}_{index}_reverse",
                seq=right_seq,
                direction=right_direction,
                start=row["right_coords"][0],
                length=row["right_coords"][1],
                tm_primer3=round(row["right_tm"], 2),
                tm=round(right_tm_bound.Tm, 2),
                bound=round(right_tm_bound.bound, 2),
                gc=round(gc_content(right_seq), 2),
                penalty=row["right_penalty"],
            )
            primer_dict[right_key] = right_primer

    return primer_dict
