# ================================================================================
# Primer classes and associated functions
#
# Author: Stefan Filges (stefan.filges@pm.me)
# Copyright (c) 2025 Stefan Filges
# ================================================================================

from dataclasses import dataclass, field


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
    snp_count: int = 0

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
    on_target_detected: bool | None = None
    dimer_score: float = None
    selected: bool = False
    snp_count: int = 0
    snp_penalty: float = 0.0

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
