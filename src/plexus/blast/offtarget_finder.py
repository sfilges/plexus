from collections import namedtuple

import pandas as pd


class AmpliconFinder:
    """
    Find all potential amplicons between a list of primer binding sites.

    This class takes a DataFrame of predicted primer binding sites (usually from BLAST)
    and identifies pairs of forward and reverse primers that could form an amplicon
    within a specified distance on the same chromosome.

    Expected bound_df schema (BLAST format 6 + annotations):
        - qseqid: Query (primer) sequence ID
        - sseqid: Subject (chromosome) sequence ID
        - sstart: Start position of alignment on subject
        - send: End position of alignment on subject
        - sstrand: Strand of alignment ('plus' or 'minus')
        - qlen: Query length
        - qend: End position of alignment on query
        - predicted_bound: Boolean indicating if the primer is predicted to bind

    Generated amplicon_df schema:
        - chrom: Chromosome ID (sseqid)
        - F_target: Name of the target for the forward primer
        - R_target: Name of the target for the reverse primer
        - F_primer: ID of the forward primer
        - R_primer: ID of the reverse primer
        - F_start: 5' genomic position of the forward primer
        - R_start: 5' genomic position of the reverse primer (on minus strand)
        - product_bp: Predicted size of the amplicon in base pairs (inclusive)
        - F_pident: Percent identity of the forward primer alignment
        - R_pident: Percent identity of the reverse primer alignment
        - F_mismatch: Number of mismatches in the forward primer alignment
        - R_mismatch: Number of mismatches in the reverse primer alignment
        - F_align_len: Alignment length of the forward primer (bp)
        - R_align_len: Alignment length of the reverse primer (bp)
        - F_evalue: E-value of the forward primer alignment
        - R_evalue: E-value of the reverse primer alignment
    """

    def __init__(self, bound_df, target_map=None):
        """
        Initialize the AmpliconFinder.

        Args:
            bound_df (pd.DataFrame): DataFrame of predicted primer binding sites.
            target_map (dict, optional): Mapping from primer IDs to target names.
        """

        # Save inputs
        self.bound_df = bound_df
        self.target_map = target_map or {}

        # Generated
        self.amplicon_df = None

    def find_amplicons(self, max_size_bp=6000):
        """
        Identify all potential amplicons within max_size_bp.

        Args:
            max_size_bp (int): Maximum predicted product size to consider.
        """

        # Storage -- can make flexible
        Amplicon = namedtuple(
            "Amplicon",
            [
                "chrom",
                "F_target",
                "R_target",
                "F_primer",
                "R_primer",
                "F_start",
                "R_start",
                "product_bp",
                "F_pident",
                "R_pident",
                "F_mismatch",
                "R_mismatch",
                "F_align_len",
                "R_align_len",
                "F_evalue",
                "R_evalue",
            ],
        )

        # Iterate and find amplicons
        amplicons = []
        for chrom, chrom_df in self.bound_df.groupby("sseqid"):
            for _, F_df in chrom_df.query("sstrand == 'plus'").iterrows():
                # Get start position of the forward primer
                # 5' position
                F_start = int(F_df["sstart"])

                # Get proximal reverse primers, if they exist
                # - Very important to get directionality of query correct
                R_df = chrom_df.query(
                    "sstrand == 'minus' and 0 < (sstart - @F_start) < @max_size_bp"
                )

                # Check for any pairs
                if R_df.shape[0] > 0:
                    F_amplicons = [
                        Amplicon(
                            chrom=chrom,
                            F_target=self.target_map.get(
                                F_df["qseqid"], F_df["qseqid"].split("_")[0]
                            ),
                            R_target=self.target_map.get(
                                row["qseqid"], row["qseqid"].split("_")[0]
                            ),
                            F_primer=F_df["qseqid"],
                            R_primer=row["qseqid"],
                            F_start=F_start,
                            R_start=row["sstart"],
                            product_bp=row["sstart"] - F_start + 1,
                            F_pident=F_df.get("pident", None),
                            R_pident=row.get("pident", None),
                            F_mismatch=F_df.get("mismatch", None),
                            R_mismatch=row.get("mismatch", None),
                            F_align_len=F_df.get("length", None),
                            R_align_len=row.get("length", None),
                            F_evalue=F_df.get("evalue", None),
                            R_evalue=row.get("evalue", None),
                        )
                        for _, row in R_df.iterrows()
                    ]
                    amplicons.extend(F_amplicons)

        # Store
        self.amplicon_df = pd.DataFrame(amplicons)
