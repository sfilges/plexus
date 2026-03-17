import numpy as np
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

    def _get_col_or_none(self, df, col):
        """Return column values as array, or None if column is missing."""
        if col in df.columns:
            return df[col].values
        return None

    def find_amplicons(self, max_size_bp=6000):
        """
        Identify all potential amplicons within max_size_bp.

        Uses sorted arrays and np.searchsorted for O(F log R) per chromosome
        instead of O(F × R) nested iteration.

        Args:
            max_size_bp (int): Maximum predicted product size to consider.
        """

        amplicon_columns = [
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
        ]
        result_chunks = []
        target_map = self.target_map

        for chrom, chrom_df in self.bound_df.groupby("sseqid"):
            fwd = chrom_df[chrom_df["sstrand"] == "plus"]
            rev = chrom_df[chrom_df["sstrand"] == "minus"]

            if fwd.empty or rev.empty:
                continue

            # Sort reverse hits by sstart for searchsorted
            rev = rev.sort_values("sstart")
            r_starts = rev["sstart"].values
            r_qseqids = rev["qseqid"].values
            r_pident = self._get_col_or_none(rev, "pident")
            r_mismatch = self._get_col_or_none(rev, "mismatch")
            r_length = self._get_col_or_none(rev, "length")
            r_evalue = self._get_col_or_none(rev, "evalue")

            f_starts = fwd["sstart"].values
            f_qseqids = fwd["qseqid"].values
            f_pident = self._get_col_or_none(fwd, "pident")
            f_mismatch = self._get_col_or_none(fwd, "mismatch")
            f_length = self._get_col_or_none(fwd, "length")
            f_evalue = self._get_col_or_none(fwd, "evalue")

            for i in range(len(f_starts)):
                f_start = f_starts[i]
                # Find reverse hits where 0 < (r_start - f_start) < max_size_bp
                lo = np.searchsorted(
                    r_starts, f_start, side="right"
                )  # r_start > f_start
                hi = np.searchsorted(
                    r_starts, f_start + max_size_bp, side="left"
                )  # r_start < f_start + max_size_bp

                if lo >= hi:
                    continue

                n_matches = hi - lo
                f_qseqid = f_qseqids[i]
                f_target = target_map.get(f_qseqid, f_qseqid.split("_")[0])

                matched_r_starts = r_starts[lo:hi]
                matched_r_qseqids = r_qseqids[lo:hi]

                r_targets = [
                    target_map.get(rq, rq.split("_")[0]) for rq in matched_r_qseqids
                ]

                chunk = {
                    "chrom": np.full(n_matches, chrom),
                    "F_target": np.full(n_matches, f_target),
                    "R_target": r_targets,
                    "F_primer": np.full(n_matches, f_qseqid),
                    "R_primer": matched_r_qseqids,
                    "F_start": np.full(n_matches, f_start, dtype=int),
                    "R_start": matched_r_starts,
                    "product_bp": matched_r_starts - f_start + 1,
                    "F_pident": np.full(
                        n_matches, f_pident[i] if f_pident is not None else None
                    ),
                    "R_pident": r_pident[lo:hi]
                    if r_pident is not None
                    else np.full(n_matches, None),
                    "F_mismatch": np.full(
                        n_matches, f_mismatch[i] if f_mismatch is not None else None
                    ),
                    "R_mismatch": r_mismatch[lo:hi]
                    if r_mismatch is not None
                    else np.full(n_matches, None),
                    "F_align_len": np.full(
                        n_matches, f_length[i] if f_length is not None else None
                    ),
                    "R_align_len": r_length[lo:hi]
                    if r_length is not None
                    else np.full(n_matches, None),
                    "F_evalue": np.full(
                        n_matches, f_evalue[i] if f_evalue is not None else None
                    ),
                    "R_evalue": r_evalue[lo:hi]
                    if r_evalue is not None
                    else np.full(n_matches, None),
                }
                result_chunks.append(pd.DataFrame(chunk, columns=amplicon_columns))

        # Store
        if result_chunks:
            self.amplicon_df = pd.concat(result_chunks, ignore_index=True)
            del result_chunks
        else:
            self.amplicon_df = pd.DataFrame(columns=amplicon_columns)
