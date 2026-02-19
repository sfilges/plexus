# Computes potential dimers betweeen pairs of primers
#
# The cross dimer or primer dimer check is an important design step to optimize
# primer performance in multiplex reactions.
#
# Examples:
# - Oli2go uses Primer3's ntthal and user-defined ΔG/Tm values
# - multiply uses Johnston et al. (2019) Sci Reports algorithm
#   https://www.nature.com/articles/s41598-018-36612-9
#
# Both use thermodynamic alignment values to estimate primer dimer alignment.
#
# Workflow: Forward and reverse primer pairs from the design task are the input.
# Starting with the sequence having fewest specific primers, these are checked
# against all other possible primers. Results involve primer pairs that don't
# exceed cross dimerization thresholds. If results contain at least one primer
# pair per sequence, each is checked against other primers in results. Finally,
# for each input sequence one primer pair forming no cross dimerization is returned.

from __future__ import annotations

import json
from dataclasses import dataclass, field
from itertools import product

from loguru import logger

from plexus.utils.root_dir import ROOT_DIR

# ================================================================================
# Define an alignment between two primers
# ================================================================================


@dataclass(order=True)
class PrimerAlignment:
    """
    Represent the alignment of two primers
    """

    primer1_name: str = field(compare=False)
    primer2_name: str = field(compare=False)
    primer1: str = field(compare=False)
    primer2: str = field(compare=False)
    score: float = field(compare=True)
    alignment: str = field(compare=False, repr=False)


# ================================================================================
# Primer alignment algorithms
# ================================================================================


class PrimerDimerPredictor:
    """
    Align two primers using an algorithm like the one described
    by Johnston et al. (2019) Sci Reports: https://www.nature.com/articles/s41598-018-36612-9

    Primary idea is to only allow:
    - Ungapped alignments
    - With 5' overhangs (i.e. extensible)

    And then to add a bonus if either 3' end is complementary in
    the highest scoring alignment
    """

    # Class-level cache: param_path -> (nn_scores, end_length, end_bonus)
    _param_cache: dict[str, tuple[dict[str, float], int, float]] = {}

    def __init__(self, param_path: str | None = None) -> None:
        """
        Initialize the PrimerDimerPredictor aligner.

        Parameters
        ----------
        param_path : str, optional
            Path to the parameters JSON file. If None, uses default path.
        """
        self.primer1: str | None = None
        self.primer2: str | None = None
        self.primer1_name: str | None = None
        self.primer2_name: str | None = None
        self.score: float | None = None
        self.best_start: int = 0
        self.best_matching: list[bool] = []
        self.shorter: str = ""
        self.longer: str = ""
        self.nn_scores: dict[str, float] = {}
        self.end_length: int = 0
        self.end_bonus: float = 0.0

        # Load parameters
        if param_path is None:
            param_path = f"{ROOT_DIR}/config/alignment_parameters.json"
        self.load_parameters(param_path)

    def set_primers(
        self, primer1: str, primer2: str, primer1_name: str, primer2_name: str
    ) -> None:
        """Set a pair of primers to align"""
        self.primer1 = primer1
        self.primer2 = primer2
        self.primer1_name = primer1_name
        self.primer2_name = primer2_name
        self.score = None

    def load_parameters(self, param_path: str) -> None:
        """
        Load parameters necessary for Primer Dimer algorithm,
        and set as attributes. Results are cached by path so that
        repeated instantiations within the same process only read
        the JSON files once.

        Parameters
        ----------
        param_path : str
            Path to the parameters JSON file
        """
        if param_path in PrimerDimerPredictor._param_cache:
            logger.debug(f"Using cached alignment parameters from: {param_path}")
            self.nn_scores, self.end_length, self.end_bonus = (
                PrimerDimerPredictor._param_cache[param_path]
            )
            return

        logger.info(f"Loading alignment parameters from: {param_path}")

        # Load parameter JSON
        with open(param_path) as f:
            params = json.load(f)

        # Load nearest neighbour model, these should all be paths
        self.nn_scores = create_nn_score_dt(
            match_json=f"{ROOT_DIR}/{params['match_scores']}",
            single_mismatch_json=f"{ROOT_DIR}/{params['single_mismatch_scores']}",
            double_mismatch_score=params["double_mismatch_score"],
        )

        # Load penalties
        self.end_length = params["end_length"]
        self.end_bonus = params["end_bonus"]

        PrimerDimerPredictor._param_cache[param_path] = (
            self.nn_scores,
            self.end_length,
            self.end_bonus,
        )

    @staticmethod
    def _calc_linear_extension_bonus(
        matching: list[bool],
        overhang_left: bool,
        overhang_right: bool,
        end_length: int,
        end_bonus: float,
    ) -> float:
        """
        Calculate a bonus score for primer-dimer alignments that would allow for *extension*

        Parameters
        ----------
        matching : list[bool]
            List of booleans indicating whether or not bases matched
            for this alignment position
        overhang_left : bool
            Is there an overhang on the left side of the matches;
            i.e. would extension be possible?
        overhang_right : bool
            As above, but right side.
        end_length : int
            Number of bases to consider for end bonus.
        end_bonus : float
            Bonus to add per aligned, extendible base.

        Returns
        -------
        float
            Bonus score in [-2*end_length*end_bonus, 0].
        """

        # Left end
        left_end = 0
        for match in matching[:end_length]:
            if not overhang_left:
                break
            if not match:
                break
            left_end += 1

        # Right end
        right_end = 0
        for match in matching[::-1][:end_length]:
            if not overhang_right:
                break
            if not match:
                break
            right_end += 1

        return (left_end + right_end) * end_bonus

    def _validate_primers_set(self) -> None:
        """
        Validate that primers have been set before alignment.

        Raises
        ------
        ValueError
            If primers have not been set via set_primers()
        """
        if self.primer1 is None or self.primer2 is None:
            raise ValueError(
                "Primers must be set before alignment. Call set_primers() first."
            )
        if not self.primer1 or not self.primer2:
            raise ValueError("Primers cannot be empty strings.")

    def align(self) -> None:
        """
        Align primers, finding highest score and its associated start position.

        Raises
        ------
        ValueError
            If primers have not been set via set_primers()
        """
        # Validate primers are set
        self._validate_primers_set()

        # Reverse complement mapping dictionary
        rc_map: dict[str, str] = {"A": "T", "T": "A", "C": "G", "G": "C"}

        # Identify longer and shorter primer
        # After validation, primer1 and primer2 are guaranteed to be str
        assert self.primer1 is not None and self.primer2 is not None
        primers = [self.primer1, self.primer2]
        primers.sort(key=len, reverse=True)
        longer, shorter = primers
        shorter = shorter[::-1]
        n_longer, n_shorter = len(longer), len(shorter)

        # Iterate over each start position
        best_start = 0
        best_score: float = 10
        best_matching: list[bool] = []
        j = 0  # Initialize j for use after inner loop
        for i in range(n_longer - 1):
            # Compute score for alignment
            matching: list[bool] = []
            current_score: float = 0
            for j in range(n_shorter - 1):
                # Compute pseudo-Gibb's
                longer_bases = longer[(i + j) : (i + j + 2)]  # 2 bases
                shorter_bases = shorter[j : (j + 2)]  # 2 bases
                nn = f"{longer_bases}/{shorter_bases}"  # Set base comparison map
                current_score += self.nn_scores[nn]

                # Compute if matching, for left-end
                matching.append(rc_map[shorter[j]] == longer[i + j])

                # Stop if reached last dinucleotide of longer primer
                if i + j == n_longer - 2:
                    break

            # Add matching state for last nucleotide
            matching.append(rc_map[shorter[j + 1]] == longer[i + j + 1])

            # Compute end bonus
            current_score += self._calc_linear_extension_bonus(
                matching,
                overhang_left=i > 0,
                overhang_right=len(matching) < n_shorter,
                end_length=self.end_length,
                end_bonus=self.end_bonus,
            )

            # Update if this is the new best score
            if current_score <= best_score:
                best_score = current_score
                best_start = i
                best_matching = matching

        # Assign to instance variables
        self.score = best_score
        self.best_start = best_start

        # Helps with getting alignment string
        self.best_matching = best_matching
        self.shorter = shorter
        self.longer = longer

    def get_alignment_string(self) -> str:
        """
        Return a string representing the best alignment between the two primers.

        Returns
        -------
        str
            ASCII representation of the alignment

        Raises
        ------
        ValueError
            If align() has not been called yet
        """
        if self.score is None:
            raise ValueError("No alignment available. Call align() first.")

        # Recover names
        if self.longer == self.primer1:
            longer_name = self.primer1_name
            shorter_name = self.primer2_name
        else:
            longer_name = self.primer2_name
            shorter_name = self.primer1_name

        # Create space, regardless of length of primer names
        assert self.primer1_name is not None and self.primer2_name is not None
        name_max = max([len(self.primer1_name), len(self.primer2_name)])
        str_template = "{:>%d}    {}\n" % name_max  # noqa: UP031

        # Create individual strings
        gap = " " * self.best_start
        longer_str = f"5-{self.longer}-3"
        match_str = f"{gap}  {''.join(['|' if m else ' ' for m in self.best_matching])}"
        shorter_str = f"{gap}3-{self.shorter}-5"

        align_str = f"Dimer Score: {self.score:.03f}\n"
        align_str += str_template.format(longer_name, longer_str)
        align_str += str_template.format("", match_str)
        align_str += str_template.format(shorter_name, shorter_str)
        align_str += "\n"

        return align_str

    def print_alignment(self) -> None:
        """Print an ASCII view of the aligned primers"""
        print(self.get_alignment_string())

    def get_primer_alignment(self) -> PrimerAlignment:
        """
        Return an alignment object.

        Returns
        -------
        PrimerAlignment
            Object containing alignment results

        Raises
        ------
        ValueError
            If align() has not been called yet
        """
        if self.score is None:
            raise ValueError("No alignment available. Call align() first.")

        assert self.primer1 is not None and self.primer2 is not None
        assert self.primer1_name is not None and self.primer2_name is not None

        return PrimerAlignment(
            primer1=self.primer1,
            primer2=self.primer2,
            primer1_name=self.primer1_name,
            primer2_name=self.primer2_name,
            score=self.score,
            alignment=self.get_alignment_string(),
        )


def create_nn_score_dt(
    match_json: str, single_mismatch_json: str, double_mismatch_score: float = 0.2
) -> dict[str, float]:
    """
    Create the Gibb's free energy nearest neighbour scoring dictionary

    Parameters
    ----------
    match_json : str
        Path to JSON file containing match scores
    single_mismatch_json : str
        Path to JSON file containing single mismatch scores
    double_mismatch_score : float, optional
        Score for double mismatches (default: 0.2)

    Returns
    -------
    dict
        Dictionary mapping dinucleotide pairs to their scores
    """
    # Load match and single mismatch .jsons
    with open(match_json) as f:
        match_dt: dict[str, float] = json.load(f)
    with open(single_mismatch_json) as f:
        single_mismatch_dt: dict[str, float] = json.load(f)

    # Set all as double mismatches; then update
    nts = ["A", "T", "C", "G"]
    nn_score_dt: dict[str, float] = {
        "".join(watson) + "/" + "".join(crick): double_mismatch_score
        for watson in product(nts, repeat=2)
        for crick in product(nts, repeat=2)
    }

    # Update
    nn_score_dt.update(match_dt)
    nn_score_dt.update(single_mismatch_dt)

    return nn_score_dt
