# Computes potential dimers betweeen pairs of primers
#
# The cross dimer or primer dimer check is an important design step to optimize primer performance in multiplex reactions. 
# What algorithms to use?
#
# Examples are:
# - Oli2go uses Primer3's ntthal and the user-defined ΔG and Tm values to check for cross dimerization.
# - multiply uses the algorithm by by Johnston et al. (2019) Sci Reports: https://www.nature.com/articles/s41598-018-36612-9
#
# Both use thermodynamic alignment values to estimate primer dimer alignment.
#
# Specific forward and reverse primer pairs resulting from the preceding design task form the input for this workflow step. It starts with 
# the input sequence that has fewest specific primers. These primers are checked against all other possible primers of the 
# other input sequences. The first results involve primer pairs which do not exceed the cross dimerization thresholds. If the 
# results contain at least one primer pair for each sequence, each one is checked against the other primers in the results. 
# Finally, for each input sequence one primer pair forming no cross dimerization with all other sequences is returned.

import json
from itertools import product
from dataclasses import dataclass, field
from multiplexdesigner.utils.root_dir import ROOT_DIR

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

    def __init__(self, param_path=None):
        """
        Initialize the PrimerDimerPredictor aligner.
        
        Parameters
        ----------
        param_path : str, optional
            Path to the parameters JSON file. If None, uses default path.
        """
        self.primer1 = None
        self.primer2 = None
        self.primer1_name = None
        self.primer2_name = None
        self.score = None
        
        # Load parameters
        if param_path is None:
            param_path = f"{ROOT_DIR}/config/alignment_parameters.json"
        self.load_parameters(param_path)

    def set_primers(self, primer1, primer2, primer1_name, primer2_name):
        """Set a pair of primers to align"""
        self.primer1 = primer1
        self.primer2 = primer2
        self.primer1_name = primer1_name
        self.primer2_name = primer2_name
        self.score = None
        
    def load_parameters(self, param_path):
        """
        Load parameters necessary for Primer Dimer algorithm,
        and set as attributes
        
        Parameters
        ----------
        param_path : str
            Path to the parameters JSON file
        """
        # Load parameter JSON
        with open(param_path, "r") as f:
            params = json.load(f)

        # Load nearest neighbour model, these should all be paths
        self.nn_scores = create_nn_score_dt(
            match_json=f"{ROOT_DIR}/{params['match_scores']}",
            single_mismatch_json=f"{ROOT_DIR}/{params['single_mismatch_scores']}",
            double_mismatch_score=params['double_mismatch_score']
        )

        # Load penalties
        self.end_length = params["end_length"]
        self.end_bonus = params["end_bonus"]
        
    @staticmethod
    def _calc_linear_extension_bonus(matching, overhang_left, overhang_right, end_length, end_bonus):
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
    
    def _validate_primers_set(self):
        """
        Validate that primers have been set before alignment.
        
        Raises
        ------
        ValueError
            If primers have not been set via set_primers()
        """
        if self.primer1 is None or self.primer2 is None:
            raise ValueError(
                "Primers must be set before alignment. "
                "Call set_primers() first."
            )
        if not self.primer1 or not self.primer2:
            raise ValueError("Primers cannot be empty strings.")
        
    def align(self):
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
        rc_map = {"A": "T", "T": "A", "C": "G", "G": "C"}
        
        # Identify longer and shorter primer
        primers = [self.primer1, self.primer2]
        primers.sort(key=len, reverse=True)
        l, s = primers
        s = s[::-1]
        nL, nS = len(l), len(s)

        # Iterate over each start position
        best_start = 0
        best_score = 10
        best_matching = []
        for i in range(nL - 1):

            # Compute score for alignment
            matching = []
            current_score = 0
            for j in range(nS - 1):

                # Compute pseudo-Gibb's
                l_bases = l[(i+j):(i+j+2)]   # 2 bases
                s_bases = s[j:(j+2)]         # 2 bases
                nn = f"{l_bases}/{s_bases}"  # Set base comparison map
                current_score += self.nn_scores[nn]

                # Compute if matching, for left-end
                matching.append(rc_map[s[j]] == l[i+j])

                # Stop if reached last dinucleotide of longer primer
                if i + j == nL - 2:
                    break

            # Add matching state for last nucleotide
            matching.append(rc_map[s[j+1]] == l[i+j+1])

            # Compute end bonus
            current_score += self._calc_linear_extension_bonus(
                matching,
                overhang_left=i > 0,
                overhang_right=len(matching) < nS, 
                end_length=self.end_length,
                end_bonus=self.end_bonus
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
        self.s = s
        self.l = l
        
    def get_alignment_string(self):
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
            raise ValueError(
                "No alignment available. Call align() first."
            )

        # Recover names, kinda ugly
        if self.l == self.primer1:
            lname = self.primer1_name
            sname = self.primer2_name
        else:
            lname = self.primer2_name
            sname = self.primer1_name

        # Create space, regardless of length of primer names
        name_max = max([len(self.primer1_name), len(self.primer2_name)])
        str_template = "{:>%d}    {}\n" % name_max

        # Create individual strings
        gap = " "*self.best_start
        lstr = f"5-{self.l}-3"
        mstr = f"{gap}  {''.join(['|' if m else ' ' for m in self.best_matching])}"
        sstr = f"{gap}3-{self.s}-5"

        align_str = f"Dimer Score: {self.score:.03f}\n"
        align_str += str_template.format(lname, lstr)
        align_str += str_template.format("", mstr)
        align_str += str_template.format(sname, sstr)
        align_str += "\n"

        return align_str
    
    def print_alignment(self):
        """Print an ASCII view of the aligned primers"""
        print(self.get_alignment_string())
    
    def get_primer_alignment(self):
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
            raise ValueError(
                "No alignment available. Call align() first."
            )
            
        return PrimerAlignment(
            primer1=self.primer1,
            primer2=self.primer2,
            primer1_name=self.primer1_name,
            primer2_name=self.primer2_name,
            score=self.score,
            alignment=self.get_alignment_string(),
        )


def create_nn_score_dt(match_json, single_mismatch_json, double_mismatch_score=0.2):
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
    with open(match_json, "r") as f:
        match_dt = json.load(f)
    with open(single_mismatch_json, "r") as f:
        single_mismatch_dt = json.load(f)

    # Set all as double mismatches; then update
    nts = ["A", "T", "C", "G"]
    nn_score_dt = {
        "".join(watson) + "/" + "".join(crick): double_mismatch_score
        for watson in product(nts, repeat=2)
        for crick in product(nts, repeat=2)
    }

    # Update
    nn_score_dt.update(match_dt)
    nn_score_dt.update(single_mismatch_dt)

    return nn_score_dt