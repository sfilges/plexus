from dataclasses import dataclass, field
from multiplexdesigner.designer.primer3_port import seqtm
from multiplexdesigner.utils.utils import gc_content, create_primer_dataframe

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
    directions = ["LEFT", "RIGHT"]

    for row in primer3_output:
        for d in directions:
            pair = {}
            pair[d] = Primer(
                seq = row['LEFT'],
                direction="F" if d == "LEFT" else "R",
                start=1,
                length=2,
                tm=1,
                gc=1
            )
        primer_pair = PrimerPair(
            F=pair["LEFT"],
            R=pair["RIGHT"]
        )
        primer_pairs.append(primer_pair)

    return primer_pairs


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
    tm_primer3: float
    tm: float
    bound: float
    gc: float
    
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
    """
    forward: Primer
    reverse: Primer
    product_length: int
    pair_id: str = field(default="", repr=False)


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
        left_seq = row['left_sequence']
        left_direction = 'forward'
        left_key = (left_seq, left_direction)
        if left_key not in primer_dict:
            left_tm_bound = seqtm(left_seq)
            left_primer = Primer(
                name = f"{junction.name}_{index}_forward",
                seq = left_seq,
                direction = left_direction,
                start = row['left_coords'][0],
                length = row['left_coords'][1],
                tm_primer3 = round(row['left_tm'], 2),
                tm = round(left_tm_bound.Tm, 2),
                bound = round(left_tm_bound.bound, 2),
                gc = round(gc_content(left_seq), 2)
            )
            primer_dict[left_key] = left_primer

        # REVERSE
        right_seq = row['right_sequence']
        right_direction = 'reverse'
        right_key = (right_seq, right_direction)
        if right_key not in primer_dict:
            right_tm_bound = seqtm(right_seq)
            right_primer = Primer(
                name=f"{junction.name}_{index}_reverse",
                seq = right_seq,
                direction = right_direction,
                start = row['right_coords'][0],
                length = row['right_coords'][1],
                tm_primer3 = round(row['right_tm'], 2),
                tm = round(right_tm_bound.Tm, 2),
                bound = round(right_tm_bound.bound, 2),
                gc = round(gc_content(right_seq), 2)
            )
            primer_dict[right_key] = right_primer

    return primer_dict


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

