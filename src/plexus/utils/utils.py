import warnings

from loguru import logger


def write_fasta_from_dict(input_dt: dict, output_fasta: str) -> None:
    """
    Write a `.fasta` file to `output_fasta` from an input dictionary
    `input_dt`

    """
    with open(output_fasta, "w") as fasta:
        for header, seq in input_dt.items():
            fasta.write(f">{header}\n")
            fasta.write(f"{seq}\n")


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


def reverse_complement(dna: str) -> str:
    """
    Returns the reverse complement of a DNA sequence.

    Args:
        dna (str): DNA sequence string (A, T, G, C)

    Returns:
        str: Reverse complement of the input DNA sequence
    """
    complement = {"A": "T", "T": "A", "G": "C", "C": "G"}

    # Handle lowercase and raise error for invalid bases
    dna = dna.upper()

    try:
        reverse_comp = "".join(complement[base] for base in reversed(dna))
        return reverse_comp
    except KeyError as e:
        raise ValueError(f"Invalid DNA base: {e.args[0]}") from e


# ===============================================
# Main function to generate candidate primers
# ===============================================


def generate_kmers(
    target_name: str,
    target_sequence: str,
    orientation: str,
    position_offset: int = 0,
    k_min: int = 18,
    k_max: int = 25,
    max_poly_X: int = 4,
    max_N: int = 0,
    min_gc: int = 30,
    max_gc: int = 70,
) -> list:
    """
    Generate k-mers as candidate primers.

    Args:
        target_name - Name of the target region, will be used for naming candidate primers
        target_sequence - Region from which to extract kmers
        orientation - Either forward or reverse
        position_offset: By how much is target_sequence shifted from the original template sequence?
        k_min - Minimum length of kmers
        k_max - Maximum length of kmers
        max_poly_X - Max number of times a base can occur in a row, e.g. AAAAA
        max_N - Max times N bases can occur anywhere in the kmer
        min_gc - Minimum GC content of kmer
        max_gc - Maximum GC content of kmer
    """
    from plexus.designer.primer import Primer

    if k_min >= k_max:
        error_msg = (
            f"Min kmer length ({k_min}) must be smaller than max kmer length ({k_max})"
        )
        logger.error(error_msg)
        raise ValueError(error_msg)

    if k_min < 10 or k_min > 20:
        warn_msg = (
            f"Provided value for kmin, {k_min}, is outside the expected range: 10 - 20"
        )
        logger.warning(warn_msg)
        warnings.warn(warn_msg, stacklevel=2)

    if k_max < 20 or k_max > 30:
        warn_msg = (
            f"Provided value for kmax, {k_max}, is outside the expected range: 20 - 30"
        )
        logger.warning(warn_msg)
        warnings.warn(warn_msg, stacklevel=2)

    kmers = []
    kmer_counter = 0
    kmer_filtered_counter = 0
    for k in range(k_min, k_max):
        # max position to search
        max_pos = len(target_sequence) + 1 - k

        for x in range(max_pos):
            kmer_counter += 1
            kmer = target_sequence[x : x + k]

            if check_kmer(
                kmer, max_poly_X=max_poly_X, max_N=max_N, min_gc=min_gc, max_gc=max_gc
            ):
                kmer_filtered_counter += 1

                # Calculate the correct start position
                if orientation == "reverse" or orientation == "right":
                    # For reverse primers, convert position to original forward coordinates
                    # The end of the kmer in the RC sequence corresponds to the start in the forward
                    start_pos = position_offset + len(target_sequence) - (x + k)
                else:
                    start_pos = position_offset + x

                # If pass: Instantiate kmer as Primer object
                good_kmer = Primer(
                    name=f"{target_name}_{kmer_filtered_counter}_{orientation}",
                    seq=kmer,
                    direction=orientation,
                    start=start_pos,
                    length=k,
                )

                kmers.append(good_kmer)

    logger.info(
        f"{orientation}: Found {kmer_filtered_counter} good kmers out of {kmer_counter} kmers checked."
    )

    return kmers


def check_kmer(
    kmer, max_poly_X: int = 4, max_N: int = 0, min_gc: int = 30, max_gc: int = 70
) -> bool:
    """
    Filter kmers (putative primers) based on GC-content, number of 'N' bases,
    and number of repeated bases (polyX).
    """
    kmer_gc = gc_content(kmer)
    return (
        min_gc <= kmer_gc <= max_gc  # Should be True
        and not check_N_in_kmers(kmer, max_N)  # Should be False
        and not find_max_poly_X(kmer, max_poly_X)  # Should be False
    )


def filter_kmers(kmers, max_poly_X=4, max_N=0):
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

    return filtered_kmers


def check_N_in_kmers(kmer, n=0):
    """
    Check if a k-mer contains more than n 'N' bases.

    Args:
        kmer (str): The k-mer sequence to check
        n (int): The maximum allowed number of N bases (default is 0)

    Returns:
        bool: True if the k-mer contains more than n 'N' bases, False otherwise
    """
    return kmer.upper().count("N") > n


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
