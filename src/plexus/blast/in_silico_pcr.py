# Run in silico PCR to evaluate primers
# https://github.com/pommevilla/ispcr

from ispcr import get_pcr_products


def ispcr_get_producr(primers, sequences):
    """
    get_pcr_products will then iterate through the sequences in sequence_file and find all products amplified by the forward and reverse primer.
    """

    # primer_file: This is currently limited to a fasta file containing two sequences, with the forward primer coming first and the reverse primer coming second
    # sequence_file: The path to the fasta file containing the sequences to test your primers against
    results = get_pcr_products(
        primer_file=primers,
        sequence_file=sequences,
        min_product_length=60,
        max_product_length=200,
    )

    return results
