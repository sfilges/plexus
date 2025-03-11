#-------------------------------------------// (2) DESIGNER //----------------------------------
#
# DESIGNER is the core engine for creating all the primer and probe candidates and computing all the
# thermodynamic aspects of design such as target unimolecular folding, primer folding, bimolecular
# hybridization, and solving the multi-state coupled equilibria for the amount bound for the desired bimolecular
# duplex. The amount bound is directly related to the amount of signal generated in a diagnostic assay.
#
# DESIGNER also analyzes each primer and probe design for a series of heuristic properties such as
# sequence complexity, polyG test, oligo length penalty, amplicon length penalty, etc. Each of the scoring effects
# are multiplied by weighting factors and combined into an overall score for each primer/ probe set.
#
# In December 2018 functionality was added to PanelPlex is support for dbSNP (for human genome
# version GRCh38 and mouse genome version mm10). PanelPlex automatically detects all the positions with
# alternative alleles with CAF > 0.01 within each of the design regions and then automatically designs the
# primers to AVOID those SNP sites.
#
# 1. Generate potential solutions TODO Generate k-mers as initial solutions
# 2. Score each solutions TODO define scoring algorithm, including SNP penalty
# 3. Select top candidates for each target (init_solutions)
# 4. BLAST against amplicons in panel and discard matches TODO implement BLAST search
# 5. Keep top N solutions (top_solutions_to_keep)
# 6. Perform multiplex picking in N^X space for N solutions and X targets TODO implement multiplex picker algorithm
#
#-----------------------------------------------------------------------------------------------

def generate_kmers(k_min, k_max, sequence):
    """
    Generate k-mers as candidate primers
    """
    kmers = []
    for k in range (k_min, k_max):
        # max position to search
        max_pos = len(sequence) + 1 - k
        
        for x in range(max_pos):
            kmer = sequence[x:x+k]
            kmers.append(kmer)

    return(kmers)

def filter_kmers(kmers, max_kmers = 100):
    print("TODO")

def filter_dbsnp(database):
    print("TODO")

def blast_search():
    print("TODO")

def main():
    region = "CACAGGGTAGAGACAGATAACAAGGGATTCCTGACACCAAAAAAAAAAATACGCTGTAGATAGCTATAACATTTCAATAGGAATCTTGGGAATC"
    k_min = min_primer_length
    k_max = max_primer_length
    
    #-----// Generate k-mer solutions for targets
    kmers = generate_kmers(k_min, k_max, region)
    print(kmers)
    print(len(kmers))

    # Test kmers against baisc criteria and rank, Retain at most max_kmers
    ## Needs to handle too few or no viable kmers
    kmers = filter_kmers(kmers, max_kmers = 100)


if __name__=="__main__":
    main()