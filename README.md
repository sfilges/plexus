# multiplexDesigner

## Install

Get the source code from GitHub:

```bash
git clone https://github.com/sfilges/multiplexDesigner
```

Then, install the software dependencies using conda:

```bash
cd multiplexDesigner
conda update conda
conda env create -f environment.yml
conda activate multiplexdesigner-run
```

Finally, install `MultiplexDesigner` itself with pip:

```bash
pip install -e .
```

Test installation by running:

```bash
multiplexdesigner
```

## Modules

### Preprocess

Import a file with user supplied mutation positions (junctions).

´´´
Name,Chrom,Five_Prime_Coordinate,Three_Prime_Coordinate
EGFR_T790M,chr7,55181378,55181378
KRAS_G12D,chr12,25245350,25245350
KRAS_G13R,chr12,25245348,25245349
BRAF_V600E,chr7,140753336,140753336
´´´

Sort and merge junctions that are close together. For example KRAS_G12D and KRAS_G13R
are merged such that the merged junction spans both junctions. It may be necessary to merge
multiple junctions.

´´´
Name,Chrom,Five_Prime_Coordinate,Three_Prime_Coordinate
EGFR_T790M,chr7,55181378,55181378
KRAS_G12D_KRAS_G13R,chr12,25245348,25245350
BRAF_V600E,chr7,140753336,140753336
´´´

Now, for each junction extract the genomic coordinates of the junction plus 200 bp padding
on either side. For example, the junction EGFR_T790M,chr7,55181378,55181378 will be
chr7:55181178-55181578.

Now search the provided reference genome (fasta file) and extract the genomic sequence for each junction. Append the
sequence to the junction table as the Design_Region, for example:

```csv
Name,Chrom,Five_Prime_Coordinate,Three_Prime_Coordinate,Design_Region
EGFR_T790M,chr7,55181378,55181378,CCAAACTCAGAGATCAGGTGACTCCGACTCCTCCTTTATCCAATGTGCTCCTCATGGCCACTGTTGCCTGGGCCTCTCTGTCATGGGGAATCCCCAGATGCACCCAGGAGGGGCCCTCTCCCACTGCATCTGTCACTTCACAGCCCTGCGTAAACGTCCCTGTGCTAGGTCTTTTGCAGGCACAGCTTTTCCTCCATGAGTACGTATTTTGAAACTCAAGATCGCATTCATGCGTCTTCACCTGGAAGGGGTCCATGTGCCCCTCCTTCTGGCCACCATGCGAAGCCACACTGACGTGCCTCTCCCTCCCTCCAGGAAGCCTACGTGATGGCCAGCGTGGACAACCCCCACGTGTGCCGCCTGCTGGGCATCTGCCTCACCTCCACCGTGCAGCTCATCACGCAGCTCATGCCCTTCGGCTGCCTCCTGGACTATGTCCGGGAACACAAAGACAATATTGGCTCCCAGTACCTGCTCAACTGGTGTGTGCAGATCGCAAAGGTAATCAGGGAAGGGAGATACGGGGAGGGGAGATAAGGAGCCAGGATCCTCACATGCGGTCTGCGCTCCTGGGATAGCAAGAGTTTGCCATGGGGATATGTGTGTGCGTGCATGCAGCACACACACATTCCTTTATTTTGGATTCAATCAAGTTGATCTTCTTGTGCACAAATCAGTGCCTGTCCCATCTGCATGTGGAAACTCTCATCAATCAGCTACCTTTGAAGAATTTTCTCTTTATTGAGTGCTCAGTGTGGTCTGATGTCTCTGTTCTTATTTCTCTGGAATTCTTTGTGAATA
```

Now we need to calculate the sequence target (junction) within the coordinates of the Design_Region. For example for the Desing_Regio, chr12:25245148-25245550 and junction KRAS_G12D_KRAS_G13R,chr12,25245348,25245350 the position of the junction is calculated as follows (including two bp padding on either side)

junction_length=25245550-25245148
jmin_coordinate=25245550-25245350-2 # End of region minus larger junction value minus 2 padding bases
jmax_coordinate=25245550-25245348+2 # End of region minus smaller junction value plus 2 padding bases

Verify that this makes sense.

Add junction_length,jmin_coordinate,jmax_coordinate as columns values to the table.

I think it might make sense to define an object class (Multplex Panel) that con tains the junction table, as well as metadata:

panel_name  = "test_panel"
genome      = "hg38"
date
panel_uuid

Furthermore we need to import design parameters from a config file, such as:

{
    "design" : {
        "primer_length_range": [18, 30],
        "amplicon_gap_range": [20, 60],
        "max_amplicon_length": 100,
        "junction_padding_bases": 3,
        "initial_solutions": 100,
        "top_solutions_to_keep": 4,
        "minimum_plexity": 5,
        "maximum_plexity": 20,
        "target_plexity": 20,
        "force_plexity": false,
        "variant_threshold": 0.01,
        "5_primer_tail": "",
        "3_prime_tail": "",
        "snp_penalty": 1.0,
        "primer_length_penalty": 1.0,
        "primer_complexity_penalty": 1.0,
        "polyA_penalty": 5,
        "polyT_penalty": 5,
        "polyC_penalty": 10,
        "polyG_penalty": 10,
        "amplicon_length_penalty": 1.0
    },
    "pcr" : {
        "annealing_temperature": 60,
        "primer_concentration": 50,
        "dntp_concentration": 0.6,
        "dna_concentration": 50,
        "mv_concentration": 50,
        "dv_concentration": 1.5,
        "dmso_concentration": 0.0,
        "dmso_fact": 0.6,
        "formamide_concentration": 0.8
    }
}

These can be added to the object as well. does an objec toriented approach make sense? Help me write the functions,

### Design

Basic primer design uses the primer3 libraries to generate individual
assays.


### Optimize


### Report


## Installation



## Usage




# Dependencies

- [primer3](https://github.com/primer3-org/primer3) through [primer3-py](https://github.com/libnano/primer3-py). Both licensed under GPL-2.0
- Blast from BioPython
https://primer3.org/manual.html#PRIMER_PRODUCT_SIZE_RANGE

https://github.com/JasonAHendry/multiply/tree/master
https://www.mfeprimer.com/

https://github.com/SemiQuant/PrimerJinn

https://github.com/matdoering/openPrimeR

https://www.sciencedirect.com/science/article/pii/S0022175920300193?via%3Dihub

https://github.com/treangenlab/Olivar


# run blast locally or through docker?

https://hub.docker.com/r/ncbi/blast

https://github.com/tamminenlab/npysearch



## Primer3 scoring

### Classic

Is based on local or global alignment scores to estimate complementarity. See [the primer3 docs](https://www.primer3plus.com/primer3plusHelp.html#PRIMER_MAX_SELF_ANY).

The scoring system gives 1.00 for complementary bases, -0.25 for a match of any base (or N) with an N, -1.00 for a mismatch, and -2.00 for a gap. Only single-base-pair gaps are allowed. For example, the alignment

   5' ATCGNA 3'
      || | |
   3' TA-CGT 5'

is allowed (and yields a score of 1.75), but the alignment

   5' ATCCGNA 3'
      ||  | |
   3' TA--CGT 5'

is not considered. Scores are non-negative, and a score of 0.00 indicates that there is no reasonable local alignment between two oligos.

### Thermodynamic

The melting temperature of the most stable structure is calculated. To calculate secondary structures nearest-neighbor parameters for perfect matches, single internal mismatches, terminal mismatches, dangling ends have been used. Also parameters for increments for length dependence of bulge and internal loops have been used. This parameter is calculated only if PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=1. The default value is 10 degrees lower than the default value of PRIMER_MIN_TM. For example, the alignment width length 15nt

  5' ATTAGATAGAGCATC 3'
  3' TAATCTATCTCGTAG 5'

is allowed (and yields a melting temperature of 32.1493 width by default Primer3 parameters), but the alignment

     T        C
  5'  GCGGCCGC GCGC 3'
  3'  CGCCGGCG CGCG 5'
     A        A

is not considered (Tm=57.0997 and the length of oligo is 14nt).