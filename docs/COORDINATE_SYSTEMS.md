# Coordinate Systems in Plexus

For a tool where off-by-one errors are safety-critical, understanding the coordinate systems used for inputs, internal processing, and outputs is essential.

## Summary

| Context | Coordinate System | Base | Inclusive |
|---------|-------------------|------|-----------|
| **Input CSV** | Genomic | 1 | Yes |
| **Output CSV** | Genomic | 1 | Yes |
| **Internal Offsets** | Design-Region Relative | 0 | Yes |
| **Reference FASTA** | Standard (pysam handles) | 1 (text) | Yes |
| **VCF/gnomAD** | Genomic | 1 | Yes |

---

## Genomic Coordinates (1-based)

Plexus follows the standard convention used by genome browsers (IGV, UCSC), VCF files, and GFF/GTF annotations.

### Input: `junctions.csv`
The columns `Five_Prime_Coordinate` and `Three_Prime_Coordinate` must be provided as **1-based inclusive** genomic coordinates.
- For a single-base target (e.g., an SNV), both coordinates should be the same.
- For a range (e.g., a deletion or a whole exon), provide the first and last base of the target.

### Output: `selected_multiplex.csv` / `candidate_pairs.csv`
All genomic coordinate columns in the output files are **1-based inclusive**:
- `Forward_Genomic_Start` / `Forward_Genomic_End`
- `Reverse_Genomic_Start` / `Reverse_Genomic_End`
- `Junction_Start` / `Junction_End`

This ensures that you can directly paste these coordinates into IGV or a genome browser to verify the primer positions.

---

## Design-Region Offsets (0-based)

Internally, Plexus extracts a padded genomic sequence around your target, known as the **Design Region**.

### `design_start`
The attribute `junction.design_start` stores the **1-based genomic coordinate** of the first base in the extracted sequence.

### Primer Offsets
The `start` attribute of a `Primer` object is a **0-based index** into the `design_region` string.
- A primer starting at the very beginning of the design region has a `start` of `0`.
- The genomic start of a primer is calculated as: `design_start + primer.start`.
- The genomic end of a primer is calculated as: `design_start + primer.start + primer.length - 1`.

---

## External Tools and Libraries

Plexus handles the conversion between its internal 1-based system and the 0-based systems used by some underlying libraries:

### pysam (VCF/FASTA access)
`pysam` uses 0-based half-open coordinates (like the BAM/BED formats).
- **Extraction**: When Plexus fetches a sequence from the FASTA, it converts the 1-based `[start, end]` range to 0-based `[start-1, end)`.
- **SNP Checking**: When querying a VCF, Plexus converts the 1-based primer range `[P_start, P_end]` to the 0-based query `fetch(chrom, P_start - 1, P_end)`.

### Primer3
`primer3-py` uses 0-based coordinates for sequence offsets. Plexus passes the `design_region` string and receives 0-based offsets, which are stored directly in the `Primer` objects.

### BLAST
The BLAST specificity check uses 1-based coordinates for its output, which Plexus parses and maps back to the original genomic coordinates for on-target/off-target classification.

---

## Developer Note: Coordinate Math

When performing arithmetic on coordinates, always verify the base:
- **Genomic to Design-Region**: `junction.start - junction.design_start` (both 1-based) yields a 0-based offset.
- **Design-Region to Genomic**: `junction.design_start + primer.offset` (1-based + 0-based) yields a 1-based genomic coordinate.
