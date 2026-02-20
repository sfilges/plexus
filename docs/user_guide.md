# Plexus Multiplex Primer Designer User Guide

## Table of Contents

1. [Getting Started](#getting-started)
2. [Operational Modes: Research vs Compliance](#operational-modes-research-vs-compliance)
3. [Core Concepts](#core-concepts)
4. [Step-by-Step Workflow](#step-by-step-workflow)
5. [Configuration Guide](#configuration-guide)
6. [Advanced Features](#advanced-features)
7. [Output Interpretation](#output-interpretation)
8. [Troubleshooting](#troubleshooting)
9. [Technical Reference](#technical-reference)
10. [Compliance and Clinical Use](#compliance-and-clinical-use)

## Getting Started

### Operational Modes: Research vs Compliance

Plexus operates in two distinct modes with different validation requirements:

**Research Mode (Default):**

- Designed for exploratory and development work
- Checksum verification is optional
- Flexible resource management
- Faster setup and iteration
- Use `--strict` flag for optional checksum verification

**Compliance Mode:**

- Designed for clinical and regulated environments
- **Mandatory checksum verification** of all resources
- Strict resource validation before each run
- Audit trail and provenance tracking
- Enforced via `--mode compliance` during initialization
- Cannot be bypassed - will error on missing checksums

**Key Differences:**

| Feature | Research Mode | Compliance Mode |
| ------- | ------------- | --------------- |
| Checksum Verification | Optional (`--strict` flag) | **Mandatory** |
| Resource Validation | Warning only | **Error on mismatch** |
| Provenance Tracking | Basic | **Comprehensive** |
| Audit Requirements | None | **Full documentation** |
| Performance Impact | None | Minimal (checksum verification) |
| Use Case | Development, research | Clinical, diagnostic, regulated |

### Setting Operational Mode

Set the mode during initialization:

```bash
# For clinical/compliance use (MANDATORY checksums)
plexus init --genome hg38 --mode compliance --checksums my_checksums.sha256

# For research use (default - checksums optional)
plexus init --genome hg38 --mode research

# Check current mode
plexus status
```

**Important Compliance Notes:**

- Once set to compliance mode, **all runs will enforce checksum verification**
- Compliance mode cannot be bypassed without explicitly changing the mode
- **All resource files must have verified checksums before use**
- Checksums are stored in `~/.plexus/data/registry.json`
- Provenance information is recorded in each output directory

### Installation

Plexus requires Python 3.10-3.13 and several bioinformatics tools:

**Prerequisites:**

- Python 3.10-3.13
- NCBI BLAST+ suite (`blastn`, `makeblastdb`)
- `bcftools` for VCF processing

**Recommended Installation (using uv):**

```bash
# Clone the repository
git clone https://github.com/sfilges/plexus
cd plexus

# Install with uv
uv pip install -e .
```

**Alternative Installation (using conda):**

```bash
# Create conda environment
conda env create -f config/environment.yml
conda activate plexus-run

# Install plexus
pip install -e .
```

**Installing dependencies on Debian/Ubuntu:**

```bash
apt install ncbi-blast+ bcftools
```

### Quick Start Example

```bash
# 1. Generate starter templates (junctions.csv, designer_config.json)
plexus template --output my_panel/
cd my_panel/

# 2. Initialize resources (downloads hg38 genome and gnomAD VCF)
#    This also verifies chromosome naming consistency.
plexus init --genome hg38 --download

# 3. Run primer design
plexus run -i junctions.csv -g hg38 -o results/
```

### System Requirements

- **Memory**: 8GB minimum, 16GB+ recommended for large panels (>50 targets)
- **Disk Space**: 20GB+ for genome resources, plus output directory space
- **CPU**: Multi-core recommended for parallel processing

## Core Concepts

### Junction

A **Junction** represents a single genomic target for primer design:

```python
Junction(
    name="EGFR_T790M",
    chrom="chr7",
    start=55181378,    # 5' coordinate
    end=55181378,      # 3' coordinate
    design_region="ATCG...",  # Extracted genomic sequence
    primer_pairs=[...]  # Designed primer pairs
)
```

**Key Attributes:**
- `name`: User-provided identifier
- `chrom`: Chromosome in genome reference format
- `start`/`end`: Genomic coordinates (1-based)
- `design_region`: Padded sequence around the target
- `primer_pairs`: List of valid PrimerPair objects

### Primer

A **Primer** represents a single oligonucleotide with thermodynamic properties:

```python
Primer(
    name="EGFR_T790M_F1",
    sequence="ATCGATCGATCG",
    tm=60.2,              # Melting temperature (°C)
    gc_content=45.5,     # GC percentage
    bound_fraction=98.7, # Fraction bound at annealing temp
    hairpin_dg=3.2,       # Hairpin ΔG (kcal/mol)
    self_dimer_dg=2.1,   # Self-dimer ΔG (kcal/mol)
    penalty=1.5           # Quality penalty score
)
```

### PrimerPair

A **PrimerPair** combines forward and reverse primers for one target:

```python
PrimerPair(
    forward=Primer(...),
    reverse=Primer(...),
    amplicon_size=120,
    tm_difference=0.8,
    pair_penalty=2.3,
    cross_dimer_dg=4.5,
    off_targets=[]
)
```

**Quality Metrics:**
- `pair_penalty`: Combined quality score (lower is better)
- `cross_dimer_dg`: Cross-dimer potential with other primers
- `off_targets`: BLAST-detected off-target sites

### MultiplexPanel

The **MultiplexPanel** is the central object containing all targets and optimization logic:

```python
MultiplexPanel(
    panel_name="cancer_hotspots",
    genome="hg38",
    junctions=[Junction(...), ...],
    config=DesignerConfig(...),
    selected_pairs=[PrimerPair(...), ...]
)
```

**Key Methods:**
- `build_selector_dataframe()`: Prepares data for optimization
- `save_candidate_pairs_to_csv()`: Exports all designed pairs
- `save_selected_multiplex_csv()`: Exports optimized selection

## Step-by-Step Workflow

### 1. Project Scaffolding

Start by generating a project template in a new directory:

```bash
plexus template --output my_new_design
cd my_new_design
```

This creates:
- `junctions.csv`: A template file for your target coordinates.
- `designer_config.json`: The default design parameters for you to customize.

### 2. Input Preparation

Edit the `junctions.csv` with your targets. The required columns are:
- **Name**: A unique identifier for the target (e.g. `BRAF_V600E`).
- **Chrom**: The chromosome name (must exactly match your FASTA, e.g. `chr7` or `7`).
- **Five_Prime_Coordinate**: The genomic position (1-based) of the target start.
- **Three_Prime_Coordinate**: The genomic position (1-based) of the target end.
- **Panel** (optional): Used for multi-panel runs.

### 3. Resource Setup

```bash
# Initialize genome resources
plexus init --genome hg38 --download

# For compliance mode (MANDATORY checksums)
plexus init --genome hg38 --mode compliance --checksums my_checksums.sha256

# Check resource status
plexus status
```

**Resource Types:**

- **FASTA**: Reference genome sequence
- **FAI**: FASTA index file
- **BLAST DB**: BLAST database for specificity checking
- **gnomAD VCF**: SNP database for overlap checking

**Naming Validation:**
During `init`, Plexus automatically compares the chromosome naming convention (e.g., `chr1` vs `1`) between your FASTA and SNP VCF. If a systematic mismatch is detected:
- In **Research Mode**, a warning is logged.
- In **Compliance Mode**, initialization fails immediately.

### 4. Running the Design Pipeline

```bash
plexus run \
  --input junctions.csv \
  --genome hg38 \
  --output results/ \
  --name my_panel \
  --preset default
```

**Pipeline Steps:**

1. **Panel Creation**: Load junctions and extract genomic regions. **Note**: Chromosome naming is re-verified here against your input CSV.
2. **Primer Design**: Generate and filter primer candidates
3. **SNP Checking**: Detect overlaps with common variants
4. **BLAST Specificity**: Check for off-target amplification
5. **Multiplex Optimization**: Select optimal primer combination
6. **Output Generation**: Save results to multiple files

### 5. Reviewing Results

Key output files in the output directory:

```text
results/
├── candidate_pairs.csv      # All designed primer pairs
├── selected_multiplex.csv   # Optimized final selection
├── top_panels.csv           # Alternative solutions
├── panel_summary.json       # Metadata and provenance
├── off_targets.csv          # BLAST specificity results
├── failed_junctions.csv    # Targets that failed design
└── provenance.json          # Tool versions and parameters
```

## Configuration Guide

### Configuration Presets

**Default Preset** (`--preset default`):

- Strict thermodynamic constraints
- Balanced Tm and fraction bound optimization
- Conservative SNP penalties
- Suitable for most applications

**Lenient Preset** (`--preset lenient`):
- Wider parameter ranges
- Higher tolerance for suboptimal primers
- Useful for difficult genomic regions
- May yield more primer pairs but lower quality

### Key Parameters and Their Impact

**Singleplex Design Parameters:**

| Parameter | Impact | Recommended Range |
| ----------- | -------- | ------------------- |
| `PRIMER_OPT_TM` | Target melting temperature | 55-65°C |
| `PRIMER_MIN_TM`/`PRIMER_MAX_TM` | Tm acceptance range | ±3-5°C from optimal |
| `PRIMER_OPT_BOUND` | Fraction bound at annealing temp | 95-99% |
| `primer_min_length`/`primer_max_length` | Primer length range | 18-25 bp typical |
| `primer_min_gc`/`primer_max_gc` | GC content constraints | 30-70% |

**Primer Pair Parameters:**

| Parameter | Impact | Recommended Range |
| ----------- | -------- | ------------------- |
| `PRIMER_PAIR_MAX_DIFF_TM` | Max Tm difference between pair | 1-5°C |
| `PRIMER_PRODUCT_OPT_SIZE` | Target amplicon size | 60-150 bp |
| `PRIMER_PRODUCT_MIN_SIZE`/`MAX_SIZE` | Amplicon size range | 50-300 bp |

**Multiplex Picker Parameters:**

| Parameter | Impact | Recommended Range |
| ----------- | -------- | ------------------- |
| `wt_cross_dimer` | Weight for cross-dimer penalties | 1-5 |
| `wt_off_target` | Weight for off-target penalties | 3-10 |
| `wt_snp_penalty` | Weight for SNP overlap penalties | 1-5 |
| `initial_solutions` | Number of solutions to generate | 100-1000 |
| `top_solutions_to_keep` | Best solutions to retain | 3-10 |

### Creating Custom Configurations

Create a JSON file with only the parameters you want to override:

```json
{
    "singleplex_design_parameters": {
        "PRIMER_OPT_TM": 62.0,
        "PRIMER_MIN_TM": 58.0,
        "PRIMER_MAX_TM": 66.0,
        "primer_min_length": 20,
        "primer_max_length": 28
    },
    "multiplex_picker_parameters": {
        "wt_cross_dimer": 2.5,
        "wt_off_target": 8.0
    }
}
```

Then use it with:

```bash
plexus run --config my_config.json ...
```

### Thermodynamic vs Tm-based Design

**Thermodynamic Design (Recommended for Multiplex):**

```json
{
    "PRIMER_WT_TM_GT": 0.0,
    "PRIMER_WT_TM_LT": 0.0,
    "PRIMER_WT_BOUND_GT": 1.0,
    "PRIMER_WT_BOUND_LT": 1.0,
    "PRIMER_OPT_BOUND": 97.0,
    "PRIMER_MIN_BOUND": 90.0,
    "PRIMER_MAX_BOUND": 110.0
}
```

**Tm-based Design (Classic Approach):**
```json
{
    "PRIMER_WT_TM_GT": 1.0,
    "PRIMER_WT_TM_LT": 1.0,
    "PRIMER_WT_BOUND_GT": 0.0,
    "PRIMER_WT_BOUND_LT": 0.0,
    "PRIMER_OPT_TM": 60.0,
    "PRIMER_MIN_TM": 57.0,
    "PRIMER_MAX_TM": 63.0
}
```

**When to Use Each:**

- **Thermodynamic**: Better for multiplex PCR where uniform binding is critical
- **Tm-based**: Simpler approach, good for singleplex or when compatibility with existing protocols is needed
- **Hybrid**: Can use intermediate weights for balanced approach

## Advanced Features

### Multi-Panel Design

Use the `Panel` column in your input CSV to create multiple independent panels:

```csv
Name,Chrom,Five_Prime_Coordinate,Three_Prime_Coordinate,Panel
EGFR_T790M,chr7,55181378,55181378,patient_A
KRAS_G12D,chr12,25245350,25245350,patient_A
TP53_R248W,chr17,7674220,7674220,patient_B
BRCA1_185delAG,chr17,41244950,41244950,patient_B
```

Run with parallel processing:

```bash
plexus run \
  --input multi_panel_junctions.csv \
  --parallel \
  --max-workers 4
```

**Output Structure:**

```text
results/
├── patient_A/
│   ├── selected_multiplex.csv
│   └── ...
└── patient_B/
    ├── selected_multiplex.csv
    └── ...
```

### Selector Algorithms

Plexus offers five multiplex optimization algorithms:

| Algorithm | Description | Best For |
| ----------- | ----------- | -------- |
| **Greedy** (default) | Iterative selection with random ordering | Most use cases, good balance |
| **Random** | Random selection from valid pairs | Quick testing, baseline comparison |
| **BruteForce** | Exhaustive search of all combinations | Small panels (<10 targets) |
| **SimulatedAnnealing** | Probabilistic optimization | Large panels, complex constraints |
| **DFS** | Depth-first search | Medium panels, thorough exploration |

**Algorithm Selection Guide:**

- **<10 targets**: `BruteForce` for optimal results
- **10-50 targets**: `Greedy` or `SimulatedAnnealing`
- **>50 targets**: `Greedy` with increased `initial_solutions`
- **Quick testing**: `Random` for baseline comparison

### Custom VCF Usage

Use your own SNP database instead of bundled gnomAD:

```bash
plexus run \
  --input junctions.csv \
  --snp-vcf /path/to/custom.vcf.gz \
  --snp-af-threshold 0.05 \
  --snp-strict
```

**VCF Requirements:**

- Must be tabix-indexed (`.vcf.gz` + `.vcf.gz.tbi`)
- Must contain AF (allele frequency) information
- Can be any genomic region or variant type

### Performance Optimization

**For Large Panels (>50 targets):**

```bash
# Use Greedy with more initial solutions
plexus run \
  --selector Greedy \
  --config large_panel_config.json
```

```json
{
    "multiplex_picker_parameters": {
        "initial_solutions": 1000,
        "top_solutions_to_keep": 5,
        "wt_cross_dimer": 1.5,
        "wt_off_target": 6.0
    }
}
```

**Memory Management:**
- Use `--skip-blast` if you've pre-validated specificity
- Use `--skip-snpcheck` if working with synthetic targets
- Process panels sequentially if memory is limited

## Output Interpretation

### candidate_pairs.csv

Contains all designed primer pairs before optimization:

| Column | Description |
|--------|-------------|
| `target_id` | Junction name |
| `pair_name` | Unique pair identifier |
| `Forward_Seq` | Forward primer sequence (binding region only) |
| `Reverse_Seq` | Reverse primer sequence (binding region only) |
| `Forward_Full_Seq` | Forward primer with adapter tail |
| `Reverse_Full_Seq` | Reverse primer with adapter tail |
| `Amplicon_Size` | Expected amplicon size (bp) |
| `Pair_Penalty` | Quality score (lower is better) |
| `Forward_Tm` | Forward primer melting temperature |
| `Reverse_Tm` | Reverse primer melting temperature |
| `Tm_Difference` | Tm difference between primers |
| `Cross_Dimer_dG` | Cross-dimer potential |
| `SNP_Penalty` | SNP overlap penalty |

### selected_multiplex.csv

The optimized final primer selection:

| Column | Description |
|--------|-------------|
| `Target` | Junction name |
| `Forward_Primer` | Selected forward primer name |
| `Reverse_Primer` | Selected reverse primer name |
| `Forward_Seq` | Forward sequence (binding region) |
| `Reverse_Seq` | Reverse sequence (binding region) |
| `Forward_Full_Seq` | Forward sequence with adapter |
| `Reverse_Full_Seq` | Reverse sequence with adapter |
| `Amplicon_Size` | Expected amplicon size |
| `Total_Cost` | Multiplex optimization cost |
| `Cross_Dimer_Risk` | Aggregate cross-dimer score |
| `Off_Target_Risk` | Aggregate off-target score |

### top_panels.csv

Alternative multiplex solutions ranked by cost:

| Column | Description |
|--------|-------------|
| `Solution_Rank` | Rank (1 = best) |
| `Total_Cost` | Optimization cost score |
| `Primer_Pairs` | List of selected pair names |
| `Cross_Dimer_Score` | Aggregate cross-dimer penalty |
| `Off_Target_Score` | Aggregate off-target penalty |
| `SNP_Score` | Aggregate SNP penalty |

### panel_summary.json

Comprehensive metadata and provenance:

```json
{
    "panel_name": "cancer_hotspots",
    "genome": "hg38",
    "num_junctions": 42,
    "num_selected_pairs": 42,
    "total_candidate_pairs": 876,
    "optimization_algorithm": "Greedy",
    "provenance": {
        "plexus_version": "0.4.6",
        "primer3_version": "2.6.1",
        "run_timestamp": "2023-11-15T14:30:00Z",
        "fasta_path": "/path/to/hg38.fa",
        "fasta_sha256": "a1b2c3...",
        "snp_vcf_path": "/path/to/gnomad.vcf.gz",
        "run_blast": true,
        "skip_snpcheck": false
    },
    "configuration": {...}
}
```

### off_targets.csv

BLAST specificity results:

| Column | Description |
|--------|-------------|
| `Primer_Name` | Primer identifier |
| `Target_Sequence` | Primer sequence |
| `Off_Target_Chrom` | Chromosome of off-target |
| `Off_Target_Start` | Start coordinate |
| `Off_Target_End` | End coordinate |
| `Alignment_Score` | BLAST alignment score |
| `E_Value` | Expect value |
| `Identity_Percent` | Sequence identity |

## Troubleshooting

### Common Error Messages

**"FASTA file not found"**
- Run `plexus init --genome hg38 --download` first
- Or provide `--fasta /path/to/genome.fa`

**"No primer pairs found for junction X"**
- Check genomic coordinates in input CSV
- Try lenient preset: `--preset lenient`
- Adjust padding: `--padding 300`
- Relax GC constraints in custom config

**"BLAST not available"**
- Install NCBI BLAST+: `conda install -c bioconda blast`
- Ensure `blastn` is in your PATH
- Or use `--skip-blast`

**"MemoryError during optimization"**
- Reduce panel size (split into multiple panels)
- Use `--selector Greedy` with lower `initial_solutions`
- Process panels sequentially instead of parallel

### Primer Design Failures

**Common Causes and Solutions:**

| Issue | Solution |
|-------|----------|
| High GC content region | Relax `primer_max_gc` constraint |
| Repetitive sequence | Increase `primer_max_poly_x` |
| No valid amplicon size | Adjust `PRIMER_PRODUCT_MIN_SIZE`/`MAX_SIZE` |
| Extreme Tm requirements | Widen `PRIMER_MIN_TM`/`MAX_TM` range |

**Example Lenient Config for Difficult Regions:**

```json
{
    "singleplex_design_parameters": {
        "primer_min_gc": 20,
        "primer_max_gc": 80,
        "primer_max_poly_x": 8,
        "PRIMER_MIN_TM": 50.0,
        "PRIMER_MAX_TM": 70.0,
        "PRIMER_MIN_BOUND": 80.0,
        "PRIMER_MAX_BOUND": 120.0
    }
}
```

### Performance Issues

**Slow BLAST Specificity Check:**
- Use `--skip-blast` for initial testing
- Create BLAST database once: `makeblastdb -in genome.fa -dbtype nucl`
- Use smaller genome subsets if possible

**Memory Constraints:**
- Process smaller batches (split input CSV)
- Use `--selector Greedy` instead of `BruteForce`
- Reduce `initial_solutions` in config
- Disable parallel processing for multi-panel runs

### Resource Verification

Check resource integrity:

```bash
plexus status
```

**Common Resource Issues:**

| Issue | Solution |
|-------|----------|
| Checksum mismatch | Re-download resources or use `--force` |
| Missing FAI index | Run `samtools faidx genome.fa` |
| Corrupt BLAST DB | Rebuild with `makeblastdb` |
| Missing tabix index | Run `tabix -p vcf snps.vcf.gz` |

### Compliance Mode Specific Issues

**"Checksum mismatch for FASTA"**
- Verify your checksums file format (sha256sum format)
- Recompute checksums: `sha256sum genome.fa > my_checksums.sha256`
- Re-register resources with correct checksums

**"No checksums stored for genome in compliance mode"**
- Must initialize with `--checksums` file in compliance mode
- Cannot use compliance mode without verified resources
- Switch to research mode or provide proper checksums

## Technical Reference

### CLI Command Reference

**Main Commands:**

```bash
# Run primer design pipeline
plexus run [OPTIONS]

# Initialize genome resources
plexus init [OPTIONS]

# Check system and resource status
plexus status

# Generate starter templates
plexus template [OPTIONS]
```

**Common Options:**

| Option | Default | Description |
|--------|---------|-------------|
| `-i, --input` | Required | Input CSV file path |
| `-g, --genome` | `hg38` | Reference genome name |
| `-o, --output` | `./output` | Output directory (for `run` and `template`) |
| `-n, --name` | `multiplex_panel` | Panel name |
| `-p, --preset` | `default` | Config preset |
| `-c, --config` | None | Custom config file |
| `-s, --selector` | `Greedy` | Optimization algorithm |
| `--skip-blast` | False | Skip BLAST specificity check |
| `--skip-snpcheck` | False | Skip SNP overlap check |
| `--snp-vcf` | Bundled | Custom VCF file path |
| `--snp-strict` | False | Discard SNP-overlapping primers |
| `--parallel` | False | Parallel multi-panel processing |

### Complete Configuration Parameters

**singleplex_design_parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `PRIMER_NUM_RETURN` | int | 10 | Number of primers to return per side |
| `PRIMER_OPT_TM` | float | 60.0 | Optimal melting temperature (°C) |
| `PRIMER_MIN_TM` | float | 57.0 | Minimum acceptable Tm (°C) |
| `PRIMER_MAX_TM` | float | 63.0 | Maximum acceptable Tm (°C) |
| `PRIMER_OPT_SIZE` | int | 22 | Optimal primer length (bp) |
| `primer_min_length` | int | 15 | Minimum primer length (bp) |
| `primer_max_length` | int | 30 | Maximum primer length (bp) |
| `PRIMER_OPT_BOUND` | float | 98.0 | Optimal fraction bound (%) |
| `PRIMER_MIN_BOUND` | float | -10.0 | Minimum fraction bound (%) |
| `PRIMER_MAX_BOUND` | float | 120.0 | Maximum fraction bound (%) |

**primer_pair_parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `PRIMER_PAIR_MAX_DIFF_TM` | float | 3.0 | Max Tm difference in pair (°C) |
| `PRIMER_PRODUCT_OPT_SIZE` | int | 60 | Optimal amplicon size (bp) |
| `PRIMER_PRODUCT_MIN_SIZE` | int | 20 | Min amplicon insert size (bp) |
| `PRIMER_PRODUCT_MAX_SIZE` | int | 120 | Max amplicon size (bp) |

**multiplex_picker_parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `initial_solutions` | int | 100 | Solutions to generate |
| `top_solutions_to_keep` | int | 4 | Best solutions to retain |
| `wt_cross_dimer` | float | 1.0 | Cross-dimer penalty weight |
| `wt_off_target` | float | 5.0 | Off-target penalty weight |
| `wt_snp_penalty` | float | 3.0 | SNP overlap penalty weight |

### Algorithm Details

**Greedy Search:**
- Iteratively selects best primer for each target
- Randomizes target order for multiple runs
- Fast and effective for most use cases
- Configurable via `initial_solutions` parameter

**Brute Force:**
- Exhaustively evaluates all possible combinations
- Guaranteed optimal solution for small panels
- Computationally expensive (O(N^K) where K = targets)
- Limited to ~10 targets due to combinatorial explosion

**Simulated Annealing:**
- Probabilistic optimization inspired by metallurgy
- Can escape local optima
- Good for complex constraint landscapes
- Requires tuning of temperature parameters

**File Format Specifications:**

**Input CSV Format:**
```csv
Name,Chrom,Five_Prime_Coordinate,Three_Prime_Coordinate,Panel
STRING,STRING,INTEGER,INTEGER,STRING(optional)
```

**Output CSV Formats:**
- All CSV files use comma delimiting
- First row contains header with column names
- String fields are unquoted unless containing commas
- Numeric fields use standard decimal notation
- Missing values represented as empty strings

**JSON Format:**
- UTF-8 encoded
- Pretty-printed with 2-space indentation
- ISO 8601 timestamps
- SHA-256 checksums in hexadecimal format

## Compliance and Clinical Use

### Checksum Requirements for Compliance Mode

Compliance mode requires SHA-256 checksums for all resource files. Checksums must be provided in `sha256sum` format:

```bash
# Create checksums file for your resources
sha256sum genome.fa > my_checksums.sha256
sha256sum snps.vcf.gz >> my_checksums.sha256

# Initialize with checksums
plexus init \
  --genome hg38 \
  --mode compliance \
  --fasta genome.fa \
  --snp-vcf snps.vcf.gz \
  --checksums my_checksums.sha256
```

**Checksums File Format:**
```
a1b2c3d4e5f6...  genome.fa
7g8h9i0j1k2l...  snps.vcf.gz
```

### Clinical Workflow Best Practices

1. **Resource Verification:**
   ```bash
   # Verify all resources before use
   plexus status

   # Check checksums explicitly
   sha256sum -c my_checksums.sha256
   ```

2. **Documented Provenance:**
   - All runs in compliance mode record comprehensive provenance
   - Provenance includes tool versions, resource checksums, timestamps
   - Output includes `provenance.json` with full audit trail

3. **Validation Requirements:**
   - **Mandatory**: SNP checking and BLAST specificity in compliance mode
   - **Recommended**: Additional wet-lab validation of selected primers
   - **Documented**: All parameters and configuration settings

4. **Reproducibility:**
   ```bash
   # Use exact configuration files
   plexus run --config validated_config.json ...

   # Version control all input files
   git add junctions.csv validated_config.json my_checksums.sha256
   ```

### Transitioning from Research to Compliance

```bash
# 1. Develop and test in research mode
plexus init --genome hg38 --mode research --download
plexus run --input test_junctions.csv ...

# 2. Generate checksums for validated resources
sha256sum ~/.plexus/data/genomes/hg38.fa > clinical_checksums.sha256
sha256sum ~/.plexus/data/snp/gnomad.vcf.gz >> clinical_checksums.sha256

# 3. Switch to compliance mode with verified resources
plexus init --genome hg38 --mode compliance --checksums clinical_checksums.sha256

# 4. Run in compliance mode (checksums enforced)
plexus run --input clinical_junctions.csv ...
```

### Compliance Mode Limitations

- **No Checksum Bypass**: Runs will fail if checksums don't match
- **Strict Validation**: All resources must pass verification
- **Performance Impact**: Additional verification steps
- **Resource Management**: More rigorous file handling requirements

### Dependency Version Handling

**Current Implementation:**
- **Version Recording**: Tool versions are captured in provenance (BLAST, bcftools, primer3-py)
- **No Version Enforcement**: Compliance mode does NOT enforce specific dependency versions
- **Version Tracking**: All versions are recorded in `provenance.json` for audit trails

**Tool Version Capture:**
```json
{
    "plexus_version": "0.4.6",
    "primer3_version": "2.6.1",
    "tool_versions": {
        "blastn": "2.12.0+",
        "bcftools": "1.15",
        "makeblastdb": "2.12.0+"
    }
}
```

**Recommendations for Clinical Use:**
1. **Document Validated Versions**: Record which tool versions were validated in your workflow
2. **Version Pinning**: Use containerization (Docker) or conda environments to ensure reproducibility
3. **Pre-run Validation**: Manually verify tool versions match your validated configuration
4. **Audit Trail**: Review provenance files to confirm versions used in each run

**Future Considerations:**
While Plexus currently records but doesn't enforce tool versions, clinical users should:
- Establish internal version requirements
- Implement version validation in their workflows
- Use the provenance data for quality assurance
- Consider containerization for complete environment control

### Audit and Documentation Requirements

For clinical/compliance use, maintain:

1. **Resource Documentation:**
   - Source of all reference files
   - Version information for genome builds
   - Date of resource acquisition
   - Checksum verification records

2. **Run Documentation:**
   - Input files (CSV, config, checksums)
   - Exact command line used
   - Plexus version and dependencies
   - Output directory with provenance files

3. **Validation Records:**
   - Wet-lab validation results
   - Any manual overrides or exceptions
   - Final primer sequences used

### Example Clinical Documentation Structure

```text
project/
├── inputs/
│   ├── clinical_junctions.csv
│   ├── validated_config.json
│   └── clinical_checksums.sha256
├── resources/
│   ├── hg38.fa (with verified checksum)
│   └── clinical_snps.vcf.gz (with verified checksum)
├── outputs/
│   ├── run_20240115/
│   │   ├── selected_multiplex.csv
│   │   ├── provenance.json
│   │   └── ...
│   └── run_20240122/
│       └── ...
└── documentation/
    ├── resource_verification.md
    ├── run_protocols.md
    └── validation_results.md
```

## Appendix: Primer Design Theory

### The Multiplex Primer Selection Problem

Given N primer pairs for each of k targets, we have 2×k×N available single primers. The goal is to pick exactly one primer pair flanking each of the k targets while minimizing:

1. **Cross-dimer potential**: Interactions between primers
2. **Off-target amplification**: Non-specific binding
3. **SNP overlaps**: Primers spanning common variants
4. **Thermodynamic mismatches**: Suboptimal binding conditions

### Optimization Cost Function

The total cost C for a multiplex M is calculated as:

```text
C(M) = Σ wt_pair_penalty × pair_penalty(p)
       + Σ wt_cross_dimer × cross_dimer_score(p_i, p_j)
       + Σ wt_off_target × off_target_score(p)
       + Σ wt_snp_penalty × snp_score(p)
```

Where the summation is over all primer pairs p in M and all primer interactions.

### Primer Binding Thermodynamics

Plexus implements the SantaLucia nearest-neighbor model for DNA duplex stability, considering:

- **Melting Temperature (Tm)**: Temperature at which 50% of primers are bound
- **Fraction Bound**: Percentage of primers bound at annealing temperature
- **ΔG Values**: Gibbs free energy for various secondary structures

The fraction bound approach is particularly advantageous for multiplex PCR as it ensures more uniform primer binding across different targets in the same reaction.

## Support and Community

### Getting Help

- **Issues**: Report bugs and request features on GitHub
- **Discussions**: Ask usage questions in GitHub Discussions
- **Documentation**: Check the latest docs in the repository

### Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Submit a pull request

### Citing Plexus

If you use Plexus in your research, please cite:

```text
Filges, S. (2026). Plexus: Automated Multiplex PCR Primer Panel Designer.
GitHub repository. https://github.com/sfilges/plexus
```

See `CITATIONS.md` for additional citation information.

## License

Plexus is licensed under the GNU General Public License v2.0 or later. See `LICENSE` for details.
