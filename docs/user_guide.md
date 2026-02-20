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
    - [Two Deployment Paths](#two-compliance-deployment-paths)
    - [The Compliance Manifest](#the-compliance-manifest)
    - [Environment Validation](#environment-validation)
    - [Compliance Provenance Record](#compliance-provenance-record)
    - [Transitioning from Research to Compliance](#transitioning-from-research-to-compliance)

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
| Checksum Verification | Optional (`--strict` flag or `--checksums`) | **Mandatory** |
| Resource Validation | Warning only | **Error on mismatch** |
| Tool Version Enforcement | None | **Exact versions via manifest** |
| Provenance Tracking | Basic | **Comprehensive + `compliance_environment`** |
| Audit Requirements | None | **Full documentation** |
| Registry Required | Yes (default) | No — `--fasta` explicit, `--checksums` for data |
| Use Case | Development, research | Clinical, diagnostic, regulated, containerized |

### Setting Operational Mode

**Priority order** — highest wins:

| Priority | Source | How to set |
| -------- | ------ | ---------- |
| 1 | `PLEXUS_MODE` env var | `export PLEXUS_MODE=compliance` or bake into Docker image |
| 2 | `~/.plexus/config.json` | `plexus init --mode compliance` |
| 3 (default) | built-in | `research` |

```bash
# Container / CI: set via env var (no workspace needed)
export PLEXUS_MODE=compliance

# Local research workstation: persist in config
plexus init --genome hg38 --mode compliance --checksums my_checksums.sha256

# Check current effective mode (respects env var)
plexus status
```

**Important Compliance Notes:**

- The `PLEXUS_MODE` env var **overrides the config file** — it is the authority for containers
- In compliance mode, `--fasta` is required (registry fallback is disabled)
- All resource files must have SHA-256 checksums for `--checksums` stateless verification
- Tool versions are validated against the bundled compliance manifest before any data is touched
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

There are two distinct resource setup paths — choose based on your deployment context:

**Path A — Local / Research (registry-based):**

```bash
# Initialize genome resources (registers FASTA + VCF, computes checksums)
plexus init --genome hg38 --download

# For compliance mode with verified resources
plexus init --genome hg38 --mode compliance \
    --fasta hg38.fa --snp-vcf gnomad.vcf.gz \
    --checksums my_checksums.sha256

# Check resource status
plexus status
```

**Path B — Container / Stateless (no registry):**

```bash
# No init needed — just create a checksums file and pass it at run time
sha256sum hg38.fa gnomad.vcf.gz > checksums.sha256

# Then at run time:
PLEXUS_MODE=compliance plexus run \
    --fasta hg38.fa \
    --snp-vcf gnomad.vcf.gz \
    --checksums checksums.sha256 \
    --input junctions.csv \
    --output results/
```

**Resource Types (Path A only):**

- **FASTA**: Reference genome sequence
- **FAI**: FASTA index file (built by `init`)
- **BLAST DB**: BLAST database for specificity checking (built by `init`)
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

### Docker & Containerized Runs

For consistent execution (especially in clinical use cases), it is recommended to run Plexus via its official Docker image. Plexus includes a `docker` command that simplifies the process of mounting host directories and translating paths.

```bash
plexus docker \
  --tag 0.5.0-strict \
  --input ./data/junctions.csv \
  --fasta ./data/hg38.fa \
  --output ./results/run_1
```

**What the `plexus docker` command does:**
- **Auto-Mounting**: Automatically creates Docker volume mounts (`-v`) for the parent directories of all your input files (input CSV, FASTA, SNP VCF, config).
- **Path Translation**: Translates host file paths (e.g. `./data/hg38.fa`) into container-relative paths (e.g. `/mnt/vol0/hg38.fa`).
- **Audit Consistency**: Ensures that the output directory is mounted correctly and writable.

#### Image Registries and Availability
The `plexus docker` command uses a **Local-First** approach:
- **Local Search**: It first checks if the specified image (registry + tag) is available locally on your machine.
- **Auto-Pull**: If the image is not found locally, it automatically attempts a `docker pull`.
- **Private Registries**: Use the `--registry` flag to specify a custom registry prefix. You must be authenticated via `docker login` before running.
- **Forced Update**: Use the `--pull` flag to force an update from the registry even if a local version is found.

```bash
# Example using a private registry and forcing a pull
plexus docker \
  --registry my-internal-registry.com/plexus \
  --tag 0.5.0-strict \
  --pull \
  --input junctions.csv ...
```

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
        "plexus_version": "0.5.0",
        "primer3_version": "2.6.1",
        "run_timestamp": "2026-02-20T14:30:00Z",
        "operational_mode": "compliance",
        "fasta_path": "/path/to/hg38.fa",
        "fasta_sha256": "a1b2c3...",
        "snp_vcf_path": "/path/to/gnomad.vcf.gz",
        "snp_vcf_sha256": "d4e5f6...",
        "run_blast": true,
        "skip_snpcheck": false,
        "compliance_environment": {
            "manifest_version": "1.0",
            "blastn": {"expected": "2.17.0", "actual": "2.17.0", "verdict": "pass"},
            "bcftools": {"expected": "1.23", "actual": "1.23", "verdict": "pass"}
        }
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
| `--strict` | False | Verify checksums via registry before running |
| `--checksums` | None | SHA-256 checksums file for stateless verification |
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

### Design Principle

> **"The code is the authority; the environment is a suspect."**

In compliance mode, Plexus does not trust the environment it runs in. It verifies both **data integrity** (checksums) and **environment integrity** (exact tool versions) before touching any file.

### Authority Hierarchy

| Priority | Authority | What it controls |
|----------|-----------|-----------------|
| 1 | Compliance manifest (bundled in package) | Exact tool versions required |
| 2 | `PLEXUS_MODE` env var | Activates compliance mode — overrides config file |
| 3 | `~/.plexus/config.json` | Convenience for local workstations |

### Two Compliance Deployment Paths

#### Path A — Local Workstation (Registry-Based)

Use `plexus init` to register resources with checksums. Compliance mode verifies against the registry automatically on every run.

```bash
# 1. Generate checksums for your validated resources
sha256sum hg38.fa gnomad.vcf.gz > clinical_checksums.sha256

# 2. Initialize in compliance mode
plexus init \
  --genome hg38 \
  --mode compliance \
  --fasta hg38.fa \
  --snp-vcf gnomad.vcf.gz \
  --checksums clinical_checksums.sha256

# 3. Run (checksums verified against registry automatically)
plexus run --input clinical_junctions.csv --output results/
```

#### Path B — Container / Stateless (No Registry)

The recommended path for automated pipelines and Docker. No `~/.plexus/` workspace is needed. `PLEXUS_MODE=compliance` is baked into the container image. Files and checksums are injected at runtime.

```bash
# Build the image (BLAST 2.17.0 + bcftools 1.23 compiled in)
docker build -f docker/DOCKERFILE -t plexus:latest .

# Run — no registry, no init
docker run \
  -v /data/hg38.fa:/mnt/hg38.fa \
  -v /data/gnomad.vcf.gz:/mnt/gnomad.vcf.gz \
  -v /data/checksums.sha256:/mnt/checksums.sha256 \
  -v /data/junctions.csv:/mnt/junctions.csv \
  -v /data/output:/mnt/output \
  plexus:latest \
  run \
    --input /mnt/junctions.csv \
    --fasta /mnt/hg38.fa \
    --snp-vcf /mnt/gnomad.vcf.gz \
    --checksums /mnt/checksums.sha256 \
    --output /mnt/output
```

**Container execution sequence:**
1. `PLEXUS_MODE=compliance` detected (baked in) → no registry consulted
2. `--fasta` is explicit → no registry lookup
3. CLI parses `checksums.sha256`, verifies FASTA + VCF → hashes stored in memory
4. `run_pipeline()` calls `validate_environment()` → exact tool versions checked against the bundled compliance manifest → `ComplianceError` in <1 s if versions don't match
5. Pipeline runs
6. `provenance.json` contains: verified checksums + `compliance_environment` verdict block

### The Compliance Manifest

The file `src/plexus/data/compliance_manifest.json` is bundled immutably inside the package — like a version number, it is part of the software's identity. It declares the exact system tool versions required for compliance-mode operation:

```json
{
  "version": "1.0",
  "description": "Compliance manifest for Plexus. Exact tool versions required for clinical operation.",
  "tools": {
    "blastn":          { "exact_version": "2.17.0", "version_regex": "(\\d+\\.\\d+\\.\\d+)" },
    "makeblastdb":     { "exact_version": "2.17.0", "version_regex": "(\\d+\\.\\d+\\.\\d+)" },
    "blast_formatter": { "exact_version": "2.17.0", "version_regex": "(\\d+\\.\\d+\\.\\d+)" },
    "bcftools":        { "exact_version": "1.23",   "version_regex": "(\\d+\\.\\d+)"        }
  }
}
```

The manifest version (`"1.0"`) is **independent of the plexus version**. It increments only when the required tool set changes — not with every plexus release.

The Dockerfile's `ARG BLAST_VERSION` / `ARG BCFTOOLS_VERSION` must match the manifest's `exact_version` values. If you re-validate against a different tool version, update both.

### Environment Validation

In compliance mode, `run_pipeline()` calls `validate_environment()` before touching any data. It:

1. Loads the bundled manifest
2. Runs each required tool with `--version` / `-version`
3. Extracts the version string using the manifest's regex
4. Compares to `exact_version`
5. Raises `ComplianceError` with a consolidated message if anything fails:

```
ComplianceError: Environment validation failed — 2 tool(s) do not match the compliance manifest:
  - blastn: expected exactly 2.17.0, found '2.14.0+'
  - bcftools: not found on PATH
```

Verdicts per tool: `"pass"` | `"fail"` | `"missing"` | `"unparseable"`

### Compliance Provenance Record

In compliance mode, `provenance.json` includes a `compliance_environment` block:

```json
{
  "plexus_version": "0.5.0",
  "operational_mode": "compliance",
  "fasta_sha256": "a1b2c3d4...",
  "snp_vcf_sha256": "e5f6g7h8...",
  "compliance_environment": {
    "manifest_version": "1.0",
    "blastn":          { "expected": "2.17.0", "actual": "2.17.0", "verdict": "pass" },
    "makeblastdb":     { "expected": "2.17.0", "actual": "2.17.0", "verdict": "pass" },
    "blast_formatter": { "expected": "2.17.0", "actual": "2.17.0", "verdict": "pass" },
    "bcftools":        { "expected": "1.23",   "actual": "1.23",   "verdict": "pass" }
  }
}
```

### Generating a Checksums File

```bash
# Generate from any set of files (sha256sum format)
sha256sum hg38.fa gnomad.vcf.gz > checksums.sha256

# Verify format
cat checksums.sha256
# a1b2c3...  hg38.fa
# 7g8h9i...  gnomad.vcf.gz
```

The `--checksums` flag on `plexus run` accepts this file. Plexus matches entries by **filename only** (not path), so the file can live anywhere.

### Transitioning from Research to Compliance

```bash
# 1. Develop and iterate in research mode
plexus init --genome hg38 --mode research --download
plexus run --input test_junctions.csv --output test_results/

# 2. Generate checksums for your validated reference files
sha256sum hg38.fa gnomad.vcf.gz > clinical_checksums.sha256

# 3a. Local workstation: switch registry to compliance mode
plexus init --genome hg38 --mode compliance \
    --fasta hg38.fa --snp-vcf gnomad.vcf.gz \
    --checksums clinical_checksums.sha256

# 3b. Container: bake PLEXUS_MODE=compliance into image (see docker/DOCKERFILE)

# 4. Run in compliance mode
plexus run --input clinical_junctions.csv --output clinical_results/
```

### Compliance Mode Troubleshooting

**"--fasta is required in compliance mode"**
- Compliance mode does not fall back to the registry — always pass `--fasta` explicitly.

**"Environment validation failed — N tool(s) do not match"**
- Tool versions on PATH do not match the compliance manifest.
- In Docker: check that the image was built with the correct `ARG BLAST_VERSION` / `ARG BCFTOOLS_VERSION`.
- On a workstation: install the exact versions listed in `src/plexus/data/compliance_manifest.json`.

**"no entry for hg38.fa in checksums file"** (compliance mode only)
- The `--checksums` file must contain an entry for the FASTA filename. Check spelling and that the file was generated from the correct path.

**"FASTA checksum mismatch"**
- The on-disk file does not match the checksums file. Verify the file has not been modified or corrupted.

### Audit and Documentation Requirements

For clinical/compliance use, maintain:

1. **Resource Documentation:**
   - Source of all reference files
   - Genome build version, acquisition date
   - Checksum verification records (`sha256sum -c checksums.sha256`)

2. **Run Documentation:**
   - Input files (CSV, config, checksums)
   - Exact command line used
   - Plexus version + compliance manifest version (from `provenance.json`)
   - Full output directory including `provenance.json`

3. **Validation Records:**
   - Wet-lab validation results for selected primers
   - Any manual overrides or exceptions
   - Final primer sequences ordered

### Example Clinical Documentation Structure

```text
project/
├── inputs/
│   ├── clinical_junctions.csv
│   ├── validated_config.json
│   └── clinical_checksums.sha256
├── resources/
│   ├── hg38.fa           (sha256 verified)
│   └── gnomad.vcf.gz     (sha256 verified)
├── outputs/
│   ├── run_2026-02-20/
│   │   ├── selected_multiplex.csv
│   │   ├── provenance.json   ← includes compliance_environment
│   │   └── ...
│   └── run_2026-03-01/
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
