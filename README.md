# plexus

[![CI](https://github.com/sfilges/plexus/actions/workflows/ci.yml/badge.svg)](https://github.com/sfilges/plexus/actions/workflows/ci.yml)
[![Python](https://img.shields.io/badge/python-3.10%E2%80%933.13-blue)](https://github.com/sfilges/plexus)
[![License: GPL v2+](https://img.shields.io/badge/License-GPL_v2+-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0)

`plexus` is a Python-based bioinformatics tool designed to automate the creation of multiplex PCR panels. It specifically targets workflows like personalized ctDNA panels, integrating genomic data processing, primer design (via `primer3`), and specificity checking (via BLAST) to generate optimized primer sets for multiple targets simultaneously.

## Features

- **Automated Primer Design**: Uses a custom k-mer enumeration algorithm (`simsen`) to generate primer candidates for each junction; `primer3-py` is used for thermodynamic filtering (hairpin/self-dimer ΔG).
- **SNP Checking**: Filters primer candidates that overlap common variants in gnomAD (or a user-supplied VCF). Supports strict mode (`--snp-strict`) to discard any SNP-overlapping pair.
- **Multiplex Optimization**: Selects the optimal primer combination minimizing cross-dimer potential. Five algorithms available: Greedy (default), Random, BruteForce, SimulatedAnnealing, DFS.
- **Specificity Checking**: Integrates BLAST to check for off-target amplification and primer specificity.
- **Multi-Panel Support**: Designs multiple independent panels from a single input CSV using a `Panel` column (e.g., for multiple patients).
- **Configuration Presets**: Includes `default` and `lenient` configuration presets for different design stringencies.
- **CLI Interface**: Easy-to-use command line interface for running the design pipeline.

## Installation

### Prerequisites

- Python 3.10–3.13
- NCBI BLAST+ (for specificity checks — `blastn` must be on `$PATH`)
- `bcftools` (for SNP checking — must be on `$PATH`)

### Using uv (Recommended)

This project uses `uv` for package management.

```bash
git clone https://github.com/sfilges/plexus
cd plexus
uv pip install -e .
```

### Using Conda

You can also set up the environment using Conda:

```bash
git clone https://github.com/sfilges/plexus
cd plexus
conda env create -f config/environment.yml
conda activate plexus-run
pip install -e .
```

## Usage

The primary interface is the `plexus` CLI.

### Basic Command

To run the complete design pipeline:

```bash
plexus run \
  --input data/junctions.csv \
  --fasta data/genome.fa \
  --output results/ \
  --name my_panel
```

### Input file format

The input file should be a CSV containing the target junctions with the following columns:

```csv
Name,Chrom,Five_Prime_Coordinate,Three_Prime_Coordinate
EGFR_T790M,chr7,55181378,55181378
KRAS_G12D,chr12,25245350,25245350
...
```

See `data/junctions.csv` for a complete example.

### Multi-Panel Input

If the input CSV contains a `Panel` column, junctions are grouped by panel value and each panel is designed independently. Results are saved to `<output>/<panel_id>/`.

```csv
Name,Chrom,Five_Prime_Coordinate,Three_Prime_Coordinate,Panel
EGFR_T790M,chr7,55181378,55181378,panel_a
KRAS_G12D,chr12,25245350,25245350,panel_a
TP53_R248W,chr17,7674220,7674220,panel_b
```

Use `--parallel` to run panels concurrently.

### SNP Checking

SNP checking requires the bundled gnomAD VCF. Download it once before the first run:

```bash
plexus init
```

Check resource status at any time:

```bash
plexus status
```

To use a custom VCF (must be tabix-indexed):

```bash
plexus run -i junctions.csv -f genome.fa --snp-vcf /path/to/custom.vcf.gz
```

Use `--snp-strict` to discard any primer pair that overlaps a SNP above the AF threshold (default 1%), rather than just penalising it.

### Key Options

| Option | Default | Description |
| --- | --- | --- |
| `-i, --input` | required | Path to input CSV file |
| `-f, --fasta` | required | Path to reference genome FASTA |
| `-o, --output` | `./output` | Output directory |
| `-n, --name` | `multiplex_panel` | Panel name |
| `-g, --genome` | `hg38` | Reference genome name |
| `-p, --preset` | `default` | Config preset (`default` or `lenient`) |
| `-c, --config` | — | Path to custom JSON config file |
| `-s, --selector` | `Greedy` | Selector algorithm: `Greedy`, `Random`, `BruteForce`, `SimulatedAnnealing`, `DFS` |
| `--skip-blast` | false | Skip BLAST specificity check |
| `--skip-snpcheck` | false | Skip SNP overlap check |
| `--snp-vcf` | bundled gnomAD | Path to tabix-indexed VCF for SNP checking |
| `--snp-af-threshold` | `0.01` | Minimum allele frequency for SNP flagging |
| `--snp-strict` | false | Discard primer pairs overlapping SNPs |
| `--padding` | `200` | Bases to extract around each junction |
| `--parallel` | false | Run multiple panels in parallel |
| `--max-workers` | auto | Max parallel workers (multi-panel mode) |

Run `plexus --help` for a full list of commands and options.

## Configuration

Design parameters are controlled via JSON config files. Use `--config` to supply your own, or `--preset` to choose a built-in starting point (`default` or `lenient`).

**Partial configs are fully supported** — you only need to include the parameters you want to change. Any omitted parameter falls back to the preset default.

### Config file structure

The config is divided into five sections:

| Section | Controls |
| --- | --- |
| `singleplex_design_parameters` | Primer length, Tm, GC%, thermodynamic thresholds, adapter tails |
| `primer_pair_parameters` | Amplicon size, Tm difference, pair penalty weights |
| `pcr_conditions` | Salt concentrations, annealing temperature, thermodynamic tables |
| `snp_check_parameters` | AF threshold, SNP penalty weight, strict mode |
| `multiplex_picker_parameters` | Plexity targets, optimization weights, selector settings |

### Minimal override example

To widen the Tm window and use a stricter cross-dimer weight, create a file such as `my_config.json`:

```json
{
    "singleplex_design_parameters": {
        "PRIMER_MIN_TM": 55.0,
        "PRIMER_MAX_TM": 66.0
    },
    "multiplex_picker_parameters": {
        "wt_cross_dimer": 3.0
    }
}
```

Then run:

```bash
plexus run -i junctions.csv -f genome.fa --config my_config.json
```

All parameters not listed in your file are inherited from the `default` preset.

### Adapter tail sequences

`forward_tail` and `reverse_tail` set the adapter sequences prepended to each primer (at the 5′ end). The defaults are SiMSen-Seq-style adapters. Override them in `singleplex_design_parameters`:

```json
{
    "singleplex_design_parameters": {
        "forward_tail": "ACACTCTTTCCCTACACGACGCTCTTCCGATCT",
        "reverse_tail": "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"
    }
}
```

The full ready-to-order sequences (`tail + binding region`) are written to the `Forward_Full_Seq` / `Reverse_Full_Seq` columns in the output CSVs. Bare binding sequences are also retained in `Forward_Seq` / `Reverse_Seq`.

### Viewing all available parameters

The built-in preset files in `config/` list every available parameter with its default value and are a useful reference when building a custom config.
