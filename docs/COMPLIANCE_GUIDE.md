# Plexus Compliance & Clinical Workflow Guide

This guide describes the **Stateless & Verified** workflow for using Plexus in clinical, diagnostic, or regulated environments. This approach ensures total reproducibility by locking tool versions in a container and enforcing data integrity via mandatory checksums.

## 1. Environment Setup

The compliance workflow requires **Docker** to ensure that the execution environment (OS, libraries, and bioinformatics tools) is identical across every run.

### System Requirements

- **Docker Engine** 20.10+
- **Disk Space**: ~30GB (for genome resources and the Docker image)
- **Memory**: 16GB+ recommended for clinical panels

### Building the Compliance Image

The Plexus Dockerfile uses a multi-stage build to compile exact tool versions (e.g., BLAST 2.17.0, bcftools 1.23) and bakes a **Compliance Manifest** into the image.

```bash
# Clone the repository
git clone https://github.com/sfilges/plexus
cd plexus

# Build the versioned compliance image
docker build -t plexus:0.5.0-strict -f docker/DOCKERFILE .
```

---

## 2. Input & Resource Preparation

In compliance mode, Plexus is **stateless**. It does not rely on a local registry (`registry.json`). Instead, it verifies every file on-the-fly against a provided checksums file.

### Organize your Data

Create a structured project directory:

```text
my_clinical_project/
├── data/
│   ├── hg38.fa                # Reference Genome
│   ├── hg38.fa.fai            # FASTA Index (required)
│   ├── hg38.nhr, .nin, ...    # BLAST Database (required)
│   └── gnomad.vcf.gz          # SNP Database (with .tbi)
├── inputs/
│   └── junctions.csv          # Your target coordinates
└── results/                   # Output directory
```

### Generate the Compliance Checksums

You must provide a `sha256sum`-formatted file containing the hashes of your reference resources. This file acts as the "Evidence of Integrity" for your audit trail.

```bash
cd my_clinical_project/data/

# Generate hashes for the genome and SNP database
sha256sum hg38.fa > clinical_resources.sha256
sha256sum gnomad.vcf.gz >> clinical_resources.sha256
```

---

## 3. Running the Pipeline

You can run the pipeline using the raw `docker run` command or the convenient `plexus docker` wrapper.

### Path A: The `plexus docker` Wrapper (Recommended)

Plexus includes a built-in command that automatically handles the complex volume mounting and path translation for you.

```bash
plexus docker \
  --tag 0.5.0-strict \
  --input ./inputs/junctions.csv \
  --fasta ./data/hg38.fa \
  --snp-vcf ./data/gnomad.vcf.gz \
  --checksums ./data/clinical_resources.sha256 \
  --output ./results/run_2026-02-20
```

**What the wrapper does for you:**

- **Host Path Resolution**: Resolves relative paths to absolute host paths.
- **Auto-Mounting**: Mounts parent directories as **read-only** volumes.
- **Path Translation**: Translates host paths to internal container paths.

### Registry & Image Availability

The `plexus docker` command uses a **Local-First** approach:

1. **Local Search**: It first checks for the image (e.g., `plexus:0.5.0-strict`) on your machine.
2. **Auto-Pull**: If missing, it attempts to pull from the registry (default: `ghcr.io/sfilges/plexus`).
3. **Private Registries**: Use the `--registry` flag to specify a custom registry prefix. You must be logged in via `docker login` before running.
4. **Forced Refresh**: Use `--pull` to force an update from the registry even if a local version exists.

```bash
plexus docker \
  --registry my-private-registry.com/clinical-tools/plexus \
  --tag 0.5.0-strict \
  --pull \
  --input junctions.csv ...
```

### Path B: The Raw `docker run` Command

For full control or integration into orchestrators (like Nextflow/Airflow), use `docker run` directly.

> **Note:** The data directory is mounted read-only (`:ro`). This means the BLAST database
> (`hg38.nhr`, `.nin`, etc.) **must be pre-built** before running — `makeblastdb` cannot write
> to a read-only mount. If you need Plexus to build the database on first run, use the
> `plexus docker` wrapper (Path A), which mounts the FASTA parent directory as `:rw`.

```bash
docker run --rm \
  -v $(pwd)/data:/data:ro \
  -v $(pwd)/inputs:/inputs:ro \
  -v $(pwd)/results:/results \
  plexus:0.5.0-strict \
  run \
    --input /inputs/junctions.csv \
    --fasta /data/hg38.fa \
    --snp-vcf /data/gnomad.vcf.gz \
    --checksums /data/clinical_resources.sha256 \
    --output /results/run_$(date +%F)
```

---

## 4. Why this is Compliant

1. **Environment Authority**: The container is "born" with `PLEXUS_MODE=compliance`. It will refuse to run if the internal `blastn` or `bcftools` versions don't match the hardcoded manifest.
2. **Data Integrity**: Plexus hashes `/data/hg38.fa` at runtime and compares it to `/data/clinical_resources.sha256`. If a single byte has changed, the pipeline halts.
3. **No Hidden State**: Because it is stateless, there is no "remembered" data in `~/.plexus`. Every run is a fresh, verified start.

---

## 5. Reviewing the Audit Trail

Every compliant run produces a `provenance.json` and a `panel_summary.json` in the output directory. These files are critical for clinical documentation.

### The Provenance Record (`provenance.json`)

This file captures the "Who, What, and How" of the run:

- **`operational_mode`**: Confirms the run was in `compliance` mode.
- **`compliance_environment`**: A report showing that every tool (BLAST, bcftools) passed the version check against the manifest.
- **`fasta_sha256`**: Records the exact hash of the genome used.
- **`tool_versions`**: Records the raw version strings of all dependencies.

### Clinical Validation Example

To verify a previous run for an auditor:

```bash
# Check the image metadata
docker inspect plexus:0.5.0-strict --format='{{json .Config.Labels}}' | jq

# Compare the hash in provenance.json with your original resource record
cat results/run_2026-02-20/provenance.json | jq .fasta_sha256
```

---

## 6. Troubleshooting Compliance Errors

### "Version mismatch for blastn"

This occurs if the Docker image was built with a different version of BLAST than what the Plexus manifest expects. Always use the provided `DOCKERFILE` which pins these versions.

### "Checksum mismatch for hg38.fa"

The genome file on disk does not match the hash in your `.sha256` file. Ensure that the genome was not compressed, decompressed, or edited after the hash was generated.

### "--fasta is required in compliance mode"

Stateless mode disables the automatic lookup of "registered" genomes. You must always provide the explicit path to the FASTA file using `--fasta`.
