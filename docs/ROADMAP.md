# Plexus Development Roadmap

## Overview

This document tracks the work required to bring Plexus to a production-ready v1.0 release and
outlines the planned improvements for subsequent minor versions. It supersedes the earlier
`ISPCR_PLAN.md`, whose content has been absorbed into the v1.1 section below.

Items are grouped by milestone. Within each milestone, items are ordered by severity/dependency,
not necessarily by implementation order.

---

## v1.0 — Production Release

The v1.0 milestone covers correctness issues, science-blocking bugs, and architecture gaps that
would compromise the reliability of a clinical-grade tool.

---

### ~~BUG-01 · `from_3prime` annotation uses alignment length instead of query length~~ ✅ Fixed in v0.4.4

**Severity: Critical · File: `src/plexus/blast/annotator.py:28`**

The `from_3prime` annotation determines whether a BLAST hit reaches the 3′ end of the primer —
the prerequisite for polymerase extension and therefore off-target amplification.

**Current (wrong):**

```python
"from_3prime": lambda row: row["qend"] == row["length"],
```

`row["length"]` is the alignment length. This causes two classes of error:

- **False positive:** A 5′-anchored partial match (qstart=1, qend=15, alignment length=15 for a
  22 bp primer) satisfies `qend == length` and is labelled as a 3′ hit. The primer cannot
  actually extend from this site but gets counted as a predicted binding event.
- **False negative:** A genuine 3′-end partial match (qstart=8, qend=22, alignment length=15)
  fails the check (`22 ≠ 15`) and is silently discarded. The primer can extend but the hit is
  not flagged.

**Fix:** Replace `row["length"]` with `row["qlen"]` (query/primer length, already present in the
BLAST column set `BLAST_COLS`).

```python
"from_3prime": lambda row: row["qend"] == row["qlen"],
```

**Impact:** Affects off-target count in every BLAST run, which feeds directly into the cost
function via `wt_off_target`. False positives inflate off-target counts and bias the selector
away from good primers; false negatives suppress real off-target warnings.

**Tests to add:** `tests/test_blast_annotator.py` — explicit cases for 5′-anchored match, 3′-end
partial match, and full-length match, asserting correct `from_3prime` values in each.

---

### ~~BUG-02 · `product_bp` off-by-one in `AmpliconFinder`~~ ✅ Fixed in v0.4.4

**Severity: Minor · File: `src/plexus/blast/offtarget_finder.py:91`**

```python
product_bp=row["sstart"] - F_start,   # wrong: off by 1
```

For a forward primer at genomic position 1000 and a reverse primer whose BLAST `sstart`
(5′ end on minus strand = rightmost genomic position) is 1200, the amplicon spans 1000–1200
inclusive = 201 bp. The current formula returns 200.

**Fix:** `product_bp = row["sstart"] - F_start + 1`

---

### ~~BUG-03 · BLAST+ v5 database detection fails~~ ✅ Fixed in v0.4.4

**Severity: Important · File: `src/plexus/blast/blast_runner.py:67`**

Database existence is checked against v4 suffixes only (`.nhr`, `.nin`, `.nsq`). BLAST+ ≥ 2.12
produces v5 databases with a different and expanded set of files (`.ndb`, `.nhr`, `.nin`, `.nsq`,
`.njs`, etc.). On systems with BLAST+ v5, an existing database may not be detected, causing
unnecessary re-indexing; or a v5 database may only have newer files and the check incorrectly
reports the database as absent.

**Fix:** Prefer checking for the `.njs` manifest file (v5) or `.nhr` (v4), falling back to
attempting `makeblastdb` only if neither is found. Alternatively, probe with a zero-hit
`blastn` dry-run to confirm the database is valid.

---

### ~~FEAT-01 · Switch BLAST to `blastn-short` task for primer-length queries~~ ✅ Fixed in v0.4.5

**Severity: Important · File: `src/plexus/blast/blast_runner.py:run()`**

For query sequences shorter than ~30 bp, NCBI recommends `-task blastn-short`, which reduces the
word size to 7 and uses scoring parameters tuned for short queries (reward 1, penalty −3, gap
costs 5/2). The current default (`word_size=11`) seeds with 11-mers, which can miss binding sites
where the 3′-terminal 11 bases contain a mismatch — exactly the scenario most relevant to PCR
off-target analysis.

**Fix:** `blastn-short` is now the default task. The hardcoded `word_size=11` in `specificity.py`
has been removed; `blastn-short` manages its own defaults (word_size=7). An explicit `word_size`
can still be passed to `BlastRunner.run()` if needed.

---

### FEAT-02 · Cross-target amplicon classification

**Severity: Important · Files: `src/plexus/blast/offtarget_finder.py`, `specificity.py`**

Currently every amplicon that does not match the intended junction coordinates of a given pair is
counted as an off-target product. This includes amplicons formed by primers from *two different
junctions in the same panel* — these are intra-panel cross-products, not genomic off-targets, and
conflating the two categories misleads users and unfairly penalises pairs that happen to sit near
another target.

**Required change:**

In `specificity.py`, after `_is_on_target()` returns False, run a second check: iterate all
*other* junctions in the panel and test `_is_on_target(prod, other_junction, other_pair)` for all
candidate pairs at that junction. If any match, classify as `cross_target` rather than
`off_target`.

Add a `cross_target_products` list field to `PrimerPair` alongside `off_target_products`.
Update `_build_enriched_pair_row()` to expose a `Cross_Target_Count` column and update
`save_off_targets_csv()` to include a `Type` column (`off_target` vs `cross_target`).

Do NOT penalise cross-target products in the cost function by default. They represent
co-amplification of other intended targets (acceptable in some applications) and can be addressed
by adjusting panel layout, not by discarding pairs.

**Tests to add:** `tests/test_blast_specificity.py` — panel with two adjacent junctions whose
primers fall within amplicon range of each other; assert cross-target classification.

---

### FEAT-03 · Remove dead code from `AmpliconFinder`

**Severity: Minor · File: `src/plexus/blast/offtarget_finder.py`**

`create_ontarget_dataframe()` and `create_offtarget_dataframe()` are never called by the
pipeline. They use an incompatible interface (a `primer_df` argument with a different schema
from what the pipeline provides) and are superseded by the coordinate-based `_is_on_target()`
logic in `specificity.py`. The TODO comment block in `__init__` and `find_amplicons()` should
be replaced with updated docstrings that reflect the current architecture.

**Fix:** Remove the two dead methods and the TODO comments. Document the expected DataFrame
schema in the class docstring.

---

### REPT-01 · Panel QC Reporting Module

**Severity: Important · Files: `src/plexus/pipeline.py`, new `src/plexus/reporting/qc.py`**

Clinical users require a high-level summary of the final panel quality and design statistics.

**Deliverables:** Create a `panel_qc.json` report containing:

- **Tm Distribution:** Mean, standard deviation, and range of Tms for all selected primers.
- **Sequence Flags:** Count of primers exceeding thresholds for high/low GC (>70% or <30%) or
  homopolymers.
- **Cross-Reactivity Matrix:** A summary matrix (junction vs junction) showing the count of
  potential dimer interactions (based on current dimer scores).

This JSON can be expanded in later versions (v1.1) to support visual plots and HTML reports.

---

### ~~AUDT-01 · Compliance & Integrity (Locked Versions & Checksums)~~ ✅ Implemented in v0.4.5

**Severity: Important · Files: `src/plexus/pipeline.py`, `utils/env.py`, `resources.py`**

To satisfy clinical audit and reproducibility requirements:

- **Tool Versions:** Capture and log the exact versions of system dependencies (`blastn`,
   `bcftools`) and library versions (`primer3-py`) in the `panel_summary.json`.
- **Data Checksums:** Calculate and store SHA-256 checksums for the input reference
   FASTA and SNP VCF files. Verify checksums at runtime in compliance mode or with `--strict`.
- **Provenance output:** `provenance.json` written to output directory with tool versions,
   file checksums, operational mode, and run timestamp. Also embedded in `panel_summary.json`.
- **Operational modes:** `research` (default) and `compliance` modes stored in global config.
   Compliance mode enforces checksum verification on every run.
- **User-provided checksums:** `--checksums` flag accepts a `sha256sum`-format file for
   verifying FASTA and VCF integrity at init time.

---

### SYS-01 · Pre-flight Disk Space Check

**Severity: Minor · Files: `src/plexus/utils/env.py`**

Before writing BLAST output or intermediate files, verify that the output directory has
sufficient free disk space. Issue a clear warning if available space is below a conservative
threshold (e.g., 2 GB). This prevents silent failures mid-pipeline on systems with full or
near-full disks.

---

### CLI-01 · `plexus init` template generation (partially addressed in v0.4.5)

**Severity: Important · File: `src/plexus/cli.py`**

Extend the existing `plexus init` command to generate starter files for a new design workspace:

- Write a template `junctions.csv` with the correct column headers and an example row.
- Write a `designer_config.json` pre-populated with the default preset values, ready to edit.

This gives new users a concrete starting point without reading the docs. Keep the command
non-interactive — parameters are passed as flags, consistent with the existing `--genome`
and `--fasta` options.

**v0.4.5 progress:** `plexus init` restructured to require explicit `--fasta`/`--snp-vcf`
by default; downloads demoted to `--download` flag. Added `--mode` and `--checksums` flags.
Template generation (junctions.csv, designer_config.json) still pending.

**Deferred to v1.x:** Interactive preset selection wizard, vendor-specific oligo ordering
formats (e.g. IDT plate layout), and workspace-level checksum enforcement profiles.

---

### ARCH-02 · Coordinate system: `design_start` 0-based vs 1-based

**Severity: Minor · Files: `src/plexus/designer/multiplexpanel.py:495`, `blast/specificity.py`**

`design_start` is assigned as `junction_start - padding` where `junction_start` comes from the
input CSV (1-based genomic coordinates). pysam's `fasta.fetch(chrom, start, end)` uses 0-based
half-open coordinates, so the fetch effectively starts one base downstream of the stored
`design_start` value. The ±5 bp tolerance in `_is_on_target()` masks this in practice, but the
coordinate semantics are inconsistent.

**Fix (preferred):** Store `design_start` as the 1-based genomic coordinate of the first extracted
base (`= junction_start - padding + 1`) and document the convention explicitly. Adjust the pysam
fetch call to use `design_start - 1` (converting to 0-based). This makes `design_start`
semantically match BLAST's 1-based `sstart` and removes the hidden off-by-one.

Alternatively, tighten the tolerance comment in `_is_on_target()` to explicitly note that 1 bp of
the tolerance budget is consumed by this coordinate offset.

---

### ~~ARCH-03 · Replace `print()` with `logger` in `blast_runner.py`~~ ✅ Fixed in v0.4.4

**Severity: Minor · File: `src/plexus/blast/blast_runner.py:67`**

```python
print(f"BLAST database '{self.db_path}' already exists.")
```

This bypasses the loguru logger and breaks log capture in tests and pipeline log files.

**Fix:** `logger.info(f"BLAST database '{self.db_path}' already exists, skipping creation.")`

---

### DOC-01 · Document input CSV format

**Severity: Important**

The required and optional column names for the junction CSV (`Name`, `Chrom`,
`Five_Prime_Coordinate`, `Three_Prime_Coordinate`, `Panel`) are inferred only from reading the
`JunctionInput` Pydantic model. No user-facing documentation exists. Add a `docs/input_format.md`
or a section in the README that specifies:

- Required columns and types
- Coordinate convention (1-based, half-open, closed?)
- Optional `Panel` column behaviour
- Chromosome naming requirements (must match FASTA header names exactly)
- Example CSV

---

### DOC-02 · Chromosome naming validation

**Severity: Important · File: `src/plexus/snpcheck/snp_data.py`**

If the FASTA uses `chr1` notation but the gnomAD VCF uses `1` (or vice versa), bcftools region
queries silently return zero variants. There is currently no check or warning for this mismatch.
Add a pre-flight validation that compares a sample of FASTA sequence names against the VCF contig
list and warns (or errors) on a naming mismatch. The existing `_get_vcf_contigs()` function
already retrieves the VCF contig set; the FASTA contig names can be obtained from pysam's
`FastaFile.references`.

---

## v1.1 — Enhanced Specificity Analysis

v1.1 focuses on improving the scientific accuracy of the BLAST-based specificity check,
incorporating thermodynamic ΔG scoring for individual binding sites and advancing toward a
genuine in-silico PCR (isPCR) model. These items are valuable but not blocking for clinical
utility; the v1.0 BLAST approach (once BUG-01 is fixed) is a reasonable heuristic.

---

### ARCH-01 · Implement plexity penalty in cost function

**Files: `src/plexus/config.py`, `src/plexus/selector/cost.py`**

`MultiplexPickerParameters` defines `target_plexity`, `minimum_plexity`, `maximum_plexity`,
`wt_lt`, and `wt_gt`, but `MultiplexCostFunction.calc_cost()` never reads these fields. The
current selectors always include every junction that has at least one valid candidate pair, so
the plexity term would currently be a no-op.

**Fix:** Implement the plexity cost term:
`cost += wt_lt × max(0, target_plexity - n_pairs) + wt_gt × max(0, n_pairs - target_plexity)`

This becomes meaningful when SPLIT-01 (automated panel splitting) is implemented, where the
selector may need to leave some targets out of a given pool to meet thermodynamic constraints.
Implementing it now keeps the config fields honest and makes the infrastructure ready.

> **Note for v1.0 users:** Currently, all junctions with at least one valid primer pair are
> always included in the multiplex. The `target_plexity`, `minimum_plexity`, and
> `maximum_plexity` config parameters are reserved for a future release that supports partial
> panel selection and pool splitting.

---

### ISPCR-01 · Thermodynamic ΔG scoring for BLAST binding sites (ntthal integration)

**Files: `src/plexus/blast/annotator.py`, new `src/plexus/blast/binding_scorer.py`**

The current binding prediction (`predicted_bound`) relies on a length/identity threshold plus
an evalue gate. Neither metric directly answers "will polymerase extend from this site?"

Replace or augment with `primer3.bindings.calc_heterodimer()` (wraps `ntthal`), which computes
a full nearest-neighbor ΔG for each primer–genomic-sequence pair. This requires extracting the
local genomic sequence at each BLAST hit position (a `pysam.FastaFile` fetch of primer-length
sequence), then calling `calc_heterodimer(primer_seq, genomic_seq)`.

**Approach:**

1. Add `binding_scorer.py` with a `BindingSiteScorer` class that accepts the BLAST hit DataFrame
   and a pysam FASTA handle, extracts local sequences at each hit, and computes ΔG per hit.
2. Add a `dg` column to the bound DataFrame.
3. Update `build_annotation_dict()` to accept a `dg_threshold` parameter (e.g., −9 kcal/mol)
   as an additional `predicted_bound` criterion alongside the existing length/evalue checks.
4. Surface per-binding-site ΔG in `off_targets.csv` for audit.

**Rationale:** The Johnston/Plexus aligner already handles intra-pair and cross-pair dimer
scoring via a 3′-end-biased heuristic. The ntthal model is more physically accurate for
individual non-specific binding events (handles mismatches in the middle of the duplex) and
is already a project dependency via `primer3-py`.

---

### ISPCR-02 · Off-target amplicon ΔG cost in the selector

**Files: `src/plexus/selector/cost.py`, `src/plexus/designer/primer.py`**

Once ISPCR-01 provides per-binding-site ΔG values, replace the current blunt off-target
count penalty (`wt_off_target × len(off_target_products)`) with a ΔG-weighted version:

```
off_target_penalty = Σ max(0, dg_threshold - ot["dg"]) for ot in off_target_products
```

This means weak off-target binding sites contribute less cost than strong ones, reducing the
frequency with which perfectly good primers are discarded for marginal BLAST hits.

Add `wt_off_target_dg` to `MultiplexPickerParameters`. Keep `wt_off_target` as a fallback
for runs where BLAST was performed without ΔG scoring (backward compatibility).

---

### ISPCR-03 · Template mispriming check (`PRIMER_MAX_TEMPLATE_MISPRIMING_TH` equivalent)

**Files: `src/plexus/designer/thal.py`, `design.py`**

Primer3's `PRIMER_MAX_TEMPLATE_MISPRIMING_TH` checks each candidate primer against the full
template sequence for internal self-priming (a primer binding elsewhere within the same design
region). For short amplicons (60–120 bp) on 400 bp design regions, the genome-wide BLAST check
already covers this. However, for repeat-rich junctions or sequences with secondary structure,
an explicit intra-template check adds confidence for clinical validation.

**Implementation:** After thermodynamic filtering in `calculate_single_primer_thermodynamics()`,
score each candidate primer against the full `design_region` string using `calc_heterodimer()`.
If the best binding score (excluding the on-target position) exceeds `PRIMER_MAX_TEMPLATE_MISPRIMING_TH`
(ΔG threshold, configurable), discard the primer.

Add the parameter to `SingleplexDesignParameters` (default: disabled/None).

**Priority note:** This is lower priority than ISPCR-01/02 because BLAST is a superset of this
check for well-characterised genomes. It matters most for novel or poorly-annotated genomes where
the BLAST database may be incomplete.

---

### ISPCR-04 · Improve AmpliconFinder: cross-target interaction matrix
**File: `src/plexus/blast/offtarget_finder.py`**

Build on FEAT-02 (v1.0 cross-target classification) to produce a structured cross-target
interaction matrix: for each ordered pair of junctions (A, B), record whether any forward primer
from A can co-amplify with any reverse primer from B. Output this as `cross_target_matrix.csv`
with junctions on both axes and a count (or minimum product size) in each cell.

This is primarily a reporting/audit feature for clinical lab documentation rather than a change
to the optimisation. It answers "do any primer pairs from different targets interact?" — a common
question during panel validation.

---

### PERF-01 · Enable parallel primer thermodynamics by default
**File: `src/plexus/designer/design.py`**

`design_multiplex_primers()` accepts a `parallel=False` parameter that enables concurrent
left/right primer thermodynamic evaluation using `ProcessPoolExecutor`. This is currently disabled
by default and not exposed in the CLI or pipeline API. For panels with many junctions, parallel
evaluation would halve thermodynamic calculation time.

**Fix:** Expose `parallel` as a pipeline/CLI option (`--parallel-design`), default to `True`
when more than N junctions are present, or always default to `True` on multi-core systems.

---

### PERF-02 · Selector scalability: pre-filter candidates before optimisation

**File: `src/plexus/selector/selectors.py`**

For large panels (50+ junctions), the Greedy and Simulated Annealing selectors evaluate the
same high-cost pairs repeatedly across iterations. A pre-filtering step that limits each junction
to the top K candidates by `pair_penalty` before running optimisation would reduce the search
space without sacrificing solution quality (pairs with very high design penalty are unlikely to
appear in the optimal solution regardless of multiplex interactions).

Add a `max_candidates_per_junction` parameter to `MultiplexPickerParameters` (default: unlimited).

---

### EXT-01 · Additional genome presets

**File: `src/plexus/resources.py`**

Add presets for:
- `hg19` / `GRCh37` (legacy clinical datasets)
- `mm39` (mouse, common model organism)
- `GRCh38_no_alt` (analysis-set FASTA without alt contigs, preferred for clinical pipelines)

Each preset requires: FASTA URL, gnomAD VCF URL (or equivalent population variant source), and
a note on expected file size.

---

### EXT-02 · Chromosome naming normalisation

**Files: `src/plexus/utils/utils.py`, pipeline input handling**

When the input CSV uses UCSC notation (`chr1`) but the FASTA or VCF uses Ensembl notation (`1`),
all downstream operations silently fail or return empty results. Add an optional
`chrom_prefix` parameter to the pipeline (or auto-detect from FASTA headers) and apply it
consistently to junction coordinates, VCF queries, and BLAST sequence IDs.

---

### REPT-02 · Visual QC Report (HTML)

**Severity: Low**

Transform the `panel_qc.json` (REPT-01) into a visual, standalone HTML report. This should
include Plotly or Seaborn charts for distributions and a searchable heatmap for the
cross-reactivity matrix, facilitating rapid review by clinical lab staff.

---

### SPLIT-01 · Automated Panel Splitting

**Severity: Future**

Implement an optimization algorithm to partition large panels into multiple pools (e.g., 2x20plex)
to minimize intra-pool dimer interactions. This addresses cases where a single multiplex is
thermodynamically impossible due to unavoidable primer-pair interactions. This is a highly
complex feature, and out of scope for the moment, but would be a great addition to the
project.

---

### TEST-01 · End-to-end integration test with real BLAST

**File: `tests/test_integration.py`**

All existing integration tests mock BLAST, the annotator, and the AmpliconFinder. Add at least
one end-to-end test that exercises the full pipeline against the fixture BLAST database in
`tests/data/` without mocking, using the synthetic FASTA already present in `conftest.py`. This
would catch regressions in the BLAST runner, annotator annotation logic (including BUG-01 above),
and AmpliconFinder pairing logic together.

---

## Summary Table

| ID | Description | Version | Severity | Status |
|---|---|---|---|---|
| ~~BUG-01~~ | ~~Fix `from_3prime` annotation (`qlen` vs `length`)~~ | ~~v1.0~~ | ~~Critical~~ | ✅ v0.4.4 |
| ~~BUG-02~~ | ~~Fix `product_bp` off-by-one in AmpliconFinder~~ | ~~v1.0~~ | ~~Minor~~ | ✅ v0.4.4 |
| ~~BUG-03~~ | ~~Fix BLAST+ v5 database detection~~ | ~~v1.0~~ | ~~Important~~ | ✅ v0.4.4 |
| ~~FEAT-01~~ | ~~Switch to `blastn-short` task for primer queries~~ | ~~v1.0~~ | ~~Important~~ | ✅ v0.4.5 |
| FEAT-02 | Cross-target vs off-target amplicon classification | v1.0 | Important | |
| FEAT-03 | Remove dead code from AmpliconFinder | v1.0 | Minor | |
| REPT-01 | Basic Panel QC Report (JSON) | v1.0 | Important | |
| ~~AUDT-01~~ | ~~Tool versions and data checksums~~ | ~~v1.0~~ | ~~Important~~ | ✅ v0.4.5 |
| SYS-01 | Pre-flight disk space check | v1.0 | Minor | |
| CLI-01 | `plexus init` template generation | v1.0 | Important | Partial v0.4.5 |
| ARCH-02 | Clarify `design_start` coordinate convention | v1.0 | Minor | |
| ~~ARCH-03~~ | ~~Replace `print()` with `logger` in blast_runner~~ | ~~v1.0~~ | ~~Minor~~ | ✅ v0.4.4 |
| DOC-01 | Document input CSV format | v1.0 | Important | |
| DOC-02 | Chromosome naming validation / pre-flight check | v1.0 | Important | |
| ARCH-01 | Implement plexity penalty in cost function | v1.1 | Important | |
| ISPCR-01 | ntthal ΔG scoring for BLAST binding sites | v1.1 | Important | |
| ISPCR-02 | ΔG-weighted off-target cost in selector | v1.1 | Important | |
| ISPCR-03 | Template mispriming check | v1.1 | Low | |
| ISPCR-04 | Cross-target interaction matrix output | v1.1 | Low | |
| PERF-01 | Parallel thermodynamics by default | v1.1 | Low | |
| PERF-02 | Pre-filter candidates per junction before optimisation | v1.1 | Low | |
| EXT-01 | Additional genome presets | v1.1 | Low | |
| EXT-02 | Chromosome naming normalisation | v1.1 | Low | |
| REPT-02 | Visual QC Report (HTML) | v1.1 | Low | |
| SPLIT-01 | Automated Panel Splitting | Future | Future | |
| TEST-01 | End-to-end integration test with real BLAST | v1.1 | Important | |
