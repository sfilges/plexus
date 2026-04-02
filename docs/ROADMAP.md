# Plexus Development Roadmap

## Overview

This document tracks the work required to bring Plexus to a production-ready v1.0 release and
outlines the planned improvements for subsequent minor versions. It supersedes the earlier
`ISPCR_PLAN.md`, whose content has been absorbed into the v1.1 section below.

Items are grouped by milestone. Within each milestone, items are ordered by severity/dependency,
not necessarily by implementation order.

---

## v1.0 — Production Release

> **Released: 2026-03-03** — All 21 items resolved (19 implemented, 1 cancelled, 1 was TEST-01).

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

### ~~FEAT-02 · Cross-target amplicon classification~~ ✅ Cancelled

**Severity: Important · Files: `src/plexus/blast/offtarget_finder.py`, `specificity.py`**

~~Currently every amplicon that does not match the intended junction coordinates of a given pair is
counted as an off-target product. This includes amplicons formed by primers from *two different
junctions in the same panel* — these are intra-panel cross-products, not genomic off-targets, and
conflating the two categories misleads users and unfairly penalises pairs that happen to sit near
another target.~~

**Cancelled — for ctDNA panels, cross-target products are equally undesirable as genomic
off-targets and are correctly treated as `off_target` by default.** Any non-intended amplicon in
the multiplex reaction consumes primers, produces spurious signal, and reduces sensitivity. The
proposed `cross_target` classification would add complexity without benefit for this use case.

---

### ~~FEAT-03 · Remove dead code from `AmpliconFinder`~~ ✅ Fixed in v0.4.5

**Severity: Minor · File: `src/plexus/blast/offtarget_finder.py`**

`create_ontarget_dataframe()` and `create_offtarget_dataframe()` are never called by the
pipeline. They use an incompatible interface (a `primer_df` argument with a different schema
from what the pipeline provides) and are superseded by the coordinate-based `_is_on_target()`
logic in `specificity.py`. The TODO comment block in `__init__` and `find_amplicons()` should
be replaced with updated docstrings that reflect the current architecture.

**Fix:** Removed the two dead methods and the TODO comments. The class docstring now documents
the expected `bound_df` input schema and the generated `amplicon_df` output schema.

---

### ~~REPT-01 · Panel QC Reporting Module~~ ✅ Implemented in v0.5.7

**Severity: Important · Files: `src/plexus/pipeline.py`, new `src/plexus/reporting/qc.py`**

Clinical users require a high-level summary of the final panel quality and design statistics.

**Deliverables:** Create a `panel_qc.json` report containing:

- **Tm Distribution:** Mean, standard deviation, and range of Tms for all selected primers.
- **Sequence Flags:** Count of primers exceeding thresholds for high/low GC (>70% or <30%) or
  homopolymers.
- **Cross-Reactivity Matrix:** A summary matrix (junction vs junction) showing the count of
  potential dimer interactions (based on current dimer scores).

**Implementation:** Added `src/plexus/reporting/qc.py` with `generate_panel_qc()`, called from
`run_pipeline()` after `panel_summary.json` is written. The call is wrapped in `try/except` so
failures are non-fatal (logged as a warning). 16 unit tests added in
`tests/test_reporting_qc.py`.

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

### ~~SYS-01 · Pre-flight Disk Space Check~~ ✅ Fixed in v0.4.5

**Severity: Minor · Files: `src/plexus/utils/env.py`**

Before writing BLAST output or intermediate files, verify that the output directory has
sufficient free disk space. Issue a clear warning if available space is below a conservative
threshold (e.g., 2 GB). This prevents silent failures mid-pipeline on systems with full or
near-full disks.

**Fix:** Implemented `check_disk_space()` in `src/plexus/utils/env.py` and integrated it into the
`run_pipeline()` pre-flight checks.

---

### ~~CLI-01 · `plexus template` generation~~ ✅ Fixed in v0.4.6

**Severity: Important · File: `src/plexus/cli.py`**

Provide a command to generate starter files for a new design workspace:

- Write a template `junctions.csv` with the correct column headers and example rows.
- Write a `designer_config.json` pre-populated with the default preset values.

This gives new users a concrete starting point without reading the docs.

**v0.4.5 progress:** `plexus init` restructured to require explicit `--fasta`/`--snp-vcf`
by default.

**v0.4.6 progress:** Implemented as a dedicated `plexus template` command.

---

### ~~ARCH-02 · Coordinate system: `design_start` 0-based vs 1-based~~ ✅ Fixed in v0.4.5

**Severity: Minor · File: `src/plexus/designer/multiplexpanel.py:533`**

`design_start` is assigned as `junction_start - padding` where `junction_start` comes from the
input CSV (1-based genomic coordinates). Investigation confirmed that `design_start` IS correctly
1-based throughout the codebase — the pysam fetch correctly converts with `design_start - 1`, and
BLAST/SNP coordinate comparisons are correct.

The actual bug was in `calculate_junction_coordinates_in_design_region()`, where a misleading
comment stated `design_start` was 0-based and an erroneous `- 1` shifted the junction-relative
coordinates by one position: `junction.start - junction.design_start - 1` instead of the correct
`junction.start - junction.design_start`. This shifted the primer design window by 1 bp, masked
by the ±3 bp `junction_padding_bases`.

**Fix:** Removed the erroneous `- 1`, corrected the comment, and documented the 1-based convention
in the `extract_design_regions_from_fasta()` docstring. Added unit tests for the coordinate
calculation.

---

### ~~ARCH-03 · Replace `print()` with `logger` in `blast_runner.py`~~ ✅ Fixed in v0.4.4

**Severity: Minor · File: `src/plexus/blast/blast_runner.py:67`**

```python
print(f"BLAST database '{self.db_path}' already exists.")
```

This bypasses the loguru logger and breaks log capture in tests and pipeline log files.

**Fix:** `logger.info(f"BLAST database '{self.db_path}' already exists, skipping creation.")`

---

### ~~DOC-01 · Document input CSV format~~ ✅ Fixed in v0.4.6

**Severity: Important**

**Completed**: Added comprehensive input CSV documentation in the updated user guide (`docs/user_guide.md`).

**Covered Topics:**
- Required columns: `Name`, `Chrom`, `Five_Prime_Coordinate`, `Three_Prime_Coordinate`
- Optional `Panel` column for multi-panel design
- Coordinate system: 1-based genomic coordinates
- Chromosome naming requirements (must match FASTA headers)
- Complete example CSV with annotations
- Multi-panel design workflow
- Coordinate system explanation and examples

**Location**: See "Step-by-Step Workflow > Input Preparation" section in the user guide.

---

### ~~DOC-02 · Chromosome naming validation~~ ✅ Fixed in v0.4.6

**Severity: Important · File: `src/plexus/snpcheck/snp_data.py`, `src/plexus/resources.py`**

If the FASTA uses `chr1` notation but the gnomAD VCF uses `1` (or vice versa), bcftools region
queries silently return zero variants. There is currently no check or warning for this mismatch.

**Fix:** Added `detect_chrom_naming_mismatch()` helper that compares contig names and detects
systematic `chr`-prefix mismatches.

- **Runtime check**: The check runs inside `intersect_vcf_with_regions()` immediately after the
  VCF contig list is obtained, emitting a descriptive `logger.warning()`.
- **Early validation**: Integrated into `plexus init` via `src/plexus/resources.py`. In `research`
  mode, a mismatch triggers a warning during resource registration. In `compliance` mode, the
  mismatch is treated as a fatal error and aborts initialization.

This ensures users are notified of naming conflicts before starting a production run.

---

### ~~TEST-01 · End-to-end integration test with real BLAST~~ ✅ Implemented in v0.4.0

**File: `tests/test_integration.py`**

~~All existing integration tests mock BLAST, the annotator, and the AmpliconFinder. Add at least
one end-to-end test that exercises the full pipeline against the fixture BLAST database in
`tests/data/` without mocking, using the synthetic FASTA already present in `conftest.py`. This
would catch regressions in the BLAST runner, annotator annotation logic (including BUG-01 above),
and AmpliconFinder pairing logic together.~~

**Implementation:** `TestFullPipelineWithBlast` in `tests/test_integration.py` runs the full
pipeline against `tests/data/fixtures.ndb` without mocking BLAST, the annotator, or
`AmpliconFinder`. Tests are conditionally skipped when BLAST+ tools are absent from `PATH`. Four
tests cover: BLAST step completion, per-pair `specificity_checked` flag, BLAST output file
creation, and off-target product list population on selected pairs.

---

### ~~REPR-01~~ · ~~Chromosome naming check at `plexus run` time in compliance mode~~ ✅ Fixed in v0.5.4

**Severity: Important · Files: `src/plexus/pipeline.py`, `src/plexus/utils/utils.py`**

`detect_chrom_naming_mismatch()` is currently called only inside `init_genome()` (i.e., at
`plexus init` time). In compliance mode, `plexus run` is stateless — it bypasses the registry
and accepts `--fasta` and `--snp-vcf` directly. A user can therefore supply a FASTA with `chr1`
notation paired with a VCF using `1` notation and the pipeline will proceed without warning:
bcftools region queries silently return zero variants, SNP penalties are never applied, and the
provenance record gives no indication that the check was skipped.

**Required change:**

1. In `run_pipeline()`, after path validation and the compliance environment guard, call
   `detect_chrom_naming_mismatch()` against the supplied `fasta_file` and `snp_vcf` when both
   are present and SNP checking is not skipped.
2. In compliance mode (`op_mode == "compliance"`), a detected mismatch must raise a descriptive
   `ValueError` (consistent with the behaviour of `init_genome()` in compliance mode).
3. In research mode, emit a `logger.warning()` — also consistent with the existing init behaviour.
4. Record the outcome in `provenance.json` under a `chrom_naming_check` key:
   `"pass"`, `"mismatch"`, `"skipped"` (SNP check disabled or VCF not provided), or
   `"unavailable"` (contigs could not be read).

**Context:** `DOC-02` (✅ v0.4.6) added this check at init time and introduced the shared
`detect_chrom_naming_mismatch()` utility. This item extends its enforcement to the stateless
run path that compliance mode uses.

**Tests to add:** `tests/test_pipeline.py` (or `test_cli.py`) — mock `detect_chrom_naming_mismatch`
to return a mismatch string; assert `ValueError` is raised in compliance mode and only a warning
is logged in research mode.

---

### ~~REPR-02 · Random seed for stochastic selectors~~ ✅ Fixed in v0.5.2

**Severity: Important · Files: `src/plexus/selector/selectors.py`, `src/plexus/pipeline.py`, `src/plexus/cli.py`**

The `GreedySearch` and `RandomSearch` selectors call `random.shuffle()` and `random.choice()`
without a fixed seed. `SimulatedAnnealing` also relies on `random`. Two runs with identical
inputs can therefore produce different selected primer pairs, which violates the reproducibility
requirement for compliant clinical use and makes debugging non-deterministic.

**Required change:**

1. Add an optional `seed: int | None = None` parameter to `MultiplexSelector.__init__()` and
   propagate it to each concrete selector's `run()` method. Before the first random call, execute
   `random.seed(seed)` if a seed is provided.
2. Add `selector_seed: int | None` to `MultiplexPickerParameters` in `config.py` (default `None`
   — existing stochastic behaviour unchanged unless a seed is set).
3. Expose `--selector-seed` on `plexus run` and thread it through `run_pipeline()` and
   `run_multi_panel()`.
4. Record the seed (or `null`) in `provenance.json` under a `selector_seed` key so that a run
   can be reproduced exactly from the audit record alone.

**Compliance note:** In compliance mode, consider whether the absence of a seed should be a
warning or a hard error. A run without a fixed seed cannot be identically reproduced, which is
relevant for regulatory documentation. At minimum, `provenance.json` should make the absence
explicit (`"selector_seed": null`).

**Tests to add:** `tests/test_selector.py` — assert that two `GreedySearch` runs with the same
seed and input produce the same solution ordering; assert that two runs without a seed are
not guaranteed to match.

---

### ~~AUDT-02 · Include `primer3-py` in the compliance manifest~~ ✅ Fixed in v0.5.5

**Severity: Important · Files: `src/plexus/data/compliance_manifest.json`, `src/plexus/utils/env.py`**

The compliance manifest previously pinned only four system tools (`blastn`, `makeblastdb`,
`blast_formatter`, `bcftools`) and said nothing about `primer3-py`, the core thermodynamic
engine for all Tm, self-dimer, hairpin, and end-stability calculations.

**Fix:** Added a `python_packages` section to `compliance_manifest.json` with `primer3-py`
pinned at `2.3.0`. Extended `validate_environment()` to check Python package versions via
`importlib.metadata.version()` using the same structured verdict format (`"pass"` | `"fail"` |
`"missing"`). Python package verdicts are surfaced in the `compliance_environment` block of
`provenance.json` alongside the tool verdicts. Manifest version bumped from `"1.0"` to `"1.1"`.
Four new tests added to `tests/test_env.py`.

---

### ~~AUDT-03 · Record run success/failure status in `provenance.json`~~ ✅ Fixed in v0.5.4

**Severity: Important · File: `src/plexus/pipeline.py`**

`provenance.json` is written at the very start of `run_pipeline()`, before any pipeline work
begins (line 288). If the pipeline fails partway through — BLAST crash, all junctions failing
design, SNP check error — the provenance file still exists and appears valid to a later auditor.
There is no field indicating whether the run completed successfully, what step it reached, or
what errors occurred.

`panel_summary.json` does record `steps_completed` and `errors`, but the compliance guide
directs clinical users to inspect `provenance.json`. The two files tell inconsistent stories
on a failed run.

**Required change:**

1. Write an initial `provenance.json` at startup as today, but with `"status": "started"` and
   `"completed_at": null`.
2. At the end of `run_pipeline()`, update (re-write) `provenance.json` with:
   - `"status": "completed"` or `"status": "failed"`
   - `"completed_at"`: ISO 8601 UTC timestamp
   - `"steps_completed"`: list from `PipelineResult.steps_completed`
   - `"errors"`: list from `PipelineResult.errors`
3. Wrap the final write in a `try/finally` so it executes even if the pipeline raises an
   unhandled exception, capturing `"status": "failed"` and the exception message.

**Tests to add:** `tests/test_pipeline.py` — mock a mid-pipeline failure; assert that
`provenance.json` exists and contains `"status": "failed"` with a non-null `completed_at`.

---

### ~~ARCH-04 · Remove genome download functionality~~ ✅ Fixed in v0.5.1

**Severity: Important · Files: `src/plexus/resources.py`, `src/plexus/snpcheck/resources.py`, `src/plexus/cli.py`**

The `--download` flag on `plexus init` and the associated download infrastructure
(`_download_fasta()`, `_download_with_progress()`, `download_gnomad_vcf()`) pull reference
files from third-party URLs (UCSC, GATK GCS bucket) at runtime. This is a convenience for
exploratory research but is actively harmful for the intended clinical use case:

- **No chain of custody.** Downloaded files have no validated provenance — they come from a
  URL, not from an institution's verified data store.
- **Fragile.** URLs rot, GCS bucket policies change, and UCSC periodically reorganizes paths.
  A broken URL in a clinical tool is disproportionately disruptive.
- **Contradicts the compliance philosophy.** Compliance mode is built around users bringing
  their own verified resources. Providing a download shortcut signals the wrong workflow.

For ctDNA panels specifically, the reference FASTA and gnomAD VCF are large, stable, and
typically already available from institutional storage. No user running a clinical panel
should be downloading a 3.1 GB genome at pipeline runtime.

**Fix:** Removed `_download_fasta()`, `_download_with_progress()` from `resources.py` and
`download_gnomad_vcf()`, `_download_with_progress()` from `snpcheck/resources.py`. Removed
`--download` flag from `plexus init` and the `download` parameter from `init_genome()`.
`GENOME_PRESETS` is now a name-validation registry containing only `description`. The VCF
fallback in `get_snp_vcf()` now looks up the registry (via `get_registered_snp_vcf()`)
instead of checking a download cache. All download tests removed; init tests cover only the
local file registration path.

---

### ~~ARCH-05 · Allow registry use in compliance mode~~ ✅ Fixed in v0.5.3

**Severity: Important · Files: `src/plexus/cli.py`, `src/plexus/resources.py`**

Compliance mode currently disables the genome registry entirely, requiring `--fasta` on every
`plexus run` invocation. This conflates two independent concerns:

- **Statelessness** — required for containers (no persistent `~/.plexus/`)
- **Checksum enforcement** — required for compliance

Containers are stateless as a natural consequence of having no persistent volume for
`~/.plexus/`, not because of anything in the compliance mode code. On a non-containerized
clinical workstation — the common deployment for a diagnostic lab — the registry is a
legitimate and auditable workflow: register the validated FASTA once, run for hundreds of
patients without repeating the path, and rely on `--strict` (auto-enabled in compliance mode)
to verify the checksum on every run.

The current restriction forces clinical workstation users to type `--fasta /path/to/hg38.fa`
on every invocation with no compliance benefit, since the registry already stores and verifies
the exact same checksum.

**Required change:**

1. Remove the `if op_mode == "compliance": ... raise typer.Exit(code=1)` guard that blocks
   registry lookup when `--fasta` is omitted in `cli.py`.
2. Compliance mode should instead guarantee that `--strict` is always active. Registry lookup
   proceeds normally; checksum verification is enforced because `should_verify` is always
   `True` in compliance mode. This is already the case when `--fasta` *is* provided — extend
   it to the registry fallback path.
3. The container use case is unaffected: without a mounted `~/.plexus/` volume, the registry
   is simply empty and the existing "genome not initialized" error fires, prompting the user
   to provide `--fasta` explicitly.
4. Update the compliance guide and `plexus status` output to reflect that the registry is
   usable in compliance mode with automatic checksum enforcement.

**Tests to update:** `tests/test_cli.py` — replace the "compliance requires --fasta" assertion
with "compliance with registry fallback verifies checksums and fails on mismatch".

---

## v1.1 — Enhanced Specificity Analysis

> **v1.1.0 released: 2026-03-13** — Vectorized BLAST annotation and off-target finding,
> HTML visual QC report (REPT-02).
>
> **v1.2.0 released: 2026-03-26** — Enriched CSV output (12 thermodynamic columns),
> per-assay cross-dimer contribution, primer-level dimer matrix + heatmap, two-stage
> candidate pair cap (partial PERF-02), target dropout for DFS selector (partial ARCH-01),
> plexity clamping, Plotly.js bundle fix.

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

**v1.2.0 progress:** Target dropout is now implemented in the DFS selector
(`allow_target_dropping`, `dropout_stringency`, `min_target_fraction`), and plexity values
are clamped to the actual input count. However, the explicit plexity cost term in
`MultiplexCostFunction.calc_cost()` is still not implemented — the config fields `wt_lt`
and `wt_gt` remain unused.

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

### SPEC-01 · Assign cross-pair off-target amplicons to primer pairs

**Severity: Important · Files: `src/plexus/blast/specificity.py`, `src/plexus/selector/cost.py`**

`AmpliconFinder.find_amplicons()` correctly discovers amplicons formed by primers from different
pairs (e.g., forward from pair A + reverse from pair B on the same chromosome). However, the
mapping step in `run_specificity_check()` (step 6, line 112+) only looks up amplicons keyed by
a single pair's own `(f_id, r_id)`. Cross-pair amplicons exist in `amplicon_map` but are never
retrieved — they silently fall through because no single `PrimerPair` owns both primer IDs.

**Impact:** For multiplex panels, cross-pair off-target products are the most likely source of
spurious amplification (many more primer combinations than same-pair). These products are
computed by BLAST and identified by `AmpliconFinder`, but never assigned to any pair. They do
not appear in `off_targets.csv`, do not affect `filter_offtarget_pairs()`, and do not contribute
to the `wt_off_target` cost term. The specificity check is therefore blind to the dominant class
of off-target risk in multiplex reactions.

**Required changes:**

1. After the per-pair mapping loop in `run_specificity_check()`, perform a second pass that
   identifies cross-pair amplicons — entries in `amplicon_map` where `F_primer` belongs to one
   pair and `R_primer` belongs to a different pair (or different junction).
2. Assign each cross-pair amplicon to **both** contributing pairs (the pair that owns the forward
   primer and the pair that owns the reverse primer), so that both are penalised. Add a
   `cross_pair: bool` flag to the product dict to distinguish these from same-pair off-targets.
3. Ensure cross-pair off-targets are surfaced in `off_targets.csv` with both contributing pair
   IDs visible.
4. Verify that `filter_offtarget_pairs()` and `wt_off_target` in the cost function correctly
   account for the increased off-target counts.

**Design consideration:** A cross-pair amplicon penalises two pairs simultaneously. If pair A's
forward primer participates in cross-products with 5 different reverse primers from other pairs,
replacing pair A eliminates all 5 — but the current per-pair count treats them independently.
Consider whether cross-pair off-targets should carry a reduced weight or whether the existing
`wt_off_target` is sufficient given that the selector can swap individual pairs.

**Relationship to other items:**
- ISPCR-04 (cross-target interaction matrix) is a reporting view of the same underlying data;
  SPEC-01 makes the data actionable in the selector.
- ISPCR-02 (ΔG-weighted off-target cost) would naturally extend to cross-pair products once
  they are tracked.

**Tests to add:** `tests/test_specificity.py` — construct a mock `bound_df` with primers from
two different pairs that form a cross-pair amplicon; assert that both pairs receive the
off-target product. Assert that the off-target CSV includes the cross-pair product with both
pair IDs.

---

### ISPCR-04 · Improve AmpliconFinder: cross-target interaction matrix
**File: `src/plexus/blast/offtarget_finder.py`**

Build on the BLAST amplicon-finding infrastructure to produce a structured cross-target
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

**Note:** The per junction parallel processing of left and right primers has been removed in v0.5.7.
However, parallel processing across junctions might still be worth implementing.

---

### PERF-02 · Selector scalability: pre-filter candidates before optimisation

**File: `src/plexus/selector/selectors.py`**

For large panels (50+ junctions), the Greedy and Simulated Annealing selectors evaluate the
same high-cost pairs repeatedly across iterations. A pre-filtering step that limits each junction
to the top K candidates by `pair_penalty` before running optimisation would reduce the search
space without sacrificing solution quality (pairs with very high design penalty are unlikely to
appear in the optimal solution regardless of multiplex interactions).

Add a `max_candidates_per_junction` parameter to `MultiplexPickerParameters` (default: unlimited).

**v1.2.0 progress:** A two-stage candidate pair cap was implemented: `max_pairs_pre_filter`
(default 500) caps immediately after pair design, and `max_pairs_per_junction` (default 100)
caps after SNP/BLAST filtering. This addresses the memory/BLAST workload concern but is a
cap-based approach rather than the quality-ranked pre-filter described above. The remaining
value of this item is implementing a `pair_penalty`-ranked top-K selection rather than an
arbitrary cap.

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

### ~~REPT-02 · Visual QC Report (HTML)~~ ✅ Implemented in v1.1.0

**Severity: Low**

Transform the `panel_qc.json` (REPT-01) into a visual, standalone HTML report with Plotly.js
charts for Tm distribution, amplicon size, GC content, and a cross-reactivity heatmap.

**Implementation:** Added `src/plexus/reporting/html_report.py` with `generate_html_report()`
and a Jinja2 template (`panel_report.html.j2`) that inlines Plotly.js basic for fully
offline-capable reports. Integrated into the pipeline after `panel_qc.json` write. Added
`plexus report` CLI command for standalone re-generation. 11 tests in
`tests/test_reporting_html.py`.

---

### AUDT-04 · Improve status command

**Severity: Low**

plexus status currently shows whether system dependencies are installed, but not
the installed version. It also does not indicate whether these versions match
the compliance_manifest.json. In mode=compliance, this mismatch should be
indicated, in mode=research just showing the versions should be enough.

---

## v2.x - Future releases

### SPLIT-01 · Automated Panel Splitting

**Severity: Important · Prerequisite: ARCH-01 (plexity cost term)**

Partition large panels into K pools (e.g., 2×20-plex) to minimize intra-pool cross-dimer
interactions. Addresses cases where a single multiplex is thermodynamically impossible due
to unavoidable primer-pair interactions.

**Approach: two-phase split using the existing DFS selector.**

#### Phase 1 — Target interaction matrix + initial partition

Before invoking the per-pool DFS, build a **target-level cross-dimer interaction matrix**:

```
For each pair of targets (A, B):
    score[A, B] = min cross-dimer cost across top candidate pairs for A and B
```

This reuses the existing `PrimerDimerPredictor`. The v1.2.0 primer-level dimer matrix in
`reporting/qc.py` already computes a similar matrix for the selected panel; this extends it
to all candidates up front.

Apply spectral clustering or greedy graph partitioning on this matrix to produce an initial
K-way target assignment. This gives Phase 2 a structured warm start.

**New file:** `src/plexus/selector/pool_assigner.py`
- `build_target_interaction_matrix(pair_lookup, primer_df) → np.ndarray`
- `partition_targets(interaction_matrix, k, strategy) → list[list[target_id]]`
- Strategies: `"spectral"` (scipy/sklearn), `"greedy"` (no extra dependency)

#### Phase 2 — Per-pool DFS + inter-pool refinement

Run the existing `DepthFirstSearch` independently on each pool's assigned targets. The DFS
runs unmodified — it already handles variable-size target sets and computes normalized
cross-dimer costs.

After per-pool DFS, run Kernighan-Lin-style local search refinement:

```
for each target T in pool_i:
    for each pool_j ≠ pool_i:
        if moving T from pool_i → pool_j reduces total cost
        and both pools still satisfy plexity constraints:
            reassign T, re-run DFS on affected pools
```

This catches cases where the Phase 1 partition was suboptimal. The refinement is cheap
(linear in N × K) because it only re-evaluates the marginal cost of moving a single target.

#### Config changes

`MultiplexPickerParameters` already stubs `allow_split_panel: bool = False` and
`max_splits: int = 2`. Add:

- `split_strategy: str = "greedy"` — partition algorithm (`"greedy"` | `"spectral"`)
- `max_plexity_per_pool: int | None = None` — per-pool ceiling (defaults to
  `maximum_plexity` if unset)
- `min_plexity_per_pool: int | None = None` — per-pool floor (defaults to
  `minimum_plexity` if unset)

#### Pipeline integration

In `run_pipeline()`, after the post-filter cap (Step 4) and before multiplex optimization
(Step 5):

1. If `allow_split_panel` is enabled and the number of targets exceeds
   `max_plexity_per_pool`, invoke the pool assigner (Phase 1).
2. Run DFS per pool (Phase 2), then local-search refinement.
3. Output becomes multiple `Multiplex` objects tagged with a `pool_id`.

#### Output changes

- `selected_multiplex.csv` and `top_panels.csv` gain a `Pool` column.
- `panel_summary.json` gains a `pools` section listing targets per pool, per-pool cost,
  and per-pool cross-dimer statistics.
- `panel_qc.json` and `panel_report.html` show per-pool and cross-pool dimer heatmaps.

#### Connection to ARCH-01

With splitting, a pool may end up with fewer targets than `target_plexity`. The plexity
cost term (`wt_lt`, `wt_gt`) is needed to discourage lopsided splits (e.g., 28 targets in
pool 1, 2 in pool 2). **Implementing ARCH-01 first is a prerequisite.**

#### Tests

- `tests/test_pool_assigner.py`: interaction matrix construction, partition balance,
  determinism with seed.
- `tests/test_selector.py`: DFS per-pool produces valid solutions, refinement improves or
  maintains cost, plexity constraints respected after reassignment.
- Integration test: end-to-end split pipeline with a synthetic 40-target panel that has two
  obvious interaction clusters.

---

## Summary Table

| ID | Description | Version | Severity | Status |
| --- | --- | --- | --- | --- |
| ~~BUG-01~~ | ~~Fix `from_3prime` annotation (`qlen` vs `length`)~~ | ~~v1.0~~ | ~~Critical~~ | ✅ v0.4.4 |
| ~~BUG-02~~ | ~~Fix `product_bp` off-by-one in AmpliconFinder~~ | ~~v1.0~~ | ~~Minor~~ | ✅ v0.4.4 |
| ~~BUG-03~~ | ~~Fix BLAST+ v5 database detection~~ | ~~v1.0~~ | ~~Important~~ | ✅ v0.4.4 |
| ~~FEAT-01~~ | ~~Switch to `blastn-short` task for primer queries~~ | ~~v1.0~~ | ~~Important~~ | ✅ v0.4.5 |
| ~~FEAT-02~~ | ~~Cross-target vs off-target amplicon classification~~ | ~~v1.0~~ | ~~Important~~ | Cancelled |
| ~~FEAT-03~~ | ~~Remove dead code from AmpliconFinder~~ | ~~v1.0~~ | ~~Minor~~ | ✅ v0.4.5 |
| ~~REPT-01~~ | ~~Basic Panel QC Report (JSON)~~ | ~~v1.0~~ | ~~Important~~ | ✅ v0.5.7 |
| ~~AUDT-01~~ | ~~Tool versions and data checksums~~ | ~~v1.0~~ | ~~Important~~ | ✅ v0.4.5 |
| ~~SYS-01~~ | ~~Pre-flight disk space check~~ | ~~v1.0~~ | ~~Minor~~ | ✅ v0.4.5 |
| ~~CLI-01~~ | ~~`plexus template` generation~~ | ~~v1.0~~ | ~~Important~~ | ✅ v0.4.6 |
| ~~ARCH-02~~ | ~~Clarify `design_start` coordinate convention~~ | ~~v1.0~~ | ~~Minor~~ | ✅ v0.4.5 |
| ~~ARCH-03~~ | ~~Replace `print()` with `logger` in blast_runner~~ | ~~v1.0~~ | ~~Minor~~ | ✅ v0.4.4 |
| ~~DOC-01~~ | ~~Document input CSV format~~ | ~~v1.0~~ | ~~Important~~ | ✅ v0.4.6 |
| ~~DOC-02~~ | ~~Chromosome naming validation / pre-flight check~~ | ~~v1.0~~ | ~~Important~~ | ✅ v0.4.6 |
| ~~REPR-01~~ | ~~Chromosome naming check at `plexus run` time in compliance mode~~ | ~~v1.0~~ | ~~Important~~ | ✅ v0.5.4 |
| ~~REPR-02~~ | ~~Random seed for stochastic selectors~~ | ~~v1.0~~ | ~~Important~~ | ✅ v0.5.2 |
| ~~AUDT-02~~ | ~~Include `primer3-py` in compliance manifest~~ | ~~v1.0~~ | ~~Important~~ | ✅ v0.5.5 |
| ~~AUDT-03~~ | ~~Record run success/failure status in `provenance.json`~~ | ~~v1.0~~ | ~~Important~~ | ✅ v0.5.4 |
| ~~SCI-01~~ | ~~Weight SNP penalties by allele frequency~~ | ~~v1.0~~ | ~~Important~~ | ✅ v0.5.6 |
| ~~ARCH-04~~ | ~~Remove genome download functionality~~ | ~~v1.0~~ | ~~Important~~ | ✅ v0.5.1 |
| ~~ARCH-05~~ | ~~Allow registry use in compliance mode~~ | ~~v1.0~~ | ~~Important~~ | ✅ v0.5.3 |
| ARCH-01 | Implement plexity penalty in cost function | v1.1 | Important | Partial (v1.2.0) |
| ISPCR-01 | ntthal ΔG scoring for BLAST binding sites | v1.1 | Important | |
| ISPCR-02 | ΔG-weighted off-target cost in selector | v1.1 | Important | |
| ISPCR-03 | Template mispriming check | v1.1 | Low | |
| SPEC-01 | Assign cross-pair off-target amplicons to primer pairs | v1.1 | Important | |
| ISPCR-04 | Cross-target interaction matrix output | v1.1 | Low | |
| PERF-01 | Parallel thermodynamics by default | v1.1 | Low | |
| PERF-02 | Pre-filter candidates per junction before optimisation | v1.1 | Low | Partial (v1.2.0) |
| EXT-01 | Additional genome presets | v1.1 | Low | |
| EXT-02 | Chromosome naming normalisation | v1.1 | Low | |
| ~~REPT-02~~ | ~~Visual QC Report (HTML)~~ | ~~v1.1~~ | ~~Low~~ | ✅ v1.1.0 |
| AUDT-04 | Improve `plexus status` command | v1.1 | Low | |
| SPLIT-01 | Automated Panel Splitting (two-phase: partition + per-pool DFS) | v2.x | Important | Planned |
| ~~TEST-01~~ | ~~End-to-end integration test with real BLAST~~ | ~~v1.0~~ | ~~Important~~ | ✅ v0.4.0 |
