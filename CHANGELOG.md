# Changelog

All notable changes to plexus will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.4.5] - 19-02-2026

### Added

- **AUDT-01 · Compliance & integrity (tool versions + data checksums)**: Pipeline runs now produce a `provenance.json` in the output directory capturing plexus version, primer3-py version, system tool versions (`blastn`, `bcftools`, `makeblastdb`), file checksums (SHA-256), operational mode, and run timestamp. Provenance is also embedded in `panel_summary.json`.
- **Two operational modes**: `research` (default) and `compliance`, set via `plexus init --mode` and stored in `~/.plexus/config.json`. In compliance mode, resource checksums are verified automatically before every pipeline run.
- **SHA-256 checksum verification**: `plexus init` computes and stores SHA-256 checksums for the reference FASTA and SNP VCF in the genome registry. Users can supply a `sha256sum`-format checksums file via `--checksums` to verify files at init time — init fails immediately on mismatch.
- **`--strict` flag on `plexus run`**: Verifies resource checksums against the registry before running the pipeline. Always enabled in compliance mode.
- **`--download` flag on `plexus init`**: Downloads from preset URLs are no longer the default. Provide local files with `--fasta`/`--snp-vcf`, or explicitly opt in to downloading with `--download`.
- **Tool version capture**: New `get_tool_version()`, `get_tool_versions()`, `get_plexus_version()`, and `get_primer3_version()` utilities in `plexus.utils.env`.
- **Enhanced `plexus status`**: Now shows operational mode and truncated SHA-256 checksums alongside resource readiness.
- **FEAT-01 · `blastn-short` task for primer queries**: `BlastRunner.run()` now uses `-task blastn-short` by default, tuned for primer-length queries (<30 bp) with word_size=7, reward 1, penalty −3, and gap costs 5/2. The hardcoded `word_size=11` in `specificity.py` has been removed. This improves sensitivity for off-target binding sites where the 3′-terminal region contains mismatches.

### Fixed

- **ARCH-02 · Off-by-one in junction coordinate calculation** (`src/plexus/designer/multiplexpanel.py:533`): `calculate_junction_coordinates_in_design_region()` computed junction-relative coordinates as `junction.start - design_start - 1`, but since both values are 1-based genomic coordinates, the correct formula is `junction.start - design_start` (no `- 1`). The erroneous subtraction shifted the primer design window by 1 bp, masked by the ±3 bp `junction_padding_bases`. A misleading comment incorrectly stated `design_start` was 0-based; corrected to document the 1-based convention. Added unit tests for the coordinate calculation.

### Changed

- **`plexus init` requires explicit files by default**: Running `plexus init` without `--fasta` (and without `--download`) now errors with an actionable message, preventing accidental multi-GB downloads.
- **`get_cache_dir()` consolidated**: Moved from `plexus.snpcheck.resources` to `plexus.resources` to serve as the single source of truth. The old import path continues to work via re-export.
- **gnomAD URLs deduplicated**: `plexus.snpcheck.resources` now derives URLs from `GENOME_PRESETS` instead of hardcoding them.

## [0.4.4] - 19-02-2026

### Fixed

- **BUG-01 · `from_3prime` annotation uses alignment length instead of query length** (`src/plexus/blast/annotator.py:28`): `row["length"]` (alignment length) replaced with `row["qlen"]` (query/primer length). The old comparison caused two classes of silent error in every BLAST run: 5′-anchored partial matches were labelled as 3′ hits (false positive — primer cannot extend but was counted as predicted_bound, inflating off-target counts); genuine 3′-end partial matches failed the check (false negative — primer can extend but the off-target was not detected). The fix uses `qlen`, which is already present in `BLAST_COLS` and loaded into the DataFrame on every run. Tests updated with physically consistent fixture data and a new `test_from_3prime_annotation_semantics` regression test covering both scenarios explicitly.
- **BUG-02 · `product_bp` off-by-one in `AmpliconFinder`** (`src/plexus/blast/offtarget_finder.py:90`): `product_bp` was computed as `row["sstart"] - F_start`, which under-counts by one base because BLAST coordinates are inclusive on both ends (a primer at 1000 paired with a reverse hit at 1200 spans 201 bp, not 200). Fixed by adding `+ 1`. A dedicated regression test `test_amplicon_size_inclusive_coordinates` and a parametrized `test_amplicon_sizes_various` suite were added to `tests/test_offtarget_finder.py`.
- **BUG-03 · BLAST+ v5 database detection fails** (`src/plexus/blast/blast_runner.py:66`): The database existence check required all three v4 files (`.nhr`, `.nin`, `.nsq`) to be present. BLAST+ ≥ 2.12 (v5) produces a `.njs` JSON manifest and may omit some v4 files, causing the check to report the database as absent and triggering an unnecessary `makeblastdb` re-run. The check now detects either a v5 database (`.njs` manifest) or a v4 database (`.nhr` header file), whichever is present. Two new tests — `test_create_database_skips_if_exists_v5` and `test_create_database_skips_if_exists_v4` — verify each path independently.
- **ARCH-03 · `print()` replaced with `logger` in `blast_runner.py`** (`src/plexus/blast/blast_runner.py:68`): The "database already exists" message was emitted via `print()`, bypassing loguru and breaking log capture in tests and pipeline log files. Replaced with `logger.info()`, consistent with the rest of the project. Fixed together with BUG-03.

## [0.4.3] - 19-02-2026

### Added

- **Independent SNP penalty weight in cost function**: `wt_snp_penalty` (default `3.0`) added to `MultiplexPickerParameters` as a first-class cost function weight, independent of `wt_pair_penalty`. Previously, `snp_penalty` was merged directly into `pair_penalty` before optimization, making it impossible to tune separately and causing SNP-affected pairs to be preferred over clean alternatives when thermodynamic scores differed only modestly. With the new default, a 3′ SNP (snp_penalty=30) now contributes 90 to cost vs. 30 before — large enough to prefer a clean pair unless its thermodynamic penalty exceeds the SNP pair's by more than 90. Set `wt_snp_penalty: 0` to disable SNP weighting entirely. Added to both `designer_default_config.json` (`3.0`) and `designer_lenient_config.json` (`1.0`).
- **Per-junction SNP observability**: `run_snp_check()` now logs a per-junction summary (`Junction X: N/M pairs overlap SNPs (K clean pair(s) available)`) immediately after processing each junction, making it trivial to diagnose whether a SNP-affected junction had clean alternatives.

### Changed

- **`pair_penalty` is now purely thermodynamic**: `run_snp_check()` no longer adds `snp_penalty` to `pair_penalty`. The `Pair_Penalty` column in output CSVs now reflects only the thermodynamic component; `SNP_Penalty` remains a separate column and is now the sole vehicle through which SNP cost enters the optimizer.

### Fixed

- **`bcftools index` contig listing incompatibility**: `_get_vcf_contigs()` previously used `bcftools index -l` (removed in bcftools 1.21+) then `--stats` (requires count metadata absent from tabix-created indices). Now uses `bcftools view --header-only` and parses `##contig=<ID=...>` lines from the VCF header, which works regardless of bcftools version or how the index was created. Previously caused 4 retries with exponential backoff (~7s wasted) before falling back to running without contig pre-filtering.
- **On-target warning noise**: `run_specificity_check()` was emitting one `WARNING` per candidate pair that had no on-target BLAST hit, producing thousands of log lines for junctions in repetitive regions (e.g. GOLGA2). Per-pair messages are now demoted to `DEBUG`; a single `WARNING` per junction summarises the count (`N/M pairs have no on-target amplicon detected`).
- **`_is_on_target()` reverse primer coordinate bug**: The expected BLAST position of the reverse primer was computed as `design_start + reverse.start`, but BLAST reports the 5′ end (rightmost base) of a minus-strand hit, so the correct value is `design_start + reverse.start + reverse.length - 1`. The old formula caused on-target amplicons to be misclassified as off-target whenever the reverse primer was longer than 1 bp.

## [0.4.2] - 19-02-2026

### Added

- **Intra-pair dimer score in cost function**: `pair.dimer_score` (F/R heterodimer, computed by the Johnston algorithm during `find_primer_pairs()`) is now included in `MultiplexCostFunction.calc_cost()`. The penalty is `wt_pair_dimer * max(0, -dimer_score)`, so only actual dimer formation (negative scores) contributes. New `wt_pair_dimer` weight (default `1.0`) added to `MultiplexPickerParameters` and `designer_default_config.json`. Set to `0.0` to disable.
- **On-target BLAST verification**: After `run_specificity_check()`, each `PrimerPair` now carries an `on_target_detected: bool | None` flag (`None` = BLAST not run; `True` = confirmed on-target amplicon found; `False` = no on-target amplicon detected, accompanied by a `logger.warning`). The field surfaces as `On_Target_Detected` in both `selected_multiplex.csv` and `top_panels.csv`.

### Fixed

- **`F_target`/`R_target` always `"SEQ"` in `AmpliconFinder`**: Synthetic primer IDs (`SEQ_N`) were split on `"_"` to recover a target name, but this always produced `"SEQ"`. `aggregate_primers()` now builds `primer_target_map` (synthetic ID → junction name) at aggregation time. `AmpliconFinder` and `BlastResultsAnnotator` accept an optional `target_map` kwarg and use it to populate `F_target`/`R_target` with real junction names (e.g. `"EXON1"`). When a sequence is shared across multiple junctions the names are joined with `"|"`. `specificity.py` passes the map through automatically; callers that omit it fall back to the old split behaviour.

## [0.4.1] - 19-02-2026

### Added

- **Pre-flight dependency checks**: The pipeline now verifies that required system tools (`blastn`, `makeblastdb`, `blast_formatter`, `bcftools`) are available on the `$PATH` at the very beginning of a run, providing clear installation instructions if any are missing.
- **System tool status**: The `plexus status` command now displays the availability of required system dependencies alongside genome resources.
- **Centralized environment utility**: New `plexus.utils.env` module for robust cross-platform executable verification.

### Changed

- Refactored `blast_runner.py` and `snp_data.py` to use centralized tool verification logic.

### Removed

- **Cleaned up dead code**: Removed unused placeholder files and methods (`archive.py`, `in_silico_pcr.py`, and `biopython_process_blast_results`) to streamline the codebase for production.
- **Removed download-resources** from CLI, as it was redundant with the addition of the init command.

## [0.4.0] - 18-02-2026

### Added

- **GC clamp filter**: New `check_gc_clamp()` function rejects primers with fewer than 1 or more than 3 G/C bases in the last 5 bases of the 3' end. Controlled by `primer_gc_clamp` config field (0=off, 1=on). Enabled by default in the `default` preset, disabled in `lenient`.
- **Position-weighted SNP penalties**: SNPs near the primer's 3' end now receive a higher penalty than those at the 5' end. New config fields `snp_3prime_window` (default 5 bp) and `snp_3prime_multiplier` (default 3.0x) in `SnpCheckParameters`. `_calc_weighted_snp_penalty()` computes orientation-aware distance from the 3' end for both forward and reverse primers.
- **Coordinate-based BLAST off-target classification**: New `_is_on_target()` helper in `specificity.py` compares BLAST hit genomic coordinates against expected primer positions derived from `junction.design_start`, replacing the product-size heuristic that could misclassify pseudogene hits as on-target.
- **Per-junction error recovery**: If primer design fails for a single junction (region too small, no k-mers found), the pipeline logs a warning and continues with the remaining junctions instead of crashing. Failed junctions are tracked in `PipelineResult.failed_junctions` and written to `failed_junctions.csv`.
- **BLAST-enabled integration tests**: New `TestFullPipelineWithBlast` test class exercises the full pipeline with `run_blast=True`. Auto-skips when BLAST+ tools are not on PATH.
- 40+ new unit tests across `test_utils.py`, `test_snpcheck.py`, `test_blast_specificity.py`, and `test_selector.py`.

### Fixed

- **BruteForce selector storage bug**: When the buffer was full and a better solution arrived, `.insert()` added to `stored_multiplexes` but never updated `stored_costs`, causing the two lists to desynchronize and the buffer to grow beyond `store_maximum`. Now uses index-based replacement consistent with the DFS selector.

### Changed

- `_count_snps_in_region()` now returns `(count, list[tuple[int, float]])` (position + allele frequency) instead of `(count, list[int])`, enabling position-aware penalty calculation.
- Default `snp_penalty_weight` increased from 5.0 to 10.0 to better reflect typical primer pair penalty magnitudes.
- `run_snp_check()` accepts new keyword arguments `snp_3prime_window` and `snp_3prime_multiplier`, passed through from config.

### Removed

- `SaltCorrectionFormula` and `ThermodynamicTable` enums from `config.py`. The code always uses SantaLucia1998 parameters; the enums suggested user-selectable alternatives that were never wired up. The `salt_correction_formula` and `thermodynamic_table` fields are removed from `PCRConditions` and both JSON preset configs.

## [0.3.4] - 18-02-2026

### Added

- `Forward_Full_Seq` / `Reverse_Full_Seq` columns in `candidate_pairs.csv` and `selected_multiplex.csv` containing the adapter tail prepended to the binding sequence, giving ready-to-order oligonucleotide sequences.
- Tail sequences are now used when computing cross-dimer scores in `Junction.find_primer_pairs()`, so dimer scores reflect the actual primer that will be ordered rather than the bare binding region.

### Changed

- Config fields renamed in `SingleplexDesignParameters`: `five_prime_tail` (alias `5_primer_tail`) → `forward_tail`; `three_prime_tail` (alias `3_prime_tail`) → `reverse_tail`. Both tails are prepended at the 5′ end of their respective primer; the new names reflect which primer each tail belongs to rather than which end of the amplicon.
- JSON preset config files (`designer_default_config.json`, `designer_lenient_config.json`): keys renamed from `5_primer_tail` / `3_prime_tail` to `forward_tail` / `reverse_tail`.

### Implementation note

`Primer.seq` is never mutated. Any `N` bases in a tail (e.g. UMI placeholders) are replaced with `A` for dimer scoring only — the nearest-neighbour thermodynamic tables used by the aligner do not handle `N`. Output CSVs receive the original tail with `N`s intact.

## [0.3.3] - 18-02-2026

### Removed

- `src/plexus/config.py`: Removed config_dict parameter and its branch from load_config(), updated docstring.
- `src/plexus/pipeline.py`: Removed config_dict parameter from run_pipeline() signature and docstring, removed from load_config() call site, and removed the now-unused from typing import Any import.
- `src/plexus/designer/multiplexpanel.py`: Removed config_dict parameter from MultiplexPanel.load_config() and its load_config() call.
- `src/plexus/orchestrator.py`: Removed config_dict from the **pipeline_kwargs docstring.
- `tests/test_config.py`: Deleted test_load_from_dict and test_priority_dict_over_file; the three DesignerConfig.from_dict() tests are retained.

## [0.3.2] - 18-02-2026

### Added

- `snp_check_parameters` section added to both JSON preset config files (`designer_default_config.json`, `designer_lenient_config.json`), making them the canonical source of truth for SNP check defaults — consistent with every other config section.
- `_check_blast_tools()` in `blast_runner.py` verifies that `blastn`, `makeblastdb`, and `blast_formatter` are on `$PATH` before any BLAST run, raising a clear `RuntimeError` with install instructions if any are missing.

### Fixed

- `snp_config.snp_strict` from the JSON config file is now correctly consulted in `pipeline.py` alongside the `--snp-strict` CLI flag (`if snp_strict or snp_config.snp_strict`). Previously the field was loaded but never read, so setting `"snp_strict": true` in a config file had no effect.

### Changed

- `run_snp_check()` — `snp_penalty_weight` is now a required keyword-only parameter (was `= 5.0`), eliminating the duplicate default that shadowed the config value.
- `filter_snp_pairs()` — return type changed from `int` to `tuple[int, list[str]]`; the second element lists junction names where the fallback (all pairs overlapped SNPs) was triggered. `pipeline.py` and all tests updated accordingly.
- BLAST runner subprocess calls converted from shell-string form (`shell=True`) to list-form, improving safety and cross-platform correctness.
- Added blast and bcftools to CI workflow.
- Updated project documentation and README.

### Removed

- 5 stale TODO comments with no actionable value (`blast_runner.py`, `thal.py`, `archive.py`, `multiplexpanel.py`, `design.py`).
- Dead `save_primer_designs()` stub method on `MultiplexPanel` (was never called by the pipeline).

## [0.3.1] - 18-02-2026

### Added

- SNP strict mode: `--snp-strict` CLI flag and `snp_strict` config parameter that discards primer pairs overlapping any SNP above the AF threshold before multiplex optimization.
- `filter_snp_pairs()` function in `checker.py` that removes SNP-overlapping pairs per junction, with a fallback to keep the least-affected pair when all pairs have SNPs.
- `snp_strict` field on `SnpCheckParameters` (default `False`).
- 5 new tests for the SNP strict filtering logic.

## [0.3.0] - 18-02-2026

### Added

- Simulated Annealing (`SimulatedAnnealing`) multiplex selector: seeds from greedy solutions, escapes local optima via probabilistic acceptance of worse solutions during cooling.
- Depth-First Search with pruning (`DFS`) multiplex selector: iterative DFS with lower-bound pruning, optional greedy seeding, and a `max_nodes` safety valve. Can find provably optimal solutions for smaller panels.
- `--selector` / `-s` CLI option to choose the multiplex selection algorithm (`Greedy`, `Random`, `BruteForce`, `SimulatedAnnealing`, `DFS`). Defaults to `Greedy`.
- `selector` parameter on `run_pipeline()` for programmatic algorithm selection.
- Tests for both new selectors, including edge cases (single-candidate targets, node limits).

## [0.2.3] - 17-02-2026

### Removed

- Removed the experimental `primer3py_design_primers` function and the `"primer3py"` / `"primer3"` design method branches from `design.py`. Only the `"simsen"` design method remains.
- Removed dead code tied to the primer3py pathway: `get_primer_dict()` from `primer.py`, and `create_primer_dataframe()` / `convert_primer_data_to_tables()` from `utils.py`.
- Removed unused `import primer3`, `import json`, `import os`, and `import pandas` from affected modules.
- Added `CITATIONS.md`
- Changed license to GPL-2.0 or later to reflect package dependencies

## [0.2.2] - 17-02-2026

### Added

- End-to-end integration tests exercising the full pipeline: load → merge → design → SNP check → optimize → output.
- 8 new integration tests in `tests/test_integration.py` marked with `@pytest.mark.integration`.
- Test fixture files in `tests/data/` (~6 KB): small FASTA contigs, remapped gnomAD VCF, and junction CSV extracted from real reference data.
- Fixture creation script `scripts/create_test_fixtures.py` for reproducing test data from full reference files.
- Session-scoped shared fixtures in `tests/conftest.py`.
- `integration` pytest marker registered in `pyproject.toml` for selective test runs (`pytest -m integration`, `pytest -m "not integration"`).

## [0.2.1] - 17-02-2026

### Added

- `download-resources` CLI command for downloading and managing gnomAD VCF resources.
- `resources.py` module to handle VCF caching, environment variables (`PLEXUS_GNOMAD_VCF`), and atomic downloads.
- Progress bar visualization for large VCF downloads using `rich`.

### Changed

- Refactored `snp_data.py` to use local gnomAD VCFs instead of Ensembl API.
- Replaced Ensembl REST API dependency with direct VCF access for improved reliability and speed.

### Fixed

- **Critical**: Fixed a file handle leak in `checker.py` by using a context manager for `pysam.VariantFile`.
- Fixed incorrect end coordinate calculation in `snp_data.py` fallback logic (now uses `max(start, end) + padding`).
- Fixed linting errors and formatting.

### Removed

- `src/plexus/snpcheck/ensembl.py` and associated Ensembl API logic.

## [0.2.0] - 13-02-2026

### Added

- SNP overlap checking for primer binding sites using local tabix-indexed
  VCF/BCF via pysam.
- New CLI options: `--snp-vcf`, `--skip-snpcheck`, `--snp-af-threshold`.
- `SnpCheckParameters` configuration model with `af_threshold` and
  `snp_penalty_weight` settings.
- SNP columns (`SNP_Count`, `SNP_Penalty`, `Forward_SNP_Count`,
  `Reverse_SNP_Count`) in output CSVs.
- `snp_count` field on `Primer`; `snp_count` and `snp_penalty` fields on
  `PrimerPair`.
- 19 new tests for the snpcheck module.
- README badges for CI, Python version, and license.

### Changed

- Updated Dockerfile to Python 3.13 with optimized layer caching.
- Capped Python version to `<3.14` due to pysam compatibility.

### Removed

- Legacy `snpcheck` module (multiply/bedtools dependencies).

## [0.1.1] - 13-02-2026

### Added

- Changes project name to `plexus`
- Multi-panel support: input CSV can now include an optional `Panel` column to
  define multiple panels in a single file. Each panel is designed and optimized
  independently with output saved to `output_dir/<panel_id>/`.
- New `orchestrator.py` module with `run_multi_panel()` entry point that detects
  panels, dispatches per-panel pipeline runs, and aggregates results.
- `MultiPanelResult` dataclass for aggregated results across panels.
- `--parallel` CLI flag to run panels concurrently via `ProcessPoolExecutor`.
- `--max-workers` CLI flag to control parallel worker count.
- `multi_panel_summary.json` output file for multi-panel runs.
- Tests for orchestrator (20 new tests) and `JunctionInput` Panel field.

### Changed

- CLI now routes through `run_multi_panel()` which delegates to `run_pipeline()`
  per panel. Single-panel CSVs (no Panel column) behave identically to before.
- Marked `primer3py` design method as experimental with warnings; it does not yet
  produce `PrimerPair` objects for downstream multiplex optimization.
- Marked `primer3` (local binary) design method as not implemented with a clear
  error message pointing users to the `simsen` method.
- `PipelineResult` is now constructed after panel creation succeeds, eliminating
  the `panel=None` / `type: ignore` workaround. The result always contains a
  valid `MultiplexPanel`.

### Removed

- Removed `run_designer()` convenience function and associated hardcoded paths
  from `designer/design.py`. Use `run_pipeline()` or the CLI instead.

## [0.1.0] - 19-01-2026

Initial release.
