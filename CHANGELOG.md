# Changelog

All notable changes to plexus will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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

## [0.1.0] - 13-02-2026

Initial release.
