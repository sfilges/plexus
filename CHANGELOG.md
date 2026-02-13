# Changelog

All notable changes to plexus will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.2] - 13-02-2026

### Changed

- Updated Dockerfile to Python 3.13 with optimized layer caching.
- Capped Python version to `<3.14` due to pysam compatibility.

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
