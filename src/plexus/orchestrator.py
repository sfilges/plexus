# ================================================================================
# Multi-panel orchestration layer
#
# Reads input CSV, groups junctions by optional Panel column,
# writes per-panel temporary CSVs, and runs run_pipeline() for each.
#
# Author: Stefan Filges (stefan@simsendiagnostics.com)
# Copyright (c) 2025 Simsen Diagnostics AB
# ================================================================================

from __future__ import annotations

import json
import shutil
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import pandas as pd
from loguru import logger

from plexus.pipeline import MultiPanelResult, PipelineResult, run_pipeline

DEFAULT_PANEL_ID = "default"


def _detect_panels(input_file: Path) -> dict[str, pd.DataFrame]:
    """
    Read the input CSV and group rows by the Panel column.

    If the Panel column is absent, all rows are assigned to a single
    panel with id DEFAULT_PANEL_ID.

    Parameters
    ----------
    input_file : Path
        Path to the input CSV file.

    Returns
    -------
    dict[str, DataFrame]
        Mapping of panel_id to DataFrame of junction rows.

    Raises
    ------
    ValueError
        If the Panel column contains empty/NaN values.
    """
    df = pd.read_csv(input_file)

    if "Panel" not in df.columns:
        return {DEFAULT_PANEL_ID: df}

    if df["Panel"].isna().any():
        raise ValueError(
            "Panel column contains empty values. "
            "Either omit the Panel column entirely or fill all rows."
        )

    grouped = {str(panel_id): group_df for panel_id, group_df in df.groupby("Panel")}
    logger.info(f"Detected {len(grouped)} panels: {list(grouped.keys())}")
    return grouped


def _write_panel_csv(df: pd.DataFrame, tmp_dir: Path, panel_id: str) -> Path:
    """
    Write a per-panel junction CSV (without the Panel column) to a temp directory.

    Parameters
    ----------
    df : DataFrame
        Junction rows for this panel.
    tmp_dir : Path
        Directory to write temporary CSV files.
    panel_id : str
        Panel identifier used in the filename.

    Returns
    -------
    Path
        Path to the written CSV file.
    """
    cols_to_write = [c for c in df.columns if c != "Panel"]
    out_path = tmp_dir / f"{panel_id}_junctions.csv"
    df[cols_to_write].to_csv(out_path, index=False)
    return out_path


def _run_single_panel(
    panel_id: str,
    panel_csv: Path,
    fasta_file: Path,
    output_dir: Path,
    **pipeline_kwargs: Any,
) -> tuple[str, PipelineResult]:
    """
    Run run_pipeline for a single panel.

    This is a top-level function (not a method or lambda) so it can be
    pickled by ProcessPoolExecutor.

    Returns
    -------
    tuple[str, PipelineResult]
        The panel_id and its PipelineResult.
    """
    panel_output_dir = output_dir / panel_id
    panel_name = pipeline_kwargs.pop("panel_name", panel_id)

    logger.info(f"Starting pipeline for panel '{panel_id}'...")

    result = run_pipeline(
        input_file=panel_csv,
        fasta_file=fasta_file,
        output_dir=panel_output_dir,
        panel_name=panel_name,
        **pipeline_kwargs,
    )
    return panel_id, result


def run_multi_panel(
    input_file: str | Path,
    fasta_file: str | Path,
    output_dir: str | Path = "./output",
    parallel: bool = False,
    max_workers: int | None = None,
    **pipeline_kwargs: Any,
) -> MultiPanelResult | PipelineResult:
    """
    Orchestrate single- or multi-panel pipeline execution.

    If the input CSV has no Panel column, delegates directly to
    run_pipeline (single panel, flat output — full backward compat).

    If a Panel column exists, groups junctions, creates per-panel
    temp CSVs, runs run_pipeline per panel, and returns a
    MultiPanelResult.

    Parameters
    ----------
    input_file : str | Path
        Path to CSV with junction coordinates and optional Panel column.
    fasta_file : str | Path
        Path to reference genome FASTA.
    output_dir : str | Path
        Root output directory.
    parallel : bool
        If True, run panels in parallel using ProcessPoolExecutor.
    max_workers : int | None
        Max parallel workers. None defaults to number of panels.
    **pipeline_kwargs
        All other arguments passed through to run_pipeline
        (genome, preset, config_path, config_dict, design_method,
        run_blast, padding).

    Returns
    -------
    MultiPanelResult | PipelineResult
        MultiPanelResult if multi-panel, PipelineResult if single panel.
    """
    input_file = Path(input_file)
    fasta_file = Path(fasta_file)
    output_dir = Path(output_dir)

    panels = _detect_panels(input_file)

    # --- Single panel, no Panel column: full backward compatibility ---
    if len(panels) == 1 and DEFAULT_PANEL_ID in panels:
        logger.info("Single panel detected, running standard pipeline.")
        return run_pipeline(
            input_file=input_file,
            fasta_file=fasta_file,
            output_dir=output_dir,
            **pipeline_kwargs,
        )

    # --- Multi-panel path ---
    logger.info(f"Multi-panel mode: {len(panels)} panels detected.")
    output_dir.mkdir(parents=True, exist_ok=True)

    # Write per-panel CSVs to a temp directory
    tmp_dir = output_dir / ".tmp_panel_csvs"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    panel_csvs: dict[str, Path] = {}
    for panel_id, df in panels.items():
        panel_csvs[panel_id] = _write_panel_csv(df, tmp_dir, panel_id)

    results: dict[str, PipelineResult] = {}

    try:
        if parallel:
            workers = max_workers or len(panels)
            logger.info(f"Running {len(panels)} panels in parallel (workers={workers})")
            with ProcessPoolExecutor(max_workers=workers) as executor:
                futures = {}
                for panel_id, csv_path in panel_csvs.items():
                    kw = dict(pipeline_kwargs)
                    kw["panel_name"] = panel_id
                    future = executor.submit(
                        _run_single_panel,
                        panel_id=panel_id,
                        panel_csv=csv_path,
                        fasta_file=fasta_file,
                        output_dir=output_dir,
                        **kw,
                    )
                    futures[future] = panel_id

                for future in as_completed(futures):
                    pid = futures[future]
                    _, result = future.result()
                    results[pid] = result
                    logger.info(f"Panel '{pid}' completed.")
        else:
            logger.info(f"Running {len(panels)} panels sequentially.")
            for panel_id, csv_path in panel_csvs.items():
                kw = dict(pipeline_kwargs)
                kw["panel_name"] = panel_id
                _, result = _run_single_panel(
                    panel_id=panel_id,
                    panel_csv=csv_path,
                    fasta_file=fasta_file,
                    output_dir=output_dir,
                    **kw,
                )
                results[panel_id] = result
    finally:
        # Clean up temp CSVs
        shutil.rmtree(tmp_dir, ignore_errors=True)

    # Build aggregated result
    multi_result = MultiPanelResult(
        panel_results=results,
        output_dir=output_dir,
        panel_ids=list(panels.keys()),
    )

    # Save aggregated summary
    summary_path = output_dir / "multi_panel_summary.json"
    with open(summary_path, "w") as f:
        json.dump(multi_result.summary_dict(), f, indent=2, default=str)
    logger.info(f"Saved multi-panel summary to {summary_path}")

    return multi_result
