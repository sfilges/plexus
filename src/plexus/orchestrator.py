# ================================================================================
# Multi-panel orchestration layer
#
# Reads input CSV, groups junctions by optional Panel column,
# writes per-panel temporary CSVs, and runs run_pipeline() for each.
#
# Author: Stefan Filges (stefan.filges@pm.me)
# Copyright (c) 2026 Stefan Filges
# ================================================================================

from __future__ import annotations

import json
import shutil
import sys
import threading
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import Manager
from pathlib import Path
from typing import Any

import pandas as pd
from loguru import logger
from rich.progress import (
    Progress,
    SpinnerColumn,
    TextColumn,
    TimeElapsedColumn,
)

from plexus.logging import restore_console_logging, suppress_console_logging
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
    _status_dict=None,
    **pipeline_kwargs: Any,
) -> tuple[str, PipelineResult]:
    """
    Run run_pipeline for a single panel.

    This is a top-level function (not a method or lambda) so it can be
    pickled by ProcessPoolExecutor.

    Parameters
    ----------
    _status_dict : dict-like or None
        Shared multiprocessing dict for reporting step progress back to
        the parent process.  Keyed by panel_id.

    Returns
    -------
    tuple[str, PipelineResult]
        The panel_id and its PipelineResult.
    """
    panel_output_dir = output_dir / panel_id
    panel_name = pipeline_kwargs.pop("panel_name", panel_id)

    # Suppress child-process console logging early (before any logger calls)
    # when the parent is showing a progress bar.
    if pipeline_kwargs.get("quiet"):
        suppress_console_logging()

    # Build a step_callback that writes to the shared dict
    step_callback = None
    if _status_dict is not None:

        def step_callback(label):
            _status_dict[panel_id] = label

        _status_dict[panel_id] = "Starting…"

    logger.info(f"Starting pipeline for panel '{panel_id}'...")

    result = run_pipeline(
        input_file=panel_csv,
        fasta_file=fasta_file,
        output_dir=panel_output_dir,
        panel_name=panel_name,
        step_callback=step_callback,
        **pipeline_kwargs,
    )

    if _status_dict is not None:
        _status_dict[panel_id] = "Done"

    return panel_id, result


def run_multi_panel(
    input_file: str | Path,
    fasta_file: str | Path,
    output_dir: str | Path = "./output",
    parallel: bool = False,
    max_workers: int | None = None,
    show_progress: bool = False,
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
        (genome, preset, config_path, design_method, run_blast, padding).

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
            show_progress=show_progress,
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

    # Panel-level progress bar for parallel mode (no per-step bars
    # because subprocesses would fight over stderr).
    _use_panel_bar = parallel and show_progress and sys.stderr.isatty()

    # Suppress parent console logging early — before any parallel INFO
    # messages — so only the progress bar and explicit console prints appear.
    if _use_panel_bar:
        suppress_console_logging()

    try:
        if parallel:
            workers = max_workers or len(panels)
            logger.info(f"Running {len(panels)} panels in parallel (workers={workers})")

            # Shared dict for child processes to report current step
            mgr = Manager() if _use_panel_bar else None
            status_dict = mgr.dict() if mgr else None

            with ProcessPoolExecutor(max_workers=workers) as executor:
                futures = {}
                for panel_id, csv_path in panel_csvs.items():
                    kw = dict(pipeline_kwargs)
                    kw["panel_name"] = panel_id
                    if _use_panel_bar:
                        kw["quiet"] = True
                    future = executor.submit(
                        _run_single_panel,
                        panel_id=panel_id,
                        panel_csv=csv_path,
                        fasta_file=fasta_file,
                        output_dir=output_dir,
                        _status_dict=status_dict,
                        **kw,
                    )
                    futures[future] = panel_id

                if _use_panel_bar:
                    progress = Progress(
                        SpinnerColumn(),
                        TextColumn("{task.description}"),
                        TimeElapsedColumn(),
                    )

                    # Route WARNING+ messages through the progress console
                    def _warning_sink(message):
                        progress.console.print(message, end="")

                    warning_handler_id = logger.add(
                        _warning_sink, level="WARNING", format="{message}"
                    )

                    # One task row per panel
                    panel_tasks = {}
                    for pid in panels:
                        panel_tasks[pid] = progress.add_task(
                            f"[bold]{pid}[/bold]: waiting…", total=None
                        )

                    # Background thread polls shared dict to update descriptions
                    _stop_poll = threading.Event()

                    def _poll_status():
                        while not _stop_poll.is_set():
                            for pid, tid in panel_tasks.items():
                                step = status_dict.get(pid)
                                if step and not progress.tasks[tid].finished:
                                    progress.update(
                                        tid,
                                        description=f"[bold]{pid}[/bold]: {step}",
                                    )
                            time.sleep(0.3)

                    poll_thread = threading.Thread(target=_poll_status, daemon=True)

                    try:
                        with progress:
                            poll_thread.start()
                            for future in as_completed(futures):
                                pid = futures[future]
                                _, result = future.result()
                                results[pid] = result
                                status = (
                                    "[green]done[/green]"
                                    if result.success
                                    else "[yellow]done (warnings)[/yellow]"
                                )
                                progress.update(
                                    panel_tasks[pid],
                                    description=f"[bold]{pid}[/bold]: {status}",
                                    total=1,
                                    completed=1,
                                )
                                logger.info(f"Panel '{pid}' completed.")
                    finally:
                        _stop_poll.set()
                        poll_thread.join(timeout=1)
                        logger.remove(warning_handler_id)
                        restore_console_logging()
                else:
                    for future in as_completed(futures):
                        pid = futures[future]
                        _, result = future.result()
                        results[pid] = result
                        logger.info(f"Panel '{pid}' completed.")

            # Clean up manager
            if mgr:
                mgr.shutdown()
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
                    show_progress=show_progress,
                    **kw,
                )
                results[panel_id] = result
    finally:
        # Ensure console logging is restored if we suppressed it
        if _use_panel_bar:
            restore_console_logging()
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
