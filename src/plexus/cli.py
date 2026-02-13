# ================================================================================
# Command-line interface for multiplex primer designer
#
# Thin wrapper around the pipeline module.
#
# Author: Stefan Filges (stefan@simsendiagnostics.com)
# Copyright (c) 2025 Simsen Diagnostics AB
# ================================================================================

from pathlib import Path
from typing import Annotated

import typer
from loguru import logger
from rich.console import Console

from plexus.version import __version__

app = typer.Typer(
    name="plexus",
    help="Design multiplex PCR primer panels.",
    no_args_is_help=True,
    rich_markup_mode="rich",
)

console = Console()


def version_callback(value: bool) -> None:
    """Show version and exit."""
    if value:
        console.print(f"[bold green]Plexus[/bold green] version {__version__}")
        raise typer.Exit()


@app.callback()
def main(
    version: Annotated[
        bool,
        typer.Option(
            "--version",
            "-v",
            callback=version_callback,
            is_eager=True,
            help="Show version and exit.",
        ),
    ] = False,
    verbose: Annotated[
        bool,
        typer.Option("--verbose", "-V", help="Enable verbose output (DEBUG level)."),
    ] = False,
) -> None:
    """Multiplex Primer Designer - Design optimized multiplex PCR primer panels."""
    log_level = "DEBUG" if verbose else "INFO"
    logger.add("plexus.log", level=log_level)


@app.command()
def run(
    input_file: Annotated[
        Path,
        typer.Option(
            "--input",
            "-i",
            help="Path to CSV file containing junction coordinates.",
        ),
    ],
    fasta_file: Annotated[
        Path,
        typer.Option(
            "--fasta",
            "-f",
            help="Path to reference genome FASTA file.",
        ),
    ],
    output_dir: Annotated[
        Path,
        typer.Option(
            "--output",
            "-o",
            help="Directory for output files.",
        ),
    ] = Path("./output"),
    panel_name: Annotated[
        str,
        typer.Option(
            "--name",
            "-n",
            help="Name for the primer panel.",
        ),
    ] = "multiplex_panel",
    genome: Annotated[
        str,
        typer.Option(
            "--genome",
            "-g",
            help="Reference genome name.",
        ),
    ] = "hg38",
    preset: Annotated[
        str,
        typer.Option(
            "--preset",
            "-p",
            help="Configuration preset (default or lenient).",
        ),
    ] = "default",
    config_file: Annotated[
        Path | None,
        typer.Option(
            "--config",
            "-c",
            help="Path to custom configuration JSON file.",
        ),
    ] = None,
    skip_blast: Annotated[
        bool,
        typer.Option(
            "--skip-blast",
            help="Skip the BLAST specificity check.",
        ),
    ] = False,
    padding: Annotated[
        int,
        typer.Option(
            "--padding",
            help="Bases to extract around each junction.",
        ),
    ] = 200,
    parallel: Annotated[
        bool,
        typer.Option(
            "--parallel",
            help="Run panels in parallel (only applies to multi-panel input).",
        ),
    ] = False,
    max_workers: Annotated[
        int | None,
        typer.Option(
            "--max-workers",
            help="Max parallel workers for multi-panel mode.",
        ),
    ] = None,
) -> None:
    """
    Run the complete multiplex primer design pipeline.

    If the input CSV contains a 'Panel' column, junctions are grouped by
    panel and each panel is designed independently. Results are saved to
    output_dir/<panel_id>/. Use --parallel to run panels concurrently.

    Example:
        plexus run -i junctions.csv -f genome.fa -o results/
    """
    from plexus.orchestrator import run_multi_panel
    from plexus.pipeline import MultiPanelResult

    console.print("[bold green]Plexus[/bold green]")
    console.print(f"  Input:  {input_file}")
    console.print(f"  FASTA:  {fasta_file}")
    console.print(f"  Output: {output_dir}")
    console.print()

    try:
        result = run_multi_panel(
            input_file=input_file,
            fasta_file=fasta_file,
            output_dir=output_dir,
            parallel=parallel,
            max_workers=max_workers,
            panel_name=panel_name,
            genome=genome,
            preset=preset,
            config_path=config_file,
            run_blast=not skip_blast,
            padding=padding,
        )

        if isinstance(result, MultiPanelResult):
            if result.success:
                console.print()
                console.print(
                    "[bold green]All panels completed successfully![/bold green]"
                )
                console.print(f"  Panels:         {len(result.panel_ids)}")
                for pid in result.panel_ids:
                    pr = result.panel_results[pid]
                    console.print(
                        f"    {pid}: {pr.num_junctions} junctions, "
                        f"{len(pr.selected_pairs)} selected pairs"
                    )
                console.print(f"  Output:         {result.output_dir}")
            else:
                console.print()
                console.print("[bold yellow]Some panels had errors:[/bold yellow]")
                for pid in result.failed_panels:
                    pr = result.panel_results[pid]
                    for error in pr.errors:
                        console.print(f"  [yellow]• {pid}: {error}[/yellow]")
        else:
            if result.success:
                console.print()
                console.print(
                    "[bold green]Pipeline completed successfully![/bold green]"
                )
                console.print(f"  Junctions:    {result.num_junctions}")
                console.print(f"  Primer pairs: {result.num_primer_pairs}")
                console.print(f"  Output:       {result.output_dir}")
            else:
                console.print()
                console.print(
                    "[bold yellow]Pipeline completed with warnings:[/bold yellow]"
                )
                for error in result.errors:
                    console.print(f"  [yellow]• {error}[/yellow]")

    except FileNotFoundError as e:
        console.print(f"[bold red]Error: {e}[/bold red]")
        raise typer.Exit(code=1) from e
    except Exception as e:
        console.print(f"[bold red]Pipeline failed: {e}[/bold red]")
        raise typer.Exit(code=1) from e


if __name__ == "__main__":
    app()
