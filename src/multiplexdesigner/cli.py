from pathlib import Path
from typing import Annotated

import typer
from loguru import logger
from rich.console import Console

from multiplexdesigner.blast.specificity import run_specificity_check
from multiplexdesigner.designer.design import run_designer
from multiplexdesigner.version import __version__

# Create main app and subcommand apps
app = typer.Typer(
    name="multiplexdesigner",
    help="Pipeline for generating multiplex PCR primers.",
    no_args_is_help=True,
    rich_markup_mode="rich",
)

DEFAULT_OUTPUT_DIR = Path("./output")

console = Console()


def version_callback(value: bool) -> None:
    """Show version and exit."""
    if value:
        console.print(
            f"[bold green]Multiplex Primer Designer[/bold green] version {__version__}"
        )
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
    """UMI Error Correct - Pipeline for analyzing barcoded amplicon sequencing data."""
    log_level = "DEBUG" if verbose else "INFO"
    logger.add("multiplexdesigner.log", level=log_level)


@app.command()
def design(
    input_file: Annotated[
        Path,
        typer.Option(
            "--input",
            "-i",
            help="Path to the input CSV file containing junctions.",
        ),
    ],
    fasta_file: Annotated[
        Path,
        typer.Option("--fasta", "-f", help="Path to the reference genome FASTA file."),
    ],
    output_dir: Annotated[
        Path,
        typer.Option("--output", "-o", help="Directory to save output files."),
    ] = DEFAULT_OUTPUT_DIR,
    skip_blast: Annotated[
        bool, typer.Option("--skip-blast", help="Skip the BLAST specificity check.")
    ] = False,
):
    """
    Run the multiplex primer design pipeline.
    """

    # Create output directory
    logger.info(f"Creating output directory: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)

    console = Console()
    console.print("[bold green]Starting design pipeline...[/bold green]")
    console.print(f"Input: {input_file}")
    console.print(f"Genome: {fasta_file}")
    console.print(f"Output: {output_dir}")

    # 1. Run Designer
    try:
        panel = run_designer(
            design_input_file=str(input_file), fasta_file=str(fasta_file)
        )
    except Exception as e:
        console.print(f"[bold red]Error during design phase: {e}[/bold red]")
        raise typer.Exit(code=1) from e

    # 2. Save Intermediate Results
    console.print("[bold blue]Saving candidate designs...[/bold blue]")
    try:
        panel.save_candidate_pairs_to_csv(str(output_dir / "candidate_pairs.csv"))
    except Exception as e:
        console.print(f"[yellow]Warning: Could not save candidate pairs: {e}[/yellow]")

    # 3. Run Specificity Check (BLAST)
    if not skip_blast:
        console.print("[bold blue]Running specificity check...[/bold blue]")
        try:
            blast_dir = output_dir / "blast"
            run_specificity_check(panel, str(blast_dir), str(fasta_file))
        except Exception as e:
            console.print(f"[bold red]Error during specificity check: {e}[/bold red]")
            # We don't exit here, as partial results might still be useful
    else:
        console.print("[yellow]Skipping BLAST specificity check.[/yellow]")

    console.print("[bold green]Pipeline completed successfully![/bold green]")


if __name__ == "__main__":
    app()
