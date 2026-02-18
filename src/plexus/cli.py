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
        Path | None,
        typer.Option(
            "--fasta",
            "-f",
            help="Path to reference genome FASTA file. If omitted, uses registered genome from `plexus init`.",
        ),
    ] = None,
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
    snp_vcf: Annotated[
        Path | None,
        typer.Option(
            "--snp-vcf",
            help="Path to tabix-indexed VCF for SNP checking. If omitted, uses bundled gnomAD AF-only VCF.",
        ),
    ] = None,
    skip_snpcheck: Annotated[
        bool,
        typer.Option(
            "--skip-snpcheck",
            help="Skip the SNP overlap check.",
        ),
    ] = False,
    snp_af_threshold: Annotated[
        float | None,
        typer.Option(
            "--snp-af-threshold",
            help="Minimum allele frequency for SNP flagging (default: 0.01).",
        ),
    ] = None,
    snp_strict: Annotated[
        bool,
        typer.Option(
            "--snp-strict",
            help="Discard primer pairs that overlap SNPs.",
        ),
    ] = False,
    selector: Annotated[
        str,
        typer.Option(
            "--selector",
            "-s",
            help="Multiplex selector algorithm: Greedy, Random, BruteForce, SimulatedAnnealing, or DFS.",
        ),
    ] = "Greedy",
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
    from plexus.resources import get_registered_fasta

    if fasta_file is None:
        fasta_file = get_registered_fasta(genome)
        if fasta_file is None:
            typer.echo(
                f"--fasta is required: genome '{genome}' is not initialized.\n"
                f"Run: plexus init --genome {genome}",
                err=True,
            )
            raise typer.Exit(code=1)

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
            snp_vcf=snp_vcf,
            skip_snpcheck=skip_snpcheck,
            snp_af_threshold=snp_af_threshold,
            snp_strict=snp_strict,
            selector=selector,
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


@app.command()
def download_resources(
    force: Annotated[
        bool,
        typer.Option(
            "--force",
            help="Re-download even if files already exist.",
        ),
    ] = False,
) -> None:
    """Download the gnomAD AF-only VCF for SNP checking."""
    from plexus.snpcheck.resources import download_gnomad_vcf, get_cache_dir

    console.print("[bold green]Plexus[/bold green] — downloading SNP resources")
    console.print(f"  Cache directory: {get_cache_dir()}")
    console.print()

    try:
        vcf_path = download_gnomad_vcf(force=force)
        console.print()
        console.print(f"[bold green]Done![/bold green] VCF saved to {vcf_path}")
    except Exception as e:
        console.print(f"[bold red]Download failed: {e}[/bold red]")
        raise typer.Exit(code=1) from e


@app.command()
def status() -> None:
    """Show Plexus version and resource status."""
    from plexus.resources import GENOME_PRESETS, genome_status
    from plexus.snpcheck.resources import resource_status_message

    console.print(f"[bold green]Plexus[/bold green] version {__version__}")
    console.print(f"  {resource_status_message()}")
    console.print()
    console.print("[bold]Genome resources:[/bold]")
    for g in GENOME_PRESETS:
        s = genome_status(g)
        marks = {True: "[green]✓[/green]", False: "[red]✗[/red]"}
        console.print(
            f"  {g}:"
            f"  FASTA {marks[s['fasta']]}"
            f"  FAI {marks[s['fai']]}"
            f"  BLAST {marks[s['blast_db']]}"
            f"  gnomAD VCF {marks[s['snp_vcf']]}"
        )


@app.command()
def init(
    genome: Annotated[
        str,
        typer.Option(
            "--genome",
            "-g",
            help="Reference genome to initialize (e.g. hg38).",
        ),
    ] = "hg38",
    fasta: Annotated[
        Path | None,
        typer.Option(
            "--fasta",
            "-f",
            help="Use a local FASTA instead of downloading.",
        ),
    ] = None,
    snp_vcf: Annotated[
        Path | None,
        typer.Option(
            "--snp-vcf",
            help="Use a local gnomAD VCF instead of downloading.",
        ),
    ] = None,
    skip_blast: Annotated[
        bool,
        typer.Option("--skip-blast", help="Skip BLAST index creation."),
    ] = False,
    skip_snp: Annotated[
        bool,
        typer.Option("--skip-snp", help="Skip SNP VCF download/registration."),
    ] = False,
    force: Annotated[
        bool,
        typer.Option(
            "--force", help="Re-download and rebuild even if resources exist."
        ),
    ] = False,
) -> None:
    """Download and index reference resources for a genome.

    Provide --fasta / --snp-vcf to use local files instead of downloading.
    """
    from plexus.resources import GENOME_PRESETS, genome_status, init_genome

    if genome not in GENOME_PRESETS:
        supported = ", ".join(GENOME_PRESETS)
        typer.echo(
            f"Unknown genome '{genome}'. Supported genomes: {supported}",
            err=True,
        )
        raise typer.Exit(code=1)

    preset = GENOME_PRESETS[genome]
    console.print(
        f"[bold green]Plexus[/bold green] — initializing resources for [bold]{genome}[/bold]"
    )
    if fasta:
        console.print(f"  FASTA:   {fasta} (local)")
    else:
        console.print(f"  FASTA:   {preset['fasta_size_note']} — downloading from UCSC")
    if snp_vcf:
        console.print(f"  SNP VCF: {snp_vcf} (local)")
    elif not skip_snp:
        console.print("  SNP VCF: downloading gnomAD AF-only VCF from GATK bucket")
    console.print()

    try:
        init_genome(
            genome=genome,
            fasta=fasta,
            snp_vcf=snp_vcf,
            force=force,
            skip_blast=skip_blast,
            skip_snp=skip_snp,
        )
    except ValueError as e:
        typer.echo(f"Error: {e}", err=True)
        raise typer.Exit(code=1) from e
    except FileNotFoundError as e:
        console.print(f"[bold red]Error: {e}[/bold red]")
        raise typer.Exit(code=1) from e
    except RuntimeError as e:
        console.print(f"[bold red]Error: {e}[/bold red]")
        raise typer.Exit(code=1) from e
    except Exception as e:
        console.print(f"[bold red]Init failed: {e}[/bold red]")
        raise typer.Exit(code=1) from e

    console.print()
    s = genome_status(genome)
    console.print(
        f"[bold green]Done![/bold green] Resources for [bold]{genome}[/bold]:"
    )
    console.print(f"  FASTA:   {'ready' if s['fasta'] else 'NOT ready'}")
    console.print(f"  FAI:     {'ready' if s['fai'] else 'NOT ready'}")
    console.print(
        f"  BLAST:   {'ready' if s['blast_db'] else ('NOT ready' if not skip_blast else 'skipped')}"
    )
    console.print(
        f"  gnomAD:  {'ready' if s['snp_vcf'] else ('NOT ready' if not skip_snp else 'skipped')}"
    )


if __name__ == "__main__":
    app()
