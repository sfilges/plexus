# ================================================================================
# Command-line interface for multiplex primer designer
#
# Thin wrapper around the pipeline module.
#
# Author: Stefan Filges (stefan.filges@pm.me)
# Copyright (c) 2026 Stefan Filges
# ================================================================================

from pathlib import Path
from typing import Annotated

import typer
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
) -> None:
    """Multiplex Primer Designer - Design optimized multiplex PCR primer panels."""


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
    selector_seed: Annotated[
        int | None,
        typer.Option(
            "--selector-seed",
            help="Random seed for stochastic selectors (Greedy, Random, SimulatedAnnealing).",
        ),
    ] = None,
    debug: Annotated[
        bool,
        typer.Option(
            "--debug",
            help="Write DEBUG-level messages to the log file (default: INFO+ only).",
        ),
    ] = False,
    strict: Annotated[
        bool,
        typer.Option(
            "--strict",
            help="Verify resource checksums before running (always on in compliance mode).",
        ),
    ] = False,
    checksums: Annotated[
        Path | None,
        typer.Option(
            "--checksums",
            help="SHA-256 checksums file (sha256sum format). Enables stateless verification "
            "without a registry. Required in compliance mode if --fasta is provided explicitly.",
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
    from plexus.resources import (
        get_operational_mode,
        get_registered_fasta,
        verify_resource_checksums,
    )

    op_mode = get_operational_mode()

    if fasta_file is None:
        # Fall back to registry for both research and compliance modes.
        # In compliance mode the checksum verification below (should_verify) is
        # automatically enabled, so this path is equally auditable.
        fasta_file = get_registered_fasta(genome)
        if fasta_file is None:
            typer.echo(
                f"--fasta is required: genome '{genome}' is not initialized.\n"
                f"Run: plexus init --genome {genome}",
                err=True,
            )
            raise typer.Exit(code=1)

    # ── Checksum verification ────────────────────────────────────────────────
    fasta_sha256: str | None = None
    snp_vcf_sha256: str | None = None

    if checksums is not None:
        # Path 1 — stateless: parse checksums file and verify on-the-fly
        from plexus.resources import _compute_sha256, _parse_checksums_file

        try:
            expected = _parse_checksums_file(checksums)
        except (FileNotFoundError, ValueError) as e:
            console.print(f"[bold red]Error reading checksums file: {e}[/bold red]")
            raise typer.Exit(code=1) from e

        # Verify FASTA
        if fasta_file.name in expected:
            actual = _compute_sha256(fasta_file)
            if actual != expected[fasta_file.name]:
                console.print(
                    f"[bold red]FASTA checksum mismatch: {fasta_file.name}[/bold red]"
                )
                raise typer.Exit(code=1)
            fasta_sha256 = actual
            console.print(f"  [green]✓[/green] FASTA verified ({fasta_sha256[:12]}…)")
        elif op_mode == "compliance":
            console.print(
                f"[bold red]Error: no entry for {fasta_file.name} in checksums file[/bold red]"
            )
            raise typer.Exit(code=1)

        # Verify SNP VCF (if applicable)
        if snp_vcf is not None and not skip_snpcheck:
            snp_vcf_path = Path(snp_vcf)
            if snp_vcf_path.name in expected:
                actual_vcf = _compute_sha256(snp_vcf_path)
                if actual_vcf != expected[snp_vcf_path.name]:
                    console.print(
                        f"[bold red]SNP VCF checksum mismatch: {snp_vcf_path.name}[/bold red]"
                    )
                    raise typer.Exit(code=1)
                snp_vcf_sha256 = actual_vcf
                console.print(
                    f"  [green]✓[/green] SNP VCF verified ({snp_vcf_sha256[:12]}…)"
                )
            elif op_mode == "compliance":
                console.print(
                    f"[bold red]Error: no entry for {snp_vcf_path.name} in checksums file[/bold red]"
                )
                raise typer.Exit(code=1)

    else:
        # Path 2 — registry-based (existing behaviour, research mode or --strict)
        should_verify = strict or (op_mode == "compliance")
        if should_verify:
            check = verify_resource_checksums(genome)
            if check.get("fasta") is False:
                console.print(
                    "[bold red]Error: FASTA checksum mismatch. "
                    "File may have been modified or corrupted since registration.[/bold red]"
                )
                raise typer.Exit(code=1)
            if check.get("snp_vcf") is False and not skip_snpcheck:
                console.print(
                    "[bold red]Error: SNP VCF checksum mismatch. "
                    "File may have been modified or corrupted since registration.[/bold red]"
                )
                raise typer.Exit(code=1)
            if check.get("fasta") is None and op_mode == "compliance":
                console.print(
                    f"[yellow]Warning: compliance mode active but no checksums stored for "
                    f"'{genome}'. Re-run `plexus init` with `--checksums` to register them.[/yellow]"
                )

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
            selector_seed=selector_seed,
            debug=debug,
            fasta_sha256=fasta_sha256,
            snp_vcf_sha256=snp_vcf_sha256,
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
def status() -> None:
    """Show Plexus version and resource status."""
    from plexus.resources import GENOME_PRESETS, genome_status, get_operational_mode
    from plexus.utils.env import check_executable

    op_mode = get_operational_mode()

    console.print(f"[bold green]Plexus[/bold green] version {__version__}")
    console.print(f"  Mode: [bold]{op_mode}[/bold]")
    if op_mode == "compliance":
        console.print(
            "  [dim]Registry lookups enforce checksum verification in compliance mode.[/dim]"
        )
    console.print()

    # Tool status
    console.print("[bold]System dependencies:[/bold]")
    marks = {True: "[green]✓[/green]", False: "[red]✗[/red]"}
    for tool in ["blastn", "makeblastdb", "blast_formatter", "bcftools"]:
        console.print(f"  {tool:<15} {marks[check_executable(tool)]}")
    console.print()

    console.print("[bold]Genome resources:[/bold]")
    for g in GENOME_PRESETS:
        s = genome_status(g)
        marks = {True: "[green]✓[/green]", False: "[red]✗[/red]"}
        fasta_cksum = (
            f" sha256:{s['fasta_sha256'][:12]}…" if s.get("fasta_sha256") else ""
        )
        vcf_cksum = (
            f" sha256:{s['snp_vcf_sha256'][:12]}…" if s.get("snp_vcf_sha256") else ""
        )
        console.print(
            f"  {g}:"
            f"  FASTA {marks[s['fasta']]}{fasta_cksum}"
            f"  FAI {marks[s['fai']]}"
            f"  BLAST {marks[s['blast_db']]}"
            f"  gnomAD VCF {marks[s['snp_vcf']]}{vcf_cksum}"
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
            help="Path to a local reference FASTA file.",
        ),
    ] = None,
    snp_vcf: Annotated[
        Path | None,
        typer.Option(
            "--snp-vcf",
            help="Path to a local tabix-indexed gnomAD VCF file.",
        ),
    ] = None,
    skip_blast: Annotated[
        bool,
        typer.Option("--skip-blast", help="Skip BLAST index creation."),
    ] = False,
    skip_snp: Annotated[
        bool,
        typer.Option("--skip-snp", help="Skip SNP VCF registration."),
    ] = False,
    force: Annotated[
        bool,
        typer.Option(
            "--force", help="Rebuild indexes even if resources already exist."
        ),
    ] = False,
    mode: Annotated[
        str | None,
        typer.Option(
            "--mode",
            help="Set operational mode: 'research' (default) or 'compliance'. "
            "Stored globally in ~/.plexus/config.json.",
        ),
    ] = None,
    checksums: Annotated[
        Path | None,
        typer.Option(
            "--checksums",
            help="SHA-256 checksums file (sha256sum format) for verifying "
            "FASTA and VCF at registration time.",
        ),
    ] = None,
) -> None:
    """Register and index reference resources for a genome.

    Provide --fasta and --snp-vcf to register local files.
    Use --checksums to verify files against known-good hashes.
    """
    from plexus.resources import (
        GENOME_PRESETS,
        genome_status,
        get_operational_mode,
        init_genome,
    )

    if genome not in GENOME_PRESETS:
        supported = ", ".join(GENOME_PRESETS)
        typer.echo(
            f"Unknown genome '{genome}'. Supported genomes: {supported}",
            err=True,
        )
        raise typer.Exit(code=1)

    if mode is not None and mode not in ("research", "compliance"):
        typer.echo(
            f"Invalid mode '{mode}'. Must be 'research' or 'compliance'.",
            err=True,
        )
        raise typer.Exit(code=1)

    if fasta is None:
        typer.echo(
            "Error: --fasta is required. Provide a local FASTA file.",
            err=True,
        )
        raise typer.Exit(code=1)
    if snp_vcf is None and not skip_snp:
        typer.echo(
            "Error: --snp-vcf is required (or use --skip-snp to skip SNP checking).",
            err=True,
        )
        raise typer.Exit(code=1)

    console.print(
        f"[bold green]Plexus[/bold green] — initializing resources for [bold]{genome}[/bold]"
    )
    console.print(f"  FASTA:     {fasta}")
    if snp_vcf:
        console.print(f"  SNP VCF:   {snp_vcf}")
    if checksums:
        console.print(f"  Checksums: {checksums}")
    if mode:
        console.print(f"  Mode:      {mode}")
    console.print()

    try:
        init_genome(
            genome=genome,
            fasta=fasta,
            snp_vcf=snp_vcf,
            force=force,
            skip_blast=skip_blast,
            skip_snp=skip_snp,
            mode=mode,
            checksums=checksums,
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
    op_mode = get_operational_mode()
    console.print(
        f"[bold green]Done![/bold green] Resources for [bold]{genome}[/bold]:"
    )
    console.print(f"  Mode:    {op_mode}")
    console.print(f"  FASTA:   {'ready' if s['fasta'] else 'NOT ready'}")
    if s.get("fasta_sha256"):
        console.print(f"           sha256:{s['fasta_sha256'][:16]}…")
    console.print(f"  FAI:     {'ready' if s['fai'] else 'NOT ready'}")
    console.print(
        f"  BLAST:   {'ready' if s['blast_db'] else ('NOT ready' if not skip_blast else 'skipped')}"
    )
    console.print(
        f"  gnomAD:  {'ready' if s['snp_vcf'] else ('NOT ready' if not skip_snp else 'skipped')}"
    )
    if s.get("snp_vcf_sha256"):
        console.print(f"           sha256:{s['snp_vcf_sha256'][:16]}…")


@app.command()
def template(
    output_dir: Annotated[
        Path,
        typer.Option(
            "--output",
            "-o",
            help="Directory to write template files.",
        ),
    ] = Path("."),
) -> None:
    """Generate starter files (junctions.csv, designer_config.json) for a new design.

    This command creates template files in the specified directory to help you
    get started with a new multiplex primer design.
    """
    from plexus.config import DesignerConfig

    output_dir.mkdir(parents=True, exist_ok=True)

    # 1. junctions.csv
    csv_path = output_dir / "junctions.csv"
    if csv_path.exists():
        console.print(f"[yellow]Warning: {csv_path} already exists, skipping.[/yellow]")
    else:
        with csv_path.open("w") as f:
            f.write("Name,Chrom,Five_Prime_Coordinate,Three_Prime_Coordinate,Panel\n")
            f.write("BRAF_V600E,chr7,140753336,140753336,panel1\n")
            f.write("EGFR_L858R,chr7,55191822,55191822,panel1\n")
        console.print(f"  [green]✓[/green] Created {csv_path}")

    # 2. designer_config.json
    json_path = output_dir / "designer_config.json"
    if json_path.exists():
        console.print(
            f"[yellow]Warning: {json_path} already exists, skipping.[/yellow]"
        )
    else:
        config = DesignerConfig()
        config.to_json_file(json_path)
        console.print(f"  [green]✓[/green] Created {json_path}")

    console.print("\n[bold green]Templates ready![/bold green]")
    console.print(f"1. Edit [bold]{csv_path.name}[/bold] with your target coordinates.")
    console.print(
        f"2. (Optional) Edit [bold]{json_path.name}[/bold] to tune design parameters."
    )
    console.print(
        f"3. Run: plexus run -i {csv_path.name} -c {json_path.name} -f /path/to/genome.fa"
    )


@app.command(
    name="docker",
    context_settings={"allow_extra_args": True, "ignore_unknown_options": True},
)
def docker_run(
    ctx: typer.Context,
    input_file: Annotated[
        Path,
        typer.Option(
            "--input", "-i", help="Path to CSV file containing junction coordinates."
        ),
    ],
    fasta_file: Annotated[
        Path,
        typer.Option(
            "--fasta",
            "-f",
            help="Path to reference genome FASTA (and adjacent BLAST DB / .fai files).",
        ),
    ],
    tag: Annotated[
        str,
        typer.Option(
            "--tag",
            help="Plexus image version tag to use (e.g. 0.5.0). Defaults to the currently installed version.",
        ),
    ] = __version__,
    registry: Annotated[
        str,
        typer.Option("--registry", help="Docker registry prefix for the plexus image."),
    ] = "ghcr.io/sfilges/plexus",
    output_dir: Annotated[
        Path,
        typer.Option("--output", "-o", help="Host directory for output files."),
    ] = Path("./output"),
    snp_vcf: Annotated[
        Path | None,
        typer.Option(
            "--snp-vcf", help="Path to tabix-indexed VCF (and adjacent .tbi index)."
        ),
    ] = None,
    checksums: Annotated[
        Path | None,
        typer.Option(
            "--checksums",
            help="SHA-256 checksums file for stateless data verification.",
        ),
    ] = None,
    config_file: Annotated[
        Path | None,
        typer.Option("--config", "-c", help="Path to custom JSON config file."),
    ] = None,
    pull: Annotated[
        bool,
        typer.Option(
            "--pull/--no-pull", help="Pull the image even if it exists locally."
        ),
    ] = False,
) -> None:
    """
    Run the plexus compliance container for a specific version.

    Wraps 'docker run': mounts parent directories of all file arguments
    (so adjacent BLAST DB, .fai, and .tbi files are included), translates
    host paths to container paths, and streams output.

    Additional 'plexus run' options (--selector, --preset, --skip-blast,
    --snp-strict, --genome, --padding, etc.) can be appended and are passed
    through unchanged.

    Example:
        plexus docker --tag 0.5.0 \\
            --fasta /data/hg38.fa \\
            --snp-vcf /data/gnomad.vcf.gz \\
            --checksums /data/checksums.sha256 \\
            --input /data/junctions.csv \\
            --output /data/results/ \\
            --skip-blast
    """
    import shutil
    import subprocess
    import sys

    # 1. Check docker is available
    if not shutil.which("docker"):
        console.print(
            "[bold red]Error: docker not found on PATH. "
            "Install Docker: https://docs.docker.com/get-docker/[/bold red]"
        )
        raise typer.Exit(code=1)

    image = f"{registry}:{tag}"

    console.print("[bold green]Plexus Docker Runner[/bold green]")
    console.print(f"  Image:  {image}")
    console.print(f"  Input:  {input_file}")
    console.print(f"  FASTA:  {fasta_file}")
    console.print(f"  Output: {output_dir}")
    console.print()

    # 2. Pull image if needed
    needs_pull = pull
    if not pull:
        inspect = subprocess.run(
            ["docker", "image", "inspect", image], capture_output=True
        )
        if inspect.returncode != 0:
            console.print(f"  Image not found locally — pulling {image} ...")
            needs_pull = True
    if needs_pull:
        result = subprocess.run(["docker", "pull", image])
        if result.returncode != 0:
            console.print(f"[bold red]Error: failed to pull {image}[/bold red]")
            raise typer.Exit(code=1)
    console.print("  [green]✓[/green] Image ready")

    # 3. Collect file args (excluding output)
    file_args: dict[str, Path | None] = {
        "input": input_file.resolve(),
        "fasta": fasta_file.resolve(),
        "snp_vcf": snp_vcf.resolve() if snp_vcf else None,
        "checksums": checksums.resolve() if checksums else None,
        "config": config_file.resolve() if config_file else None,
    }

    # 4. Build volume mounts: unique parent dirs → /mnt/vol0, /mnt/vol1, ...
    dir_map: dict[str, str] = {}
    counter = 0
    for path in file_args.values():
        if path is None:
            continue
        parent = str(path.parent)
        if parent not in dir_map:
            dir_map[parent] = f"/mnt/vol{counter}"
            counter += 1

    fasta_parent = str(file_args["fasta"].parent)

    volume_flags: list[str] = []
    for host_dir, mount_pt in dir_map.items():
        mode = "rw" if host_dir == fasta_parent else "ro"
        volume_flags += ["-v", f"{host_dir}:{mount_pt}:{mode}"]

    # Output dir — writable
    abs_output = output_dir.resolve()
    abs_output.mkdir(parents=True, exist_ok=True)
    volume_flags += ["-v", f"{abs_output}:/mnt/output"]

    # 5. Translate host paths → container paths
    def to_container(path: Path) -> str:
        return f"{dir_map[str(path.parent)]}/{path.name}"

    # 6. Build plexus run args
    run_args = [
        "--input",
        to_container(file_args["input"]),
        "--fasta",
        to_container(file_args["fasta"]),
        "--output",
        "/mnt/output",
    ]
    if file_args["snp_vcf"]:
        run_args += ["--snp-vcf", to_container(file_args["snp_vcf"])]
    if file_args["checksums"]:
        run_args += ["--checksums", to_container(file_args["checksums"])]
    if file_args["config"]:
        run_args += ["--config", to_container(file_args["config"])]

    # 7. Extra args passed through verbatim (non-file plexus run flags)
    extra_args: list[str] = ctx.args

    # 8. TTY flag for rich output
    tty_flag = ["-t"] if sys.stdout.isatty() else []

    # 9. Assemble and run
    cmd = [
        "docker",
        "run",
        "--rm",
        *tty_flag,
        *volume_flags,
        image,
        "run",
        *run_args,
        *extra_args,
    ]
    console.print(f"  Running: {' '.join(cmd)}\n")
    result = subprocess.run(cmd)
    raise typer.Exit(code=result.returncode)


if __name__ == "__main__":
    app()
