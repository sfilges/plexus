# ================================================================================
# Docker wrapper command for plexus
#
# Extracted from cli.py to keep the main CLI module focused.
# Registers the 'docker' command on the shared Typer app.
# ================================================================================

from pathlib import Path
from typing import Annotated

import typer
from rich.console import Console

from plexus.cli import app
from plexus.version import __version__

console = Console()


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
