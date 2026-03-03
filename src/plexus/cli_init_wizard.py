# ================================================================================
# Interactive init wizard for `plexus init`
#
# Extracted from cli.py to keep the main CLI module focused.
# These are plain helper functions — the init command itself stays in cli.py.
# ================================================================================

from pathlib import Path

import typer
from rich.console import Console

console = Console()


def _is_interactive() -> bool:
    """Return True if stdin is an interactive TTY."""
    import sys

    return sys.stdin.isatty()


def _prompt_path(label: str, *, must_exist: bool = True) -> Path:
    """Prompt for a file path, re-asking until a valid path is given."""
    from rich.prompt import Prompt

    while True:
        raw = Prompt.ask(f"  {label}")
        if not raw or not raw.strip():
            console.print("    [red]Path cannot be empty.[/red]")
            continue
        p = Path(raw.strip()).expanduser().resolve()
        if must_exist and not p.is_file():
            console.print(f"    [red]File not found: {p}[/red]")
            continue
        return p


def _run_init_wizard(
    genome_presets: list[str],
) -> dict:
    """Interactive wizard for `plexus init`. Returns a dict of collected parameters."""

    from rich.panel import Panel
    from rich.prompt import Confirm, Prompt

    console.print()
    console.print(
        Panel(
            "[bold]Plexus — Resource Initialization Wizard[/bold]\n"
            "[dim]Answer the prompts below to register reference resources.\n"
            "Press Ctrl-C at any time to abort.[/dim]",
            border_style="green",
        )
    )

    # 1. Genome
    default_genome = genome_presets[0] if genome_presets else "hg38"
    genome = Prompt.ask(
        "  [bold]Genome[/bold]",
        choices=genome_presets,
        default=default_genome,
    )

    # 2. FASTA
    console.print()
    fasta = _prompt_path("[bold]Reference FASTA file[/bold]")

    # 3. SNP VCF
    console.print()
    register_vcf = Confirm.ask(
        "  [bold]Register a SNP VCF[/bold] (gnomAD)?", default=True
    )
    snp_vcf: Path | None = None
    skip_snp = True
    if register_vcf:
        snp_vcf = _prompt_path("[bold]SNP VCF file[/bold] (tabix-indexed .vcf.gz)")
        skip_snp = False

    # 4. Operational mode
    console.print()
    mode = Prompt.ask(
        "  [bold]Operational mode[/bold]",
        choices=["research", "compliance"],
        default="research",
    )

    # 5. Checksums
    console.print()
    use_checksums = Confirm.ask(
        "  [bold]Verify files with a checksums file[/bold]?", default=False
    )
    checksums: Path | None = None
    if use_checksums:
        checksums = _prompt_path("[bold]Checksums file[/bold] (sha256sum format)")

    # 6. BLAST index
    console.print()
    build_blast = Confirm.ask("  [bold]Build BLAST index[/bold]?", default=True)

    # 7. Force rebuild
    force = False
    if build_blast:
        force = Confirm.ask(
            "  [bold]Force rebuild[/bold] existing indexes?", default=False
        )

    # ── Summary ──────────────────────────────────────────────────────────────
    console.print()
    summary_lines = [
        f"  Genome:     [bold]{genome}[/bold]",
        f"  FASTA:      {fasta}",
    ]
    if snp_vcf:
        summary_lines.append(f"  SNP VCF:    {snp_vcf}")
    else:
        summary_lines.append("  SNP VCF:    [dim]skipped[/dim]")
    summary_lines.append(f"  Mode:       {mode}")
    if checksums:
        summary_lines.append(f"  Checksums:  {checksums}")
    summary_lines.append(
        f"  BLAST:      {'build' if build_blast else '[dim]skip[/dim]'}"
    )
    if force:
        summary_lines.append("  Force:      yes")

    console.print(
        Panel(
            "\n".join(summary_lines),
            title="[bold]Summary[/bold]",
            border_style="cyan",
        )
    )

    if not Confirm.ask("  [bold]Proceed with initialization?[/bold]", default=True):
        console.print("[yellow]Aborted.[/yellow]")
        raise typer.Exit(code=0)

    return {
        "genome": genome,
        "fasta": fasta,
        "snp_vcf": snp_vcf,
        "skip_snp": skip_snp,
        "skip_blast": not build_blast,
        "force": force,
        "mode": mode,
        "checksums": checksums,
    }
