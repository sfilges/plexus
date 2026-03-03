# ================================================================================
# Main pipeline for multiplex primer panel design
#
# This module orchestrates the complete primer design workflow:
#   1. Load configuration and create panel
#   2. Design primers for each junction
#   3. Check specificity via BLAST (optional)
#   4. Multiplex optimization — select best primer combination
#   5. Save outputs
#
# Author: Stefan Filges (stefan.filges@pm.me)
# Copyright (c) 2026 Stefan Filges
# ================================================================================

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from loguru import logger

from plexus.config import DesignerConfig, load_config
from plexus.designer.design import design_primers
from plexus.designer.multiplexpanel import MultiplexPanel, panel_factory
from plexus.logging import configure_file_logging
from plexus.utils.env import (
    check_disk_space,
    get_missing_tools,
    get_plexus_version,
    get_primer3_version,
    get_tool_versions,
)


@dataclass
class PipelineResult:
    """Result of running the primer design pipeline."""

    panel: MultiplexPanel
    output_dir: Path
    config: DesignerConfig
    steps_completed: list[str] = field(default_factory=list)
    errors: list[str] = field(default_factory=list)
    multiplex_solutions: list = field(default_factory=list)
    selected_pairs: list = field(default_factory=list)
    failed_junctions: list = field(default_factory=list)

    @property
    def success(self) -> bool:
        """Returns True if pipeline completed without errors."""
        return len(self.errors) == 0

    @property
    def num_junctions(self) -> int:
        """Number of junctions in the panel."""
        return len(self.panel.junctions)

    @property
    def num_primer_pairs(self) -> int:
        """Total number of primer pairs across all junctions."""
        return sum(len(j.primer_pairs) for j in self.panel.junctions if j.primer_pairs)


@dataclass
class MultiPanelResult:
    """Aggregated result across multiple panel pipeline runs."""

    panel_results: dict[str, PipelineResult]
    output_dir: Path
    panel_ids: list[str]

    @property
    def success(self) -> bool:
        """True if all panels completed without errors."""
        return all(r.success for r in self.panel_results.values())

    @property
    def failed_panels(self) -> list[str]:
        """Panel IDs that had errors."""
        return [pid for pid, r in self.panel_results.items() if not r.success]

    @property
    def total_junctions(self) -> int:
        """Total junctions across all panels."""
        return sum(r.num_junctions for r in self.panel_results.values())

    @property
    def total_selected_pairs(self) -> int:
        """Total selected primer pairs across all panels."""
        return sum(len(r.selected_pairs) for r in self.panel_results.values())

    def summary_dict(self) -> dict:
        """Return a JSON-serializable summary across all panels."""
        return {
            "num_panels": len(self.panel_ids),
            "panel_ids": self.panel_ids,
            "total_junctions": self.total_junctions,
            "total_selected_pairs": self.total_selected_pairs,
            "all_successful": self.success,
            "failed_panels": self.failed_panels,
            "per_panel": {
                pid: {
                    "junctions": r.num_junctions,
                    "primer_pairs": r.num_primer_pairs,
                    "selected_pairs": len(r.selected_pairs),
                    "success": r.success,
                    "errors": r.errors,
                }
                for pid, r in self.panel_results.items()
            },
        }


def _collect_provenance(
    fasta_file: Path,
    snp_vcf: str | Path | None,
    genome: str,
    run_blast: bool,
    skip_snpcheck: bool,
    fasta_sha256: str | None = None,
    snp_vcf_sha256: str | None = None,
    compliance_report: dict | None = None,
    selector_seed: int | None = None,
    chrom_naming_check: str = "skipped",
) -> dict:
    """Collect tool versions, resource paths, and checksums for provenance."""
    import datetime

    from plexus.resources import _load_registry, get_operational_mode

    tool_versions = get_tool_versions(["blastn", "bcftools", "makeblastdb"])

    # Fall back to registry lookup only when caller did not supply hashes
    if fasta_sha256 is None or snp_vcf_sha256 is None:
        try:
            registry = _load_registry()
            for entry in registry.values():
                if fasta_sha256 is None and entry.get("fasta") == str(fasta_file):
                    fasta_sha256 = entry.get("fasta_sha256")
                vcf_str = str(snp_vcf) if snp_vcf else None
                if (
                    snp_vcf_sha256 is None
                    and vcf_str
                    and entry.get("snp_vcf") == vcf_str
                ):
                    snp_vcf_sha256 = entry.get("snp_vcf_sha256")
        except Exception:
            pass

    record: dict = {
        "plexus_version": get_plexus_version(),
        "primer3_version": get_primer3_version(),
        "tool_versions": tool_versions,
        "run_timestamp": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "operational_mode": get_operational_mode(),
        "genome": genome,
        "fasta_path": str(fasta_file),
        "fasta_sha256": fasta_sha256,
        "snp_vcf_path": str(snp_vcf) if snp_vcf else None,
        "snp_vcf_sha256": snp_vcf_sha256,
        "run_blast": run_blast,
        "skip_snpcheck": skip_snpcheck,
        "selector_seed": selector_seed,
        "chrom_naming_check": chrom_naming_check,
        "status": "started",
        "completed_at": None,
    }

    if compliance_report is not None:
        record["compliance_environment"] = compliance_report

    return record


def _resolve_source_vcf(
    snp_vcf: str | Path | None,
    genome: str,
) -> Path | None:
    """Return the effective SNP VCF source path without intersecting.

    Priority: user arg → $PLEXUS_SNP_VCF → registry.
    Returns None if nothing is available.
    """
    import os

    from plexus.snpcheck.resources import ENV_SNP_VCF

    if snp_vcf is not None:
        return Path(snp_vcf)
    env_vcf = os.environ.get(ENV_SNP_VCF)
    if env_vcf:
        return Path(env_vcf)
    from plexus.resources import get_registered_snp_vcf

    return get_registered_snp_vcf(genome)


def run_pipeline(
    input_file: str | Path,
    fasta_file: str | Path,
    output_dir: str | Path = "./output",
    panel_name: str = "multiplex_panel",
    genome: str = "hg38",
    preset: str = "default",
    config_path: str | Path | None = None,
    design_method: str = "plexus",
    run_blast: bool = True,
    padding: int = 200,
    snp_vcf: str | Path | None = None,
    skip_snpcheck: bool = False,
    snp_af_threshold: float | None = None,
    snp_strict: bool = False,
    selector: str = "Greedy",
    selector_seed: int | None = None,
    blast_num_threads: int = 4,
    debug: bool = False,
    fasta_sha256: str | None = None,
    snp_vcf_sha256: str | None = None,
) -> PipelineResult:
    """
    Run the complete multiplex primer design pipeline.

    Parameters
    ----------
    input_file : str | Path
        Path to CSV file containing junction coordinates.
    fasta_file : str | Path
        Path to reference genome FASTA file.
    output_dir : str | Path
        Directory for output files (default: "./output").
    panel_name : str
        Name for the primer panel (default: "multiplex_panel").
    genome : str
        Reference genome name (default: "hg38").
    preset : str
        Configuration preset ("default" or "lenient").
    config_path : str | Path | None
        Path to custom configuration JSON file.
    design_method : str
        Primer design algorithm (default: "plexus").
    run_blast : bool
        Whether to run BLAST specificity check (default: True).
    padding : int
        Bases to extract around each junction (default: 200).
    selector : str
        Multiplex selector algorithm. One of "Greedy", "Random",
        "BruteForce", "SimulatedAnnealing", or "DFS" (default: "Greedy").

    Returns
    -------
    PipelineResult
        Object containing the panel, outputs, and status information.

    Raises
    ------
    FileNotFoundError
        If input_file or fasta_file doesn't exist.
    ValidationError
        If configuration is invalid.
    """
    # Convert paths
    input_file = Path(input_file)
    fasta_file = Path(fasta_file)
    output_dir = Path(output_dir)

    # Validate inputs exist
    if not input_file.exists():
        raise FileNotFoundError(f"Input file not found: {input_file}")
    if not fasta_file.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_file}")

    # --- Compliance environment guard ---
    from plexus.resources import get_operational_mode
    from plexus.utils.env import ComplianceError, validate_environment

    op_mode = get_operational_mode()
    compliance_report: dict | None = None
    if op_mode == "compliance":
        try:
            compliance_report = validate_environment(
                need_blast=run_blast,
                need_snp=not skip_snpcheck,
            )
        except ComplianceError as e:
            logger.error(f"Compliance environment check failed: {e}")
            raise

    # --- Pre-flight dependency check ---
    missing = get_missing_tools(need_blast=run_blast, need_snp=not skip_snpcheck)
    if missing:
        msg = (
            f"Missing required system dependencies: {', '.join(missing)}. "
            "Install them via conda: 'conda install -c bioconda blast bcftools' "
            "or your system package manager."
        )
        logger.error(msg)
        raise RuntimeError(msg)

    # Check for sufficient disk space in output directory
    check_disk_space(output_dir)

    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load configuration
    logger.info("Loading configuration...")
    config = load_config(
        preset=preset,
        config_path=str(config_path) if config_path else None,
    )

    # Override config seed if provided as argument
    if selector_seed is not None:
        config.multiplex_picker_parameters.selector_seed = selector_seed

    # ── Chromosome naming consistency check ──────────────────────────────────
    effective_snp_vcf = (
        _resolve_source_vcf(snp_vcf, genome) if not skip_snpcheck else None
    )
    chrom_naming_check = "skipped"
    if not skip_snpcheck and effective_snp_vcf is not None:
        from plexus.utils.utils import (
            detect_chrom_naming_mismatch,
            get_fasta_contigs,
            get_vcf_contigs,
        )

        try:
            fasta_contigs = get_fasta_contigs(fasta_file)
            vcf_contigs = get_vcf_contigs(effective_snp_vcf)
            mismatch = detect_chrom_naming_mismatch(fasta_contigs, vcf_contigs)
            if mismatch:
                chrom_naming_check = "mismatch"
                logger.warning(
                    f"Chromosome naming mismatch between FASTA and VCF: {mismatch}"
                )
                if op_mode == "compliance":
                    raise ValueError(
                        f"Chromosome naming mismatch detected: {mismatch} "
                        "Ensure FASTA and VCF use the same chromosome naming convention "
                        "(e.g., both 'chr1' or both '1'). Aborting in compliance mode."
                    )
            else:
                chrom_naming_check = "pass"
        except ValueError:
            raise  # re-raise the compliance ValueError above
        except Exception as e:
            chrom_naming_check = "unavailable"
            logger.warning(f"Could not verify chromosome naming consistency: {e}")

    # ── Write provenance record ──────────────────────────────────────────────
    provenance = _collect_provenance(
        fasta_file=fasta_file,
        snp_vcf=effective_snp_vcf,
        genome=genome,
        run_blast=run_blast,
        skip_snpcheck=skip_snpcheck,
        fasta_sha256=fasta_sha256,
        snp_vcf_sha256=snp_vcf_sha256,
        compliance_report=compliance_report,
        selector_seed=config.multiplex_picker_parameters.selector_seed,
        chrom_naming_check=chrom_naming_check,
    )
    import json as _json

    provenance_path = output_dir / "provenance.json"
    with provenance_path.open("w") as _pf:
        _json.dump(provenance, _pf, indent=2, default=str)

    # Sentinels for the finally handler — result/panel may not exist if an early step fails
    result: PipelineResult | None = None
    panel: MultiplexPanel | None = None
    _exc: BaseException | None = None

    try:
        # Enable file logging to output directory
        log_file = configure_file_logging(str(output_dir), debug=debug)
        logger.info(f"Log file: {log_file}")
        logger.info(f"Provenance written to {provenance_path}")

        # =========================================================================
        # Step 1: Create panel and load junctions
        # =========================================================================
        logger.info(f"Creating panel '{panel_name}' with {genome} reference...")

        panel = panel_factory(
            name=panel_name,
            genome=genome,
            design_input_file=str(input_file),
            fasta_file=str(fasta_file),
            preset=preset,
            config_file=str(config_path) if config_path else None,
            padding=padding,
        )
        # Update panel config to use our loaded config
        panel.config = config
        logger.info(f"Loaded {len(panel.junctions)} junctions")

        # Initialize result now that we have a valid panel
        result = PipelineResult(
            panel=panel,
            output_dir=output_dir,
            config=config,
            steps_completed=["panel_created"],
        )

        # =========================================================================
        # Step 2: Design primers
        # =========================================================================
        logger.info(f"Designing primers using '{design_method}' method...")

        try:
            panel = design_primers(panel, method=design_method)
            result.steps_completed.append("primers_designed")

            # Capture any junctions that failed during design
            if hasattr(panel, "failed_junctions") and panel.failed_junctions:
                result.failed_junctions = panel.failed_junctions
                for fj in panel.failed_junctions:
                    err_msg = getattr(fj, "_design_error", "no valid primer pairs")
                    result.errors.append(
                        f"Junction '{fj.name}' failed primer design: {err_msg}"
                    )

            if not panel.junctions:
                logger.error("All junctions failed primer design. Cannot continue.")
                result.errors.append("All junctions failed primer design.")
                return result

            total_pairs = sum(
                len(j.primer_pairs) for j in panel.junctions if j.primer_pairs
            )
            logger.info(
                f"Designed {total_pairs} primer pairs across {len(panel.junctions)} junctions"
            )
        except Exception as e:
            logger.error(f"Primer design failed: {e}")
            result.errors.append(f"Primer design failed: {e}")
            raise

        # =========================================================================
        # Step 3: Save intermediate results
        # =========================================================================
        logger.info("Saving candidate primer pairs...")

        try:
            pairs_file = output_dir / "candidate_pairs.csv"
            panel.save_candidate_pairs_to_csv(str(pairs_file))
            result.steps_completed.append("candidates_saved")
            logger.info(f"Saved candidate pairs to {pairs_file}")
        except Exception as e:
            logger.warning(f"Could not save candidate pairs: {e}")
            result.errors.append(f"Save candidates failed: {e}")

        # =========================================================================
        # Step 3.5: SNP overlap check (optional)
        # =========================================================================
        snp_config = config.snp_check_parameters
        _run_snpcheck = not skip_snpcheck

        if _run_snpcheck:
            af_thresh = (
                snp_af_threshold
                if snp_af_threshold is not None
                else snp_config.af_threshold
            )

            try:
                from plexus.snpcheck.snp_data import get_snp_vcf

                resolved_vcf = get_snp_vcf(
                    panel=panel,
                    output_dir=output_dir,
                    user_vcf=snp_vcf,
                    padding=padding,
                    genome=genome,
                )

                from plexus.snpcheck.checker import run_snp_check

                run_snp_check(
                    panel=panel,
                    vcf_path=str(resolved_vcf),
                    af_threshold=af_thresh,
                    snp_penalty_weight=snp_config.snp_penalty_weight,
                    snp_3prime_window=snp_config.snp_3prime_window,
                    snp_3prime_multiplier=snp_config.snp_3prime_multiplier,
                    snp_af_weight=snp_config.snp_af_weight,
                )
                result.steps_completed.append("snp_checked")

                if snp_strict or snp_config.snp_strict:
                    from plexus.snpcheck.checker import filter_snp_pairs

                    n_removed, fallback_junctions = filter_snp_pairs(panel)
                    logger.info(
                        f"SNP strict mode: removed {n_removed} primer pairs overlapping SNPs"
                    )
                    for name in fallback_junctions:
                        result.errors.append(
                            f"SNP strict: '{name}' — no SNP-free pairs found; least-affected pair kept"
                        )
                    result.steps_completed.append("snp_strict_filtered")
            except Exception as e:
                logger.error(f"SNP check failed: {e}")
                result.errors.append(f"SNP check failed: {e}")
        else:
            logger.info("Skipping SNP check")
            result.steps_completed.append("snp_check_skipped")

        # =========================================================================
        # Step 4: Run BLAST specificity check (optional)
        # =========================================================================
        if run_blast:
            logger.info(
                f"Running BLAST specificity check (num_threads={blast_num_threads})..."
            )

            try:
                from plexus.blast.specificity import (
                    filter_offtarget_pairs,
                    run_specificity_check,
                )

                blast_dir = output_dir / "blast"
                blast_config = config.blast_parameters
                run_specificity_check(
                    panel,
                    str(blast_dir),
                    str(fasta_file),
                    num_threads=blast_num_threads,
                    length_threshold=blast_config.length_threshold,
                    evalue_threshold=blast_config.evalue_threshold,
                    max_mismatches=blast_config.max_mismatches,
                    three_prime_tolerance=blast_config.three_prime_tolerance,
                    max_amplicon_size=blast_config.max_amplicon_size,
                    ontarget_tolerance=blast_config.ontarget_tolerance,
                    blast_evalue=blast_config.blast_evalue,
                    blast_word_size=blast_config.blast_word_size,
                    blast_reward=blast_config.blast_reward,
                    blast_penalty=blast_config.blast_penalty,
                    blast_max_hsps=blast_config.blast_max_hsps,
                    blast_dust=blast_config.blast_dust,
                )
                result.steps_completed.append("specificity_checked")
                logger.info("Specificity check complete")

                n_removed, fallback_junctions = filter_offtarget_pairs(panel)
                if n_removed > 0:
                    logger.info(
                        f"Off-target filter: removed {n_removed} primer pair(s) "
                        "with off-target products"
                    )
                for name in fallback_junctions:
                    result.errors.append(
                        f"Off-target filter: '{name}' — no clean pairs; "
                        "least-affected pair kept"
                    )
                result.steps_completed.append("offtarget_filtered")
            except ImportError as e:
                logger.warning(f"BLAST module not available: {e}")
                result.errors.append(f"BLAST not available: {e}")
            except Exception as e:
                logger.error(f"Specificity check failed: {e}")
                result.errors.append(f"Specificity check failed: {e}")
        else:
            logger.info("Skipping BLAST specificity check")
            result.steps_completed.append("specificity_skipped")

        # =========================================================================
        # Step 5: Multiplex optimization
        # =========================================================================
        logger.info("Running multiplex optimization...")

        try:
            from plexus.selector.cost import MultiplexCostFunction
            from plexus.selector.selectors import selector_collection

            selector_df = panel.build_selector_dataframe()
            pair_lookup = panel.build_pair_lookup()

            if selector_df.empty:
                logger.warning("No primer pairs available for multiplex optimization.")
            elif selector not in selector_collection:
                raise ValueError(
                    f"Unknown selector '{selector}'. "
                    f"Available: {', '.join(selector_collection)}"
                )
            else:
                cost_fn = MultiplexCostFunction(
                    pair_lookup, config.multiplex_picker_parameters
                )

                # Warn if stochastic selector is used without a seed
                stochastic_selectors = ["Greedy", "Random", "SimulatedAnnealing"]
                if (
                    selector in stochastic_selectors
                    and config.multiplex_picker_parameters.selector_seed is None
                ):
                    if op_mode == "compliance":
                        logger.warning(
                            f"Compliance mode: stochastic selector '{selector}' used without a seed. "
                            "This run cannot be identically reproduced."
                        )
                    else:
                        logger.info(
                            f"Stochastic selector '{selector}' used without a seed. "
                            "Results will vary between runs."
                        )

                selector_cls = selector_collection[selector]
                selector_obj = selector_cls(
                    selector_df,
                    cost_fn,
                    seed=config.multiplex_picker_parameters.selector_seed,
                )
                logger.info(f"Using '{selector}' selector algorithm.")
                if selector in ("Greedy", "Random"):
                    solutions = selector_obj.run(
                        N=config.multiplex_picker_parameters.initial_solutions
                    )
                else:
                    solutions = selector_obj.run()

                # Sort by cost and keep top solutions
                solutions.sort(key=lambda m: m.cost)
                top_n = config.multiplex_picker_parameters.top_solutions_to_keep
                result.multiplex_solutions = solutions[:top_n]

                # Auto-apply the best solution
                if solutions:
                    best = solutions[0]
                    selected = []
                    for pair_id in best.primer_pairs:
                        pair = pair_lookup.get(pair_id)
                        if pair:
                            pair.selected = True
                            selected.append(pair)
                    result.selected_pairs = selected
                    logger.info(
                        f"Selected {len(selected)} primer pairs (best cost: {best.cost:.2f})"
                    )

                result.steps_completed.append("multiplex_optimized")
        except Exception as e:
            logger.error(f"Multiplex optimization failed: {e}")
            result.errors.append(f"Multiplex optimization failed: {e}")

        # =========================================================================
        # Step 6: Save final results
        # =========================================================================
        logger.info("Saving final results...")

        try:
            # Selected multiplex (best solution)
            if result.selected_pairs:
                panel.save_selected_multiplex_csv(
                    str(output_dir / "selected_multiplex.csv"),
                    result.selected_pairs,
                )

            # Top N panel solutions
            if result.multiplex_solutions:
                panel.save_top_panels_csv(
                    str(output_dir / "top_panels.csv"),
                    result.multiplex_solutions,
                )

            # Off-target details for selected pairs
            if result.selected_pairs:
                panel.save_off_targets_csv(
                    str(output_dir / "off_targets.csv"),
                    result.selected_pairs,
                )

            # Panel QC report (REPT-01)
            if result.selected_pairs:
                try:
                    from plexus.reporting.qc import generate_panel_qc

                    qc_data = generate_panel_qc(panel.junctions)
                    qc_path = output_dir / "panel_qc.json"
                    with qc_path.open("w") as f:
                        _json.dump(qc_data, f, indent=2)
                    logger.info(f"Wrote panel QC report to {qc_path.name}")
                except Exception as e:
                    logger.warning(f"Could not write panel QC report: {e}")
                    result.errors.append(f"Panel QC report failed: {e}")

            # Failed junctions report
            if result.failed_junctions:
                import pandas as pd

                rows = [
                    {
                        "Junction": fj.name,
                        "Chrom": fj.chrom,
                        "Start": fj.start,
                        "End": fj.end,
                        "Error": getattr(fj, "_design_error", "unknown"),
                    }
                    for fj in result.failed_junctions
                ]
                pd.DataFrame(rows).to_csv(
                    output_dir / "failed_junctions.csv", index=False
                )
                logger.info(
                    f"Wrote {len(rows)} failed junction(s) to failed_junctions.csv"
                )

            result.steps_completed.append("final_results_saved")
        except Exception as e:
            logger.warning(f"Could not save final results: {e}")
            result.errors.append(f"Save final results failed: {e}")

        # =========================================================================
        # Summary
        # =========================================================================
        logger.info("=" * 60)
        logger.info("Pipeline Summary")
        logger.info("=" * 60)
        logger.info(f"Panel: {panel.panel_name}")
        logger.info(f"Junctions: {len(panel.junctions)}")
        logger.info(f"Total candidate primer pairs: {result.num_primer_pairs}")
        logger.info(f"Selected primer pairs: {len(result.selected_pairs)}")
        if result.multiplex_solutions:
            logger.info(
                f"Top {len(result.multiplex_solutions)} multiplex solutions "
                f"(best cost: {result.multiplex_solutions[0].cost:.2f})"
            )
        logger.info(f"Steps completed: {', '.join(result.steps_completed)}")
        if result.errors:
            logger.warning(f"Errors encountered: {len(result.errors)}")
            for err in result.errors:
                logger.warning(f"  - {err}")
        logger.info(f"Output directory: {output_dir}")
        logger.info("=" * 60)

        return result

    except BaseException as exc:
        _exc = exc
        raise

    finally:
        import datetime as _dt

        _completed_at = _dt.datetime.now(_dt.timezone.utc).isoformat()

        # Re-read the file to preserve all original fields written at startup
        try:
            with provenance_path.open() as _f:
                _final_prov = _json.load(_f)
        except Exception:
            _final_prov = dict(provenance)

        if _exc is not None:
            _final_prov["status"] = "failed"
            _final_prov["errors"] = result.errors if result else [str(_exc)]
        else:
            _final_prov["status"] = "completed"
            _final_prov["errors"] = result.errors if result else []

        _final_prov["completed_at"] = _completed_at
        _final_prov["steps_completed"] = result.steps_completed if result else []

        try:
            with provenance_path.open("w") as _f:
                _json.dump(_final_prov, _f, indent=2, default=str)
        except Exception:
            pass  # Never mask the original pipeline exception

        # Rewrite panel_summary.json with finalized provenance
        try:
            if panel is not None and result is not None:
                panel.save_panel_summary_json(
                    str(output_dir / "panel_summary.json"),
                    result,
                    provenance=_final_prov,
                )
        except Exception:
            pass  # Never mask the original pipeline exception
