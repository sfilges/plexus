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
# Author: Stefan Filges (stefan@simsendiagnostics.com)
# Copyright (c) 2025 Simsen Diagnostics AB
# ================================================================================

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from loguru import logger

from multiplexdesigner.config import DesignerConfig, load_config
from multiplexdesigner.designer.design import design_primers
from multiplexdesigner.designer.multiplexpanel import MultiplexPanel, panel_factory
from multiplexdesigner.logging import configure_file_logging


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


def run_pipeline(
    input_file: str | Path,
    fasta_file: str | Path,
    output_dir: str | Path = "./output",
    panel_name: str = "multiplex_panel",
    genome: str = "hg38",
    preset: str = "default",
    config_path: str | Path | None = None,
    config_dict: dict[str, Any] | None = None,
    design_method: str = "simsen",
    run_blast: bool = True,
    padding: int = 200,
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
    config_dict : dict | None
        Configuration dictionary (highest priority).
    design_method : str
        Primer design algorithm ("simsen", "primer3py", or "primer3").
    run_blast : bool
        Whether to run BLAST specificity check (default: True).
    padding : int
        Bases to extract around each junction (default: 200).

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

    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    # Enable file logging to output directory
    log_file = configure_file_logging(str(output_dir))
    logger.info(f"Log file: {log_file}")

    # Load configuration
    logger.info("Loading configuration...")
    config = load_config(
        preset=preset,
        config_path=str(config_path) if config_path else None,
        config_dict=config_dict,
    )

    # Initialize result
    result = PipelineResult(
        panel=None,  # type: ignore[arg-type]
        output_dir=output_dir,
        config=config,
    )

    # =========================================================================
    # Step 1: Create panel and load junctions
    # =========================================================================
    logger.info(f"Creating panel '{panel_name}' with {genome} reference...")

    try:
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
        result.panel = panel
        result.steps_completed.append("panel_created")
        logger.info(f"Loaded {len(panel.junctions)} junctions")
    except Exception as e:
        logger.error(f"Failed to create panel: {e}")
        result.errors.append(f"Panel creation failed: {e}")
        raise

    # =========================================================================
    # Step 2: Design primers
    # =========================================================================
    logger.info(f"Designing primers using '{design_method}' method...")

    try:
        panel = design_primers(panel, method=design_method)
        result.steps_completed.append("primers_designed")

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
    # Step 4: Run BLAST specificity check (optional)
    # =========================================================================
    if run_blast:
        logger.info("Running BLAST specificity check...")

        try:
            from multiplexdesigner.blast.specificity import run_specificity_check

            blast_dir = output_dir / "blast"
            run_specificity_check(panel, str(blast_dir), str(fasta_file))
            result.steps_completed.append("specificity_checked")
            logger.info("Specificity check complete")
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
        from multiplexdesigner.selector.cost import MultiplexCostFunction
        from multiplexdesigner.selector.selectors import selector_collection

        selector_df = panel.build_selector_dataframe()
        pair_lookup = panel.build_pair_lookup()

        if selector_df.empty:
            logger.warning("No primer pairs available for multiplex optimization.")
        else:
            cost_fn = MultiplexCostFunction(
                pair_lookup, config.multiplex_picker_parameters
            )
            selector_cls = selector_collection["Greedy"]
            selector = selector_cls(selector_df, cost_fn)
            solutions = selector.run(
                N=config.multiplex_picker_parameters.initial_solutions
            )

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

        # Panel summary JSON
        panel.save_panel_summary_json(
            str(output_dir / "panel_summary.json"),
            result,
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
