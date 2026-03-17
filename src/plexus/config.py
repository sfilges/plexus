# ================================================================================
# Configuration models for multiplex primer design
#
# Uses Pydantic for runtime validation of configuration parameters.
#
# Author: Stefan Filges (stefan.filges@pm.me)
# Copyright (c) 2026 Stefan Filges
# ================================================================================

from __future__ import annotations

import json
from importlib.resources import files
from pathlib import Path
from typing import Any, Literal

from loguru import logger
from pydantic import BaseModel, Field, model_validator


class SingleplexDesignParameters(BaseModel):
    """Parameters for individual primer design."""

    # Number of primers to return
    PRIMER_NUM_RETURN: int = Field(default=10, ge=1, le=1000)

    # Temperature parameters (°C)
    PRIMER_OPT_TM: float = Field(default=60.0, ge=50.0, le=72.0)
    PRIMER_MIN_TM: float = Field(default=57.0, ge=40.0, le=70.0)
    PRIMER_MAX_TM: float = Field(default=63.0, ge=50.0, le=80.0)

    # Size parameters (bp)
    PRIMER_OPT_SIZE: int = Field(default=22, ge=15, le=36)
    primer_min_length: int = Field(default=18, ge=10, le=30)
    primer_max_length: int = Field(default=28, ge=15, le=40)

    # Bound parameters (fraction bound at annealing temp)
    PRIMER_OPT_BOUND: float = Field(default=95.0, ge=0.0, le=100.0)
    PRIMER_MIN_BOUND: float = Field(default=-10.0, ge=-100.0, le=100.0)
    PRIMER_MAX_BOUND: float = Field(default=120.0, ge=0.0, le=200.0)

    # Junction parameters
    junction_padding_bases: int = Field(default=3, ge=0, le=50)

    # Tail sequences (optional adapter sequences prepended to each primer)
    forward_tail: str = Field(
        default="GGACACTCTTTCCCTACACGACGCTCTTCCGATCTAAAAAAAAAAAAAAAAAAAATGGGAAAGAGTGTCC",
        alias="forward_tail",
    )
    reverse_tail: str = Field(
        default="GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT",
        alias="reverse_tail",
    )

    # Penalty weights
    primer_length_penalty: float = Field(default=1.0, ge=0.0)
    primer_complexity_penalty: float = Field(default=1.0, ge=0.0)
    amplicon_length_penalty: float = Field(default=1.0, ge=0.0)

    # GC content parameters (%)
    PRIMER_OPT_GC_PERCENT: float = Field(default=50.0, ge=20.0, le=80.0)
    primer_min_gc: int = Field(default=30, ge=0, le=100)
    primer_max_gc: int = Field(default=70, ge=0, le=100)
    primer_gc_clamp: int = Field(
        default=0,
        ge=0,
        le=5,
        description="Enable 3' GC clamp filter (0=off, 1=on). Requires 1-3 G/C in last 5 bases.",
    )

    # Sequence composition constraints
    primer_max_poly_x: int = Field(default=4, ge=1, le=10)
    primer_max_poly_gc: int = Field(default=3, ge=1, le=10)
    primer_max_n: int = Field(default=0, ge=0, le=5)

    # Thermodynamic thresholds (°C)
    PRIMER_MAX_SELF_ANY_TH: float = Field(default=45.0, ge=0.0, le=100.0)
    PRIMER_MAX_SELF_END_TH: float = Field(default=35.0, ge=0.0, le=100.0)
    PRIMER_MAX_HAIRPIN_TH: float = Field(default=24.0, ge=0.0, le=100.0)
    PRIMER_MAX_END_STABILITY: float = Field(default=4.5, ge=0.0, le=15.0)
    PRIMER_MAX_TEMPLATE_MISPRIMING_TH: float = Field(default=35.0, ge=0.0, le=100.0)

    # Penalty weights for primer scoring
    PRIMER_WT_SIZE_LT: float = Field(default=1.0, ge=0.0)
    PRIMER_WT_SIZE_GT: float = Field(default=1.0, ge=0.0)
    PRIMER_WT_TM_GT: float = Field(default=1.0, ge=0.0)
    PRIMER_WT_TM_LT: float = Field(default=1.0, ge=0.0)
    PRIMER_WT_BOUND_GT: float = Field(default=1.0, ge=0.0)
    PRIMER_WT_BOUND_LT: float = Field(default=1.0, ge=0.0)
    PRIMER_WT_GC_PERCENT_GT: float = Field(default=0.0, ge=0.0)
    PRIMER_WT_GC_PERCENT_LT: float = Field(default=0.0, ge=0.0)
    PRIMER_WT_SELF_ANY_TH: float = Field(default=0.0, ge=0.0)
    PRIMER_WT_SELF_END_TH: float = Field(default=0.0, ge=0.0)
    PRIMER_WT_HAIRPIN_TH: float = Field(default=0.0, ge=0.0)
    PRIMER_WT_END_STABILITY: float = Field(default=0.0, ge=0.0)

    model_config = {"populate_by_name": True}

    @model_validator(mode="after")
    def validate_tm_range(self) -> SingleplexDesignParameters:
        """Validate that min_tm <= opt_tm <= max_tm."""
        if not (self.PRIMER_MIN_TM <= self.PRIMER_OPT_TM <= self.PRIMER_MAX_TM):
            raise ValueError(
                f"Tm values must satisfy: min ({self.PRIMER_MIN_TM}) <= opt ({self.PRIMER_OPT_TM}) <= max ({self.PRIMER_MAX_TM})"
            )
        return self

    @model_validator(mode="after")
    def validate_length_range(self) -> SingleplexDesignParameters:
        """Validate that min_length <= opt_size <= max_length."""
        if not (
            self.primer_min_length <= self.PRIMER_OPT_SIZE <= self.primer_max_length
        ):
            raise ValueError(
                f"Length values must satisfy: min ({self.primer_min_length}) <= opt ({self.PRIMER_OPT_SIZE}) <= max ({self.primer_max_length})"
            )
        return self

    @model_validator(mode="after")
    def validate_gc_range(self) -> SingleplexDesignParameters:
        """Validate that min_gc <= opt_gc <= max_gc."""
        if not (self.primer_min_gc <= self.PRIMER_OPT_GC_PERCENT <= self.primer_max_gc):
            raise ValueError(
                f"GC values must satisfy: min ({self.primer_min_gc}) <= opt ({self.PRIMER_OPT_GC_PERCENT}) <= max ({self.primer_max_gc})"
            )
        return self


class PrimerPairParameters(BaseModel):
    """Parameters for primer pair selection."""

    PRIMER_PAIR_MAX_DIFF_TM: float = Field(default=3.0, ge=0.0, le=10.0)
    PRIMER_PRODUCT_OPT_SIZE: int = Field(default=60, ge=30, le=500)
    PRIMER_PRODUCT_MIN_INSERT_SIZE: int = Field(default=20, ge=0, le=200)
    PRIMER_PRODUCT_MAX_INSERT_SIZE: int = Field(default=60, ge=10, le=500)
    PRIMER_PRODUCT_MAX_SIZE: int = Field(default=120, ge=40, le=1000)

    # Penalty weights for pair scoring
    PRIMER_PAIR_WT_PR_PENALTY: float = Field(default=1.0, ge=0.0)
    PRIMER_PAIR_WT_DIFF_TM: float = Field(default=0.0, ge=0.0)
    PRIMER_PAIR_WT_PRODUCT_SIZE_LT: float = Field(default=0.5, ge=0.0)
    PRIMER_PAIR_WT_PRODUCT_SIZE_GT: float = Field(default=2.0, ge=0.0)

    @model_validator(mode="after")
    def validate_insert_size_range(self) -> PrimerPairParameters:
        """Validate that min_insert <= max_insert <= max_product."""
        if self.PRIMER_PRODUCT_MIN_INSERT_SIZE > self.PRIMER_PRODUCT_MAX_INSERT_SIZE:
            raise ValueError(
                f"Insert size: min ({self.PRIMER_PRODUCT_MIN_INSERT_SIZE}) must be <= max ({self.PRIMER_PRODUCT_MAX_INSERT_SIZE})"
            )
        return self


class PCRConditions(BaseModel):
    """PCR reaction conditions for thermodynamic calculations."""

    annealing_temperature: float = Field(default=60.0, ge=45.0, le=72.0)
    primer_concentration: float = Field(default=50.0, ge=1.0, le=1000.0)  # nM
    dntp_concentration: float = Field(default=0.6, ge=0.0, le=10.0)  # mM
    dna_concentration: float = Field(default=50.0, ge=0.0, le=1000.0)  # nM
    mv_concentration: float = Field(default=50.0, ge=0.0, le=500.0)  # mM (monovalent)
    dv_concentration: float = Field(default=1.5, ge=0.0, le=50.0)  # mM (divalent)
    dmso_concentration: float = Field(default=0.0, ge=0.0, le=20.0)  # %
    dmso_fact: float = Field(default=0.6, ge=0.0, le=1.0)
    formamide_concentration: float = Field(default=0.8, ge=0.0, le=10.0)  # M


class MultiplexPickerParameters(BaseModel):
    """Parameters for multiplex panel optimization."""

    initial_solutions: int = Field(default=100, ge=1, le=10000)
    top_solutions_to_keep: int = Field(default=4, ge=1, le=100)
    selector_seed: int | None = Field(
        default=None,
        description="Random seed for stochastic selectors (Greedy, Random, SimulatedAnnealing).",
    )

    target_plexity: int = Field(default=20, ge=2, le=500)
    minimum_plexity: int = Field(default=10, ge=1, le=500)
    maximum_plexity: int = Field(default=50, ge=2, le=1000)

    plexity_wt_lt: float = Field(default=1.0, ge=0.0)
    plexity_wt_gt: float = Field(default=1.0, ge=0.0)

    force_plexity: bool = Field(default=False)
    allow_split_panel: bool = Field(default=False)
    max_splits: int = Field(default=2, ge=1, le=10)

    # Cost function weights for multiplex optimization
    wt_pair_penalty: float = Field(default=1.0, ge=0.0)
    wt_off_target: float = Field(default=5.0, ge=0.0)
    wt_cross_dimer: float = Field(default=1.0, ge=0.0)
    wt_pair_dimer: float = Field(
        default=1.0,
        ge=0.0,
        description="Weight for intra-pair F/R dimer score in cost function.",
    )
    wt_snp_penalty: float = Field(
        default=3.0,
        ge=0.0,
        description="Weight for SNP penalty in cost function. "
        "Independent of wt_pair_penalty. Set to 0 to disable.",
    )

    @model_validator(mode="after")
    def validate_plexity_range(self) -> MultiplexPickerParameters:
        """Validate that min <= target <= max plexity."""
        if not (self.minimum_plexity <= self.target_plexity <= self.maximum_plexity):
            raise ValueError(
                f"Plexity values must satisfy: min ({self.minimum_plexity}) <= target ({self.target_plexity}) <= max ({self.maximum_plexity})"
            )
        return self


class SnpCheckParameters(BaseModel):
    """Parameters for SNP overlap checking."""

    af_threshold: float = Field(
        default=0.01, ge=0.0, le=1.0, description="Minimum allele frequency to flag"
    )
    snp_penalty_weight: float = Field(
        default=10.0, ge=0.0, description="Base penalty per SNP overlapping a primer"
    )
    snp_3prime_window: int = Field(
        default=5,
        ge=1,
        le=15,
        description="Number of bases from 3' end considered high-impact",
    )
    snp_3prime_multiplier: float = Field(
        default=3.0,
        ge=1.0,
        description="Penalty multiplier for SNPs within the 3' window",
    )
    snp_strict: bool = Field(
        default=False,
        description="Discard primer pairs overlapping any SNP above af_threshold",
    )
    snp_af_weight: float = Field(
        default=0.0,
        ge=0.0,
        description=(
            "Exponent for AF-based penalty scaling, normalised to af_threshold. "
            "0.0 = no scaling (all passing SNPs penalised equally). "
            "1.0 = linear scaling (SNP at 10× threshold → 10× penalty). "
            "0.5 = sqrt scaling (recommended for most panels)."
        ),
    )


class BlastParameters(BaseModel):
    """Parameters for BLAST specificity checking."""

    length_threshold: int = Field(
        default=15,
        ge=5,
        le=30,
        description=(
            "Minimum 3'-anchored alignment length (bp) to predict primer binding. "
            "A BLAST hit with at least this many bases aligned from the 3' end "
            "is classified as 'predicted bound'."
        ),
    )
    evalue_threshold: float = Field(
        default=10.0,
        gt=0.0,
        description=(
            "E-value cutoff for predicted binding. Hits with e-value below this "
            "threshold (and anchored at the 3' end) are classified as 'predicted bound'. "
            "High default (10) is appropriate for short primer queries where "
            "e-values are naturally large."
        ),
    )
    max_mismatches: int = Field(
        default=3,
        ge=0,
        le=5,
        description=(
            "Maximum mismatches allowed in a 3'-anchored alignment for it to be "
            "classified as 'predicted bound'. A primer with 1-3 mismatches in a "
            "15+ bp 3' stretch will still extend in PCR."
        ),
    )
    three_prime_tolerance: int = Field(
        default=3,
        ge=0,
        le=5,
        description=(
            "Maximum number of unaligned bases at the 3' end of the primer for "
            "a BLAST hit to still be considered '3'-anchored'. BLAST's local "
            "alignment clips terminal mismatches, so a primer with mismatches "
            "at or near the 3' end will have qend < qlen. A tolerance of 3 "
            "catches these real binding events (e.g. DCAF12L1 paralogs) that "
            "Primer-BLAST detects via its semi-global alignment."
        ),
    )
    blast_evalue: float = Field(
        default=30000.0,
        gt=0.0,
        description=(
            "E-value passed to blastn via -evalue. Matches the Primer-BLAST "
            "default (30000) to ensure short, gapped primer alignments are "
            "not silently dropped by BLAST before plexus can evaluate them."
        ),
    )
    blast_word_size: int = Field(
        default=7,
        ge=4,
        le=28,
        description=(
            "Word size passed to blastn via -word_size. Matches the "
            "Primer-BLAST / blastn-short default of 7, which is appropriate "
            "for primer-length queries (<30 bp)."
        ),
    )
    blast_reward: int = Field(
        default=1,
        ge=1,
        le=5,
        description="Match reward passed to blastn via -reward.",
    )
    blast_penalty: int = Field(
        default=-1,
        ge=-5,
        le=-1,
        description=(
            "Mismatch penalty passed to blastn via -penalty. "
            "The blastn-short default (-3) causes 3-mismatch paralog alignments "
            "to score too low (E > 100,000) for BLAST to report them. "
            "A gentler penalty (-1) keeps these alignments reportable "
            "(E ~ 70 for a 21bp primer with 3 mismatches) while the annotator's "
            "3'-anchored filters still reject true noise."
        ),
    )
    blast_max_hsps: int = Field(
        default=100,
        ge=1,
        le=10000,
        description=(
            "Maximum HSPs per query-subject pair passed to blastn via -max_hsps. "
            "Each chromosome is one subject; 100 is generous for real binding sites "
            "and avoids wasting time on thousands of low-quality HSPs."
        ),
    )
    blast_dust: str = Field(
        default="yes",
        pattern=r"^(yes|no|\d+ \d+ \d+)$",
        description=(
            "Low-complexity filter passed to blastn via -dust. "
            "blastn-short defaults to 'no'; enabling it ('yes') avoids extending "
            "seeds in Alu/LINE regions, improving runtime without affecting "
            "off-target detection."
        ),
    )
    max_bound_per_primer: int | None = Field(
        default=10000,
        ge=100,
        description=(
            "Maximum predicted binding sites to retain per primer, ranked by "
            "lowest e-value. Acts as a safety net against primers in highly "
            "repetitive regions (e.g. Alu elements) that can cause combinatorial "
            "explosion in amplicon finding. Normal primers have 700-7000 bound "
            "sites; only extreme outliers (>10000) are affected. "
            "Set to None to disable."
        ),
    )
    max_amplicon_size: int = Field(
        default=2000,
        ge=100,
        le=50000,
        description=(
            "Maximum distance (bp) between two predicted-bound primers for them "
            "to form an amplicon. Products larger than this are unlikely to "
            "amplify efficiently under standard PCR conditions."
        ),
    )
    ontarget_tolerance: int = Field(
        default=5,
        ge=0,
        le=50,
        description=(
            "Coordinate tolerance (bp) for classifying a BLAST amplicon as on-target. "
            "A hit is on-target if both primer positions are within this distance "
            "of the expected genomic coordinates. Absorbs minor BLAST alignment "
            "shifts while still distinguishing on-target from pseudogene hits."
        ),
    )


class DesignerConfig(BaseModel):
    """Complete configuration for multiplex primer panel design."""

    singleplex_design_parameters: SingleplexDesignParameters = Field(
        default_factory=SingleplexDesignParameters
    )
    primer_pair_parameters: PrimerPairParameters = Field(
        default_factory=PrimerPairParameters
    )
    pcr_conditions: PCRConditions = Field(default_factory=PCRConditions)
    multiplex_picker_parameters: MultiplexPickerParameters = Field(
        default_factory=MultiplexPickerParameters
    )
    snp_check_parameters: SnpCheckParameters = Field(default_factory=SnpCheckParameters)
    blast_parameters: BlastParameters = Field(default_factory=BlastParameters)

    @classmethod
    def from_json_file(cls, file_path: str | Path) -> DesignerConfig:
        """
        Load configuration from a JSON file.

        Parameters
        ----------
        file_path : str | Path
            Path to the JSON configuration file.

        Returns
        -------
        DesignerConfig
            Validated configuration object.

        Raises
        ------
        FileNotFoundError
            If the configuration file does not exist.
        ValidationError
            If the configuration fails validation.
        """
        path = Path(file_path)
        if not path.exists():
            raise FileNotFoundError(f"Configuration file not found: {path}")

        with open(path) as f:
            data = json.load(f)

        return cls.model_validate(data)

    @classmethod
    def from_dict(cls, config_dict: dict[str, Any]) -> DesignerConfig:
        """
        Load configuration from a dictionary.

        Parameters
        ----------
        config_dict : dict
            Dictionary containing configuration parameters.

        Returns
        -------
        DesignerConfig
            Validated configuration object.

        Raises
        ------
        ValidationError
            If the configuration fails validation.
        """
        return cls.model_validate(config_dict)

    @classmethod
    def from_preset(
        cls, preset: Literal["default", "lenient"] = "default"
    ) -> DesignerConfig:
        """
        Load a preset configuration.

        Parameters
        ----------
        preset : str
            Preset name, either "default" or "lenient".

        Returns
        -------
        DesignerConfig
            Validated configuration object.

        Raises
        ------
        ValueError
            If the preset name is not recognized.
        """
        preset_files = {
            "default": "designer_default_config.json",
            "lenient": "designer_lenient_config.json",
        }
        if preset not in preset_files:
            raise ValueError(f"Unknown preset: {preset}. Use 'default' or 'lenient'.")

        data = json.loads(
            files("plexus.data").joinpath(preset_files[preset]).read_text()
        )
        return cls.model_validate(data)

    def to_dict(self) -> dict[str, Any]:
        """
        Export configuration to a dictionary.

        Returns
        -------
        dict
            Configuration as a dictionary.
        """
        return self.model_dump(by_alias=True)

    def to_json_file(self, file_path: str | Path) -> None:
        """
        Save configuration to a JSON file.

        Parameters
        ----------
        file_path : str | Path
            Path to the output JSON file.
        """
        path = Path(file_path)
        with open(path, "w") as f:
            json.dump(self.to_dict(), f, indent=4)
        logger.info(f"Configuration saved to: {path}")


def load_config(
    preset: str = "default",
    config_path: str | None = None,
) -> DesignerConfig:
    """
    Load and validate design configuration.

    Priority order: config_path > preset

    Parameters
    ----------
    preset : str
        Preset configuration name ("default" or "lenient").
    config_path : str | None
        Path to a custom configuration JSON file.

    Returns
    -------
    DesignerConfig
        Validated configuration object.

    Raises
    ------
    ValidationError
        If the configuration fails validation.
    FileNotFoundError
        If the specified config file does not exist.
    """
    if config_path is not None:
        logger.info(f"Loading config from: {config_path}")
        return DesignerConfig.from_json_file(config_path)

    if preset not in ("default", "lenient"):
        logger.warning(
            f"Preset value `{preset}` must be either 'default' or 'lenient'. Using default instead."
        )
        preset = "default"

    logger.info(f"Loading preset configuration: {preset}")
    return DesignerConfig.from_preset(preset)  # type: ignore[arg-type]
