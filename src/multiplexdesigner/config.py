# ================================================================================
# Configuration models for multiplex primer design
#
# Uses Pydantic for runtime validation of configuration parameters.
#
# Author: Stefan Filges (stefan@simsendiagnostics.com)
# Copyright (c) 2025 Simsen Diagnostics AB
# ================================================================================

from __future__ import annotations

import json
from enum import Enum
from pathlib import Path
from typing import Any, Literal

from loguru import logger
from pydantic import BaseModel, Field, model_validator

from multiplexdesigner.utils.root_dir import ROOT_DIR


class SaltCorrectionFormula(str, Enum):
    """Supported salt correction formulas for Tm calculation."""

    SANTALUCIA1998 = "santalucia1998"
    OWCZARZY2004 = "owczarzy2004"
    SCHILDKRAUT1965 = "schildkraut1965"


class ThermodynamicTable(str, Enum):
    """Supported thermodynamic tables for Tm calculation."""

    SANTALUCIA1998 = "santalucia1998"
    BRESLAUER1986 = "breslauer1986"


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

    # Tail sequences (optional adapter sequences)
    five_prime_tail: str = Field(
        default="GGACACTCTTTCCCTACACGACGCTCTTCCGATCTAAAAAAAAAAAAAAAAAAAATGGGAAAGAGTGTCC",
        alias="5_primer_tail",
    )
    three_prime_tail: str = Field(
        default="GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT",
        alias="3_prime_tail",
    )

    # Penalty weights
    primer_length_penalty: float = Field(default=1.0, ge=0.0)
    primer_complexity_penalty: float = Field(default=1.0, ge=0.0)
    amplicon_length_penalty: float = Field(default=1.0, ge=0.0)

    # GC content parameters (%)
    PRIMER_OPT_GC_PERCENT: float = Field(default=50.0, ge=20.0, le=80.0)
    primer_min_gc: int = Field(default=30, ge=0, le=100)
    primer_max_gc: int = Field(default=70, ge=0, le=100)
    primer_gc_clamp: int = Field(default=0, ge=0, le=5)

    # Sequence composition constraints
    primer_max_poly_x: int = Field(default=5, ge=1, le=10)
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

    salt_correction_formula: SaltCorrectionFormula = Field(
        default=SaltCorrectionFormula.SANTALUCIA1998
    )
    thermodynamic_table: ThermodynamicTable = Field(
        default=ThermodynamicTable.SANTALUCIA1998
    )


class MultiplexPickerParameters(BaseModel):
    """Parameters for multiplex panel optimization."""

    initial_solutions: int = Field(default=100, ge=1, le=10000)
    top_solutions_to_keep: int = Field(default=4, ge=1, le=100)

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

    @model_validator(mode="after")
    def validate_plexity_range(self) -> MultiplexPickerParameters:
        """Validate that min <= target <= max plexity."""
        if not (self.minimum_plexity <= self.target_plexity <= self.maximum_plexity):
            raise ValueError(
                f"Plexity values must satisfy: min ({self.minimum_plexity}) <= target ({self.target_plexity}) <= max ({self.maximum_plexity})"
            )
        return self


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
        if preset == "default":
            config_path = Path(ROOT_DIR) / "config" / "designer_default_config.json"
        elif preset == "lenient":
            config_path = Path(ROOT_DIR) / "config" / "designer_lenient_config.json"
        else:
            raise ValueError(f"Unknown preset: {preset}. Use 'default' or 'lenient'.")

        return cls.from_json_file(config_path)

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
    config_dict: dict[str, Any] | None = None,
) -> DesignerConfig:
    """
    Load and validate design configuration.

    Priority order: config_dict > config_path > preset

    Parameters
    ----------
    preset : str
        Preset configuration name ("default" or "lenient").
    config_path : str | None
        Path to a custom configuration JSON file.
    config_dict : dict | None
        Configuration dictionary (highest priority).

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
    if config_dict is not None:
        logger.info("Loading config from dictionary.")
        return DesignerConfig.from_dict(config_dict)

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
