# ================================================================================
# Tests for configuration loading and validation
# ================================================================================

import json
import tempfile
from pathlib import Path

import pytest
from pydantic import ValidationError

from multiplexdesigner.config import (
    DesignerConfig,
    MultiplexPickerParameters,
    PCRConditions,
    PrimerPairParameters,
    SaltCorrectionFormula,
    SingleplexDesignParameters,
    ThermodynamicTable,
    load_config,
)


class TestSingleplexDesignParameters:
    """Tests for SingleplexDesignParameters model."""

    def test_default_values(self):
        """Test that default values are set correctly."""
        params = SingleplexDesignParameters()
        assert params.PRIMER_OPT_TM == 60.0
        assert params.PRIMER_MIN_TM == 57.0
        assert params.PRIMER_MAX_TM == 63.0
        assert params.PRIMER_OPT_SIZE == 22
        assert params.primer_min_length == 18
        assert params.primer_max_length == 28
        assert params.junction_padding_bases == 3

    def test_valid_custom_values(self):
        """Test creating params with valid custom values."""
        params = SingleplexDesignParameters(
            PRIMER_OPT_TM=62.0,
            PRIMER_MIN_TM=59.0,
            PRIMER_MAX_TM=65.0,
            PRIMER_OPT_SIZE=20,
            primer_min_length=18,
            primer_max_length=25,
        )
        assert params.PRIMER_OPT_TM == 62.0
        assert params.PRIMER_MIN_TM == 59.0

    def test_tm_range_validation_min_greater_than_opt(self):
        """Test that min_tm > opt_tm raises ValidationError."""
        with pytest.raises(ValidationError) as exc_info:
            SingleplexDesignParameters(
                PRIMER_MIN_TM=65.0,
                PRIMER_OPT_TM=60.0,
                PRIMER_MAX_TM=68.0,
            )
        assert "Tm values must satisfy" in str(exc_info.value)

    def test_tm_range_validation_opt_greater_than_max(self):
        """Test that opt_tm > max_tm raises ValidationError."""
        with pytest.raises(ValidationError) as exc_info:
            SingleplexDesignParameters(
                PRIMER_MIN_TM=55.0,
                PRIMER_OPT_TM=65.0,
                PRIMER_MAX_TM=62.0,
            )
        assert "Tm values must satisfy" in str(exc_info.value)

    def test_length_range_validation(self):
        """Test that invalid length range raises ValidationError."""
        with pytest.raises(ValidationError) as exc_info:
            SingleplexDesignParameters(
                primer_min_length=25,
                PRIMER_OPT_SIZE=20,
                primer_max_length=28,
            )
        assert "Length values must satisfy" in str(exc_info.value)

    def test_gc_range_validation(self):
        """Test that invalid GC range raises ValidationError."""
        with pytest.raises(ValidationError) as exc_info:
            SingleplexDesignParameters(
                primer_min_gc=60,
                PRIMER_OPT_GC_PERCENT=50.0,
                primer_max_gc=70,
            )
        assert "GC values must satisfy" in str(exc_info.value)

    def test_tm_out_of_bounds(self):
        """Test that Tm values outside allowed bounds raise ValidationError."""
        with pytest.raises(ValidationError):
            SingleplexDesignParameters(PRIMER_OPT_TM=80.0)  # max is 72

    def test_negative_penalty_rejected(self):
        """Test that negative penalty weights are rejected."""
        with pytest.raises(ValidationError):
            SingleplexDesignParameters(primer_length_penalty=-1.0)

    def test_alias_for_tail_sequences(self):
        """Test that alias names work for tail sequences."""
        params = SingleplexDesignParameters.model_validate(
            {"5_primer_tail": "ATCG", "3_prime_tail": "GCTA"}
        )
        assert params.five_prime_tail == "ATCG"
        assert params.three_prime_tail == "GCTA"


class TestPrimerPairParameters:
    """Tests for PrimerPairParameters model."""

    def test_default_values(self):
        """Test that default values are set correctly."""
        params = PrimerPairParameters()
        assert params.PRIMER_PAIR_MAX_DIFF_TM == 3.0
        assert params.PRIMER_PRODUCT_OPT_SIZE == 60
        assert params.PRIMER_PRODUCT_MAX_SIZE == 120

    def test_insert_size_validation(self):
        """Test that min_insert > max_insert raises ValidationError."""
        with pytest.raises(ValidationError) as exc_info:
            PrimerPairParameters(
                PRIMER_PRODUCT_MIN_INSERT_SIZE=100,
                PRIMER_PRODUCT_MAX_INSERT_SIZE=50,
            )
        assert "Insert size" in str(exc_info.value)

    def test_valid_insert_size_range(self):
        """Test valid insert size range."""
        params = PrimerPairParameters(
            PRIMER_PRODUCT_MIN_INSERT_SIZE=30,
            PRIMER_PRODUCT_MAX_INSERT_SIZE=80,
        )
        assert params.PRIMER_PRODUCT_MIN_INSERT_SIZE == 30
        assert params.PRIMER_PRODUCT_MAX_INSERT_SIZE == 80


class TestPCRConditions:
    """Tests for PCRConditions model."""

    def test_default_values(self):
        """Test that default values are set correctly."""
        conditions = PCRConditions()
        assert conditions.annealing_temperature == 60.0
        assert conditions.primer_concentration == 50.0
        assert (
            conditions.salt_correction_formula == SaltCorrectionFormula.SANTALUCIA1998
        )
        assert conditions.thermodynamic_table == ThermodynamicTable.SANTALUCIA1998

    def test_enum_values(self):
        """Test that enum values work correctly."""
        conditions = PCRConditions(
            salt_correction_formula="owczarzy2004",
            thermodynamic_table="breslauer1986",
        )
        assert conditions.salt_correction_formula == SaltCorrectionFormula.OWCZARZY2004
        assert conditions.thermodynamic_table == ThermodynamicTable.BRESLAUER1986

    def test_invalid_enum_value(self):
        """Test that invalid enum values raise ValidationError."""
        with pytest.raises(ValidationError):
            PCRConditions(salt_correction_formula="invalid_formula")

    def test_concentration_bounds(self):
        """Test that concentrations outside bounds raise ValidationError."""
        with pytest.raises(ValidationError):
            PCRConditions(primer_concentration=0.5)  # min is 1.0


class TestMultiplexPickerParameters:
    """Tests for MultiplexPickerParameters model."""

    def test_default_values(self):
        """Test that default values are set correctly."""
        params = MultiplexPickerParameters()
        assert params.target_plexity == 20
        assert params.minimum_plexity == 10
        assert params.maximum_plexity == 50
        assert params.force_plexity is False

    def test_plexity_range_validation(self):
        """Test that invalid plexity range raises ValidationError."""
        with pytest.raises(ValidationError) as exc_info:
            MultiplexPickerParameters(
                minimum_plexity=30,
                target_plexity=20,
                maximum_plexity=50,
            )
        assert "Plexity values must satisfy" in str(exc_info.value)

    def test_valid_plexity_range(self):
        """Test valid plexity range."""
        params = MultiplexPickerParameters(
            minimum_plexity=5,
            target_plexity=15,
            maximum_plexity=25,
        )
        assert params.minimum_plexity == 5
        assert params.target_plexity == 15
        assert params.maximum_plexity == 25


class TestDesignerConfig:
    """Tests for DesignerConfig model."""

    def test_default_config(self):
        """Test creating config with all defaults."""
        config = DesignerConfig()
        assert config.singleplex_design_parameters.PRIMER_OPT_TM == 60.0
        assert config.primer_pair_parameters.PRIMER_PRODUCT_MAX_SIZE == 120
        assert config.pcr_conditions.annealing_temperature == 60.0
        assert config.multiplex_picker_parameters.target_plexity == 20

    def test_from_dict_partial(self):
        """Test loading config from dict with partial values."""
        config = DesignerConfig.from_dict(
            {
                "singleplex_design_parameters": {
                    "PRIMER_OPT_TM": 62.0,
                },
                "pcr_conditions": {
                    "annealing_temperature": 58.0,
                },
            }
        )
        # Custom values
        assert config.singleplex_design_parameters.PRIMER_OPT_TM == 62.0
        assert config.pcr_conditions.annealing_temperature == 58.0
        # Defaults preserved
        assert config.primer_pair_parameters.PRIMER_PRODUCT_MAX_SIZE == 120

    def test_from_dict_invalid(self):
        """Test that invalid dict raises ValidationError."""
        with pytest.raises(ValidationError):
            DesignerConfig.from_dict(
                {
                    "singleplex_design_parameters": {
                        "PRIMER_MIN_TM": 70.0,  # Greater than opt
                        "PRIMER_OPT_TM": 60.0,
                    }
                }
            )

    def test_from_preset_default(self):
        """Test loading default preset."""
        config = DesignerConfig.from_preset("default")
        assert config.singleplex_design_parameters.PRIMER_OPT_TM == 60.0

    def test_from_preset_lenient(self):
        """Test loading lenient preset."""
        config = DesignerConfig.from_preset("lenient")
        assert config is not None

    def test_from_preset_invalid(self):
        """Test that invalid preset raises ValueError."""
        with pytest.raises(ValueError) as exc_info:
            DesignerConfig.from_preset("invalid_preset")
        assert "Unknown preset" in str(exc_info.value)

    def test_from_json_file(self):
        """Test loading config from JSON file."""
        config_data = {
            "singleplex_design_parameters": {
                "PRIMER_OPT_TM": 61.0,
                "PRIMER_MIN_TM": 58.0,
                "PRIMER_MAX_TM": 64.0,
            },
            "primer_pair_parameters": {},
            "pcr_conditions": {},
            "multiplex_picker_parameters": {},
        }

        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            json.dump(config_data, f)
            temp_path = f.name

        try:
            config = DesignerConfig.from_json_file(temp_path)
            assert config.singleplex_design_parameters.PRIMER_OPT_TM == 61.0
        finally:
            Path(temp_path).unlink()

    def test_from_json_file_not_found(self):
        """Test that missing file raises FileNotFoundError."""
        with pytest.raises(FileNotFoundError):
            DesignerConfig.from_json_file("/nonexistent/path/config.json")

    def test_to_dict(self):
        """Test exporting config to dict."""
        config = DesignerConfig()
        data = config.to_dict()

        assert "singleplex_design_parameters" in data
        assert "primer_pair_parameters" in data
        assert "pcr_conditions" in data
        assert "multiplex_picker_parameters" in data
        assert data["singleplex_design_parameters"]["PRIMER_OPT_TM"] == 60.0

    def test_to_dict_uses_aliases(self):
        """Test that to_dict uses field aliases."""
        config = DesignerConfig()
        data = config.to_dict()

        # Should use alias "5_primer_tail" not "five_prime_tail"
        assert "5_primer_tail" in data["singleplex_design_parameters"]

    def test_to_json_file(self):
        """Test saving config to JSON file."""
        config = DesignerConfig()

        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            temp_path = f.name

        try:
            config.to_json_file(temp_path)

            with open(temp_path) as f:
                data = json.load(f)

            assert data["singleplex_design_parameters"]["PRIMER_OPT_TM"] == 60.0
        finally:
            Path(temp_path).unlink()

    def test_roundtrip(self):
        """Test that config survives roundtrip through dict."""
        original = DesignerConfig()
        data = original.to_dict()
        restored = DesignerConfig.from_dict(data)

        assert (
            original.singleplex_design_parameters.PRIMER_OPT_TM
            == restored.singleplex_design_parameters.PRIMER_OPT_TM
        )
        assert (
            original.primer_pair_parameters.PRIMER_PRODUCT_MAX_SIZE
            == restored.primer_pair_parameters.PRIMER_PRODUCT_MAX_SIZE
        )


class TestLoadConfig:
    """Tests for the load_config convenience function."""

    def test_load_preset_default(self):
        """Test loading default preset via load_config."""
        config = load_config(preset="default")
        assert isinstance(config, DesignerConfig)
        assert config.singleplex_design_parameters.PRIMER_OPT_TM == 60.0

    def test_load_preset_lenient(self):
        """Test loading lenient preset via load_config."""
        config = load_config(preset="lenient")
        assert isinstance(config, DesignerConfig)

    def test_load_invalid_preset_falls_back(self):
        """Test that invalid preset falls back to default."""
        config = load_config(preset="invalid")
        assert isinstance(config, DesignerConfig)

    def test_load_from_dict(self):
        """Test loading from dict via load_config."""
        config = load_config(
            config_dict={
                "singleplex_design_parameters": {
                    "PRIMER_OPT_TM": 62.0,
                }
            }
        )
        assert config.singleplex_design_parameters.PRIMER_OPT_TM == 62.0

    def test_load_from_file(self):
        """Test loading from file via load_config."""
        config_data = {
            "singleplex_design_parameters": {"PRIMER_OPT_TM": 63.0},
        }

        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            json.dump(config_data, f)
            temp_path = f.name

        try:
            config = load_config(config_path=temp_path)
            assert config.singleplex_design_parameters.PRIMER_OPT_TM == 63.0
        finally:
            Path(temp_path).unlink()

    def test_priority_dict_over_file(self):
        """Test that config_dict takes priority over config_path."""
        config_data = {"singleplex_design_parameters": {"PRIMER_OPT_TM": 59.0}}

        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            json.dump(config_data, f)
            temp_path = f.name

        try:
            config = load_config(
                config_path=temp_path,
                config_dict={"singleplex_design_parameters": {"PRIMER_OPT_TM": 61.0}},
            )
            # Dict value should win
            assert config.singleplex_design_parameters.PRIMER_OPT_TM == 61.0
        finally:
            Path(temp_path).unlink()

    def test_priority_file_over_preset(self):
        """Test that config_path takes priority over preset."""
        config_data = {"singleplex_design_parameters": {"PRIMER_OPT_TM": 63.0}}

        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            json.dump(config_data, f)
            temp_path = f.name

        try:
            config = load_config(preset="default", config_path=temp_path)
            # File value should win over preset
            assert config.singleplex_design_parameters.PRIMER_OPT_TM == 63.0
        finally:
            Path(temp_path).unlink()
