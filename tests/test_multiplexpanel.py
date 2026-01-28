# ================================================================================
# Tests for MultiplexPanel and Junction validation
# ================================================================================

import tempfile
from pathlib import Path

import pytest

from multiplexdesigner.designer.multiplexpanel import MultiplexPanel


class TestMultiplexPanelImport:
    """Tests for importing and validating junctions in MultiplexPanel."""

    def test_import_valid_csv(self):
        """Test importing a valid CSV file."""
        csv_content = (
            "Name,Chrom,Five_Prime_Coordinate,Three_Prime_Coordinate\n"
            "Junction1,chr1,100,200\n"
            "Junction2,chr2,300,400\n"
        )
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write(csv_content)
            temp_path = f.name

        try:
            panel = MultiplexPanel("test_panel", "hg38")
            panel.import_junctions_csv(temp_path)

            assert len(panel.junctions) == 2
            assert panel.junctions[0].name == "Junction1"
            assert panel.junctions[0].chrom == "chr1"
            assert panel.junctions[0].start == 100
            assert panel.junctions[0].end == 200
        finally:
            Path(temp_path).unlink()

    def test_import_missing_columns(self):
        """Test importing a CSV with missing required columns."""
        # Missing Three_Prime_Coordinate
        csv_content = "Name,Chrom,Five_Prime_Coordinate\n" "Junction1,chr1,100\n"
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write(csv_content)
            temp_path = f.name

        try:
            panel = MultiplexPanel("test_panel", "hg38")
            # Should raise ValueError because validation fails for the missing field
            with pytest.raises(ValueError) as exc_info:
                panel.import_junctions_csv(temp_path)
            assert "Invalid junction data" in str(exc_info.value)
        finally:
            Path(temp_path).unlink()

    def test_import_invalid_types(self):
        """Test importing a CSV with invalid data types (e.g., string for coordinate)."""
        csv_content = (
            "Name,Chrom,Five_Prime_Coordinate,Three_Prime_Coordinate\n"
            "Junction1,chr1,not_a_number,200\n"
        )
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write(csv_content)
            temp_path = f.name

        try:
            panel = MultiplexPanel("test_panel", "hg38")
            with pytest.raises(ValueError) as exc_info:
                panel.import_junctions_csv(temp_path)
            assert "Invalid junction data" in str(exc_info.value)
        finally:
            Path(temp_path).unlink()

    def test_import_extra_columns_ignored(self):
        """Test that extra columns are ignored and do not cause failure."""
        csv_content = (
            "Name,Chrom,Five_Prime_Coordinate,Three_Prime_Coordinate,Extra_Col\n"
            "Junction1,chr1,100,200,some_info\n"
        )
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write(csv_content)
            temp_path = f.name

        try:
            panel = MultiplexPanel("test_panel", "hg38")
            panel.import_junctions_csv(temp_path)
            assert len(panel.junctions) == 1
            assert panel.junctions[0].name == "Junction1"
        finally:
            Path(temp_path).unlink()
