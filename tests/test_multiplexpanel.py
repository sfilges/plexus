# ================================================================================
# Tests for MultiplexPanel and Junction validation
# ================================================================================

import tempfile
from pathlib import Path

import pytest

from plexus.designer.multiplexpanel import Junction, JunctionInput, MultiplexPanel


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
        csv_content = "Name,Chrom,Five_Prime_Coordinate\nJunction1,chr1,100\n"
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

    def test_import_csv_with_panel_column(self):
        """CSV with Panel column imports without error."""
        csv_content = (
            "Name,Chrom,Five_Prime_Coordinate,Three_Prime_Coordinate,Panel\n"
            "J1,chr1,100,200,PanelA\n"
            "J2,chr2,300,400,PanelB\n"
        )
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write(csv_content)
            temp_path = f.name

        try:
            panel = MultiplexPanel("test_panel", "hg38")
            panel.import_junctions_csv(temp_path)
            assert len(panel.junctions) == 2
        finally:
            Path(temp_path).unlink()


class TestJunctionCoordinateCalculation:
    """Tests for calculate_junction_coordinates_in_design_region()."""

    def test_junction_relative_coordinates(self):
        """Junction position is correctly converted to 0-based design region index."""
        # Setup: junction at genomic 1-based position 1000, design_start at 800
        # design_region is 401 bases (genomic 800–1200 inclusive)
        panel = MultiplexPanel("test", "hg38")
        junction = Junction(
            name="J1",
            chrom="chr1",
            start=1000,
            end=1000,
            design_region="A" * 401,
            design_start=800,
            design_end=1200,
        )
        panel.junctions = [junction]
        # No config loaded → default padding of 3bp
        panel.calculate_junction_coordinates_in_design_region()

        # junction.start(1000) - design_start(800) = 200 (0-based index)
        # jmin = 200 - 3 = 197, jmax = 200 + 3 = 203
        assert junction.jmin_coordinate == 197
        assert junction.jmax_coordinate == 203

    def test_junction_relative_coordinates_asymmetric(self):
        """Junction spanning multiple bases produces correct min/max."""
        panel = MultiplexPanel("test", "hg38")
        junction = Junction(
            name="J1",
            chrom="chr1",
            start=1000,
            end=1010,
            design_region="A" * 421,
            design_start=800,
            design_end=1220,
        )
        panel.junctions = [junction]
        panel.calculate_junction_coordinates_in_design_region()

        # five_rel = 1000 - 800 = 200, three_rel = 1010 - 800 = 210
        # jmin = 200 - 3 = 197, jmax = 210 + 3 = 213
        assert junction.jmin_coordinate == 197
        assert junction.jmax_coordinate == 213

    def test_junction_coordinates_clamped_to_bounds(self):
        """Coordinates are clamped to design region boundaries."""
        panel = MultiplexPanel("test", "hg38")
        junction = Junction(
            name="J1",
            chrom="chr1",
            start=802,
            end=802,
            design_region="A" * 10,
            design_start=800,
            design_end=809,
        )
        panel.junctions = [junction]
        panel.calculate_junction_coordinates_in_design_region()

        # five_rel = 802 - 800 = 2, padding 3 → jmin = max(0, 2-3) = 0
        # jmax = min(9, 2+3) = 5
        assert junction.jmin_coordinate == 0
        assert junction.jmax_coordinate == 5


class TestJunctionInputModel:
    """Tests for the JunctionInput Pydantic model."""

    def test_panel_field_optional_absent(self):
        """JunctionInput works without Panel field."""
        j = JunctionInput(
            **{
                "Name": "J1",
                "Chrom": "chr1",
                "Five_Prime_Coordinate": 100,
                "Three_Prime_Coordinate": 200,
            }
        )
        assert j.panel is None

    def test_panel_field_optional_present(self):
        """JunctionInput accepts Panel field."""
        j = JunctionInput(
            **{
                "Name": "J1",
                "Chrom": "chr1",
                "Five_Prime_Coordinate": 100,
                "Three_Prime_Coordinate": 200,
                "Panel": "MyPanel",
            }
        )
        assert j.panel == "MyPanel"
