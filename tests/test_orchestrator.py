# ================================================================================
# Tests for multi-panel orchestrator
# ================================================================================

import json
from pathlib import Path
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest

from plexus.orchestrator import (
    DEFAULT_PANEL_ID,
    _detect_panels,
    _write_panel_csv,
    run_multi_panel,
)
from plexus.pipeline import MultiPanelResult, PipelineResult

# ================================================================================
# Tests for _detect_panels
# ================================================================================


class TestDetectPanels:
    """Tests for _detect_panels()."""

    def test_no_panel_column_returns_single_default(self, tmp_path):
        """CSV without Panel column returns single default panel."""
        csv_path = tmp_path / "junctions.csv"
        csv_path.write_text(
            "Name,Chrom,Five_Prime_Coordinate,Three_Prime_Coordinate\n"
            "J1,chr1,100,200\n"
            "J2,chr2,300,400\n"
        )

        panels = _detect_panels(csv_path)
        assert list(panels.keys()) == [DEFAULT_PANEL_ID]
        assert len(panels[DEFAULT_PANEL_ID]) == 2

    def test_panel_column_groups_correctly(self, tmp_path):
        """CSV with Panel column groups junctions by panel."""
        csv_path = tmp_path / "junctions.csv"
        csv_path.write_text(
            "Name,Chrom,Five_Prime_Coordinate,Three_Prime_Coordinate,Panel\n"
            "J1,chr1,100,200,PanelA\n"
            "J2,chr2,300,400,PanelA\n"
            "J3,chr3,500,600,PanelB\n"
        )

        panels = _detect_panels(csv_path)
        assert sorted(panels.keys()) == ["PanelA", "PanelB"]
        assert len(panels["PanelA"]) == 2
        assert len(panels["PanelB"]) == 1

    def test_empty_panel_values_raises(self, tmp_path):
        """Panel column with NaN raises ValueError."""
        csv_path = tmp_path / "junctions.csv"
        csv_path.write_text(
            "Name,Chrom,Five_Prime_Coordinate,Three_Prime_Coordinate,Panel\n"
            "J1,chr1,100,200,PanelA\n"
            "J2,chr2,300,400,\n"
        )

        with pytest.raises(ValueError, match="empty values"):
            _detect_panels(csv_path)

    def test_single_panel_value(self, tmp_path):
        """All rows with same Panel value returns one group (not default)."""
        csv_path = tmp_path / "junctions.csv"
        csv_path.write_text(
            "Name,Chrom,Five_Prime_Coordinate,Three_Prime_Coordinate,Panel\n"
            "J1,chr1,100,200,MyPanel\n"
            "J2,chr2,300,400,MyPanel\n"
        )

        panels = _detect_panels(csv_path)
        assert list(panels.keys()) == ["MyPanel"]
        assert len(panels["MyPanel"]) == 2


# ================================================================================
# Tests for _write_panel_csv
# ================================================================================


class TestWritePanelCsv:
    """Tests for _write_panel_csv()."""

    def test_panel_column_stripped(self, tmp_path):
        """Written CSV should not contain the Panel column."""
        df = pd.DataFrame(
            {
                "Name": ["J1"],
                "Chrom": ["chr1"],
                "Five_Prime_Coordinate": [100],
                "Three_Prime_Coordinate": [200],
                "Panel": ["PanelA"],
            }
        )

        out_path = _write_panel_csv(df, tmp_path, "PanelA")
        result_df = pd.read_csv(out_path)

        assert "Panel" not in result_df.columns
        assert list(result_df.columns) == [
            "Name",
            "Chrom",
            "Five_Prime_Coordinate",
            "Three_Prime_Coordinate",
        ]

    def test_junction_columns_preserved(self, tmp_path):
        """All junction data is preserved in the output CSV."""
        df = pd.DataFrame(
            {
                "Name": ["J1", "J2"],
                "Chrom": ["chr1", "chr2"],
                "Five_Prime_Coordinate": [100, 300],
                "Three_Prime_Coordinate": [200, 400],
                "Panel": ["A", "A"],
            }
        )

        out_path = _write_panel_csv(df, tmp_path, "A")
        result_df = pd.read_csv(out_path)

        assert len(result_df) == 2
        assert result_df.iloc[0]["Name"] == "J1"
        assert result_df.iloc[1]["Five_Prime_Coordinate"] == 300

    def test_filename_uses_panel_id(self, tmp_path):
        """Output filename includes the panel ID."""
        df = pd.DataFrame(
            {
                "Name": ["J1"],
                "Chrom": ["chr1"],
                "Five_Prime_Coordinate": [100],
                "Three_Prime_Coordinate": [200],
            }
        )

        out_path = _write_panel_csv(df, tmp_path, "my_panel")
        assert out_path.name == "my_panel_junctions.csv"


# ================================================================================
# Tests for run_multi_panel
# ================================================================================


class TestRunMultiPanel:
    """Integration tests for run_multi_panel()."""

    def _make_mock_result(self, output_dir="/tmp/output"):
        """Create a mock PipelineResult."""
        mock_panel = MagicMock()
        mock_panel.junctions = [MagicMock()]
        mock_panel.junctions[0].primer_pairs = [MagicMock()]

        return PipelineResult(
            panel=mock_panel,
            output_dir=Path(output_dir),
            config=MagicMock(),
            steps_completed=["panel_created", "primers_designed"],
            errors=[],
        )

    @patch("plexus.orchestrator.run_pipeline")
    def test_single_panel_returns_pipeline_result(self, mock_run, tmp_path):
        """No Panel column returns PipelineResult (not MultiPanelResult)."""
        csv_path = tmp_path / "junctions.csv"
        csv_path.write_text(
            "Name,Chrom,Five_Prime_Coordinate,Three_Prime_Coordinate\nJ1,chr1,100,200\n"
        )
        fasta_path = tmp_path / "genome.fa"
        fasta_path.write_text(">chr1\nACGT\n")

        mock_run.return_value = self._make_mock_result()

        result = run_multi_panel(
            input_file=csv_path,
            fasta_file=fasta_path,
            output_dir=tmp_path / "output",
        )

        assert isinstance(result, PipelineResult)
        mock_run.assert_called_once()

    @patch("plexus.orchestrator.run_pipeline")
    def test_multi_panel_returns_multi_panel_result(self, mock_run, tmp_path):
        """Panel column returns MultiPanelResult."""
        csv_path = tmp_path / "junctions.csv"
        csv_path.write_text(
            "Name,Chrom,Five_Prime_Coordinate,Three_Prime_Coordinate,Panel\n"
            "J1,chr1,100,200,PanelA\n"
            "J2,chr2,300,400,PanelB\n"
        )
        fasta_path = tmp_path / "genome.fa"
        fasta_path.write_text(">chr1\nACGT\n")

        mock_run.return_value = self._make_mock_result()

        result = run_multi_panel(
            input_file=csv_path,
            fasta_file=fasta_path,
            output_dir=tmp_path / "output",
        )

        assert isinstance(result, MultiPanelResult)
        assert sorted(result.panel_ids) == ["PanelA", "PanelB"]
        assert len(result.panel_results) == 2
        assert mock_run.call_count == 2

    @patch("plexus.orchestrator.run_pipeline")
    def test_output_dirs_created_per_panel(self, mock_run, tmp_path):
        """Each panel gets output_dir/panel_id/ passed to run_pipeline."""
        csv_path = tmp_path / "junctions.csv"
        csv_path.write_text(
            "Name,Chrom,Five_Prime_Coordinate,Three_Prime_Coordinate,Panel\n"
            "J1,chr1,100,200,A\n"
            "J2,chr2,300,400,B\n"
        )
        fasta_path = tmp_path / "genome.fa"
        fasta_path.write_text(">chr1\nACGT\n")

        mock_run.return_value = self._make_mock_result()
        output_dir = tmp_path / "output"

        run_multi_panel(
            input_file=csv_path,
            fasta_file=fasta_path,
            output_dir=output_dir,
        )

        # Check run_pipeline was called with per-panel output dirs
        call_output_dirs = sorted(
            str(call[1]["output_dir"]) for call in mock_run.call_args_list
        )
        assert str(output_dir / "A") in call_output_dirs
        assert str(output_dir / "B") in call_output_dirs

    @patch("plexus.orchestrator.run_pipeline")
    def test_multi_panel_summary_json_written(self, mock_run, tmp_path):
        """Multi-panel mode writes multi_panel_summary.json."""
        csv_path = tmp_path / "junctions.csv"
        csv_path.write_text(
            "Name,Chrom,Five_Prime_Coordinate,Three_Prime_Coordinate,Panel\n"
            "J1,chr1,100,200,A\n"
            "J2,chr2,300,400,B\n"
        )
        fasta_path = tmp_path / "genome.fa"
        fasta_path.write_text(">chr1\nACGT\n")

        mock_run.return_value = self._make_mock_result()
        output_dir = tmp_path / "output"

        run_multi_panel(
            input_file=csv_path,
            fasta_file=fasta_path,
            output_dir=output_dir,
        )

        summary_path = output_dir / "multi_panel_summary.json"
        assert summary_path.exists()

        with open(summary_path) as f:
            summary = json.load(f)

        assert summary["num_panels"] == 2
        assert sorted(summary["panel_ids"]) == ["A", "B"]
        assert summary["all_successful"] is True

    @patch("plexus.orchestrator.run_pipeline")
    def test_backward_compat_single_panel_flat_output(self, mock_run, tmp_path):
        """Single-panel CSV produces flat output_dir/ (no subfolder)."""
        csv_path = tmp_path / "junctions.csv"
        csv_path.write_text(
            "Name,Chrom,Five_Prime_Coordinate,Three_Prime_Coordinate\nJ1,chr1,100,200\n"
        )
        fasta_path = tmp_path / "genome.fa"
        fasta_path.write_text(">chr1\nACGT\n")

        mock_run.return_value = self._make_mock_result()
        output_dir = tmp_path / "output"

        run_multi_panel(
            input_file=csv_path,
            fasta_file=fasta_path,
            output_dir=output_dir,
        )

        # run_pipeline should be called with the flat output_dir, not a subfolder
        call_kwargs = mock_run.call_args[1]
        assert call_kwargs["output_dir"] == output_dir

    @patch("plexus.orchestrator.run_pipeline")
    def test_temp_csvs_cleaned_up(self, mock_run, tmp_path):
        """Temp panel CSVs are cleaned up after execution."""
        csv_path = tmp_path / "junctions.csv"
        csv_path.write_text(
            "Name,Chrom,Five_Prime_Coordinate,Three_Prime_Coordinate,Panel\n"
            "J1,chr1,100,200,A\n"
            "J2,chr2,300,400,B\n"
        )
        fasta_path = tmp_path / "genome.fa"
        fasta_path.write_text(">chr1\nACGT\n")

        mock_run.return_value = self._make_mock_result()
        output_dir = tmp_path / "output"

        run_multi_panel(
            input_file=csv_path,
            fasta_file=fasta_path,
            output_dir=output_dir,
        )

        tmp_dir = output_dir / ".tmp_panel_csvs"
        assert not tmp_dir.exists()


# ================================================================================
# Tests for MultiPanelResult
# ================================================================================


class TestMultiPanelResult:
    """Tests for the MultiPanelResult dataclass."""

    def _make_result(self, errors=None):
        mock_panel = MagicMock()
        mock_panel.junctions = [MagicMock()]
        mock_panel.junctions[0].primer_pairs = [MagicMock()]

        return PipelineResult(
            panel=mock_panel,
            output_dir=Path("/tmp/output"),
            config=MagicMock(),
            steps_completed=["panel_created"],
            errors=errors or [],
            selected_pairs=[MagicMock()],
        )

    def test_success_all_panels_ok(self):
        """success is True when all panels have no errors."""
        result = MultiPanelResult(
            panel_results={"A": self._make_result(), "B": self._make_result()},
            output_dir=Path("/tmp"),
            panel_ids=["A", "B"],
        )
        assert result.success is True
        assert result.failed_panels == []

    def test_success_false_on_error(self):
        """success is False when any panel has errors."""
        result = MultiPanelResult(
            panel_results={
                "A": self._make_result(),
                "B": self._make_result(errors=["something failed"]),
            },
            output_dir=Path("/tmp"),
            panel_ids=["A", "B"],
        )
        assert result.success is False
        assert result.failed_panels == ["B"]

    def test_total_junctions(self):
        """total_junctions sums across all panels."""
        result = MultiPanelResult(
            panel_results={"A": self._make_result(), "B": self._make_result()},
            output_dir=Path("/tmp"),
            panel_ids=["A", "B"],
        )
        assert result.total_junctions == 2  # 1 junction per panel

    def test_summary_dict(self):
        """summary_dict returns expected structure."""
        result = MultiPanelResult(
            panel_results={"A": self._make_result()},
            output_dir=Path("/tmp"),
            panel_ids=["A"],
        )
        summary = result.summary_dict()

        assert summary["num_panels"] == 1
        assert summary["panel_ids"] == ["A"]
        assert summary["all_successful"] is True
        assert "A" in summary["per_panel"]
        assert summary["per_panel"]["A"]["success"] is True
