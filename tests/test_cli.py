# ================================================================================
# Tests for the command-line interface
# ================================================================================

import json
import re
import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

from typer.testing import CliRunner

from plexus.cli import app
from plexus.pipeline import PipelineResult
from plexus.version import __version__

runner = CliRunner()


def strip_ansi(text: str) -> str:
    """Strip ANSI escape sequences from a string."""
    ansi_escape = re.compile(r"\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])")
    return ansi_escape.sub("", text)


class TestCLIBasics:
    """Tests for basic CLI functionality."""

    def test_help(self):
        """Test that --help works."""
        result = runner.invoke(app, ["--help"], env={"COLUMNS": "120"})
        output = strip_ansi(result.output)
        assert result.exit_code == 0
        assert "Design multiplex PCR primer panels" in output

    def test_version(self):
        """Test that --version shows version."""
        result = runner.invoke(app, ["--version"])
        assert result.exit_code == 0
        assert __version__ in result.output

    def test_no_args_shows_help(self):
        """Test that no arguments shows help/usage."""
        result = runner.invoke(app, [], env={"COLUMNS": "120"})
        output = strip_ansi(result.output)
        # Typer shows usage when no args provided
        assert "Usage:" in output
        assert "COMMAND" in output


class TestRunCommand:
    """Tests for the run command."""

    def test_run_help(self):
        """Test that run --help works."""
        result = runner.invoke(app, ["run", "--help"], env={"COLUMNS": "120"})
        output = strip_ansi(result.output)
        assert result.exit_code == 0
        assert "Run the complete multiplex primer design pipeline" in output
        assert "--input" in output
        assert "--fasta" in output
        assert "--output" in output

    def test_run_missing_input(self):
        """Test that run fails without --input."""
        result = runner.invoke(app, ["run", "--fasta", "genome.fa"])
        assert (
            result.exit_code == 2
        )  # Typer uses exit code 2 for missing required options
        assert "Missing option" in result.output or "--input" in result.output

    def test_run_missing_fasta(self):
        """Test that run fails without --fasta."""
        result = runner.invoke(app, ["run", "--input", "junctions.csv"])
        assert result.exit_code == 2
        assert "Missing option" in result.output or "--fasta" in result.output

    def test_run_input_file_not_found(self):
        """Test that run fails with non-existent input file."""
        with tempfile.NamedTemporaryFile(suffix=".fa", delete=False) as f:
            fasta_path = f.name
            f.write(b">chr1\nACGT\n")

        try:
            result = runner.invoke(
                app,
                [
                    "run",
                    "--input",
                    "/nonexistent/junctions.csv",
                    "--fasta",
                    fasta_path,
                ],
            )
            assert result.exit_code == 1
            assert (
                "not found" in result.output.lower() or "error" in result.output.lower()
            )
        finally:
            Path(fasta_path).unlink()

    def test_run_fasta_file_not_found(self):
        """Test that run fails with non-existent FASTA file."""
        with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as f:
            csv_path = f.name
            f.write(b"Name,Chrom,Position\n")

        try:
            result = runner.invoke(
                app,
                [
                    "run",
                    "--input",
                    csv_path,
                    "--fasta",
                    "/nonexistent/genome.fa",
                ],
            )
            assert result.exit_code == 1
            assert (
                "not found" in result.output.lower() or "error" in result.output.lower()
            )
        finally:
            Path(csv_path).unlink()

    @patch("plexus.orchestrator.run_pipeline")
    def test_run_success(self, mock_run_pipeline):
        """Test successful run with mocked pipeline."""
        # Create mock result
        mock_panel = MagicMock()
        mock_panel.junctions = [MagicMock(), MagicMock()]
        for j in mock_panel.junctions:
            j.primer_pairs = [MagicMock(), MagicMock(), MagicMock()]

        mock_result = PipelineResult(
            panel=mock_panel,
            output_dir=Path("/tmp/output"),
            config=MagicMock(),
            steps_completed=["panel_created", "primers_designed", "candidates_saved"],
            errors=[],
        )
        mock_run_pipeline.return_value = mock_result

        # Create temp files
        with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as csv_f:
            csv_path = csv_f.name
            csv_f.write(b"Name,Chrom,Position\ntest,chr1,1000\n")

        with tempfile.NamedTemporaryFile(suffix=".fa", delete=False) as fa_f:
            fasta_path = fa_f.name
            fa_f.write(b">chr1\nACGTACGTACGT\n")

        try:
            result = runner.invoke(
                app,
                [
                    "run",
                    "--input",
                    csv_path,
                    "--fasta",
                    fasta_path,
                    "--output",
                    "/tmp/test_output",
                ],
            )

            assert result.exit_code == 0
            assert "completed successfully" in result.output.lower()
            assert "Junctions:" in result.output
            assert "Primer pairs:" in result.output

            # Verify pipeline was called with correct arguments
            mock_run_pipeline.assert_called_once()
            call_kwargs = mock_run_pipeline.call_args[1]
            assert str(call_kwargs["input_file"]) == csv_path
            assert str(call_kwargs["fasta_file"]) == fasta_path
            assert str(call_kwargs["output_dir"]) == "/tmp/test_output"
        finally:
            Path(csv_path).unlink()
            Path(fasta_path).unlink()

    @patch("plexus.orchestrator.run_pipeline")
    def test_run_with_warnings(self, mock_run_pipeline):
        """Test run that completes with warnings."""
        mock_panel = MagicMock()
        mock_panel.junctions = [MagicMock()]
        mock_panel.junctions[0].primer_pairs = []

        mock_result = PipelineResult(
            panel=mock_panel,
            output_dir=Path("/tmp/output"),
            config=MagicMock(),
            steps_completed=["panel_created", "primers_designed"],
            errors=["Save candidates failed: test error"],
        )
        mock_run_pipeline.return_value = mock_result

        with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as csv_f:
            csv_path = csv_f.name
            csv_f.write(b"Name,Chrom,Position\n")

        with tempfile.NamedTemporaryFile(suffix=".fa", delete=False) as fa_f:
            fasta_path = fa_f.name
            fa_f.write(b">chr1\nACGT\n")

        try:
            result = runner.invoke(
                app,
                ["run", "--input", csv_path, "--fasta", fasta_path],
            )

            assert result.exit_code == 0
            assert "warnings" in result.output.lower()
            assert "Save candidates failed" in result.output
        finally:
            Path(csv_path).unlink()
            Path(fasta_path).unlink()

    @patch("plexus.orchestrator.run_pipeline")
    def test_run_with_all_options(self, mock_run_pipeline):
        """Test run with all CLI options specified."""
        mock_panel = MagicMock()
        mock_panel.junctions = []

        mock_result = PipelineResult(
            panel=mock_panel,
            output_dir=Path("/tmp/output"),
            config=MagicMock(),
            steps_completed=["panel_created"],
            errors=[],
        )
        mock_run_pipeline.return_value = mock_result

        with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as csv_f:
            csv_path = csv_f.name
            csv_f.write(b"Name,Chrom,Position\n")

        with tempfile.NamedTemporaryFile(suffix=".fa", delete=False) as fa_f:
            fasta_path = fa_f.name
            fa_f.write(b">chr1\nACGT\n")

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".json", delete=False
        ) as cfg_f:
            config_path = cfg_f.name
            json.dump({"singleplex_design_parameters": {}}, cfg_f)

        try:
            result = runner.invoke(
                app,
                [
                    "run",
                    "--input",
                    csv_path,
                    "--fasta",
                    fasta_path,
                    "--output",
                    "/tmp/custom_output",
                    "--name",
                    "my_panel",
                    "--genome",
                    "hg19",
                    "--preset",
                    "lenient",
                    "--config",
                    config_path,
                    "--skip-blast",
                    "--padding",
                    "300",
                ],
            )

            assert result.exit_code == 0

            # Verify all options were passed correctly
            call_kwargs = mock_run_pipeline.call_args[1]
            assert call_kwargs["panel_name"] == "my_panel"
            assert call_kwargs["genome"] == "hg19"
            assert call_kwargs["preset"] == "lenient"
            assert str(call_kwargs["config_path"]) == config_path
            assert call_kwargs["run_blast"] is False
            assert call_kwargs["padding"] == 300
        finally:
            Path(csv_path).unlink()
            Path(fasta_path).unlink()
            Path(config_path).unlink()

    @patch("plexus.orchestrator.run_pipeline")
    def test_run_skip_blast(self, mock_run_pipeline):
        """Test that --skip-blast flag is passed correctly."""
        mock_panel = MagicMock()
        mock_panel.junctions = []

        mock_result = PipelineResult(
            panel=mock_panel,
            output_dir=Path("/tmp/output"),
            config=MagicMock(),
            steps_completed=[],
            errors=[],
        )
        mock_run_pipeline.return_value = mock_result

        with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as csv_f:
            csv_path = csv_f.name
            csv_f.write(b"Name,Chrom,Position\n")

        with tempfile.NamedTemporaryFile(suffix=".fa", delete=False) as fa_f:
            fasta_path = fa_f.name
            fa_f.write(b">chr1\nACGT\n")

        try:
            # Without --skip-blast
            runner.invoke(
                app,
                ["run", "--input", csv_path, "--fasta", fasta_path],
            )
            assert mock_run_pipeline.call_args[1]["run_blast"] is True

            mock_run_pipeline.reset_mock()

            # With --skip-blast
            runner.invoke(
                app,
                ["run", "--input", csv_path, "--fasta", fasta_path, "--skip-blast"],
            )
            assert mock_run_pipeline.call_args[1]["run_blast"] is False
        finally:
            Path(csv_path).unlink()
            Path(fasta_path).unlink()

    @patch("plexus.orchestrator.run_pipeline")
    def test_run_pipeline_exception(self, mock_run_pipeline):
        """Test that pipeline exceptions are handled gracefully."""
        mock_run_pipeline.side_effect = Exception("Pipeline crashed!")

        with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as csv_f:
            csv_path = csv_f.name
            csv_f.write(b"Name,Chrom,Position\n")

        with tempfile.NamedTemporaryFile(suffix=".fa", delete=False) as fa_f:
            fasta_path = fa_f.name
            fa_f.write(b">chr1\nACGT\n")

        try:
            result = runner.invoke(
                app,
                ["run", "--input", csv_path, "--fasta", fasta_path],
            )

            assert result.exit_code == 1
            assert (
                "failed" in result.output.lower() or "crashed" in result.output.lower()
            )
        finally:
            Path(csv_path).unlink()
            Path(fasta_path).unlink()


class TestCLIShortFlags:
    """Tests for short flag aliases."""

    @patch("plexus.orchestrator.run_pipeline")
    def test_short_flags(self, mock_run_pipeline):
        """Test that short flags work correctly."""
        mock_panel = MagicMock()
        mock_panel.junctions = []

        mock_result = PipelineResult(
            panel=mock_panel,
            output_dir=Path("/tmp/output"),
            config=MagicMock(),
            steps_completed=[],
            errors=[],
        )
        mock_run_pipeline.return_value = mock_result

        with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as csv_f:
            csv_path = csv_f.name
            csv_f.write(b"Name,Chrom,Position\n")

        with tempfile.NamedTemporaryFile(suffix=".fa", delete=False) as fa_f:
            fasta_path = fa_f.name
            fa_f.write(b">chr1\nACGT\n")

        try:
            result = runner.invoke(
                app,
                [
                    "run",
                    "-i",
                    csv_path,
                    "-f",
                    fasta_path,
                    "-o",
                    "/tmp/out",
                    "-n",
                    "test",
                    "-g",
                    "hg38",
                    "-p",
                    "default",
                ],
            )

            assert result.exit_code == 0

            call_kwargs = mock_run_pipeline.call_args[1]
            assert str(call_kwargs["input_file"]) == csv_path
            assert str(call_kwargs["fasta_file"]) == fasta_path
            assert str(call_kwargs["output_dir"]) == "/tmp/out"
            assert call_kwargs["panel_name"] == "test"
            assert call_kwargs["genome"] == "hg38"
            assert call_kwargs["preset"] == "default"
        finally:
            Path(csv_path).unlink()
            Path(fasta_path).unlink()
