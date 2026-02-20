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

    def test_run_missing_fasta_no_registry(self):
        """Test that run fails without --fasta when genome is not in registry."""
        from unittest.mock import patch

        with patch("plexus.resources.get_registered_fasta", return_value=None):
            result = runner.invoke(app, ["run", "--input", "junctions.csv"])
        assert result.exit_code == 1
        assert "not initialized" in result.output or "plexus init" in result.output

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


class TestInitCommand:
    """Tests for the init command."""

    def test_init_help(self):
        """Test that init --help works."""
        result = runner.invoke(app, ["init", "--help"], env={"COLUMNS": "120"})
        output = strip_ansi(result.output)
        assert result.exit_code == 0
        assert "--genome" in output
        assert "--fasta" in output
        assert "--skip-blast" in output

    def test_init_unknown_genome(self):
        """Test that init rejects an unknown genome."""
        result = runner.invoke(app, ["init", "--genome", "mm10_nonexistent"])
        assert result.exit_code == 1
        assert "Unknown genome" in (result.output + (result.stderr or ""))

    @patch("plexus.resources.init_genome")
    @patch("plexus.resources.genome_status")
    @patch("plexus.resources.get_operational_mode", return_value="research")
    def test_init_with_local_fasta_and_vcf(self, _mock_mode, mock_status, mock_init):
        """Test init with local FASTA and VCF skips downloads."""
        mock_status.return_value = {
            "fasta": True,
            "fai": True,
            "blast_db": True,
            "snp_vcf": True,
            "fasta_sha256": "a" * 64,
            "snp_vcf_sha256": "b" * 64,
        }

        with tempfile.NamedTemporaryFile(suffix=".fa", delete=False) as fa_f:
            fasta_path = fa_f.name
            fa_f.write(b">chr1\nACGT\n")

        with tempfile.NamedTemporaryFile(suffix=".vcf.gz", delete=False) as vcf_f:
            vcf_path = vcf_f.name

        try:
            result = runner.invoke(
                app,
                [
                    "init",
                    "--genome",
                    "hg38",
                    "--fasta",
                    fasta_path,
                    "--snp-vcf",
                    vcf_path,
                ],
            )
            assert result.exit_code == 0
            assert "Done!" in strip_ansi(result.output)
            mock_init.assert_called_once()
            call_kwargs = mock_init.call_args[1]
            assert call_kwargs["genome"] == "hg38"
            assert str(call_kwargs["fasta"]) == fasta_path
        finally:
            Path(fasta_path).unlink()
            Path(vcf_path).unlink()

    @patch("plexus.resources.init_genome")
    @patch("plexus.resources.genome_status")
    @patch("plexus.resources.get_operational_mode", return_value="research")
    def test_init_skip_flags(self, _mock_mode, mock_status, mock_init):
        """Test that --skip-blast and --skip-snp are passed through."""
        mock_status.return_value = {
            "fasta": True,
            "fai": True,
            "blast_db": False,
            "snp_vcf": False,
            "fasta_sha256": None,
            "snp_vcf_sha256": None,
        }

        with tempfile.NamedTemporaryFile(suffix=".fa", delete=False) as fa_f:
            fasta_path = fa_f.name
            fa_f.write(b">chr1\nACGT\n")

        try:
            result = runner.invoke(
                app,
                [
                    "init",
                    "--genome",
                    "hg38",
                    "--fasta",
                    fasta_path,
                    "--skip-blast",
                    "--skip-snp",
                ],
            )
            assert result.exit_code == 0
            call_kwargs = mock_init.call_args[1]
            assert call_kwargs["skip_blast"] is True
            assert call_kwargs["skip_snp"] is True
        finally:
            Path(fasta_path).unlink()


class TestRunRegistryFallback:
    """Tests for the --fasta registry fallback in run."""

    @patch("plexus.orchestrator.run_pipeline")
    @patch("plexus.resources.get_registered_fasta")
    def test_run_uses_registered_fasta(self, mock_get_fasta, mock_run_pipeline):
        """Test that run uses the registered FASTA when --fasta is omitted."""
        import tempfile

        with tempfile.NamedTemporaryFile(suffix=".fa", delete=False) as fa_f:
            fasta_path = Path(fa_f.name)
            fa_f.write(b">chr1\nACGT\n")

        mock_get_fasta.return_value = fasta_path

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

        try:
            result = runner.invoke(app, ["run", "--input", csv_path])
            assert result.exit_code == 0
            mock_get_fasta.assert_called_once_with("hg38")
            call_kwargs = mock_run_pipeline.call_args[1]
            assert call_kwargs["fasta_file"] == fasta_path
        finally:
            fasta_path.unlink()
            Path(csv_path).unlink()

    @patch("plexus.resources.get_registered_fasta", return_value=None)
    def test_run_no_fasta_no_registry_exits(self, _mock):
        """Test run exits with error when --fasta omitted and genome not registered."""
        with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as csv_f:
            csv_path = csv_f.name
            csv_f.write(b"Name,Chrom,Position\n")
        try:
            result = runner.invoke(app, ["run", "--input", csv_path])
            assert result.exit_code == 1
            combined = strip_ansi(result.output)
            assert "not initialized" in combined or "plexus init" in combined
        finally:
            Path(csv_path).unlink()


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


class TestInitNewFlags:
    """Tests for --mode, --checksums flags on init."""

    def test_init_requires_fasta(self):
        """init without --fasta should error."""
        result = runner.invoke(app, ["init", "--genome", "hg38", "--skip-snp"])
        assert result.exit_code == 1
        assert "--fasta" in (result.output + (result.stderr or ""))

    def test_init_requires_snp_vcf(self):
        """init without --snp-vcf should error (unless --skip-snp)."""
        with tempfile.NamedTemporaryFile(suffix=".fa", delete=False) as fa_f:
            fasta_path = fa_f.name
            fa_f.write(b">chr1\nACGT\n")
        try:
            result = runner.invoke(
                app, ["init", "--genome", "hg38", "--fasta", fasta_path]
            )
            assert result.exit_code == 1
            assert "--snp-vcf" in (result.output + (result.stderr or ""))
        finally:
            Path(fasta_path).unlink()

    @patch("plexus.resources.init_genome")
    @patch("plexus.resources.genome_status")
    @patch("plexus.resources.get_operational_mode", return_value="compliance")
    def test_init_with_mode_compliance(self, _mock_mode, mock_status, mock_init):
        """init --mode compliance passes through to init_genome."""
        mock_status.return_value = {
            "fasta": True,
            "fai": True,
            "blast_db": True,
            "snp_vcf": True,
            "fasta_sha256": "a" * 64,
            "snp_vcf_sha256": None,
        }
        with tempfile.NamedTemporaryFile(suffix=".fa", delete=False) as fa_f:
            fasta_path = fa_f.name
            fa_f.write(b">chr1\nACGT\n")
        try:
            result = runner.invoke(
                app,
                ["init", "--fasta", fasta_path, "--skip-snp", "--mode", "compliance"],
            )
            assert result.exit_code == 0
            assert mock_init.call_args[1]["mode"] == "compliance"
        finally:
            Path(fasta_path).unlink()

    def test_init_invalid_mode(self):
        """init with invalid mode exits with error."""
        result = runner.invoke(app, ["init", "--mode", "invalid_mode"])
        assert result.exit_code == 1
        assert "Invalid mode" in (result.output + (result.stderr or ""))

    @patch("plexus.resources.init_genome")
    @patch("plexus.resources.genome_status")
    @patch("plexus.resources.get_operational_mode", return_value="research")
    def test_init_with_checksums(self, _mock_mode, mock_status, mock_init):
        """init --checksums passes the path to init_genome."""
        mock_status.return_value = {
            "fasta": True,
            "fai": True,
            "blast_db": True,
            "snp_vcf": True,
            "fasta_sha256": "a" * 64,
            "snp_vcf_sha256": None,
        }
        with tempfile.NamedTemporaryFile(suffix=".fa", delete=False) as fa_f:
            fasta_path = fa_f.name
            fa_f.write(b">chr1\nACGT\n")
        with tempfile.NamedTemporaryFile(suffix=".txt", delete=False) as ck_f:
            checksums_path = ck_f.name
            ck_f.write(b"abc123  genome.fa\n")
        try:
            result = runner.invoke(
                app,
                [
                    "init",
                    "--fasta",
                    fasta_path,
                    "--skip-snp",
                    "--checksums",
                    checksums_path,
                ],
            )
            assert result.exit_code == 0
            assert str(mock_init.call_args[1]["checksums"]) == checksums_path
        finally:
            Path(fasta_path).unlink()
            Path(checksums_path).unlink()


class TestRunStrictFlag:
    """Tests for --strict flag on run command."""

    @patch("plexus.resources.verify_resource_checksums")
    @patch("plexus.resources.get_operational_mode", return_value="research")
    @patch("plexus.orchestrator.run_pipeline")
    def test_run_strict_passes_with_valid_checksums(
        self, mock_pipeline, _mock_mode, mock_verify
    ):
        """--strict with matching checksums allows pipeline to run."""
        mock_verify.return_value = {"fasta": True, "snp_vcf": True}
        mock_panel = MagicMock()
        mock_panel.junctions = []
        mock_pipeline.return_value = PipelineResult(
            panel=mock_panel,
            output_dir=Path("/tmp/output"),
            config=MagicMock(),
            steps_completed=[],
            errors=[],
        )

        with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as csv_f:
            csv_path = csv_f.name
            csv_f.write(b"Name,Chrom,Position\n")
        with tempfile.NamedTemporaryFile(suffix=".fa", delete=False) as fa_f:
            fasta_path = fa_f.name
            fa_f.write(b">chr1\nACGT\n")
        try:
            result = runner.invoke(
                app,
                ["run", "--input", csv_path, "--fasta", fasta_path, "--strict"],
            )
            assert result.exit_code == 0
            mock_verify.assert_called_once_with("hg38")
        finally:
            Path(csv_path).unlink()
            Path(fasta_path).unlink()

    @patch("plexus.resources.verify_resource_checksums")
    @patch("plexus.resources.get_operational_mode", return_value="research")
    def test_run_strict_fails_on_fasta_mismatch(self, _mock_mode, mock_verify):
        """--strict with FASTA checksum mismatch should exit with error."""
        mock_verify.return_value = {"fasta": False, "snp_vcf": True}

        with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as csv_f:
            csv_path = csv_f.name
            csv_f.write(b"Name,Chrom,Position\n")
        with tempfile.NamedTemporaryFile(suffix=".fa", delete=False) as fa_f:
            fasta_path = fa_f.name
            fa_f.write(b">chr1\nACGT\n")
        try:
            result = runner.invoke(
                app,
                ["run", "--input", csv_path, "--fasta", fasta_path, "--strict"],
            )
            assert result.exit_code == 1
            assert "checksum mismatch" in strip_ansi(result.output).lower()
        finally:
            Path(csv_path).unlink()
            Path(fasta_path).unlink()

    @patch("plexus.resources.verify_resource_checksums")
    @patch("plexus.resources.get_operational_mode", return_value="compliance")
    def test_run_compliance_mode_auto_verifies(self, _mock_mode, mock_verify):
        """Compliance mode triggers checksum verification without --strict."""
        mock_verify.return_value = {"fasta": False, "snp_vcf": None}

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
            mock_verify.assert_called_once()
        finally:
            Path(csv_path).unlink()
            Path(fasta_path).unlink()


class TestStatelessRun:
    """Tests for the stateless --checksums run path."""

    @patch("plexus.orchestrator.run_multi_panel")
    @patch("plexus.resources._compute_sha256")
    @patch("plexus.resources.get_operational_mode", return_value="compliance")
    def test_run_stateless_success(self, _mock_mode, mock_sha, mock_run_multi):
        """--checksums allows run to skip registry in compliance mode."""
        from plexus.pipeline import MultiPanelResult

        # Setup mock hashes
        mock_sha.side_effect = lambda p: f"hash_of_{p.name}"

        # Setup mock result
        mock_result = MagicMock(spec=MultiPanelResult)
        mock_result.success = True
        mock_result.panel_ids = ["panel1"]
        mock_result.panel_results = {
            "panel1": MagicMock(num_junctions=1, selected_pairs=[])
        }
        mock_result.output_dir = Path("/tmp")
        mock_run_multi.return_value = mock_result

        with (
            tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as csv_f,
            tempfile.NamedTemporaryFile(suffix=".fa", delete=False) as fa_f,
            tempfile.NamedTemporaryFile(suffix=".txt", delete=False) as ck_f,
        ):
            csv_path = Path(csv_f.name)
            fa_path = Path(fa_f.name)
            ck_path = Path(ck_f.name)

            # Write a valid sha256sum file
            ck_path.write_text(f"hash_of_{fa_path.name}  {fa_path.name}\n")

            try:
                result = runner.invoke(
                    app,
                    [
                        "run",
                        "-i",
                        str(csv_path),
                        "-f",
                        str(fa_path),
                        "--checksums",
                        str(ck_path),
                    ],
                )

                assert result.exit_code == 0
                assert "FASTA verified" in result.output

                # Verify orchestrator was called with the computed hash
                call_kwargs = mock_run_multi.call_args[1]
                assert call_kwargs["fasta_sha256"] == f"hash_of_{fa_path.name}"
            finally:
                csv_path.unlink()
                fa_path.unlink()
                ck_path.unlink()

    @patch("plexus.resources.get_operational_mode", return_value="compliance")
    def test_run_compliance_requires_fasta(self, _mock_mode):
        """In compliance mode, run must have --fasta explicitly (stateless)."""
        with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as csv_f:
            try:
                result = runner.invoke(app, ["run", "-i", csv_f.name])
                assert result.exit_code == 1
                assert "--fasta is required in compliance mode" in strip_ansi(
                    result.output
                )
            finally:
                Path(csv_f.name).unlink()
