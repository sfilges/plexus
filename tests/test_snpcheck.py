# ================================================================================
# Tests for SNP overlap checking module
# ================================================================================

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from typer.testing import CliRunner

from plexus.cli import app
from plexus.config import DesignerConfig, SnpCheckParameters
from plexus.designer.primer import Primer, PrimerPair
from plexus.snpcheck.checker import (
    _count_snps_in_region,
    _get_allele_frequency,
    _primer_genomic_coords,
    run_snp_check,
)
from plexus.snpcheck.resources import (
    DEFAULT_CACHE_DIR,
    ENV_DATA_DIR,
    ENV_SNP_VCF,
    GNOMAD_TBI_FILENAME,
    GNOMAD_VCF_FILENAME,
    get_cache_dir,
    is_resource_available,
)
from plexus.snpcheck.snp_data import (
    _write_regions_bed,
    get_snp_vcf,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

runner = CliRunner()


def _make_primer(name="test_fwd", direction="forward", start=10, length=20):
    return Primer(
        name=name, seq="A" * length, direction=direction, start=start, length=length
    )


def _make_pair(fwd_start=10, fwd_len=20, rev_start=50, rev_len=20):
    fwd = _make_primer("fwd", "forward", fwd_start, fwd_len)
    rev = _make_primer("rev", "reverse", rev_start, rev_len)
    return PrimerPair(
        forward=fwd,
        reverse=rev,
        insert_size=rev_start - (fwd_start + fwd_len),
        amplicon_sequence="A" * 80,
        amplicon_length=80,
        pair_penalty=1.0,
        pair_id="test_pair_0",
    )


def _make_junction(design_start=1000, chrom="chr7"):
    junction = MagicMock()
    junction.name = "TEST_JUNCTION"
    junction.chrom = chrom
    junction.start = design_start + 200
    junction.end = design_start + 200
    junction.design_start = design_start
    junction.design_end = design_start + 400
    junction.primer_pairs = [_make_pair()]
    return junction


def _make_panel(junctions=None):
    panel = MagicMock()
    panel.junctions = junctions or [_make_junction()]
    return panel


# ---------------------------------------------------------------------------
# Coordinate mapping
# ---------------------------------------------------------------------------


class TestPrimerGenomicCoords:
    def test_basic_mapping(self):
        """Primer at relative start=10 with design_start=1000 -> genomic 1010."""
        primer = _make_primer(start=10, length=20)
        junction = _make_junction(design_start=1000)
        chrom, gstart, gend = _primer_genomic_coords(primer, junction)
        assert chrom == "chr7"
        assert gstart == 1010
        assert gend == 1030

    def test_zero_offset(self):
        """Primer at the start of the design region."""
        primer = _make_primer(start=0, length=25)
        junction = _make_junction(design_start=5000)
        _, gstart, gend = _primer_genomic_coords(primer, junction)
        assert gstart == 5000
        assert gend == 5025


# ---------------------------------------------------------------------------
# Allele frequency extraction
# ---------------------------------------------------------------------------


class TestGetAlleleFrequency:
    def test_af_field_single(self):
        record = MagicMock()
        record.info = {"AF": 0.05}
        assert _get_allele_frequency(record) == 0.05

    def test_af_field_list(self):
        record = MagicMock()
        record.info = {"AF": [0.01, 0.03]}
        assert _get_allele_frequency(record) == 0.03

    def test_caf_field(self):
        record = MagicMock()
        record.info = {"CAF": ["0.95", "0.03", "0.02"]}
        assert _get_allele_frequency(record) == 0.03

    def test_caf_with_missing(self):
        record = MagicMock()
        record.info = {"CAF": ["0.95", ".", "0.05"]}
        assert _get_allele_frequency(record) == 0.05

    def test_no_frequency_info(self):
        record = MagicMock()
        record.info = {}
        assert _get_allele_frequency(record) is None


# ---------------------------------------------------------------------------
# SNP counting
# ---------------------------------------------------------------------------


class TestCountSnpsInRegion:
    def test_no_variants(self):
        vcf = MagicMock()
        vcf.fetch.return_value = []
        count, positions = _count_snps_in_region(vcf, "chr7", 1000, 1020, 0.01)
        assert count == 0
        assert positions == []

    def test_variants_above_threshold(self):
        record1 = MagicMock()
        record1.info = {"AF": 0.05}
        record1.pos = 1005
        record2 = MagicMock()
        record2.info = {"AF": 0.001}  # Below threshold
        record2.pos = 1010

        vcf = MagicMock()
        vcf.fetch.return_value = [record1, record2]

        count, positions = _count_snps_in_region(vcf, "chr7", 1000, 1020, 0.01)
        assert count == 1
        assert positions == [1005]

    def test_contig_not_in_vcf(self):
        """ValueError from pysam when contig doesn't exist should be handled."""
        vcf = MagicMock()
        vcf.fetch.side_effect = ValueError("contig not found")
        count, positions = _count_snps_in_region(vcf, "chrUn", 100, 200, 0.01)
        assert count == 0
        assert positions == []


# ---------------------------------------------------------------------------
# run_snp_check integration
# ---------------------------------------------------------------------------


class TestRunSnpCheck:
    def test_applies_penalty(self):
        """SNP overlaps should increase snp_count and pair_penalty."""
        panel = _make_panel()
        pair = panel.junctions[0].primer_pairs[0]
        original_penalty = pair.pair_penalty

        # Mock pysam.VariantFile
        mock_record = MagicMock()
        mock_record.info = {"AF": 0.05}
        mock_record.pos = 1015

        with patch("pysam.VariantFile") as mock_variant_file:
            mock_vcf = MagicMock()
            mock_vcf.fetch.return_value = [mock_record]
            mock_variant_file.return_value = mock_vcf

            run_snp_check(
                panel=panel,
                vcf_path="/fake/dbsnp.vcf.gz",
                af_threshold=0.01,
                snp_penalty_weight=5.0,
            )

        # Each primer gets 1 SNP -> pair gets 2 SNPs
        assert pair.snp_count == 2
        assert pair.snp_penalty == 10.0
        assert pair.pair_penalty == original_penalty + 10.0

    def test_no_snps_no_penalty(self):
        """No SNPs means penalty unchanged."""
        panel = _make_panel()
        pair = panel.junctions[0].primer_pairs[0]
        original_penalty = pair.pair_penalty

        with patch("pysam.VariantFile") as mock_variant_file:
            mock_vcf = MagicMock()
            mock_vcf.fetch.return_value = []
            mock_variant_file.return_value = mock_vcf

            run_snp_check(
                panel=panel,
                vcf_path="/fake/dbsnp.vcf.gz",
                af_threshold=0.01,
                snp_penalty_weight=5.0,
            )

        assert pair.snp_count == 0
        assert pair.snp_penalty == 0.0
        assert pair.pair_penalty == original_penalty

    def test_skips_junction_without_design_start(self):
        """Junctions without design_start are skipped gracefully."""
        junction = _make_junction()
        junction.design_start = None
        panel = _make_panel([junction])

        with patch("pysam.VariantFile") as mock_variant_file:
            mock_vcf = MagicMock()
            mock_variant_file.return_value = mock_vcf

            # Should not raise
            run_snp_check(panel=panel, vcf_path="/fake/dbsnp.vcf.gz")

        # fetch should never be called
        mock_vcf.fetch.assert_not_called()


# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------


class TestSnpCheckConfig:
    def test_default_values(self):
        params = SnpCheckParameters()
        assert params.af_threshold == 0.01
        assert params.snp_penalty_weight == 5.0

    def test_custom_values(self):
        params = SnpCheckParameters(af_threshold=0.05, snp_penalty_weight=10.0)
        assert params.af_threshold == 0.05
        assert params.snp_penalty_weight == 10.0

    def test_included_in_designer_config(self):
        config = DesignerConfig()
        assert hasattr(config, "snp_check_parameters")
        assert config.snp_check_parameters.af_threshold == 0.01

    def test_af_threshold_validation(self):
        from pydantic import ValidationError

        with pytest.raises(ValidationError):
            SnpCheckParameters(af_threshold=-0.1)
        with pytest.raises(ValidationError):
            SnpCheckParameters(af_threshold=1.5)


# ---------------------------------------------------------------------------
# PrimerPair SNP fields
# ---------------------------------------------------------------------------


class TestPrimerPairSnpFields:
    def test_default_zero(self):
        pair = _make_pair()
        assert pair.snp_count == 0
        assert pair.snp_penalty == 0.0

    def test_primer_default_zero(self):
        primer = _make_primer()
        assert primer.snp_count == 0


# ---------------------------------------------------------------------------
# SNP data module (snp_data.py)
# ---------------------------------------------------------------------------


class TestWriteRegionsBed:
    def test_generates_bed(self, tmp_path):
        """BED file has correct chr-style regions."""
        panel = _make_panel([_make_junction(design_start=1000, chrom="chr7")])
        bed = _write_regions_bed(panel, tmp_path, padding=200)

        lines = bed.read_text().strip().split("\n")
        assert len(lines) == 1
        parts = lines[0].split("\t")
        assert parts[0] == "chr7"
        assert int(parts[1]) == 800  # 1000 - 200
        assert int(parts[2]) > 800

    def test_skips_junction_without_design_start(self, tmp_path):
        junction = _make_junction()
        junction.design_start = None
        panel = _make_panel([junction])
        bed = _write_regions_bed(panel, tmp_path, padding=200)

        content = bed.read_text().strip()
        assert content == ""

    def test_multiple_junctions(self, tmp_path):
        junctions = [
            _make_junction(design_start=1000, chrom="chr7"),
            _make_junction(design_start=5000, chrom="chr1"),
        ]
        panel = _make_panel(junctions)
        bed = _write_regions_bed(panel, tmp_path, padding=200)

        lines = bed.read_text().strip().split("\n")
        assert len(lines) == 2


class TestGetSnpVcf:
    def test_user_vcf_returned_directly(self, tmp_path):
        """User-provided VCF is returned without intersection."""
        fake_vcf = tmp_path / "user.vcf.gz"
        fake_vcf.touch()

        result = get_snp_vcf(
            panel=_make_panel(),
            output_dir=tmp_path,
            user_vcf=str(fake_vcf),
        )
        assert result == fake_vcf

    def test_user_vcf_not_found_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError, match="SNP VCF not found"):
            get_snp_vcf(
                panel=_make_panel(),
                output_dir=tmp_path,
                user_vcf="/nonexistent/file.vcf.gz",
            )

    def test_env_var_resolution(self, tmp_path):
        """$PLEXUS_SNP_VCF should be used when no user VCF is provided."""
        fake_vcf = tmp_path / "env.vcf.gz"
        fake_vcf.touch()

        with patch.dict("os.environ", {ENV_SNP_VCF: str(fake_vcf)}):
            result = get_snp_vcf(
                panel=_make_panel(),
                output_dir=tmp_path,
            )
        assert result == fake_vcf

    def test_env_var_not_found_raises(self, tmp_path):
        """$PLEXUS_SNP_VCF pointing to missing file should raise."""
        with patch.dict("os.environ", {ENV_SNP_VCF: "/nonexistent/env.vcf.gz"}):
            with pytest.raises(FileNotFoundError, match="PLEXUS_SNP_VCF"):
                get_snp_vcf(
                    panel=_make_panel(),
                    output_dir=tmp_path,
                )

    def test_user_vcf_takes_priority_over_env(self, tmp_path):
        """--snp-vcf should win over $PLEXUS_SNP_VCF."""
        user_vcf = tmp_path / "user.vcf.gz"
        user_vcf.touch()
        env_vcf = tmp_path / "env.vcf.gz"
        env_vcf.touch()

        with patch.dict("os.environ", {ENV_SNP_VCF: str(env_vcf)}):
            result = get_snp_vcf(
                panel=_make_panel(),
                output_dir=tmp_path,
                user_vcf=str(user_vcf),
            )
        assert result == user_vcf

    def test_cached_gnomad_fallback(self, tmp_path):
        """When no user VCF or env var, uses cached gnomAD with bcftools intersection."""
        panel = _make_panel()

        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stderr = ""

        with (
            patch(
                "plexus.snpcheck.snp_data.is_resource_available",
                return_value=True,
            ),
            patch(
                "plexus.snpcheck.snp_data.get_cached_vcf_path",
                return_value=tmp_path / "gnomad.vcf.gz",
            ),
            patch("plexus.snpcheck.snp_data._check_bcftools", return_value="bcftools"),
            patch("subprocess.run", return_value=mock_result) as mock_run,
        ):
            # Create fake cached file
            (tmp_path / "gnomad.vcf.gz").touch()

            result = get_snp_vcf(
                panel=panel,
                output_dir=tmp_path,
            )

            # bcftools should be called (view + index)
            assert mock_run.call_count == 2
            assert result == tmp_path / "snpcheck_regions.vcf.gz"

    def test_no_vcf_available_raises(self, tmp_path):
        """When nothing is available, should raise with actionable instructions."""
        with patch(
            "plexus.snpcheck.snp_data.is_resource_available",
            return_value=False,
        ):
            with pytest.raises(FileNotFoundError, match="plexus download-resources"):
                get_snp_vcf(
                    panel=_make_panel(),
                    output_dir=tmp_path,
                )


# ---------------------------------------------------------------------------
# Resources module (resources.py)
# ---------------------------------------------------------------------------


class TestCacheDir:
    def test_default_cache_dir(self):
        with patch.dict("os.environ", {}, clear=True):
            result = get_cache_dir()
        assert result == DEFAULT_CACHE_DIR

    def test_env_override(self, tmp_path):
        with patch.dict("os.environ", {ENV_DATA_DIR: str(tmp_path)}):
            result = get_cache_dir()
        assert result == tmp_path


class TestIsResourceAvailable:
    def test_both_files_present(self, tmp_path):
        with patch("plexus.snpcheck.resources.get_cache_dir", return_value=tmp_path):
            (tmp_path / GNOMAD_VCF_FILENAME).touch()
            (tmp_path / GNOMAD_TBI_FILENAME).touch()
            assert is_resource_available() is True

    def test_vcf_missing(self, tmp_path):
        with patch("plexus.snpcheck.resources.get_cache_dir", return_value=tmp_path):
            (tmp_path / GNOMAD_TBI_FILENAME).touch()
            assert is_resource_available() is False

    def test_tbi_missing(self, tmp_path):
        with patch("plexus.snpcheck.resources.get_cache_dir", return_value=tmp_path):
            (tmp_path / GNOMAD_VCF_FILENAME).touch()
            assert is_resource_available() is False

    def test_neither_present(self, tmp_path):
        with patch("plexus.snpcheck.resources.get_cache_dir", return_value=tmp_path):
            assert is_resource_available() is False


class TestDownloadGnomadVcf:
    def test_skips_if_already_cached(self, tmp_path):
        """download_gnomad_vcf should skip download when files exist."""
        from plexus.snpcheck.resources import download_gnomad_vcf

        with patch("plexus.snpcheck.resources.get_cache_dir", return_value=tmp_path):
            (tmp_path / GNOMAD_VCF_FILENAME).touch()
            (tmp_path / GNOMAD_TBI_FILENAME).touch()

            with patch("plexus.snpcheck.resources._download_with_progress") as mock_dl:
                result = download_gnomad_vcf()
                mock_dl.assert_not_called()
            assert result == tmp_path / GNOMAD_VCF_FILENAME

    def test_downloads_when_missing(self, tmp_path):
        """download_gnomad_vcf should call _download_with_progress for both files."""
        from plexus.snpcheck.resources import download_gnomad_vcf

        with (
            patch("plexus.snpcheck.resources.get_cache_dir", return_value=tmp_path),
            patch("plexus.snpcheck.resources._download_with_progress") as mock_dl,
        ):
            download_gnomad_vcf()
            assert mock_dl.call_count == 2

    def test_force_redownloads(self, tmp_path):
        """download_gnomad_vcf(force=True) should download even when cached."""
        from plexus.snpcheck.resources import download_gnomad_vcf

        with patch("plexus.snpcheck.resources.get_cache_dir", return_value=tmp_path):
            (tmp_path / GNOMAD_VCF_FILENAME).touch()
            (tmp_path / GNOMAD_TBI_FILENAME).touch()

            with patch("plexus.snpcheck.resources._download_with_progress") as mock_dl:
                download_gnomad_vcf(force=True)
                assert mock_dl.call_count == 2

    def test_cleanup_on_failure(self, tmp_path):
        """Partial .part file should be removed on download failure."""
        from plexus.snpcheck.resources import _download_with_progress

        dest = tmp_path / "test.vcf.gz"

        with patch("urllib.request.urlopen", side_effect=OSError("network error")):
            with pytest.raises(OSError, match="network error"):
                _download_with_progress("http://example.com/file", dest)

        # .part file should not remain
        assert not dest.with_suffix(".gz.part").exists()
        assert not dest.exists()


# ---------------------------------------------------------------------------
# CLI commands
# ---------------------------------------------------------------------------


class TestDownloadResourcesCli:
    def test_success(self):
        with patch(
            "plexus.snpcheck.resources.download_gnomad_vcf",
            return_value=Path("/fake/gnomad.vcf.gz"),
        ):
            result = runner.invoke(app, ["download-resources"])
        assert result.exit_code == 0
        assert "Done" in result.output

    def test_failure(self):
        with patch(
            "plexus.snpcheck.resources.download_gnomad_vcf",
            side_effect=OSError("network error"),
        ):
            result = runner.invoke(app, ["download-resources"])
        assert result.exit_code == 1
        assert "Download failed" in result.output


class TestStatusCli:
    def test_shows_version_and_resource_status(self):
        with patch(
            "plexus.snpcheck.resources.resource_status_message",
            return_value="gnomAD VCF: not downloaded",
        ):
            result = runner.invoke(app, ["status"])
        assert result.exit_code == 0
        assert "Plexus" in result.output
        assert "gnomAD VCF" in result.output
