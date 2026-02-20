from unittest.mock import patch

from plexus.utils.env import (
    check_disk_space,
    check_executable,
    get_missing_tools,
    get_plexus_version,
    get_primer3_version,
    get_tool_version,
    get_tool_versions,
)


def test_check_executable():
    # Test with something that definitely exists (like 'ls' or 'sh' on unix)
    assert check_executable("sh") is True
    # Test with something that shouldn't exist
    assert check_executable("non_existent_tool_12345") is False


def test_get_missing_tools_all_present():
    with patch("plexus.utils.env.check_executable", return_value=True):
        missing = get_missing_tools(need_blast=True, need_snp=True)
        assert len(missing) == 0


def test_get_missing_tools_some_missing():
    def mock_check(name):
        return name not in ["blastn", "bcftools"]

    with patch("plexus.utils.env.check_executable", side_effect=mock_check):
        missing = get_missing_tools(need_blast=True, need_snp=True)
        assert "blastn" in missing
        assert "bcftools" in missing
        assert "makeblastdb" not in missing
        assert "blast_formatter" not in missing


def test_get_missing_tools_conditional():
    def mock_check(name):
        return False  # All missing

    with patch("plexus.utils.env.check_executable", side_effect=mock_check):
        # Only check SNP
        missing = get_missing_tools(need_blast=False, need_snp=True)
        assert missing == ["bcftools"]

        # Only check BLAST
        missing = get_missing_tools(need_blast=True, need_snp=False)
        assert sorted(missing) == sorted(["blastn", "makeblastdb", "blast_formatter"])


# ---------------------------------------------------------------------------
# Tool version capture
# ---------------------------------------------------------------------------


def test_get_tool_version_existing():
    """Tool on PATH returns a non-empty version string."""
    # 'sh' is always available and --version typically produces output
    with patch("plexus.utils.env.check_executable", return_value=True):
        with patch("subprocess.run") as mock_run:
            mock_run.return_value.stdout = "toolname 1.2.3\nExtra info"
            mock_run.return_value.stderr = ""
            result = get_tool_version("toolname")
            assert result == "toolname 1.2.3"


def test_get_tool_version_missing():
    """Tool not on PATH returns None."""
    with patch("plexus.utils.env.check_executable", return_value=False):
        assert get_tool_version("nonexistent_tool_xyz") is None


def test_get_tool_version_falls_back_to_stderr():
    """Some tools write version info to stderr."""
    with patch("plexus.utils.env.check_executable", return_value=True):
        with patch("subprocess.run") as mock_run:
            mock_run.return_value.stdout = ""
            mock_run.return_value.stderr = "tool v3.0"
            result = get_tool_version("sometool")
            assert result == "tool v3.0"


def test_get_tool_versions_dict():
    """Batch version lookup returns correct dict shape."""
    with patch("plexus.utils.env.get_tool_version") as mock_ver:
        mock_ver.side_effect = lambda name: f"{name}-1.0" if name != "missing" else None
        result = get_tool_versions(["blastn", "bcftools", "missing"])
        assert result == {
            "blastn": "blastn-1.0",
            "bcftools": "bcftools-1.0",
            "missing": None,
        }


def test_get_plexus_version():
    """Returns a non-empty version string."""
    version = get_plexus_version()
    assert isinstance(version, str)
    assert len(version) > 0


def test_get_primer3_version():
    """Returns a string or None."""
    version = get_primer3_version()
    assert version is None or isinstance(version, str)


def test_check_disk_space():
    """Verify that check_disk_space correctly reports status based on thresholds."""
    with patch("shutil.disk_usage") as mock_usage:
        # Sufficient space: 10 GB free
        mock_usage.return_value.free = 10 * 1024**3
        assert check_disk_space("/tmp", threshold_gb=2.0) is True

        # Insufficient space: 1 GB free
        mock_usage.return_value.free = 1 * 1024**3
        assert check_disk_space("/tmp", threshold_gb=2.0) is False


class TestComplianceManifest:
    """Tests for the compliance manifest and environment validation."""

    def test_load_manifest(self):
        """Manifest can be loaded as a dict."""
        from plexus.utils.env import load_compliance_manifest

        manifest = load_compliance_manifest()
        assert isinstance(manifest, dict)
        assert "tools" in manifest
        assert "blastn" in manifest["tools"]

    @patch("plexus.utils.env.get_tool_version")
    def test_validate_environment_pass(self, mock_ver):
        """Returns passing verdict when versions match exactly."""
        from plexus.utils.env import validate_environment

        mock_ver.side_effect = lambda name: {
            "blastn": "blastn: 2.17.0+",
            "makeblastdb": "makeblastdb: 2.17.0+",
            "blast_formatter": "blast_formatter: 2.17.0+",
            "bcftools": "bcftools 1.23",
        }.get(name)

        verdict = validate_environment(need_blast=True, need_snp=True)
        assert verdict["blastn"]["verdict"] == "pass"
        assert verdict["bcftools"]["verdict"] == "pass"

    @patch("plexus.utils.env.get_tool_version")
    def test_validate_environment_fail(self, mock_ver):
        """Raises ComplianceError when versions mismatch."""
        from plexus.utils.env import ComplianceError, validate_environment

        mock_ver.side_effect = lambda name: {
            "blastn": "blastn: 2.12.0+",  # Mismatch
            "bcftools": "bcftools 1.23",
        }.get(name)

        import pytest

        with pytest.raises(ComplianceError) as excinfo:
            validate_environment(need_blast=True, need_snp=True)

        assert "blastn: expected exactly 2.17.0, found '2.12.0'" in str(excinfo.value)

    @patch("plexus.utils.env.get_tool_version")
    def test_validate_environment_missing(self, mock_ver):
        """Raises ComplianceError when tool is missing."""
        from plexus.utils.env import ComplianceError, validate_environment

        mock_ver.return_value = None

        import pytest

        with pytest.raises(ComplianceError) as excinfo:
            validate_environment(need_blast=True)

        assert "blastn: not found on PATH" in str(excinfo.value)
