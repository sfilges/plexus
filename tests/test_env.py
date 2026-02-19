from unittest.mock import patch

from plexus.utils.env import check_executable, get_missing_tools


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
