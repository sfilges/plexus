import subprocess
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from plexus.designer.multiplexpanel import MultiplexPanel
from plexus.snpcheck.snp_data import (
    _get_vcf_contigs,
    _write_regions_bed,
)
from plexus.utils.utils import run_command


def test_run_command_success():
    result = run_command(["echo", "hello"], retries=1)
    assert result.stdout.strip() == "hello"
    assert result.returncode == 0


def test_run_command_retry_success():
    # Mock subprocess.run to fail once and then succeed
    mock_run = MagicMock()
    mock_run.side_effect = [
        subprocess.CalledProcessError(1, ["cmd"], stderr="fail"),
        subprocess.CompletedProcess(["cmd"], 0, stdout="success", stderr=""),
    ]

    with patch("subprocess.run", mock_run):
        # We need to mock time.sleep so it doesn't actually wait
        with patch("time.sleep"):
            result = run_command(["cmd"], retries=2)
            assert result.stdout == "success"
            assert mock_run.call_count == 2


def test_run_command_failure_after_retries():
    mock_run = MagicMock()
    mock_run.side_effect = subprocess.CalledProcessError(1, ["cmd"], stderr="fail")

    with patch("subprocess.run", mock_run):
        with patch("time.sleep"):
            with pytest.raises(subprocess.CalledProcessError):
                run_command(["cmd"], retries=2)
            assert mock_run.call_count == 3


def test_write_regions_bed_filtering(tmp_path):
    # Mock a panel with junctions on different contigs
    panel = MagicMock(spec=MultiplexPanel)
    j1 = MagicMock()
    j1.chrom = "chr1"
    j1.design_start = 100
    j1.design_end = 200

    j2 = MagicMock()
    j2.chrom = "chr2"
    j2.design_start = 300
    j2.design_end = 400

    panel.junctions = [j1, j2]

    valid_contigs = {"chr1"}

    bed_path = _write_regions_bed(
        panel, tmp_path, padding=0, valid_contigs=valid_contigs
    )

    assert bed_path.exists()
    content = bed_path.read_text()
    assert "chr1\t100\t200" in content
    assert "chr2" not in content


def test_get_vcf_contigs_mock():
    mock_result = MagicMock()
    mock_result.stdout = (
        "##fileformat=VCFv4.1\n"
        "##contig=<ID=chr1,length=248956422>\n"
        "##contig=<ID=chr2,length=242193529>\n"
        "##contig=<ID=chr3,length=198295559>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    )

    with patch("plexus.snpcheck.snp_data.run_command", return_value=mock_result):
        with patch("plexus.snpcheck.snp_data._check_bcftools", return_value="bcftools"):
            contigs = _get_vcf_contigs(Path("fake.vcf.gz"))
            assert contigs == {"chr1", "chr2", "chr3"}
