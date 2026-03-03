from unittest.mock import patch

import pandas as pd
import pytest

from plexus.blast.blast_runner import BlastRunner


@pytest.fixture
def runner(tmp_path):
    fasta = tmp_path / "ref.fasta"
    fasta.touch()
    query = tmp_path / "query.fasta"
    query.touch()
    return BlastRunner(str(query), str(fasta))


def test_create_database_command(runner):
    with patch("subprocess.run") as mock_run:
        with patch("os.path.isfile", return_value=False):
            runner.create_database()

            # Check if makeblastdb was called with the right list arguments
            args, kwargs = mock_run.call_args
            command = args[0]
            assert "makeblastdb" in command
            assert "-in" in command
            assert runner.reference_fasta in command
            assert "-out" in command
            assert runner.db_path in command


def test_create_database_skips_if_exists(runner):
    with patch("subprocess.run") as mock_run:
        with patch("os.path.isfile", return_value=True):
            runner.create_database()
            mock_run.assert_not_called()


def test_create_database_skips_if_exists_v5(runner):
    """v5 database: .njs manifest present, .nhr absent — must skip makeblastdb."""

    def isfile_v5(path):
        return path.endswith(".njs")

    with patch("subprocess.run") as mock_run:
        with patch("os.path.isfile", side_effect=isfile_v5):
            runner.create_database()
            mock_run.assert_not_called()


def test_create_database_skips_if_exists_v4(runner):
    """v4 database: .nhr present, .njs absent — must skip makeblastdb."""

    def isfile_v4(path):
        return path.endswith(".nhr")

    with patch("subprocess.run") as mock_run:
        with patch("os.path.isfile", side_effect=isfile_v4):
            runner.create_database()
            mock_run.assert_not_called()


def test_run_command_defaults(runner, tmp_path):
    """Default run uses -task blastn-short and omits -word_size."""
    output_archive = tmp_path / "output.asn"
    runner.db_path = "/fake/db"

    with patch("subprocess.run") as mock_run:
        runner.run(str(output_archive))

        args, kwargs = mock_run.call_args
        command = args[0]
        assert "blastn" in command
        assert "-db" in command
        assert runner.db_path in command
        assert "-query" in command
        assert runner.input_fasta in command
        assert "-task" in command
        assert "blastn-short" in command
        assert "-word_size" not in command
        assert "-outfmt" in command
        assert "-out" in command
        assert str(output_archive) in command


def test_run_command_explicit_word_size(runner, tmp_path):
    """Explicit word_size is included in the command."""
    output_archive = tmp_path / "output.asn"
    runner.db_path = "/fake/db"

    with patch("subprocess.run") as mock_run:
        runner.run(str(output_archive), word_size=11)

        args, kwargs = mock_run.call_args
        command = args[0]
        assert "-task" in command
        assert "blastn-short" in command
        assert "-word_size" in command
        assert "11" in command


def test_run_command_custom_task(runner, tmp_path):
    """Custom task parameter is passed through."""
    output_archive = tmp_path / "output.asn"
    runner.db_path = "/fake/db"

    with patch("subprocess.run") as mock_run:
        runner.run(str(output_archive), task="blastn")

        args, kwargs = mock_run.call_args
        command = args[0]
        assert "-task" in command
        assert "blastn" in command
        assert "blastn-short" not in command


def test_reformat_output_command(runner, tmp_path):
    runner.output_archive = str(tmp_path / "output.asn")
    output_table = tmp_path / "output.txt"

    with patch("subprocess.run") as mock_run:
        runner.reformat_output_as_table(str(output_table))

        args, kwargs = mock_run.call_args
        command = args[0]
        assert "blast_formatter" in command
        assert "-archive" in command
        assert runner.output_archive in command
        assert "-out" in command
        assert str(output_table) in command
        assert any("6 qseqid sseqid" in arg for arg in command)


def test_run_emits_num_threads_when_gt_1(runner, tmp_path):
    output_archive = str(tmp_path / "out.asn")
    runner.db_path = "/fake/db"
    with patch("subprocess.run") as mock_run:
        runner.run(output_archive, num_threads=4)
        cmd = mock_run.call_args[0][0]
        assert "-num_threads" in cmd
        assert "4" in cmd


def test_run_omits_num_threads_when_1(runner, tmp_path):
    output_archive = str(tmp_path / "out.asn")
    runner.db_path = "/fake/db"
    with patch("subprocess.run") as mock_run:
        runner.run(output_archive, num_threads=1)
        cmd = mock_run.call_args[0][0]
        assert "-num_threads" not in cmd


def test_run_includes_evalue_when_provided(runner, tmp_path):
    """evalue parameter is included in the blastn command when set."""
    output_archive = str(tmp_path / "out.asn")
    runner.db_path = "/fake/db"
    with patch("subprocess.run") as mock_run:
        runner.run(output_archive, evalue=1000.0)
        cmd = mock_run.call_args[0][0]
        assert "-evalue" in cmd
        assert "1000.0" in cmd


def test_run_omits_evalue_when_none(runner, tmp_path):
    """evalue parameter is omitted from the blastn command when None."""
    output_archive = str(tmp_path / "out.asn")
    runner.db_path = "/fake/db"
    with patch("subprocess.run") as mock_run:
        runner.run(output_archive, evalue=None)
        cmd = mock_run.call_args[0][0]
        assert "-evalue" not in cmd


def test_run_includes_reward_penalty_when_provided(runner, tmp_path):
    """reward and penalty parameters are included in the blastn command."""
    output_archive = str(tmp_path / "out.asn")
    runner.db_path = "/fake/db"
    with patch("subprocess.run") as mock_run:
        runner.run(output_archive, reward=1, penalty=-1)
        cmd = mock_run.call_args[0][0]
        assert "-reward" in cmd
        assert "1" in cmd
        assert "-penalty" in cmd
        assert "-1" in cmd


def test_run_omits_reward_penalty_when_none(runner, tmp_path):
    """reward and penalty are omitted from the blastn command when None."""
    output_archive = str(tmp_path / "out.asn")
    runner.db_path = "/fake/db"
    with patch("subprocess.run") as mock_run:
        runner.run(output_archive, reward=None, penalty=None)
        cmd = mock_run.call_args[0][0]
        assert "-reward" not in cmd
        assert "-penalty" not in cmd


def test_get_dataframe(runner, tmp_path):
    output_table = tmp_path / "output.txt"
    runner.output_table = str(output_table)

    # Create a dummy table
    dummy_data = ["Q1", "S1", 100, 20, 0, 0, 1, 20, 100, 119, 0.0, 40.0, "plus", 20]
    df_content = "	".join(map(str, dummy_data))
    output_table.write_text(df_content)

    df = runner.get_dataframe()
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 1
    assert df.iloc[0]["qseqid"] == "Q1"
