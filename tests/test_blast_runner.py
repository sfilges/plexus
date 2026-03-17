import warnings
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
    """Default run with output_table uses -task blastn-short, outfmt 6, and omits -word_size."""
    output_table = tmp_path / "output.txt"
    runner.db_path = "/fake/db"

    with patch("subprocess.run") as mock_run:
        runner.run(output_table=str(output_table))

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
        assert any("6 qseqid sseqid" in arg for arg in command)
        assert "-out" in command
        assert str(output_table) in command


def test_run_command_explicit_word_size(runner, tmp_path):
    """Explicit word_size is included in the command."""
    output_table = tmp_path / "output.txt"
    runner.db_path = "/fake/db"

    with patch("subprocess.run") as mock_run:
        runner.run(output_table=str(output_table), word_size=11)

        args, kwargs = mock_run.call_args
        command = args[0]
        assert "-task" in command
        assert "blastn-short" in command
        assert "-word_size" in command
        assert "11" in command


def test_run_command_custom_task(runner, tmp_path):
    """Custom task parameter is passed through."""
    output_table = tmp_path / "output.txt"
    runner.db_path = "/fake/db"

    with patch("subprocess.run") as mock_run:
        runner.run(output_table=str(output_table), task="blastn")

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
    output_table = str(tmp_path / "out.txt")
    runner.db_path = "/fake/db"
    with patch("subprocess.run") as mock_run:
        runner.run(output_table=output_table, num_threads=4)
        cmd = mock_run.call_args[0][0]
        assert "-num_threads" in cmd
        assert "4" in cmd


def test_run_omits_num_threads_when_1(runner, tmp_path):
    output_table = str(tmp_path / "out.txt")
    runner.db_path = "/fake/db"
    with patch("subprocess.run") as mock_run:
        runner.run(output_table=output_table, num_threads=1)
        cmd = mock_run.call_args[0][0]
        assert "-num_threads" not in cmd


def test_run_includes_evalue_when_provided(runner, tmp_path):
    """evalue parameter is included in the blastn command when set."""
    output_table = str(tmp_path / "out.txt")
    runner.db_path = "/fake/db"
    with patch("subprocess.run") as mock_run:
        runner.run(output_table=output_table, evalue=1000.0)
        cmd = mock_run.call_args[0][0]
        assert "-evalue" in cmd
        assert "1000.0" in cmd


def test_run_omits_evalue_when_none(runner, tmp_path):
    """evalue parameter is omitted from the blastn command when None."""
    output_table = str(tmp_path / "out.txt")
    runner.db_path = "/fake/db"
    with patch("subprocess.run") as mock_run:
        runner.run(output_table=output_table, evalue=None)
        cmd = mock_run.call_args[0][0]
        assert "-evalue" not in cmd


def test_run_includes_reward_penalty_when_provided(runner, tmp_path):
    """reward and penalty parameters are included in the blastn command."""
    output_table = str(tmp_path / "out.txt")
    runner.db_path = "/fake/db"
    with patch("subprocess.run") as mock_run:
        runner.run(output_table=output_table, reward=1, penalty=-1)
        cmd = mock_run.call_args[0][0]
        assert "-reward" in cmd
        assert "1" in cmd
        assert "-penalty" in cmd
        assert "-1" in cmd


def test_run_omits_reward_penalty_when_none(runner, tmp_path):
    """reward and penalty are omitted from the blastn command when None."""
    output_table = str(tmp_path / "out.txt")
    runner.db_path = "/fake/db"
    with patch("subprocess.run") as mock_run:
        runner.run(output_table=output_table, reward=None, penalty=None)
        cmd = mock_run.call_args[0][0]
        assert "-reward" not in cmd
        assert "-penalty" not in cmd


def test_run_includes_max_hsps_when_provided(runner, tmp_path):
    """max_hsps parameter is included in the blastn command when set."""
    output_table = str(tmp_path / "out.txt")
    runner.db_path = "/fake/db"
    with patch("subprocess.run") as mock_run:
        runner.run(output_table=output_table, max_hsps=100)
        cmd = mock_run.call_args[0][0]
        assert "-max_hsps" in cmd
        assert "100" in cmd


def test_run_omits_max_hsps_when_none(runner, tmp_path):
    """max_hsps parameter is omitted from the blastn command when None."""
    output_table = str(tmp_path / "out.txt")
    runner.db_path = "/fake/db"
    with patch("subprocess.run") as mock_run:
        runner.run(output_table=output_table, max_hsps=None)
        cmd = mock_run.call_args[0][0]
        assert "-max_hsps" not in cmd


def test_run_includes_dust_when_provided(runner, tmp_path):
    """dust parameter is included in the blastn command when set."""
    output_table = str(tmp_path / "out.txt")
    runner.db_path = "/fake/db"
    with patch("subprocess.run") as mock_run:
        runner.run(output_table=output_table, dust="yes")
        cmd = mock_run.call_args[0][0]
        assert "-dust" in cmd
        assert "yes" in cmd


def test_run_omits_dust_when_none(runner, tmp_path):
    """dust parameter is omitted from the blastn command when None."""
    output_table = str(tmp_path / "out.txt")
    runner.db_path = "/fake/db"
    with patch("subprocess.run") as mock_run:
        runner.run(output_table=output_table, dust=None)
        cmd = mock_run.call_args[0][0]
        assert "-dust" not in cmd


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


def test_get_dataframe_uses_compact_dtypes(runner, tmp_path):
    """DataFrame uses categorical strings and downcasted numerics."""
    output_table = tmp_path / "output.txt"
    runner.output_table = str(output_table)

    dummy_data = ["Q1", "S1", 100, 20, 0, 0, 1, 20, 100, 119, 0.0, 40.0, "plus", 20]
    df_content = "	".join(map(str, dummy_data))
    output_table.write_text(df_content)

    df = runner.get_dataframe()

    # String columns should be categorical
    for col in ("qseqid", "sseqid", "sstrand"):
        assert df[col].dtype.name == "category", f"{col} should be categorical"

    # Numeric columns should be downcasted
    assert df["sstart"].dtype == "int32"
    assert df["send"].dtype == "int32"
    assert df["length"].dtype == "int32"
    assert df["mismatch"].dtype == "int32"
    assert df["evalue"].dtype == "float32"
    assert df["pident"].dtype == "float32"
    assert df["bitscore"].dtype == "float32"
    assert df["qlen"].dtype == "int16"


def test_run_direct_tabular_output(runner, tmp_path):
    """output_table= produces outfmt 6 and sets self.output_table."""
    output_table = str(tmp_path / "table.txt")
    runner.db_path = "/fake/db"

    with patch("subprocess.run"):
        runner.run(output_table=output_table)

    assert runner.output_table == output_table
    assert not hasattr(runner, "output_archive")


def test_run_legacy_archive_emits_deprecation_warning(runner, tmp_path):
    """output_archive= still works but emits DeprecationWarning."""
    output_archive = str(tmp_path / "out.asn")
    runner.db_path = "/fake/db"

    with patch("subprocess.run") as mock_run:
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            runner.run(output_archive)

            assert len(w) == 1
            assert issubclass(w[0].category, DeprecationWarning)
            assert "output_archive is deprecated" in str(w[0].message)

        cmd = mock_run.call_args[0][0]
        assert "-outfmt" in cmd
        idx = cmd.index("-outfmt")
        assert cmd[idx + 1] == "11"
    assert runner.output_archive == output_archive


def test_run_raises_when_no_output_specified(runner):
    """ValueError when neither output_table nor output_archive is given."""
    runner.db_path = "/fake/db"

    with pytest.raises(ValueError, match="Either output_table or output_archive"):
        runner.run()
