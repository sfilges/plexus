from unittest.mock import MagicMock, patch

import pandas as pd
import pytest

from plexus.blast.specificity import (
    _is_on_target,
    filter_offtarget_pairs,
    run_specificity_check,
)
from plexus.designer.multiplexpanel import MultiplexPanel


@pytest.fixture
def mock_panel():
    panel = MagicMock(spec=MultiplexPanel)
    panel.unique_primer_map = {"ACTG": "P1_F", "CAGT": "P1_R"}
    panel.primer_target_map = {"P1_F": "JUNCTION1", "P1_R": "JUNCTION1"}

    # Mock Junction
    junction = MagicMock()
    junction.chrom = "chr7"
    junction.name = "J1"
    junction.design_start = 1000

    # Mock PrimerPair
    # junction.design_start=1000, fwd_start=10 -> expected F_start=1010
    # junction.design_start=1000, rev_start=180, rev_length=22
    #   -> expected R_start = 1000 + 180 + 22 - 1 = 1201
    pair = MagicMock()
    pair.forward.seq = "ACTG"
    pair.forward.start = 10
    pair.forward.length = 4
    pair.reverse.seq = "CAGT"
    pair.reverse.start = 180
    pair.reverse.length = 22
    pair.amplicon_length = 200
    pair.pair_id = "J1_P1"

    junction.primer_pairs = [pair]
    panel.junctions = [junction]

    return panel


def test_blast_archive_removed_after_run(mock_panel, tmp_path):
    """The BLAST archive file is cleaned up after a successful specificity check."""
    blast_archive = tmp_path / "blast_archive"

    with (
        patch("plexus.blast.specificity.BlastRunner") as MockRunner,
        patch("plexus.blast.specificity.BlastResultsAnnotator") as MockAnnotator,
        patch("plexus.blast.specificity.AmpliconFinder") as MockFinder,
    ):
        runner_instance = MockRunner.return_value

        def fake_run(output_archive, **kwargs):
            # Simulate BLAST creating the archive file
            open(output_archive, "w").close()

        runner_instance.run.side_effect = fake_run
        runner_instance.get_dataframe.return_value = pd.DataFrame({"dummy": [1]})

        annotator_instance = MockAnnotator.return_value
        annotator_instance.get_predicted_bound.return_value = pd.DataFrame(
            {"dummy_bound": [1]}
        )

        finder_instance = MockFinder.return_value
        finder_instance.amplicon_df = pd.DataFrame()

        run_specificity_check(mock_panel, str(tmp_path), "fake_genome.fa")

        assert not blast_archive.exists(), "BLAST archive should be removed after run"


def test_run_specificity_check_integration(mock_panel, tmp_path):
    # Setup mocks for internal classes
    with (
        patch("plexus.blast.specificity.BlastRunner") as MockRunner,
        patch("plexus.blast.specificity.BlastResultsAnnotator") as MockAnnotator,
        patch("plexus.blast.specificity.AmpliconFinder") as MockFinder,
        patch("os.makedirs"),
    ):
        # 1. Mock Runner behavior
        runner_instance = MockRunner.return_value
        runner_instance.get_dataframe.return_value = pd.DataFrame({"dummy": [1]})

        # 2. Mock Annotator behavior
        annotator_instance = MockAnnotator.return_value
        annotator_instance.get_predicted_bound.return_value = pd.DataFrame(
            {"dummy_bound": [1]}
        )

        # 3. Mock Finder behavior — off-target at wrong coordinates
        finder_instance = MockFinder.return_value
        finder_instance.amplicon_df = pd.DataFrame(
            [
                {
                    "chrom": "chr7",
                    "F_primer": "P1_F",
                    "R_primer": "P1_R",
                    "product_bp": 500,
                    "F_start": 5000,  # Wrong coordinates (pseudogene)
                    "R_start": 5500,
                }
            ]
        )

        # Run
        run_specificity_check(mock_panel, str(tmp_path), "fake_genome.fa")

        # Verify AmpliconFinder was called with the target_map
        import unittest.mock

        MockFinder.assert_called_once_with(
            unittest.mock.ANY,  # bound_df
            target_map={"P1_F": "JUNCTION1", "P1_R": "JUNCTION1"},
        )

        # Verify
        pair = mock_panel.junctions[0].primer_pairs[0]
        assert pair.specificity_checked is True
        assert len(pair.off_target_products) == 1
        assert pair.off_target_products[0]["product_bp"] == 500
        # Coordinates are off-target (5000 vs expected ~1010) so on-target not detected
        assert pair.on_target_detected is False


def test_run_specificity_check_forwards_num_threads(mock_panel, tmp_path):
    with (
        patch("plexus.blast.specificity.BlastRunner") as MockRunner,
        patch("plexus.blast.specificity.BlastResultsAnnotator") as MockAnnotator,
        patch("plexus.blast.specificity.AmpliconFinder") as MockFinder,
        patch("os.makedirs"),
    ):
        runner_instance = MockRunner.return_value
        runner_instance.get_dataframe.return_value = pd.DataFrame({"dummy": [1]})

        annotator_instance = MockAnnotator.return_value
        annotator_instance.get_predicted_bound.return_value = pd.DataFrame(
            {"dummy_bound": [1]}
        )

        finder_instance = MockFinder.return_value
        finder_instance.amplicon_df = pd.DataFrame()

        run_specificity_check(mock_panel, str(tmp_path), "genome.fa", num_threads=6)

        _, kwargs = runner_instance.run.call_args
        assert kwargs.get("num_threads") == 6


def test_run_specificity_check_forwards_blast_parameters(mock_panel, tmp_path):
    """Custom BLAST parameters are threaded through to annotator and finder."""
    with (
        patch("plexus.blast.specificity.BlastRunner") as MockRunner,
        patch("plexus.blast.specificity.BlastResultsAnnotator") as MockAnnotator,
        patch("plexus.blast.specificity.AmpliconFinder") as MockFinder,
        patch("os.makedirs"),
    ):
        runner_instance = MockRunner.return_value
        runner_instance.get_dataframe.return_value = pd.DataFrame({"dummy": [1]})

        annotator_instance = MockAnnotator.return_value
        annotator_instance.get_predicted_bound.return_value = pd.DataFrame(
            {"dummy_bound": [1]}
        )

        finder_instance = MockFinder.return_value
        finder_instance.amplicon_df = pd.DataFrame()

        run_specificity_check(
            mock_panel,
            str(tmp_path),
            "genome.fa",
            length_threshold=20,
            evalue_threshold=5.0,
            max_mismatches=1,
            max_amplicon_size=5000,
            blast_evalue=500.0,
            blast_word_size=11,
            blast_reward=2,
            blast_penalty=-3,
        )

        # Verify annotator received custom thresholds
        annotator_instance.build_annotation_dict.assert_called_once_with(
            length_threshold=20, evalue_threshold=5.0, max_mismatches=1
        )

        # Verify finder received custom max amplicon size
        finder_instance.find_amplicons.assert_called_once_with(max_size_bp=5000)

        # Verify blast parameters were forwarded to runner.run()
        _, run_kwargs = runner_instance.run.call_args
        assert run_kwargs.get("evalue") == 500.0
        assert run_kwargs.get("word_size") == 11
        assert run_kwargs.get("reward") == 2
        assert run_kwargs.get("penalty") == -3


def test_run_specificity_check_no_hits(mock_panel, tmp_path):
    with (
        patch("plexus.blast.specificity.BlastRunner") as MockRunner,
        patch("os.makedirs"),
    ):
        runner_instance = MockRunner.return_value
        runner_instance.get_dataframe.return_value = pd.DataFrame()  # Empty

        # Ensure it's False or unset initially
        pair = mock_panel.junctions[0].primer_pairs[0]
        pair.specificity_checked = False

        # Should return early without error
        run_specificity_check(mock_panel, str(tmp_path), "fake_genome.fa")

        assert pair.specificity_checked is False
        # on_target_detected should remain at its initial mock value (not set by the function)
        assert pair.on_target_detected != True  # noqa: E712  — mock hasn't been set to True


def test_run_specificity_check_on_target_detected(mock_panel, tmp_path):
    """When the amplicon is at the correct coordinates, on_target_detected is True."""
    with (
        patch("plexus.blast.specificity.BlastRunner") as MockRunner,
        patch("plexus.blast.specificity.BlastResultsAnnotator") as MockAnnotator,
        patch("plexus.blast.specificity.AmpliconFinder") as MockFinder,
        patch("os.makedirs"),
    ):
        runner_instance = MockRunner.return_value
        runner_instance.get_dataframe.return_value = pd.DataFrame({"dummy": [1]})

        annotator_instance = MockAnnotator.return_value
        annotator_instance.get_predicted_bound.return_value = pd.DataFrame(
            {"dummy_bound": [1]}
        )

        # Amplicon at the CORRECT coordinates:
        # F_start = design_start + fwd_start = 1000 + 10 = 1010
        # R_start = design_start + rev_start + rev_length - 1 = 1000 + 180 + 22 - 1 = 1201
        finder_instance = MockFinder.return_value
        finder_instance.amplicon_df = pd.DataFrame(
            [
                {
                    "chrom": "chr7",
                    "F_primer": "P1_F",
                    "R_primer": "P1_R",
                    "product_bp": 191,
                    "F_start": 1010,
                    "R_start": 1201,
                }
            ]
        )

        run_specificity_check(mock_panel, str(tmp_path), "fake_genome.fa")

        pair = mock_panel.junctions[0].primer_pairs[0]
        assert pair.specificity_checked is True
        assert pair.on_target_detected is True
        assert len(pair.off_target_products) == 0


def test_run_specificity_check_on_target_not_detected(mock_panel, tmp_path):
    """When no amplicon matches target coordinates, on_target_detected is False."""
    with (
        patch("plexus.blast.specificity.BlastRunner") as MockRunner,
        patch("plexus.blast.specificity.BlastResultsAnnotator") as MockAnnotator,
        patch("plexus.blast.specificity.AmpliconFinder") as MockFinder,
        patch("os.makedirs"),
    ):
        runner_instance = MockRunner.return_value
        runner_instance.get_dataframe.return_value = pd.DataFrame({"dummy": [1]})

        annotator_instance = MockAnnotator.return_value
        annotator_instance.get_predicted_bound.return_value = pd.DataFrame(
            {"dummy_bound": [1]}
        )

        # BLAST only found an off-target hit (pair maps somewhere, but not the intended locus)
        finder_instance = MockFinder.return_value
        finder_instance.amplicon_df = pd.DataFrame(
            [
                {
                    "chrom": "chr7",
                    "F_primer": "P1_F",
                    "R_primer": "P1_R",
                    "product_bp": 500,
                    "F_start": 9000,  # Wrong coordinates
                    "R_start": 9500,
                }
            ]
        )

        run_specificity_check(mock_panel, str(tmp_path), "fake_genome.fa")

        pair = mock_panel.junctions[0].primer_pairs[0]
        assert pair.specificity_checked is True
        assert pair.on_target_detected is False
        assert len(pair.off_target_products) == 1


# ---------------------------------------------------------------------------
# _is_on_target coordinate-based classification
# ---------------------------------------------------------------------------


class TestIsOnTarget:
    """Tests for coordinate-based on-target classification."""

    def _make_junction(self, chrom="chr7", design_start=1000):
        j = MagicMock()
        j.chrom = chrom
        j.design_start = design_start
        return j

    def _make_pair(self, fwd_start=10, fwd_length=22, rev_start=180, rev_length=22):
        p = MagicMock()
        p.forward.start = fwd_start
        p.forward.length = fwd_length
        p.reverse.start = rev_start
        p.reverse.length = rev_length
        return p

    def test_same_chrom_correct_coords_is_on_target(self):
        """Same chrom + correct coordinates -> on-target."""
        junction = self._make_junction(chrom="chr7", design_start=1000)
        pair = self._make_pair(
            fwd_start=10, fwd_length=22, rev_start=180, rev_length=22
        )
        prod = {
            "chrom": "chr7",
            # F_start = design_start + fwd_start = 1000 + 10 = 1010
            "F_start": 1010,
            # R_start = design_start + rev_start + rev_length - 1 = 1000 + 180 + 22 - 1 = 1201
            "R_start": 1201,
        }
        assert _is_on_target(prod, junction, pair) is True

    def test_same_chrom_wrong_coords_is_off_target(self):
        """Same chrom + wrong coordinates (pseudogene) -> off-target."""
        junction = self._make_junction(chrom="chr7", design_start=1000)
        pair = self._make_pair(
            fwd_start=10, fwd_length=22, rev_start=180, rev_length=22
        )
        prod = {
            "chrom": "chr7",
            "F_start": 5000,
            "R_start": 5191,
        }
        assert _is_on_target(prod, junction, pair) is False

    def test_same_chrom_similar_product_size_different_coords(self):
        """Same chrom + similar product size but different coords -> off-target (the old bug)."""
        junction = self._make_junction(chrom="chr7", design_start=1000)
        pair = self._make_pair(
            fwd_start=10, fwd_length=22, rev_start=180, rev_length=22
        )
        prod = {
            "chrom": "chr7",
            "product_bp": 191,
            "F_start": 50000,
            "R_start": 50191,
        }
        assert _is_on_target(prod, junction, pair) is False

    def test_different_chrom_is_off_target(self):
        """Different chromosome -> always off-target."""
        junction = self._make_junction(chrom="chr7", design_start=1000)
        pair = self._make_pair(
            fwd_start=10, fwd_length=22, rev_start=180, rev_length=22
        )
        prod = {
            "chrom": "chr1",
            "F_start": 1010,
            "R_start": 1201,
        }
        assert _is_on_target(prod, junction, pair) is False

    def test_tolerance_allows_small_offset(self):
        """Coordinates within tolerance (5bp) -> on-target."""
        junction = self._make_junction(chrom="chr7", design_start=1000)
        pair = self._make_pair(
            fwd_start=10, fwd_length=22, rev_start=180, rev_length=22
        )
        prod = {
            "chrom": "chr7",
            "F_start": 1013,  # 3bp off from expected 1010
            "R_start": 1204,  # 3bp off from expected 1201
        }
        assert _is_on_target(prod, junction, pair) is True

    def test_beyond_tolerance_is_off_target(self):
        """Coordinates beyond tolerance (>5bp) -> off-target."""
        junction = self._make_junction(chrom="chr7", design_start=1000)
        pair = self._make_pair(
            fwd_start=10, fwd_length=22, rev_start=180, rev_length=22
        )
        prod = {
            "chrom": "chr7",
            "F_start": 1016,  # 6bp off from expected 1010
            "R_start": 1201,  # correct reverse
        }
        assert _is_on_target(prod, junction, pair) is False

    def test_custom_tolerance(self):
        """A larger tolerance classifies previously off-target hits as on-target."""
        junction = self._make_junction(chrom="chr7", design_start=1000)
        pair = self._make_pair(
            fwd_start=10, fwd_length=22, rev_start=180, rev_length=22
        )
        prod = {
            "chrom": "chr7",
            "F_start": 1020,  # 10bp off from expected 1010
            "R_start": 1211,  # 10bp off from expected 1201
        }
        # Default tolerance=5 -> off-target
        assert _is_on_target(prod, junction, pair) is False
        # Custom tolerance=10 -> on-target
        assert _is_on_target(prod, junction, pair, tolerance=10) is True

    def test_missing_design_start_defaults_to_zero(self):
        """Junction without design_start uses 0 as default."""
        junction = self._make_junction(chrom="chr7", design_start=None)
        pair = self._make_pair(
            fwd_start=10, fwd_length=22, rev_start=180, rev_length=22
        )
        prod = {
            "chrom": "chr7",
            # F_start = 0 + 10 = 10
            "F_start": 10,
            # R_start = 0 + 180 + 22 - 1 = 201
            "R_start": 201,
        }
        assert _is_on_target(prod, junction, pair) is True


# ---------------------------------------------------------------------------
# filter_offtarget_pairs
# ---------------------------------------------------------------------------


class TestFilterOfftargetPairs:
    """Tests for filter_offtarget_pairs()."""

    @staticmethod
    def _make_pair(pair_id, off_targets=0):
        p = MagicMock()
        p.pair_id = pair_id
        p.off_target_products = [{"dummy": i} for i in range(off_targets)]
        return p

    @staticmethod
    def _make_junction(name, pairs):
        j = MagicMock()
        j.name = name
        j.primer_pairs = list(pairs)
        return j

    def test_all_clean_no_removal(self):
        """When no pairs have off-targets, nothing is removed."""
        panel = MagicMock(spec=MultiplexPanel)
        panel.junctions = [
            self._make_junction("J1", [self._make_pair("P1"), self._make_pair("P2")]),
        ]

        removed, fallbacks = filter_offtarget_pairs(panel)

        assert removed == 0
        assert fallbacks == []
        assert len(panel.junctions[0].primer_pairs) == 2

    def test_mixed_dirty_removed(self):
        """Dirty pairs are removed when clean alternatives exist."""
        clean1 = self._make_pair("P1", off_targets=0)
        clean2 = self._make_pair("P2", off_targets=0)
        dirty = self._make_pair("P3", off_targets=2)

        panel = MagicMock(spec=MultiplexPanel)
        panel.junctions = [
            self._make_junction("J1", [clean1, dirty, clean2]),
        ]

        removed, fallbacks = filter_offtarget_pairs(panel)

        assert removed == 1
        assert fallbacks == []
        assert panel.junctions[0].primer_pairs == [clean1, clean2]

    def test_all_dirty_fallback_keeps_least_affected(self):
        """When all pairs have off-targets, the one with fewest is kept."""
        worst = self._make_pair("P1", off_targets=5)
        bad = self._make_pair("P2", off_targets=3)
        least = self._make_pair("P3", off_targets=1)

        panel = MagicMock(spec=MultiplexPanel)
        panel.junctions = [
            self._make_junction("J1", [worst, bad, least]),
        ]

        removed, fallbacks = filter_offtarget_pairs(panel)

        assert removed == 2
        assert fallbacks == ["J1"]
        assert panel.junctions[0].primer_pairs == [least]

    def test_empty_junction_no_crash(self):
        """A junction with no primer pairs is skipped gracefully."""
        panel = MagicMock(spec=MultiplexPanel)
        panel.junctions = [self._make_junction("J1", [])]

        removed, fallbacks = filter_offtarget_pairs(panel)

        assert removed == 0
        assert fallbacks == []

    def test_multiple_junctions(self):
        """Filter operates independently per junction."""
        # J1: mixed — dirty removed
        j1_clean = self._make_pair("J1_P1", off_targets=0)
        j1_dirty = self._make_pair("J1_P2", off_targets=3)

        # J2: all dirty — fallback
        j2_bad = self._make_pair("J2_P1", off_targets=4)
        j2_least = self._make_pair("J2_P2", off_targets=1)

        panel = MagicMock(spec=MultiplexPanel)
        panel.junctions = [
            self._make_junction("J1", [j1_clean, j1_dirty]),
            self._make_junction("J2", [j2_bad, j2_least]),
        ]

        removed, fallbacks = filter_offtarget_pairs(panel)

        assert removed == 2  # 1 from J1, 1 from J2
        assert fallbacks == ["J2"]
        assert panel.junctions[0].primer_pairs == [j1_clean]
        assert panel.junctions[1].primer_pairs == [j2_least]
