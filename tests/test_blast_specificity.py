from unittest.mock import MagicMock, patch

import pandas as pd
import pytest

from plexus.blast.specificity import _is_on_target, run_specificity_check
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
