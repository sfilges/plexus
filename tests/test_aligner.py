"""Tests for plexus.aligner.align — PrimerDimerPredictor and helpers."""

from __future__ import annotations

from importlib.resources import as_file, files

import pytest

from plexus.aligner.align import (
    PrimerAlignment,
    PrimerDimerPredictor,
    create_nn_score_dt,
)

# ---------------------------------------------------------------------------
# Shared paths (resolved once at module level)
# ---------------------------------------------------------------------------

_DATA_PKG = files("plexus.data")
_DOUBLE_MM_SCORE = 0.2

# All 16 Watson-Crick dinucleotide pairs present in match.json
_WC_PAIRS = {
    "AT/TA",
    "TA/AT",
    "AA/TT",
    "TT/AA",
    "AC/TG",
    "GT/CA",
    "CA/GT",
    "TG/AC",
    "TC/AG",
    "GA/CT",
    "AG/TC",
    "CT/GA",
    "CG/GC",
    "GC/CG",
    "GG/CC",
    "CC/GG",
}


# ============================================================================
# TestCreateNnScoreDict
# ============================================================================


class TestCreateNnScoreDict:
    """Tests for the create_nn_score_dt() helper."""

    @pytest.fixture(scope="class")
    def nn_scores(self) -> dict[str, float]:
        with (
            as_file(_DATA_PKG / "nn_model" / "match.json") as match_path,
            as_file(_DATA_PKG / "nn_model" / "single_mismatch.json") as mm_path,
        ):
            return create_nn_score_dt(
                match_json=str(match_path),
                single_mismatch_json=str(mm_path),
                double_mismatch_score=_DOUBLE_MM_SCORE,
            )

    def test_all_16_match_pairs_present(self, nn_scores):
        """All 16 Watson-Crick dinucleotide pairs appear as keys."""
        for pair in _WC_PAIRS:
            assert pair in nn_scores, f"Missing match pair: {pair}"

    def test_match_scores_are_negative(self, nn_scores):
        """All 16 Watson-Crick match-pair scores are < 0 (thermodynamically stable)."""
        for pair in _WC_PAIRS:
            assert (
                nn_scores[pair] < 0
            ), f"Expected negative score for {pair}, got {nn_scores[pair]}"

    def test_double_mismatch_score_applied(self, nn_scores):
        """A key not in either JSON (e.g. 'AA/CC') gets double_mismatch_score."""
        assert nn_scores["AA/CC"] == pytest.approx(_DOUBLE_MM_SCORE)

    def test_single_mismatch_overrides_double(self, nn_scores):
        """A key in single_mismatch JSON (e.g. 'AG/TT') does NOT equal double_mismatch_score."""
        # 'AG/TT' is present in single_mismatch.json with score 0.71
        assert nn_scores["AG/TT"] != pytest.approx(_DOUBLE_MM_SCORE)
        assert "AG/TT" in nn_scores


# ============================================================================
# TestLinearExtensionBonus
# ============================================================================


class TestLinearExtensionBonus:
    """Tests for PrimerDimerPredictor._calc_linear_extension_bonus() static method."""

    def test_no_overhang_gives_zero_bonus(self):
        bonus = PrimerDimerPredictor._calc_linear_extension_bonus(
            matching=[True, True, True, True],
            overhang_left=False,
            overhang_right=False,
            end_length=4,
            end_bonus=-0.5,
        )
        assert bonus == pytest.approx(0.0)

    def test_left_overhang_all_matching_gives_bonus(self):
        end_length = 4
        end_bonus = -0.5
        bonus = PrimerDimerPredictor._calc_linear_extension_bonus(
            matching=[True] * end_length,
            overhang_left=True,
            overhang_right=False,
            end_length=end_length,
            end_bonus=end_bonus,
        )
        assert bonus == pytest.approx(end_length * end_bonus)

    def test_right_overhang_all_matching_gives_bonus(self):
        end_length = 4
        end_bonus = -0.5
        bonus = PrimerDimerPredictor._calc_linear_extension_bonus(
            matching=[True] * end_length,
            overhang_left=False,
            overhang_right=True,
            end_length=end_length,
            end_bonus=end_bonus,
        )
        assert bonus == pytest.approx(end_length * end_bonus)

    def test_partial_match_stops_at_first_mismatch(self):
        """For left overhang: matching = [True, False, True, True] → only first counted."""
        bonus = PrimerDimerPredictor._calc_linear_extension_bonus(
            matching=[True, False, True, True],
            overhang_left=True,
            overhang_right=False,
            end_length=4,
            end_bonus=-0.5,
        )
        # Only the first True before the False contributes: 1 * (-0.5) = -0.5
        assert bonus == pytest.approx(1 * -0.5)

    def test_no_matches_gives_zero_bonus(self):
        """All matching=False → returns 0.0 regardless of overhang flags."""
        bonus = PrimerDimerPredictor._calc_linear_extension_bonus(
            matching=[False, False, False],
            overhang_left=True,
            overhang_right=True,
            end_length=3,
            end_bonus=-0.5,
        )
        assert bonus == pytest.approx(0.0)


# ============================================================================
# TestPrimerDimerPredictor
# ============================================================================


class TestPrimerDimerPredictor:
    """Tests for PrimerDimerPredictor using real parameter JSON files."""

    @pytest.fixture
    def predictor(self) -> PrimerDimerPredictor:
        return PrimerDimerPredictor()

    # ------------------------------------------------------------------
    # State / validation guards
    # ------------------------------------------------------------------

    def test_set_primers_stores_sequences(self, predictor):
        predictor.set_primers("ATCGATCG", "CGATCGAT", "fwd", "rev")
        assert predictor.primer1 == "ATCGATCG"
        assert predictor.primer2 == "CGATCGAT"
        assert predictor.primer1_name == "fwd"
        assert predictor.primer2_name == "rev"
        assert predictor.score is None  # reset on set_primers

    def test_align_without_set_primers_raises(self, predictor):
        with pytest.raises(ValueError):
            predictor.align()

    def test_get_alignment_string_before_align_raises(self, predictor):
        predictor.set_primers("ATCGATCG", "CGATCGAT", "fwd", "rev")
        with pytest.raises(ValueError):
            predictor.get_alignment_string()

    def test_get_primer_alignment_before_align_raises(self, predictor):
        predictor.set_primers("ATCGATCG", "CGATCGAT", "fwd", "rev")
        with pytest.raises(ValueError):
            predictor.get_primer_alignment()

    # ------------------------------------------------------------------
    # Scoring behaviour
    # ------------------------------------------------------------------

    def test_complementary_pair_scores_negative(self, predictor):
        """primer1 and its reverse complement → strongly negative score."""
        # RC("ATCGATCG") = "CGATCGAT"
        predictor.set_primers("ATCGATCG", "CGATCGAT", "p1", "p2")
        predictor.align()
        assert predictor.score < 0

    def test_unrelated_pair_scores_positive(self, predictor):
        """Poly-A vs poly-C: no Watson-Crick pairing → positive score."""
        predictor.set_primers("AAAAAAAAAA", "CCCCCCCCCC", "p1", "p2")
        predictor.align()
        assert predictor.score > 0

    def test_score_is_float_after_align(self, predictor):
        predictor.set_primers("ATCGATCG", "TTTTTTTT", "p1", "p2")
        predictor.align()
        assert isinstance(predictor.score, float)

    def test_equal_length_primers_no_crash(self, predictor):
        """Equal-length (8-base) primers should complete without error."""
        predictor.set_primers("ATCGATCG", "GCTAGCTA", "p1", "p2")
        predictor.align()  # should not raise

    def test_self_primer_dimer(self, predictor):
        """Self-complementary sequence (GCGCGCGC) should give score < 0."""
        # RC("GCGCGCGC") == "GCGCGCGC" (palindrome)
        predictor.set_primers("GCGCGCGC", "GCGCGCGC", "self", "self")
        predictor.align()
        assert predictor.score < 0

    # ------------------------------------------------------------------
    # Return-value types
    # ------------------------------------------------------------------

    def test_get_primer_alignment_fields(self, predictor):
        predictor.set_primers("ATCGATCG", "CGATCGAT", "fwd", "rev")
        predictor.align()
        pa = predictor.get_primer_alignment()
        assert isinstance(pa, PrimerAlignment)
        assert pa.primer1_name == "fwd"
        assert pa.primer2_name == "rev"
        assert pa.score == pytest.approx(predictor.score)

    def test_get_alignment_string_contains_sequence(self, predictor):
        predictor.set_primers("ATCGATCG", "CGATCGAT", "fwd", "rev")
        predictor.align()
        aln_str = predictor.get_alignment_string()
        # The longer sequence appears in 5'-...-3' notation
        assert "ATCGATCG" in aln_str or "CGATCGAT" in aln_str

    # ------------------------------------------------------------------
    # Caching
    # ------------------------------------------------------------------

    def test_parameter_cache_populated(self):
        """After construction the class-level cache should have at least one entry."""
        PrimerDimerPredictor._param_cache.clear()
        PrimerDimerPredictor()
        assert len(PrimerDimerPredictor._param_cache) >= 1

    # ------------------------------------------------------------------
    # Extension bonus integration
    # ------------------------------------------------------------------

    def test_extension_bonus_improves_score_for_3prime_overlap(self, predictor):
        """A primer pair where the 3′ end of the shorter primer is in the aligned
        region (overhang_left=True, consecutive matches) should score more negatively
        than a pair where the complementary region is at the 5′ end (no overhang)."""
        # primer2 is the reverse complement of ATCGATCG
        primer2 = "CGATCGAT"

        # primer1_with_3prime: non-complementary GCGC at 5' end, ATCGATCG at 3' end.
        # Best alignment: i=4 (overhang_left=True) → 3′ end of shorter is in region
        # → extension bonus applied → more negative score.
        predictor.set_primers("GCGCATCGATCG", primer2, "p1_3prime", "p2")
        predictor.align()
        score_with_3prime = predictor.score

        # primer1_no_3prime: ATCGATCG at 5' end, non-complementary GCGC at 3' end.
        # Best alignment: i=0 (overhang_left=False) → no extension bonus.
        predictor.set_primers("ATCGATCGGCGC", primer2, "p1_no3prime", "p2")
        predictor.align()
        score_without_3prime = predictor.score

        assert score_with_3prime < score_without_3prime
