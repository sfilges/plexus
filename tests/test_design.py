"""Unit tests for plexus.designer.design module."""

from unittest.mock import MagicMock, patch

import pytest

from plexus.config import RescueTier
from plexus.designer.design import design_multiplex_primers, design_primers
from plexus.designer.multiplexpanel import Junction, MultiplexPanel, PrimerDesigns
from plexus.designer.primer import Primer, PrimerPair

# ================================================================================
# Helpers
# ================================================================================


def make_primer(
    name="p1",
    seq="ACGTACGTACGTACGTAC",
    direction="forward",
    start=10,
    length=18,
    tm=60.0,
    penalty=0.5,
):
    return Primer(
        name=name,
        seq=seq,
        direction=direction,
        start=start,
        length=length,
        tm=tm,
        penalty=penalty,
    )


def make_pair():
    fwd = make_primer("fwd", start=10, length=18)
    rev = make_primer("rev", direction="reverse", start=50, length=18)
    return PrimerPair(
        forward=fwd,
        reverse=rev,
        insert_size=22,
        amplicon_sequence="ACGT" * 20,
        amplicon_length=80,
    )


def _many_primers(n=150, direction="forward"):
    return [
        make_primer(name=f"p{i}", direction=direction, start=i + 1, length=18)
        for i in range(n)
    ]


# ================================================================================
# Fixtures
# ================================================================================


@pytest.fixture
def mock_config():
    config = MagicMock()
    s = config.singleplex_design_parameters
    s.primer_min_length = 18
    s.primer_max_length = 25
    s.primer_max_poly_x = 4
    s.primer_max_n = 0
    s.primer_min_gc = 30
    s.primer_max_gc = 70
    s.primer_gc_clamp = 0
    s.forward_tail = ""
    s.reverse_tail = ""
    s.primer_max_poly_gc = 3

    p = config.primer_pair_parameters
    p.PRIMER_PRODUCT_MIN_INSERT_SIZE = 40
    p.PRIMER_PRODUCT_MAX_SIZE = 300
    p.PRIMER_PAIR_MAX_DIFF_TM = 3.0
    p.PRIMER_PRODUCT_OPT_SIZE = 150
    p.PRIMER_PAIR_WT_PR_PENALTY = 1.0
    p.PRIMER_PAIR_WT_PRODUCT_SIZE_GT = 1.0
    p.PRIMER_PAIR_WT_PRODUCT_SIZE_LT = 1.0
    p.PRIMER_PAIR_WT_DIFF_TM = 1.0
    p.max_pairs_pre_filter = None  # No cap in tests by default
    p.max_pairs_per_junction = None

    # Rescue defaults
    config.enable_rescue = True
    config.rescue_tiers = [
        RescueTier(
            PRIMER_MIN_TM=56.0,
            PRIMER_MAX_TM=64.0,
            PRIMER_PAIR_MAX_DIFF_TM=4.0,
            PRIMER_MAX_HAIRPIN_TH=30.0,
            PRIMER_MAX_END_STABILITY=5.0,
            description="Tier 1",
        ),
        RescueTier(
            PRIMER_MIN_TM=55.0,
            PRIMER_MAX_TM=65.0,
            PRIMER_PAIR_MAX_DIFF_TM=5.0,
            PRIMER_MAX_HAIRPIN_TH=35.0,
            PRIMER_MAX_END_STABILITY=5.5,
            PRIMER_PRODUCT_MAX_SIZE=130,
            description="Tier 2",
        ),
    ]
    config.apply_rescue_tier = MagicMock(side_effect=lambda idx: config)
    return config


@pytest.fixture
def minimal_junction():
    """Junction with a 401-char design region and symmetric flanks around position 201."""
    design_region = ("ACGT" * 101)[:401]
    return Junction(
        name="J1",
        chrom="chr1",
        start=100,
        end=101,
        design_region=design_region,
        junction_length=401,
        jmin_coordinate=201,
        jmax_coordinate=202,
    )


@pytest.fixture
def minimal_panel(minimal_junction, mock_config):
    panel = MultiplexPanel("test_panel")
    panel.junctions = [minimal_junction]
    panel.config = mock_config
    return panel


# ================================================================================
# Tests
# ================================================================================


class TestDesignPrimersRouter:
    def test_plexus_method_delegates(self, minimal_panel):
        with patch("plexus.designer.design.design_multiplex_primers") as mock_dmp:
            mock_dmp.return_value = minimal_panel
            design_primers(minimal_panel, method="plexus")
            mock_dmp.assert_called_once_with(minimal_panel, on_junction_done=None)

    def test_unknown_method_raises_value_error(self, minimal_panel):
        with pytest.raises(ValueError, match="Unknown design method"):
            design_primers(minimal_panel, method="xyz")


class TestRegionExtraction:
    """Verify generate_kmers is called with correctly-sliced regions."""

    def _run_happy_path(self, mock_kmers, mock_thermo, minimal_junction, minimal_panel):
        primers = _many_primers(150)
        mock_kmers.return_value = primers
        mock_thermo.return_value = (primers, "eval string")
        minimal_junction.find_primer_pairs = MagicMock(return_value=[make_pair()])
        design_multiplex_primers(minimal_panel)

    @patch("plexus.designer.design.calculate_single_primer_thermodynamics")
    @patch("plexus.designer.design.generate_kmers")
    def test_forward_sequence_sliced_correctly(
        self, mock_kmers, mock_thermo, minimal_junction, minimal_panel
    ):
        self._run_happy_path(mock_kmers, mock_thermo, minimal_junction, minimal_panel)

        fwd_call = mock_kmers.call_args_list[0]
        expected = minimal_junction.design_region[0 : minimal_junction.jmin_coordinate]
        assert fwd_call.kwargs["target_sequence"] == expected

    @patch("plexus.designer.design.calculate_single_primer_thermodynamics")
    @patch("plexus.designer.design.generate_kmers")
    def test_reverse_sequence_is_rc_of_right_region(
        self, mock_kmers, mock_thermo, minimal_junction, minimal_panel
    ):
        from plexus.utils.utils import reverse_complement

        self._run_happy_path(mock_kmers, mock_thermo, minimal_junction, minimal_panel)

        rev_call = mock_kmers.call_args_list[1]
        j = minimal_junction
        expected = reverse_complement(
            j.design_region[j.jmax_coordinate : j.junction_length]
        )
        assert rev_call.kwargs["target_sequence"] == expected

    @patch("plexus.designer.design.calculate_single_primer_thermodynamics")
    @patch("plexus.designer.design.generate_kmers")
    def test_forward_offset_is_1(
        self, mock_kmers, mock_thermo, minimal_junction, minimal_panel
    ):
        self._run_happy_path(mock_kmers, mock_thermo, minimal_junction, minimal_panel)
        assert mock_kmers.call_args_list[0].kwargs["position_offset"] == 0

    @patch("plexus.designer.design.calculate_single_primer_thermodynamics")
    @patch("plexus.designer.design.generate_kmers")
    def test_reverse_offset_is_jmax_coordinate(
        self, mock_kmers, mock_thermo, minimal_junction, minimal_panel
    ):
        self._run_happy_path(mock_kmers, mock_thermo, minimal_junction, minimal_panel)
        assert (
            mock_kmers.call_args_list[1].kwargs["position_offset"]
            == minimal_junction.jmax_coordinate
        )


class TestKmerValidation:
    """Verify kmer count checks produce the correct errors and warnings."""

    @patch("plexus.designer.design.generate_kmers")
    def test_empty_left_kmers_fails_junction(
        self, mock_kmers, minimal_junction, minimal_panel
    ):
        mock_kmers.return_value = []
        design_multiplex_primers(minimal_panel)

        assert minimal_junction in minimal_panel.failed_junctions
        assert minimal_junction._design_error

    @patch("plexus.designer.design.generate_kmers")
    def test_empty_right_kmers_fails_junction(
        self, mock_kmers, minimal_junction, minimal_panel
    ):
        many = _many_primers(150)
        mock_kmers.side_effect = [many, []]
        design_multiplex_primers(minimal_panel)

        assert minimal_junction in minimal_panel.failed_junctions
        assert minimal_junction._design_error

    def test_short_left_region_raises(self, mock_config):
        design_region = ("ACGT" * 101)[:401]
        junction = Junction(
            name="J_short_left",
            chrom="chr1",
            start=100,
            end=101,
            design_region=design_region,
            junction_length=401,
            jmin_coordinate=5,  # left_region = design_region[1:5] = 4 chars < 2*25
            jmax_coordinate=202,
        )
        panel = MultiplexPanel("test_panel")
        panel.junctions = [junction]
        panel.config = mock_config

        design_multiplex_primers(panel)

        assert junction in panel.failed_junctions
        assert junction._design_error

    def test_short_right_region_raises(self, mock_config):
        design_region = ("ACGT" * 101)[:401]
        junction = Junction(
            name="J_short_right",
            chrom="chr1",
            start=100,
            end=101,
            design_region=design_region,
            junction_length=401,
            jmin_coordinate=201,
            jmax_coordinate=396,  # right_region = design_region[396:401] = 5 chars < 2*25
        )
        panel = MultiplexPanel("test_panel")
        panel.junctions = [junction]
        panel.config = mock_config

        with patch("plexus.designer.design.generate_kmers") as mock_kmers:
            mock_kmers.return_value = _many_primers(150)
            design_multiplex_primers(panel)

        assert junction in panel.failed_junctions
        assert junction._design_error

    @patch("plexus.designer.design.calculate_single_primer_thermodynamics")
    @patch("plexus.designer.design.generate_kmers")
    def test_fewer_than_100_kmers_left_warns(
        self, mock_kmers, mock_thermo, minimal_junction, minimal_panel
    ):
        few = _many_primers(50)
        many = _many_primers(150)
        mock_kmers.side_effect = [few, many]
        mock_thermo.return_value = (many, "eval")
        minimal_junction.find_primer_pairs = MagicMock(return_value=[make_pair()])

        with pytest.warns(UserWarning, match="Fewer than 100 left kmers"):
            design_multiplex_primers(minimal_panel)

    @patch("plexus.designer.design.calculate_single_primer_thermodynamics")
    @patch("plexus.designer.design.generate_kmers")
    def test_fewer_than_100_kmers_right_warns(
        self, mock_kmers, mock_thermo, minimal_junction, minimal_panel
    ):
        many = _many_primers(150)
        few = _many_primers(50)
        mock_kmers.side_effect = [many, few]
        mock_thermo.return_value = (many, "eval")
        minimal_junction.find_primer_pairs = MagicMock(return_value=[make_pair()])

        with pytest.warns(UserWarning, match="Fewer than 100 right kmers"):
            design_multiplex_primers(minimal_panel)


class TestJunctionOutcomes:
    """Verify junction partitioning into panel.junctions / panel.failed_junctions."""

    def _setup_success(self, mock_kmers, mock_thermo, junction):
        primers = _many_primers(150)
        mock_kmers.return_value = primers
        mock_thermo.return_value = (primers, "eval")
        junction.find_primer_pairs = MagicMock(return_value=[make_pair()])

    @patch("plexus.designer.design.calculate_single_primer_thermodynamics")
    @patch("plexus.designer.design.generate_kmers")
    def test_successful_junction_retained_in_panel(
        self, mock_kmers, mock_thermo, minimal_junction, minimal_panel
    ):
        self._setup_success(mock_kmers, mock_thermo, minimal_junction)
        design_multiplex_primers(minimal_panel)
        assert minimal_junction in minimal_panel.junctions

    @patch("plexus.designer.design.generate_kmers")
    def test_failed_junction_excluded_from_panel(
        self, mock_kmers, minimal_junction, minimal_panel
    ):
        mock_kmers.return_value = []
        design_multiplex_primers(minimal_panel)
        assert minimal_junction not in minimal_panel.junctions

    @patch("plexus.designer.design.generate_kmers")
    def test_failed_junction_in_failed_junctions_list(
        self, mock_kmers, minimal_junction, minimal_panel
    ):
        mock_kmers.return_value = []
        design_multiplex_primers(minimal_panel)
        assert minimal_junction in minimal_panel.failed_junctions

    @patch("plexus.designer.design.generate_kmers")
    def test_failed_junction_has_design_error(
        self, mock_kmers, minimal_junction, minimal_panel
    ):
        mock_kmers.return_value = []
        design_multiplex_primers(minimal_panel)
        assert minimal_junction._design_error

    @patch("plexus.designer.design.calculate_single_primer_thermodynamics")
    @patch("plexus.designer.design.generate_kmers")
    def test_mixed_one_success_one_failure(self, mock_kmers, mock_thermo, mock_config):
        design_region = ("ACGT" * 101)[:401]

        good_junction = Junction(
            name="J_good",
            chrom="chr1",
            start=100,
            end=101,
            design_region=design_region,
            junction_length=401,
            jmin_coordinate=201,
            jmax_coordinate=202,
        )
        bad_junction = Junction(
            name="J_bad",
            chrom="chr1",
            start=200,
            end=201,
            design_region=design_region,
            junction_length=401,
            jmin_coordinate=5,  # short left region → raises before generate_kmers
            jmax_coordinate=202,
        )

        panel = MultiplexPanel("test_panel")
        panel.junctions = [good_junction, bad_junction]
        panel.config = mock_config

        primers = _many_primers(150)
        mock_kmers.return_value = primers
        mock_thermo.return_value = (primers, "eval")
        good_junction.find_primer_pairs = MagicMock(return_value=[make_pair()])

        design_multiplex_primers(panel)

        assert len(panel.junctions) == 1
        assert good_junction in panel.junctions
        assert len(panel.failed_junctions) == 1
        assert bad_junction in panel.failed_junctions


class TestPrimerDesignsAttached:
    """Verify primer_designs is correctly set on successfully designed junctions."""

    def _run_success(self, mock_kmers, mock_thermo, minimal_junction, minimal_panel):
        primers = _many_primers(150)
        mock_kmers.return_value = primers
        mock_thermo.return_value = (primers, "eval_detail")
        minimal_junction.find_primer_pairs = MagicMock(return_value=[make_pair()])
        design_multiplex_primers(minimal_panel)

    @patch("plexus.designer.design.calculate_single_primer_thermodynamics")
    @patch("plexus.designer.design.generate_kmers")
    def test_primer_designs_set_on_junction(
        self, mock_kmers, mock_thermo, minimal_junction, minimal_panel
    ):
        self._run_success(mock_kmers, mock_thermo, minimal_junction, minimal_panel)
        pd = minimal_junction.primer_designs
        assert isinstance(pd, PrimerDesigns)
        assert pd.name == f"{minimal_junction.name}_designs"
        assert pd.target == minimal_junction.name

    @patch("plexus.designer.design.calculate_single_primer_thermodynamics")
    @patch("plexus.designer.design.generate_kmers")
    def test_primer_designs_eval_string_contains_both_sides(
        self, mock_kmers, mock_thermo, minimal_junction, minimal_panel
    ):
        self._run_success(mock_kmers, mock_thermo, minimal_junction, minimal_panel)
        eval_string = minimal_junction.primer_designs.eval_string
        assert "LEFT PRIMERS" in eval_string
        assert "RIGHT PRIMERS" in eval_string

    @patch("plexus.designer.design.calculate_single_primer_thermodynamics")
    @patch("plexus.designer.design.generate_kmers")
    def test_primer_table_has_two_entries(
        self, mock_kmers, mock_thermo, minimal_junction, minimal_panel
    ):
        self._run_success(mock_kmers, mock_thermo, minimal_junction, minimal_panel)
        primer_table = minimal_junction.primer_designs.primer_table
        assert isinstance(primer_table, list)
        assert len(primer_table) == 2


class TestRescueDesign:
    """Tests for tiered rescue when default primer design fails."""

    @patch("plexus.designer.design.calculate_single_primer_thermodynamics")
    @patch("plexus.designer.design.generate_kmers")
    def test_rescue_triggered_when_default_fails(
        self, mock_kmers, mock_thermo, minimal_junction, minimal_panel
    ):
        """Default attempt yields 0 pairs, tier 1 rescue succeeds."""
        primers = _many_primers(150)
        mock_kmers.return_value = primers
        mock_thermo.return_value = (primers, "eval")

        # Default: 0 pairs; rescue: has pairs
        rescue_pair = make_pair()
        call_count = [0]

        def find_pairs_side_effect(**kwargs):
            call_count[0] += 1
            if call_count[0] == 1:
                return []  # Default attempt fails
            return [rescue_pair]  # Rescue succeeds

        minimal_junction.find_primer_pairs = MagicMock(
            side_effect=find_pairs_side_effect
        )
        design_multiplex_primers(minimal_panel)

        assert minimal_junction in minimal_panel.junctions
        assert minimal_junction._rescue_tier == 1
        assert rescue_pair.rescue_tier == 1

    @patch("plexus.designer.design.calculate_single_primer_thermodynamics")
    @patch("plexus.designer.design.generate_kmers")
    def test_rescue_tier2_triggered_when_tier1_fails(
        self, mock_kmers, mock_thermo, minimal_junction, minimal_panel
    ):
        """Default and tier 1 fail, tier 2 succeeds."""
        primers = _many_primers(150)
        mock_kmers.return_value = primers
        mock_thermo.return_value = (primers, "eval")

        rescue_pair = make_pair()
        call_count = [0]

        def find_pairs_side_effect(**kwargs):
            call_count[0] += 1
            if call_count[0] <= 2:
                return []  # Default + tier 1 fail
            return [rescue_pair]  # Tier 2 succeeds

        minimal_junction.find_primer_pairs = MagicMock(
            side_effect=find_pairs_side_effect
        )
        design_multiplex_primers(minimal_panel)

        assert minimal_junction in minimal_panel.junctions
        assert minimal_junction._rescue_tier == 2
        assert rescue_pair.rescue_tier == 2

    @patch("plexus.designer.design.calculate_single_primer_thermodynamics")
    @patch("plexus.designer.design.generate_kmers")
    def test_rescue_disabled_skips_retry(
        self, mock_kmers, mock_thermo, minimal_junction, minimal_panel
    ):
        """With rescue disabled, no retry occurs and junction fails."""
        minimal_panel.config.enable_rescue = False

        primers = _many_primers(150)
        mock_kmers.return_value = primers
        mock_thermo.return_value = (primers, "eval")
        minimal_junction.find_primer_pairs = MagicMock(return_value=[])

        design_multiplex_primers(minimal_panel)

        assert minimal_junction in minimal_panel.failed_junctions
        assert minimal_junction not in minimal_panel.junctions
        # find_primer_pairs should be called exactly once (no rescue)
        assert minimal_junction.find_primer_pairs.call_count == 1

    @patch("plexus.designer.design.calculate_single_primer_thermodynamics")
    @patch("plexus.designer.design.generate_kmers")
    def test_rescue_tags_all_pairs(
        self, mock_kmers, mock_thermo, minimal_junction, minimal_panel
    ):
        """All pairs from a rescued junction have the correct rescue_tier."""
        primers = _many_primers(150)
        mock_kmers.return_value = primers
        mock_thermo.return_value = (primers, "eval")

        pair1, pair2 = make_pair(), make_pair()
        pair2.pair_id = "pair2"
        call_count = [0]

        def find_pairs_side_effect(**kwargs):
            call_count[0] += 1
            if call_count[0] == 1:
                return []
            return [pair1, pair2]

        minimal_junction.find_primer_pairs = MagicMock(
            side_effect=find_pairs_side_effect
        )
        design_multiplex_primers(minimal_panel)

        assert pair1.rescue_tier == 1
        assert pair2.rescue_tier == 1

    @patch("plexus.designer.design.calculate_single_primer_thermodynamics")
    @patch("plexus.designer.design.generate_kmers")
    def test_rescue_exhausted_reports_design_error(
        self, mock_kmers, mock_thermo, minimal_junction, minimal_panel
    ):
        """All tiers fail — junction._design_error set with diagnostic info."""
        primers = _many_primers(150)
        mock_kmers.return_value = primers
        mock_thermo.return_value = (primers, "eval")
        minimal_junction.find_primer_pairs = MagicMock(return_value=[])

        design_multiplex_primers(minimal_panel)

        assert minimal_junction in minimal_panel.failed_junctions
        assert "no valid primer pairs" in minimal_junction._design_error

    @patch("plexus.designer.design.calculate_single_primer_thermodynamics")
    @patch("plexus.designer.design.generate_kmers")
    def test_kmers_not_regenerated_during_rescue(
        self, mock_kmers, mock_thermo, minimal_junction, minimal_panel
    ):
        """generate_kmers is called exactly twice (fwd+rev), regardless of rescue."""
        primers = _many_primers(150)
        mock_kmers.return_value = primers
        mock_thermo.return_value = (primers, "eval")
        minimal_junction.find_primer_pairs = MagicMock(return_value=[])

        design_multiplex_primers(minimal_panel)

        # 2 calls: forward + reverse — no extra calls during rescue
        assert mock_kmers.call_count == 2

    @patch("plexus.designer.design.calculate_single_primer_thermodynamics")
    @patch("plexus.designer.design.generate_kmers")
    def test_standard_design_has_rescue_tier_zero(
        self, mock_kmers, mock_thermo, minimal_junction, minimal_panel
    ):
        """Junctions that succeed on default design get rescue_tier=0."""
        primers = _many_primers(150)
        mock_kmers.return_value = primers
        mock_thermo.return_value = (primers, "eval")
        pair = make_pair()
        minimal_junction.find_primer_pairs = MagicMock(return_value=[pair])

        design_multiplex_primers(minimal_panel)

        assert minimal_junction._rescue_tier == 0
        assert pair.rescue_tier == 0
