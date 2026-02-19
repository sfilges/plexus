# ================================================================================
# Tests for MultiplexCostFunction
# ================================================================================

import pytest

from plexus.config import MultiplexPickerParameters
from plexus.designer.primer import Primer, PrimerPair
from plexus.selector.cost import MultiplexCostFunction


def _make_primer(name="fwd", direction="forward", seq="ACGTACGTACGTACGTACGT"):
    return Primer(name=name, seq=seq, direction=direction, start=0, length=len(seq))


def _make_pair(pair_penalty=100.0, snp_penalty=0.0, pair_id="pair_0"):
    fwd = _make_primer("fwd", "forward")
    rev = _make_primer("rev", "reverse", seq="TGCATGCATGCATGCATGCA")
    pair = PrimerPair(
        forward=fwd,
        reverse=rev,
        insert_size=20,
        amplicon_sequence="A" * 60,
        amplicon_length=60,
        pair_penalty=pair_penalty,
        pair_id=pair_id,
    )
    pair.snp_penalty = snp_penalty
    return pair


def _make_config(**kwargs):
    defaults = dict(
        target_plexity=2,
        minimum_plexity=1,
        maximum_plexity=10,
        wt_pair_penalty=1.0,
        wt_off_target=0.0,
        wt_cross_dimer=0.0,
        wt_pair_dimer=0.0,
        wt_snp_penalty=3.0,
    )
    defaults.update(kwargs)
    return MultiplexPickerParameters(**defaults)


class TestSnpPenaltyCostTerm:
    def test_snp_penalty_used_independently(self):
        """wt_snp_penalty scales pair.snp_penalty as a separate cost term."""
        pair = _make_pair(pair_penalty=100.0, snp_penalty=30.0)
        config = _make_config(wt_pair_penalty=1.0, wt_snp_penalty=3.0)
        cf = MultiplexCostFunction({pair.pair_id: pair}, config)
        cost = cf.calc_cost([pair.pair_id])
        # 1.0 * 100 (pair) + 3.0 * 30 (snp) = 190
        assert cost == pytest.approx(190.0)

    def test_wt_snp_penalty_zero_excludes_snp(self):
        """Setting wt_snp_penalty=0 ignores SNP penalty entirely."""
        pair = _make_pair(pair_penalty=100.0, snp_penalty=30.0)
        config = _make_config(wt_pair_penalty=1.0, wt_snp_penalty=0.0)
        cf = MultiplexCostFunction({pair.pair_id: pair}, config)
        cost = cf.calc_cost([pair.pair_id])
        assert cost == pytest.approx(100.0)

    def test_no_snp_penalty_when_zero(self):
        """Pair with snp_penalty=0 contributes nothing from SNP term."""
        pair = _make_pair(pair_penalty=50.0, snp_penalty=0.0)
        config = _make_config(wt_snp_penalty=3.0)
        cf = MultiplexCostFunction({pair.pair_id: pair}, config)
        assert cf.calc_cost([pair.pair_id]) == pytest.approx(50.0)

    def test_snp_penalty_additive_across_pairs(self):
        """SNP cost accumulates correctly across multiple pairs."""
        p1 = _make_pair(pair_penalty=100.0, snp_penalty=30.0, pair_id="p1")
        p2 = _make_pair(pair_penalty=80.0, snp_penalty=0.0, pair_id="p2")
        config = _make_config(wt_pair_penalty=1.0, wt_snp_penalty=2.0)
        cf = MultiplexCostFunction({"p1": p1, "p2": p2}, config)
        cost = cf.calc_cost(["p1", "p2"])
        # 1.0 * 100 + 2.0 * 30 + 1.0 * 80 + 2.0 * 0 = 240
        assert cost == pytest.approx(240.0)
