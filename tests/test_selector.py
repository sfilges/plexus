# ================================================================================
# Tests for selector modules (Multiplex, GreedySearch, RandomSearch, BruteForce)
# ================================================================================

from unittest.mock import MagicMock, patch

import pandas as pd
import pytest

from plexus.config import MultiplexPickerParameters
from plexus.selector.cost import MultiplexCostFunction
from plexus.selector.multiplex import Multiplex
from plexus.selector.selectors import (
    BruteForce,
    DepthFirstSearch,
    GreedySearch,
    RandomSearch,
    SimulatedAnnealing,
)


class TestMultiplex:
    def test_defaults(self):
        m = Multiplex()
        assert m.primer_pairs == []
        assert m.cost == 0.0


@pytest.fixture
def selector_inputs():
    df = pd.DataFrame(
        {
            "target_id": ["T1", "T1", "T2", "T2"],
            "pair_name": ["P1a", "P1b", "P2a", "P2b"],
        }
    )
    cost_fn = MagicMock()
    cost_fn.calc_cost = MagicMock(return_value=1.0)
    return df, cost_fn


class TestGreedySearch:
    def test_returns_multiplexes(self, selector_inputs):
        df, cost_fn = selector_inputs
        results = GreedySearch(df, cost_fn).run(N=5)
        assert len(results) == 5
        assert all(isinstance(m, Multiplex) for m in results)
        assert all(len(m.primer_pairs) == 2 for m in results)


class TestRandomSearch:
    def test_returns_multiplexes(self, selector_inputs):
        df, cost_fn = selector_inputs
        results = RandomSearch(df, cost_fn).run(N=5)
        assert len(results) == 5
        assert all(isinstance(m, Multiplex) for m in results)


class TestBruteForce:
    def test_returns_multiplexes(self, selector_inputs):
        df, cost_fn = selector_inputs
        results = BruteForce(df, cost_fn).run(store_maximum=10)
        assert len(results) == 4  # 2 pairs x 2 pairs = 4 combos
        assert all(isinstance(m, Multiplex) for m in results)

    def test_store_maximum_respected(self, selector_inputs):
        """Verify buffer doesn't grow beyond store_maximum."""
        df, cost_fn = selector_inputs
        results = BruteForce(df, cost_fn).run(store_maximum=2)
        assert len(results) <= 2  # 4 combos, only 2 stored

    def test_store_maximum_replacement_path(self):
        """Larger fixture to exercise the replacement path properly."""
        # 3 targets x 3 pairs = 27 combos
        df = pd.DataFrame(
            {
                "target_id": ["T1"] * 3 + ["T2"] * 3 + ["T3"] * 3,
                "pair_name": [
                    "P1a",
                    "P1b",
                    "P1c",
                    "P2a",
                    "P2b",
                    "P2c",
                    "P3a",
                    "P3b",
                    "P3c",
                ],
            }
        )
        # Assign varying costs so the replacement path is exercised
        call_count = [0]

        def varying_cost(pairs):
            call_count[0] += 1
            return float(call_count[0] % 7)  # Cycles 1-6 then 0

        cost_fn = MagicMock()
        cost_fn.calc_cost = varying_cost

        results = BruteForce(df, cost_fn).run(store_maximum=5)
        assert len(results) <= 5
        assert all(isinstance(m, Multiplex) for m in results)


class TestSimulatedAnnealing:
    def test_returns_correct_number(self, selector_inputs):
        df, cost_fn = selector_inputs
        results = SimulatedAnnealing(df, cost_fn).run(
            N_restarts=3, steps_per_restart=50, greedy_seed_iterations=5
        )
        assert len(results) == 3
        assert all(isinstance(m, Multiplex) for m in results)
        assert all(len(m.primer_pairs) == 2 for m in results)

    def test_single_candidate_per_target(self):
        df = pd.DataFrame(
            {
                "target_id": ["T1", "T2"],
                "pair_name": ["P1", "P2"],
            }
        )
        cost_fn = MagicMock()
        cost_fn.calc_cost = MagicMock(return_value=0.5)
        results = SimulatedAnnealing(df, cost_fn).run(N_restarts=3)
        assert len(results) == 1
        assert results[0].primer_pairs == ["P1", "P2"]


class TestDepthFirstSearch:
    def test_returns_all_combinations(self, selector_inputs):
        df, cost_fn = selector_inputs
        results = DepthFirstSearch(df, cost_fn).run(
            seed_with_greedy=False, store_maximum=200
        )
        assert len(results) == 4  # 2 x 2 = 4 combinations
        assert all(isinstance(m, Multiplex) for m in results)

    def test_with_greedy_seed(self, selector_inputs):
        df, cost_fn = selector_inputs
        results = DepthFirstSearch(df, cost_fn).run(
            seed_with_greedy=True, greedy_seed_iterations=5, store_maximum=200
        )
        assert len(results) == 4
        assert all(isinstance(m, Multiplex) for m in results)

    def test_max_nodes_limit(self, selector_inputs):
        df, cost_fn = selector_inputs
        results = DepthFirstSearch(df, cost_fn).run(seed_with_greedy=False, max_nodes=2)
        # With max_nodes=2 we can only visit 2 nodes, so fewer solutions
        assert len(results) < 4


class TestSelectorReproducibility:
    @pytest.fixture
    def large_selector_inputs(self):
        # Larger fixture to make stochastic differences likely
        n_targets = 10
        df = pd.DataFrame(
            {
                "target_id": [f"T{i}" for i in range(n_targets) for _ in range(5)],
                "pair_name": [f"P{i}_{j}" for i in range(n_targets) for j in range(5)],
            }
        )
        cost_fn = MagicMock()
        # Random but consistent costs for pairs to ensure different paths
        costs = {p: (hash(p) % 100) / 10.0 for p in df["pair_name"]}
        cost_fn.calc_cost = lambda pairs: sum(costs[p] for p in pairs)
        return df, cost_fn

    def test_greedy_reproducibility(self, large_selector_inputs):
        df, cost_fn = large_selector_inputs
        seed = 42
        results1 = GreedySearch(df, cost_fn, seed=seed).run(N=5)
        results2 = GreedySearch(df, cost_fn, seed=seed).run(N=5)

        for m1, m2 in zip(results1, results2, strict=False):
            assert m1.primer_pairs == m2.primer_pairs
            assert m1.cost == pytest.approx(m2.cost)

    def test_random_reproducibility(self, large_selector_inputs):
        df, cost_fn = large_selector_inputs
        seed = 123
        results1 = RandomSearch(df, cost_fn, seed=seed).run(N=5)
        results2 = RandomSearch(df, cost_fn, seed=seed).run(N=5)

        for m1, m2 in zip(results1, results2, strict=False):
            assert m1.primer_pairs == m2.primer_pairs
            assert m1.cost == pytest.approx(m2.cost)

    def test_simulated_annealing_reproducibility(self, large_selector_inputs):
        df, cost_fn = large_selector_inputs
        seed = 777
        results1 = SimulatedAnnealing(df, cost_fn, seed=seed).run(
            N_restarts=2, steps_per_restart=10
        )
        results2 = SimulatedAnnealing(df, cost_fn, seed=seed).run(
            N_restarts=2, steps_per_restart=10
        )

        for m1, m2 in zip(results1, results2, strict=False):
            assert m1.primer_pairs == m2.primer_pairs
            assert m1.cost == pytest.approx(m2.cost)

    def test_stochastic_without_seed_differs(self, large_selector_inputs):
        """Verify that without a seed, runs are likely different."""
        import random

        df, cost_fn = large_selector_inputs

        # Run 1
        random.seed(1)
        results1 = GreedySearch(df, cost_fn).run(N=1)

        # Run 2
        random.seed(2)
        results2 = GreedySearch(df, cost_fn).run(N=1)

        assert results1[0].primer_pairs != results2[0].primer_pairs


# ================================================================================
# Tests for MultiplexCostFunction
# ================================================================================


def _make_mock_pair(
    pair_id, pair_penalty=0.0, off_targets=0, dimer_score=None, snp_penalty=0.0
):
    """Create a mock PrimerPair with specified attributes."""
    pair = MagicMock()
    pair.pair_id = pair_id
    pair.pair_penalty = pair_penalty
    pair.snp_penalty = snp_penalty
    pair.off_target_products = [object()] * off_targets
    pair.dimer_score = dimer_score
    pair.forward.seq = f"ACGT_{pair_id}_F"
    pair.forward.name = f"{pair_id}_F"
    pair.reverse.seq = f"ACGT_{pair_id}_R"
    pair.reverse.name = f"{pair_id}_R"
    return pair


class TestMultiplexCostFunction:
    def _make_cost_fn(self, pairs, **config_overrides):
        config = MultiplexPickerParameters(**config_overrides)
        pair_lookup = {p.pair_id: p for p in pairs}
        with patch.object(
            MultiplexCostFunction, "__init__", wraps=MultiplexCostFunction.__init__
        ):
            return MultiplexCostFunction(pair_lookup, config)

    def test_pair_dimer_strong_dimer_increases_cost(self):
        """A pair with strong dimer_score=-5.0 costs more than one with dimer_score=0.0."""
        pair_clean = _make_mock_pair("clean", pair_penalty=1.0, dimer_score=0.0)
        pair_dimer = _make_mock_pair("dimer", pair_penalty=1.0, dimer_score=-5.0)

        config = MultiplexPickerParameters(
            target_plexity=2,
            minimum_plexity=1,
            maximum_plexity=5,
            wt_pair_dimer=1.0,
            wt_cross_dimer=0.0,  # isolate intra-pair dimer contribution
        )
        pair_lookup = {
            "clean": pair_clean,
            "dimer": pair_dimer,
        }

        with patch("plexus.selector.cost.PrimerDimerPredictor"):
            cost_fn = MultiplexCostFunction(pair_lookup, config)

        cost_clean = cost_fn.calc_cost(["clean"])
        cost_dimer = cost_fn.calc_cost(["dimer"])

        assert cost_dimer > cost_clean
        assert cost_dimer - cost_clean == pytest.approx(5.0)

    def test_pair_dimer_positive_score_no_penalty(self):
        """A pair with dimer_score > 0 should NOT add to cost (max(0, ...) clamps it)."""
        pair = _make_mock_pair("pos", pair_penalty=0.0, dimer_score=3.0)

        config = MultiplexPickerParameters(
            target_plexity=2,
            minimum_plexity=1,
            maximum_plexity=5,
            wt_pair_dimer=1.0,
            wt_cross_dimer=0.0,
        )
        pair_lookup = {"pos": pair}

        with patch("plexus.selector.cost.PrimerDimerPredictor"):
            cost_fn = MultiplexCostFunction(pair_lookup, config)

        assert cost_fn.calc_cost(["pos"]) == pytest.approx(0.0)

    def test_pair_dimer_none_score_treated_as_zero(self):
        """A pair with dimer_score=None should not raise and contributes 0 cost."""
        pair = _make_mock_pair("none_score", pair_penalty=0.0, dimer_score=None)

        config = MultiplexPickerParameters(
            target_plexity=2,
            minimum_plexity=1,
            maximum_plexity=5,
            wt_pair_dimer=1.0,
            wt_cross_dimer=0.0,
        )
        pair_lookup = {"none_score": pair}

        with patch("plexus.selector.cost.PrimerDimerPredictor"):
            cost_fn = MultiplexCostFunction(pair_lookup, config)

        assert cost_fn.calc_cost(["none_score"]) == pytest.approx(0.0)

    def test_pair_dimer_weight_zero_disables_contribution(self):
        """Setting wt_pair_dimer=0.0 should suppress all intra-pair dimer cost."""
        pair = _make_mock_pair("weighted", pair_penalty=0.0, dimer_score=-10.0)

        config = MultiplexPickerParameters(
            target_plexity=2,
            minimum_plexity=1,
            maximum_plexity=5,
            wt_pair_dimer=0.0,
            wt_cross_dimer=0.0,
        )
        pair_lookup = {"weighted": pair}

        with patch("plexus.selector.cost.PrimerDimerPredictor"):
            cost_fn = MultiplexCostFunction(pair_lookup, config)

        assert cost_fn.calc_cost(["weighted"]) == pytest.approx(0.0)


# ================================================================================
# Tests for DFS target dropout
# ================================================================================


class TestDFSTargetDropout:
    def test_dropout_disabled_unchanged(self, selector_inputs):
        """allow_target_dropping=False produces identical results to default."""
        df, cost_fn = selector_inputs
        results_default = DepthFirstSearch(df, cost_fn).run(seed_with_greedy=False)
        results_nodrop = DepthFirstSearch(df, cost_fn).run(
            seed_with_greedy=False, allow_target_dropping=False
        )
        assert len(results_default) == len(results_nodrop)
        for m_def, m_nod in zip(
            sorted(results_default, key=lambda m: str(m.primer_pairs)),
            sorted(results_nodrop, key=lambda m: str(m.primer_pairs)),
            strict=False,
        ):
            assert m_def.primer_pairs == m_nod.primer_pairs
            assert m_def.dropped_targets == []
            assert m_nod.dropped_targets == []

    def test_drops_toxic_target(self):
        """DFS drops a target whose inclusion dramatically increases cost."""
        df = pd.DataFrame(
            {
                "target_id": ["T1", "T1", "T2", "T2", "T3", "T3"],
                "pair_name": ["P1a", "P1b", "P2a", "P2b", "P3a", "P3b"],
            }
        )

        def toxic_cost(pairs):
            base = len(pairs) * 0.1
            if any(p.startswith("P3") for p in pairs):
                base += 100.0
            return base

        cost_fn = MagicMock()
        cost_fn.calc_cost = toxic_cost

        results = DepthFirstSearch(df, cost_fn).run(
            seed_with_greedy=False,
            allow_target_dropping=True,
            dropout_penalty=5.0,
            min_target_floor=2,
        )
        best = min(results, key=lambda m: m.cost)
        # Best solution should drop T3 (cost 100) for penalty of 5
        assert len(best.primer_pairs) == 2
        assert "T3" in best.dropped_targets

    def test_min_target_floor_enforced(self):
        """Cannot drop below min_target_floor."""
        df = pd.DataFrame(
            {
                "target_id": ["T1", "T2", "T3"],
                "pair_name": ["P1", "P2", "P3"],
            }
        )
        cost_fn = MagicMock()
        cost_fn.calc_cost = MagicMock(return_value=1.0)

        results = DepthFirstSearch(df, cost_fn).run(
            seed_with_greedy=False,
            allow_target_dropping=True,
            dropout_penalty=0.01,
            min_target_floor=3,
        )
        for m in results:
            assert len(m.dropped_targets) == 0
            assert len(m.primer_pairs) == 3

    def test_multiplex_dropped_targets_default(self):
        """Multiplex backward compat: dropped_targets defaults to []."""
        m = Multiplex()
        assert m.dropped_targets == []
        m2 = Multiplex(cost=1.0, primer_pairs=["P1"])
        assert m2.dropped_targets == []

    def test_solutions_include_both_dropped_and_full(self):
        """Solutions buffer contains both full and reduced panels."""
        df = pd.DataFrame(
            {
                "target_id": ["T1", "T2", "T3"],
                "pair_name": ["P1", "P2", "P3"],
            }
        )

        def mild_cost(pairs):
            base = len(pairs) * 1.0
            if any(p == "P3" for p in pairs):
                base += 2.0
            return base

        cost_fn = MagicMock()
        cost_fn.calc_cost = mild_cost

        results = DepthFirstSearch(df, cost_fn).run(
            seed_with_greedy=False,
            allow_target_dropping=True,
            dropout_penalty=1.5,
            min_target_floor=2,
            store_maximum=200,
        )
        # Should have solutions with 3 targets AND solutions with 2 targets
        full_solutions = [m for m in results if len(m.dropped_targets) == 0]
        dropped_solutions = [m for m in results if len(m.dropped_targets) > 0]
        assert len(full_solutions) > 0
        assert len(dropped_solutions) > 0

    def test_adaptive_penalty_computation(self):
        """Verify _compute_dropout_penalty returns reasonable value."""
        df = pd.DataFrame(
            {
                "target_id": ["T1", "T2", "T3"],
                "pair_name": ["P1", "P2", "P3"],
            }
        )

        # Cost function where T3 contributes most of the cost
        def cost_with_t3_toxic(pairs):
            base = len(pairs) * 0.5
            if "P3" in pairs:
                base += 10.0
            return base

        cost_fn = MagicMock()
        cost_fn.calc_cost = cost_with_t3_toxic

        dfs = DepthFirstSearch(df, cost_fn)
        pair_to_target = {"P1": "T1", "P2": "T2", "P3": "T3"}
        ordered_targets = ["T1", "T2", "T3"]
        greedy_sol = Multiplex(primer_pairs=["P1", "P2", "P3"], cost=11.5)

        penalty = dfs._compute_dropout_penalty(
            pair_to_target, ordered_targets, greedy_sol, stringency=0.5
        )
        # T1 and T2 marginal costs are ~0.5 each (just their individual contribution)
        # T3 marginal cost is ~10.5 (the toxic contribution)
        # Sorted: [~0.5, ~0.5, ~10.5]. At stringency=0.5, idx=1 -> ~0.5
        assert penalty > 0
        assert penalty < 10.0  # should not be at the toxic level with stringency=0.5
