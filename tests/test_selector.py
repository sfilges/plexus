# ================================================================================
# Tests for selector modules (Multiplex, GreedySearch, RandomSearch, BruteForce)
# ================================================================================

from unittest.mock import MagicMock

import pandas as pd
import pytest

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
