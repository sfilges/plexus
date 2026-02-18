import math
import random
from abc import ABC, abstractmethod
from functools import reduce
from itertools import product

from loguru import logger

from .multiplex import Multiplex

# ================================================================================
# Abstract class for selection algorithm
#
# ================================================================================


class MultiplexSelector(ABC):
    def __init__(self, primer_df, cost_function):
        self.primer_df = primer_df
        self.cost_function = cost_function

    @abstractmethod
    def run(self):
        """
        Run the selection method

        """
        pass


# ================================================================================
# Concrete selection algorithms
#
# ================================================================================


class GreedySearch(MultiplexSelector):
    """
    Try to find the optimal multiplex using a greedy search algorithm

    Note, that it would be possible to compute exhaustively
    the number of possible permutations through the greedy algorithm;

    Alternatively; I could create a unique set of orders *instead*
    of shuffling; depending on the ratio of the number of iterations
    to the number of permutations, this would remove some redundant
    calculation

    """

    def run(self, N=10_000):
        """
        Run a greedy search algorithm for the lowest cost multiplex

        """

        # Get every UNIQUE primer pair, for each target
        # NB: from `primer_df` these are doubled, must use set()
        target_pairs = {
            target_id: list(set(target_df["pair_name"]))
            for target_id, target_df in self.primer_df.groupby("target_id")
        }

        # IDs
        target_ids = list(target_pairs)

        # Iterate
        multiplexes = []
        logger.info(
            f"Running greedy search with {N} iterations over {len(target_ids)} targets..."
        )
        for ix in range(N):
            # Prepare empty new multiplex
            multiplex = []

            # Shuffle the target IDs in place
            random.shuffle(target_ids)

            # Compute scores of each possible pair
            for target_id in target_ids:
                costs = [
                    self.cost_function.calc_cost(multiplex + [primer_pair])
                    for primer_pair in target_pairs[target_id]
                ]

                # Add max scoring from this step
                idxmax = costs.index(min(costs))
                multiplex.append(target_pairs[target_id][idxmax])

            # Add to list of all multiplexes
            multiplexes.append(Multiplex(cost=min(costs), primer_pairs=multiplex))

            if (ix + 1) % max(1, N // 10) == 0:
                logger.debug(f"Greedy search: {ix + 1}/{N} iterations complete")

        logger.info(f"Greedy search complete. Generated {len(multiplexes)} solutions.")
        return multiplexes


class BruteForce(MultiplexSelector):
    def run(self, store_maximum=200):
        """
        Run a brute force search for the highest scoring multiplex

        Note, we control the maximum number of multiplexes stored
        using the `store_maximum` argument; otherwise we can get
        exceptionally long lists of multiplexes

        """
        # Split target pairs into list of sets
        target_pairs = [
            set(target_df["pair_name"])
            for _, target_df in self.primer_df.groupby("target_id")
        ]

        # Compute number of iterations required
        total_N = reduce(lambda a, b: a * b, [len(t) for t in target_pairs])
        logger.info(
            f"Found {int(self.primer_df.shape[0] / 2)} primer pairs across {len(target_pairs)} targets."
        )
        logger.info(f"A total of {total_N} possible multiplexes exist.")

        # Iterate over all possible multiplexes
        stored_multiplexes = []
        stored_costs = []
        for ix, primer_pairs in enumerate(product(*target_pairs)):
            # Create the multiplex
            multiplex = Multiplex(
                primer_pairs=primer_pairs,
                cost=self.cost_function.calc_cost(primer_pairs),
            )

            # Store
            if len(stored_multiplexes) < store_maximum:
                stored_multiplexes.append(multiplex)
                stored_costs.append(multiplex.cost)
                highest_stored_cost = max(stored_costs)
            elif multiplex.cost < highest_stored_cost:
                stored_multiplexes.insert(
                    stored_costs.index(highest_stored_cost), multiplex
                )

            if (ix + 1) % max(1, total_N // 10) == 0:
                logger.debug(f"Brute force: {ix + 1}/{total_N} iterations complete")

        logger.info(
            f"Brute force complete. Stored {len(stored_multiplexes)} solutions."
        )
        return stored_multiplexes


class RandomSearch(MultiplexSelector):
    def run(self, N=10_000):
        """Run the random  selection algorithm"""
        # Get target pairs
        target_pairs = {
            target_id: list(set(target_df["pair_name"]))
            for target_id, target_df in self.primer_df.groupby("target_id")
        }

        # Iterate
        multiplexes = []
        logger.info(f"Running random search with {N} iterations...")
        for ix in range(N):
            # Randomly generate a multiplex
            multiplex = [random.choice(pairs) for _, pairs in target_pairs.items()]

            # Compute the cost
            cost = self.cost_function.calc_cost(multiplex)

            # Store
            multiplexes.append(Multiplex(cost=cost, primer_pairs=multiplex))

            if (ix + 1) % max(1, N // 10) == 0:
                logger.debug(f"Random search: {ix + 1}/{N} iterations complete")

        logger.info(f"Random search complete. Generated {len(multiplexes)} solutions.")
        return multiplexes


class SimulatedAnnealing(MultiplexSelector):
    """
    Simulated annealing selector that escapes local optima by
    probabilistically accepting worse solutions during cooling.

    Seeds each restart with a greedy solution, then iteratively
    swaps one target's primer pair per step.

    """

    def run(
        self,
        N_restarts=10,
        steps_per_restart=1000,
        T_initial=10.0,
        cooling_rate=0.995,
        T_min=0.01,
        greedy_seed_iterations=100,
    ):
        target_pairs = {
            target_id: list(set(target_df["pair_name"]))
            for target_id, target_df in self.primer_df.groupby("target_id")
        }
        target_ids = list(target_pairs)

        # Targets with more than one candidate can be swapped
        swappable = [tid for tid in target_ids if len(target_pairs[tid]) > 1]

        # Edge case: no swappable targets — only one possible solution
        if not swappable:
            only_solution = [target_pairs[tid][0] for tid in target_ids]
            cost = self.cost_function.calc_cost(only_solution)
            return [Multiplex(cost=cost, primer_pairs=only_solution)]

        multiplexes = []
        logger.info(
            f"Running simulated annealing with {N_restarts} restarts, "
            f"{steps_per_restart} steps each..."
        )
        for restart in range(N_restarts):
            # Seed with greedy solution
            current = self._greedy_seed(
                target_pairs, target_ids, greedy_seed_iterations
            )
            current_cost = self.cost_function.calc_cost(current)
            best = list(current)
            best_cost = current_cost

            T = T_initial
            for _ in range(steps_per_restart):
                if T < T_min:
                    break

                # Pick a random swappable target and swap to a different pair
                tid = random.choice(swappable)
                tid_idx = target_ids.index(tid)
                old_pair = current[tid_idx]
                alternatives = [p for p in target_pairs[tid] if p != old_pair]
                new_pair = random.choice(alternatives)

                # Evaluate neighbour
                neighbour = list(current)
                neighbour[tid_idx] = new_pair
                neighbour_cost = self.cost_function.calc_cost(neighbour)

                delta = neighbour_cost - current_cost
                if delta <= 0 or random.random() < math.exp(-delta / T):
                    current = neighbour
                    current_cost = neighbour_cost

                if current_cost < best_cost:
                    best = list(current)
                    best_cost = current_cost

                T *= cooling_rate

            multiplexes.append(Multiplex(cost=best_cost, primer_pairs=best))
            logger.debug(
                f"Simulated annealing restart {restart + 1}/{N_restarts}: "
                f"best cost = {best_cost:.4f}"
            )

        logger.info(
            f"Simulated annealing complete. Generated {len(multiplexes)} solutions."
        )
        return multiplexes

    def _greedy_seed(self, target_pairs, target_ids, N):
        """Run a small greedy search and return the best solution."""
        best = None
        best_cost = float("inf")
        for _ in range(N):
            ids = list(target_ids)
            random.shuffle(ids)
            multiplex = []
            for tid in ids:
                costs = [
                    self.cost_function.calc_cost(multiplex + [p])
                    for p in target_pairs[tid]
                ]
                multiplex.append(target_pairs[tid][costs.index(min(costs))])
            cost = self.cost_function.calc_cost(multiplex)
            if cost < best_cost:
                best = multiplex
                best_cost = cost
        return best


class DepthFirstSearch(MultiplexSelector):
    """
    Depth-first search with pruning for finding optimal multiplex solutions.

    Orders targets by ascending number of candidates and prunes branches
    whose partial cost plus a lower-bound on remaining cost exceeds the
    best known solution.

    """

    def run(
        self,
        store_maximum=200,
        seed_with_greedy=True,
        greedy_seed_iterations=100,
        max_nodes=10_000_000,
        target_ordering="ascending_candidates",
    ):
        target_pairs = {
            target_id: list(set(target_df["pair_name"]))
            for target_id, target_df in self.primer_df.groupby("target_id")
        }

        # Order targets
        if target_ordering == "ascending_candidates":
            ordered_targets = sorted(target_pairs, key=lambda t: len(target_pairs[t]))
        else:
            ordered_targets = list(target_pairs)

        n_targets = len(ordered_targets)

        # Sort candidates within each target by individual cost (cheapest first)
        sorted_candidates = {}
        min_individual_cost = {}
        for tid in ordered_targets:
            pairs = target_pairs[tid]
            pair_costs = [(p, self.cost_function.calc_cost([p])) for p in pairs]
            pair_costs.sort(key=lambda x: x[1])
            sorted_candidates[tid] = [p for p, _ in pair_costs]
            min_individual_cost[tid] = pair_costs[0][1]

        # Precompute suffix sums of minimum individual costs for lower bounds
        suffix_min = [0.0] * (n_targets + 1)
        for i in range(n_targets - 1, -1, -1):
            suffix_min[i] = suffix_min[i + 1] + min_individual_cost[ordered_targets[i]]

        # Optionally seed best_known_cost from greedy
        best_known_cost = float("inf")
        if seed_with_greedy:
            greedy = GreedySearch(self.primer_df, self.cost_function)
            greedy_results = greedy.run(N=greedy_seed_iterations)
            if greedy_results:
                best_known_cost = min(m.cost for m in greedy_results)

        # Iterative DFS with explicit stack
        # Stack entries: (depth, assignment_list, partial_cost)
        stored_multiplexes = []
        stored_costs = []
        nodes_visited = 0

        total_combos = reduce(
            lambda a, b: a * b, [len(target_pairs[t]) for t in ordered_targets]
        )
        logger.info(
            f"Running DFS over {n_targets} targets ({total_combos} total combinations), "
            f"max_nodes={max_nodes}..."
        )

        # Push initial candidates for depth 0 in reverse order (cheapest first via LIFO)
        stack = []
        tid0 = ordered_targets[0]
        for candidate in reversed(sorted_candidates[tid0]):
            stack.append((0, [candidate]))

        while stack:
            if nodes_visited >= max_nodes:
                logger.warning(
                    f"DFS hit max_nodes limit ({max_nodes}). Stopping search."
                )
                break

            depth, assignment = stack.pop()
            nodes_visited += 1

            partial_cost = self.cost_function.calc_cost(assignment)

            # Prune: if buffer is full, check lower bound
            if len(stored_multiplexes) >= store_maximum:
                lower_bound = partial_cost + suffix_min[depth + 1]
                if lower_bound >= best_known_cost:
                    continue

            # Complete solution
            if depth + 1 == n_targets:
                if len(stored_multiplexes) < store_maximum:
                    stored_multiplexes.append(
                        Multiplex(cost=partial_cost, primer_pairs=list(assignment))
                    )
                    stored_costs.append(partial_cost)
                    if partial_cost < best_known_cost:
                        best_known_cost = partial_cost
                elif partial_cost < max(stored_costs):
                    worst_idx = stored_costs.index(max(stored_costs))
                    stored_multiplexes[worst_idx] = Multiplex(
                        cost=partial_cost, primer_pairs=list(assignment)
                    )
                    stored_costs[worst_idx] = partial_cost
                    best_known_cost = min(best_known_cost, partial_cost)
                continue

            # Expand next level — push in reverse so cheapest is popped first
            next_tid = ordered_targets[depth + 1]
            for candidate in reversed(sorted_candidates[next_tid]):
                stack.append((depth + 1, assignment + [candidate]))

        logger.info(
            f"DFS complete. Visited {nodes_visited} nodes, "
            f"stored {len(stored_multiplexes)} solutions."
        )
        return stored_multiplexes


# ================================================================================
# Collection of selection algorithms
#
# ================================================================================


selector_collection = {
    "Greedy": GreedySearch,
    "Random": RandomSearch,
    "BruteForce": BruteForce,
    "SimulatedAnnealing": SimulatedAnnealing,
    "DFS": DepthFirstSearch,
}
