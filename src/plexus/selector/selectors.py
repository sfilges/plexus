# -------------------------------------------// (3) OPTIMISER //----------------------------------
#
# ThermoBLAST is the algorithm that is used to scan groups of primer candidates against the human
# genome to deduce background false amplicons penalties (the penalties are a function of the percent
# bound for each of the primers and the length of the false amplicon). ThermoBLAST is also used to
# exclude primer candidates that cross-hybridize with desired amplicons. These two ThermoBLAST runs
# are shown in solid green boxes in the flowchart figure above. The ThermoBLAST step maximizes
# specificity by penalizing primers that have strong off-target hybridizations. See the white paper about
# ThermoBLAST to learn more about how it works (download from www.dnasoftware.com).
#
# Lastly, MultiPick performs a mixing and matching of FP and RP candidates from different targets
# to make “multiplex solution sets”. MultiPick is run in 2 phases: first MultiPick analyzes all 100
# solutions from all the targets to determine all the exclusions and to then form a shorter list of
# 4 primer pair candidates for each of the N-plex targets. This shorter list of primer candidates
# is submitted to ThermoBLAST against the human genome to deduce all the possible false amplicons.
# For a 20-plex example, ThermoBLAST would consider 160 primers (i.e. 4 FP and 4 RP for 20 targets)
# against the human genome – such massive runs are only possible due to the high-performance
# cloud computing of Amazon Web Services and also clever algorithms written by DNA Software.
# In phase 2, MultiPick considers all 4N possible multiplex combinations that are possible and
# finds the combination(s) with the best score(s) consisting of Forward Primer and Reverse Primer
# (and optionally Probe and/or Reverse Transcription primer). For large N-plexes, the exponential
# explosion requires a sophisticated algorithm to solve the problem. DNA Software has developed s
# uch an algorithm – that uses a depth-first search with pruning algorithm, an approach used in
# 21st century artificial intelligence applications.
# -----------------------------------------------------------------------------------------------

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


# ================================================================================
# Collection of selection algorithms
#
# ================================================================================


selector_collection = {
    "Greedy": GreedySearch,
    "Random": RandomSearch,
    "BruteForce": BruteForce,
}
