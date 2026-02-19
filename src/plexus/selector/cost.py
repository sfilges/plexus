from __future__ import annotations

from itertools import combinations

from plexus.aligner import PrimerDimerPredictor
from plexus.config import MultiplexPickerParameters
from plexus.designer.primer import PrimerPair


class MultiplexCostFunction:
    """Score a multiplex combination based on design penalties, off-targets, and cross-dimers.

    The selectors call ``calc_cost(primer_pair_ids)`` where *primer_pair_ids*
    is a list of ``pair_id`` strings (one per target).

    Parameters
    ----------
    pair_lookup : dict[str, PrimerPair]
        Mapping of pair_id -> PrimerPair object.
    config : MultiplexPickerParameters
        Weights and parameters for the cost function.
    """

    def __init__(
        self,
        pair_lookup: dict[str, PrimerPair],
        config: MultiplexPickerParameters,
    ) -> None:
        self.pair_lookup = pair_lookup
        self.wt_pair_penalty = config.wt_pair_penalty
        self.wt_off_target = config.wt_off_target
        self.wt_cross_dimer = config.wt_cross_dimer
        self.wt_pair_dimer = config.wt_pair_dimer
        self.wt_snp_penalty = config.wt_snp_penalty

        # Cache for pairwise dimer scores: (seq_a, seq_b) -> score
        self._dimer_cache: dict[tuple[str, str], float] = {}
        self._dimer_predictor = PrimerDimerPredictor()

    def calc_cost(self, primer_pair_ids: list[str]) -> float:
        """Calculate the total cost for a candidate multiplex.

        Parameters
        ----------
        primer_pair_ids : list[str]
            List of pair_id strings representing the candidate multiplex.

        Returns
        -------
        float
            Total cost (lower is better).
        """
        pairs = [self.pair_lookup[pid] for pid in primer_pair_ids]

        cost = 0.0

        # 1. Sum individual pair penalties and off-target penalties
        for pair in pairs:
            cost += self.wt_pair_penalty * (pair.pair_penalty or 0.0)
            cost += self.wt_snp_penalty * (pair.snp_penalty or 0.0)
            cost += self.wt_off_target * len(pair.off_target_products)
            cost += self.wt_pair_dimer * max(0, -(pair.dimer_score or 0.0))

        # 2. All-pairwise cross-dimer interactions
        if self.wt_cross_dimer > 0:
            cost += self.wt_cross_dimer * self._calc_cross_dimer_penalty(pairs)

        return cost

    def _calc_cross_dimer_penalty(self, pairs: list[PrimerPair]) -> float:
        """Compute sum of cross-dimer scores for all primer pairs in the multiplex.

        Every forward and reverse primer is checked against every other primer
        in the multiplex. Scores are cached by sequence pair to avoid
        redundant thermodynamic calculations.

        The dimer predictor returns a *dG-like* score where **lower** (more
        negative) values indicate stronger (worse) dimer formation. We negate
        the score so that worse dimers contribute a *higher* cost.
        """
        # Collect all unique primer sequences in this multiplex
        all_primers: list[tuple[str, str]] = []  # (seq, name)
        for pair in pairs:
            all_primers.append((pair.forward.seq, pair.forward.name))
            all_primers.append((pair.reverse.seq, pair.reverse.name))

        penalty = 0.0

        for (seq_a, name_a), (seq_b, name_b) in combinations(all_primers, 2):
            # Canonical key: sort to avoid (a,b) vs (b,a) duplicates
            cache_key = tuple(sorted((seq_a, seq_b)))

            if cache_key not in self._dimer_cache:
                self._dimer_predictor.set_primers(seq_a, seq_b, name_a, name_b)
                self._dimer_predictor.align()
                score = self._dimer_predictor.score or 0.0
                self._dimer_cache[cache_key] = score

            dimer_score = self._dimer_cache[cache_key]

            # Lower dimer_score = stronger dimer = worse.
            # Convert to positive penalty: negate so more negative -> higher cost.
            if dimer_score < 0:
                penalty += abs(dimer_score)

        return penalty
