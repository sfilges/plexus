from dataclasses import dataclass, field


@dataclass
class Multiplex:
    """A multiplex solution — one primer pair selected per target."""

    primer_pairs: list = field(default_factory=list)  # list of pair_id strings
    cost: float = 0.0
