from dataclasses import dataclass, field
from itertools import product
from typing import Tuple


@dataclass(unsafe_hash=True, order=True)
class Multiplex:
    """
    Encapsulate information about a multiplex in a manner
    that allows for sorting (primarily by cost) and hashing,
    such that a unique set of multiplexes can be found.

    To make hashable, the `.primer_pairs` attribute is
    co-erced into a tuple during `__post_init__`.

    """

    cost: float = field(compare=False)
    primer_pairs: Tuple[str]
    method: str = ""

    def __post_init__(self):
        if not isinstance(self.primer_pairs, Tuple):
            self.primer_pairs.sort()
            self.primer_pairs = tuple(self.primer_pairs)

    def get_primer_names(self):
        directions = ["F", "R"]
        return [f"{p}_{d}" for p, d in product(self.primer_pairs, directions)]
