"""Generator interface: the Evo 2-shaped "designer" of the pipeline."""

from __future__ import annotations

import abc
from typing import List, Sequence as Seq

from ..types import Sequence


class Generator(abc.ABC):
    """Proposes new candidate sequences, optionally conditioned on parents.

    Concrete implementations wrap Evo 2 (:class:`Evo2Generator`) or a
    dependency-free evolutionary stand-in (:class:`MockGenerator`).
    """

    @abc.abstractmethod
    def generate(
        self,
        parents: Seq[Sequence],
        n: int,
        *,
        temperature: float = 0.8,
        mutation_rate: float = 0.03,
    ) -> List[Sequence]:
        """Return ``n`` new sequences seeded from ``parents``.

        ``parents`` may be empty for unconditional generation. Implementations
        must return exactly ``n`` sequences with unique ids.
        """
