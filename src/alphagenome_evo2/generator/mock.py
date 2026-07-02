"""A deterministic, dependency-free stand-in for Evo 2.

Rather than a 40B-parameter genomic language model, this performs classic
*directed evolution* on the parent pool:

* **point mutation** at ``mutation_rate`` per base (``temperature`` scales it);
* **recombination** (single-point crossover between two parents);
* **motif transplantation** — occasionally copy a short window from one parent
  into another, mimicking Evo 2's "in-context" reuse of functional grammar;
* **unconditional seeding** when the parent pool is empty.

It is fully seeded, so runs are reproducible, and it is good enough to let the
AlphaGenome↔Evo2 loop demonstrably climb the fitness surface. Swap in
:class:`~alphagenome_evo2.generator.evo2.Evo2Generator` for real generation.
"""

from __future__ import annotations

import random
from typing import List, Sequence as Seq

from ..types import Sequence
from .base import Generator

_BASES = "ACGT"


class MockGenerator(Generator):
    def __init__(self, seed: int = 0, default_length: int = 200):
        self._rng = random.Random(seed)
        self.default_length = default_length
        self._counter = 0

    def _new_id(self, prefix: str) -> str:
        self._counter += 1
        return f"{prefix}-{self._counter:06d}"

    def generate(
        self,
        parents: Seq[Sequence],
        n: int,
        *,
        temperature: float = 0.8,
        mutation_rate: float = 0.03,
    ) -> List[Sequence]:
        if n <= 0:
            return []
        parents = list(parents)
        rate = max(0.0, min(1.0, mutation_rate * (0.5 + temperature)))
        out: List[Sequence] = []
        for _ in range(n):
            if not parents:
                out.append(self._random_sequence(self.default_length))
                continue
            roll = self._rng.random()
            if len(parents) >= 2 and roll < 0.35:
                child = self._recombine(parents)
            elif len(parents) >= 2 and roll < 0.5:
                child = self._transplant(parents)
            else:
                child = self._mutate(self._rng.choice(parents), rate)
            out.append(child)
        return out

    def _random_sequence(self, length: int) -> Sequence:
        s = "".join(self._rng.choice(_BASES) for _ in range(length))
        return Sequence(id=self._new_id("denovo"), seq=s, metadata={"origin": "denovo"})

    def _mutate(self, parent: Sequence, rate: float) -> Sequence:
        chars = list(parent.seq)
        n_mut = 0
        for i, base in enumerate(chars):
            if self._rng.random() < rate:
                alt = self._rng.choice([b for b in _BASES if b != base])
                chars[i] = alt
                n_mut += 1
        # Rare indels keep length from being perfectly conserved.
        if self._rng.random() < 0.1 and len(chars) > 10:
            pos = self._rng.randrange(len(chars))
            if self._rng.random() < 0.5:
                chars.insert(pos, self._rng.choice(_BASES))
            else:
                del chars[pos]
        return Sequence(
            id=self._new_id("mut"),
            seq="".join(chars),
            metadata={"origin": "mutation", "parent": parent.id, "n_mutations": n_mut},
        )

    def _recombine(self, parents: Seq[Sequence]) -> Sequence:
        a, b = self._rng.sample(list(parents), 2)
        cut_a = self._rng.randrange(1, len(a.seq))
        cut_b = self._rng.randrange(1, len(b.seq))
        child = a.seq[:cut_a] + b.seq[cut_b:]
        return Sequence(
            id=self._new_id("xover"),
            seq=child,
            metadata={"origin": "recombination", "parents": [a.id, b.id]},
        )

    def _transplant(self, parents: Seq[Sequence]) -> Sequence:
        donor, host = self._rng.sample(list(parents), 2)
        win = min(len(donor.seq), max(6, len(donor.seq) // 5))
        d_start = self._rng.randrange(0, max(1, len(donor.seq) - win))
        motif = donor.seq[d_start:d_start + win]
        h_start = self._rng.randrange(0, max(1, len(host.seq) - win))
        child = host.seq[:h_start] + motif + host.seq[h_start + win:]
        return Sequence(
            id=self._new_id("graft"),
            seq=child,
            metadata={"origin": "transplant", "host": host.id, "donor": donor.id},
        )
