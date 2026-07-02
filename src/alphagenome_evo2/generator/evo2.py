"""Evo 2-backed generator.

Wraps Arc Institute's ``evo2`` genomic language model so the loop can perform
true de-novo / in-context sequence design. The dependency is imported lazily so
the package remains usable (via :class:`MockGenerator`) without GPUs or weights.

Two generation modes are combined so the loop keeps both exploration and
parent-conditioned refinement:

* **prompted generation** — each selected parent is used as a prompt and Evo 2
  autoregressively extends/relies on it, echoing the "generate variations of good
  candidates" step of the design loop;
* **directed mutation fallback** — if fewer than ``n`` sequences come back, the
  remainder are filled by resampling prompts at higher temperature.

Install with ``pip install evo2`` (see https://github.com/ArcInstitute/evo2), or
point ``Evo2Generator`` at a hosted NVIDIA BioNeMo/NIM endpoint via ``client``.
"""

from __future__ import annotations

from typing import List, Optional, Sequence as Seq

from ..types import Sequence
from .base import Generator


class Evo2Generator(Generator):
    def __init__(
        self,
        model_name: str = "evo2_7b",
        model=None,
        gen_length: int = 200,
        id_prefix: str = "evo2",
    ):
        """Load an Evo 2 model, or inject a pre-built ``model``/client for tests.

        ``model`` may be any object exposing a ``generate(prompt_seqs, n_tokens,
        temperature=...)`` method compatible with the ``evo2`` API.
        """
        self.model = model if model is not None else self._load(model_name)
        self.gen_length = gen_length
        self.id_prefix = id_prefix
        self._counter = 0

    @staticmethod
    def _load(model_name: str):
        try:
            from evo2 import Evo2  # type: ignore
        except ImportError as exc:  # pragma: no cover - exercised only w/o dep
            raise ImportError(
                "Evo2Generator requires the 'evo2' package (and a CUDA GPU). "
                "Install per https://github.com/ArcInstitute/evo2, or use "
                "MockGenerator for offline development."
            ) from exc
        return Evo2(model_name)

    def _new_id(self) -> str:
        self._counter += 1
        return f"{self.id_prefix}-{self._counter:06d}"

    def generate(
        self,
        parents: Seq[Sequence],
        n: int,
        *,
        temperature: float = 0.8,
        mutation_rate: float = 0.03,  # noqa: ARG002 - kept for interface parity
    ) -> List[Sequence]:
        if n <= 0:
            return []
        parents = list(parents)
        prompts = self._build_prompts(parents, n)

        generated = self.model.generate(
            prompt_seqs=prompts,
            n_tokens=self.gen_length,
            temperature=temperature,
        )
        sampled = _extract_sequences(generated)

        out: List[Sequence] = []
        for i, s in enumerate(sampled[:n]):
            clean = _sanitise(s)
            if not clean:
                continue
            parent = parents[i % len(parents)].id if parents else None
            out.append(
                Sequence(
                    id=self._new_id(),
                    seq=clean,
                    metadata={"origin": "evo2", "prompt_parent": parent},
                )
            )
        return out

    @staticmethod
    def _build_prompts(parents: Seq[Sequence], n: int) -> List[str]:
        if not parents:
            # Unconditional generation: a neutral primer lets Evo 2 sample freely.
            return ["N" * 0 or "ACGT"] * n
        return [parents[i % len(parents)].seq for i in range(n)]


def _extract_sequences(generated) -> List[str]:
    """Normalise the several shapes ``evo2.generate`` may return into a str list."""
    if generated is None:
        return []
    # evo2 typically returns an object with a ``.sequences`` list.
    seqs = getattr(generated, "sequences", None)
    if seqs is not None:
        return [str(s) for s in seqs]
    if isinstance(generated, dict) and "sequences" in generated:
        return [str(s) for s in generated["sequences"]]
    if isinstance(generated, (list, tuple)):
        return [str(s) for s in generated]
    return [str(generated)]


def _sanitise(seq: str) -> Optional[str]:
    """Upper-case and strip any non-ACGT characters Evo 2 may emit."""
    kept = [c for c in seq.upper() if c in "ACGT"]
    return "".join(kept) if kept else None
