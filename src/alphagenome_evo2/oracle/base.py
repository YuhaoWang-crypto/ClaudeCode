"""Oracle interface: the AlphaGenome-shaped "judge" of the pipeline."""

from __future__ import annotations

import abc
from typing import Dict, Iterable, List, Sequence as Seq

from ..types import Candidate, CellContext, OraclePrediction, Sequence


class Oracle(abc.ABC):
    """Predicts cell-context regulatory read-outs for DNA sequences.

    Concrete implementations wrap AlphaGenome (:class:`AlphaGenomeOracle`) or a
    deterministic stand-in (:class:`MockOracle`). Both must be safe to call on
    batches so the loop can amortise model/API overhead.
    """

    @abc.abstractmethod
    def predict(
        self,
        sequence: Sequence,
        contexts: Seq[CellContext],
    ) -> Dict[str, OraclePrediction]:
        """Return a ``{context.name: OraclePrediction}`` map for one sequence."""

    def predict_batch(
        self,
        sequences: Iterable[Sequence],
        contexts: Seq[CellContext],
    ) -> List[Candidate]:
        """Score a batch, returning candidates carrying their predictions.

        Override for backends with a native batch endpoint; the default loops.
        """
        out: List[Candidate] = []
        for s in sequences:
            preds = self.predict(s, contexts)
            out.append(Candidate(sequence=s, predictions=preds))
        return out
