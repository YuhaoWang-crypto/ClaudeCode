"""Gradient-conflict-aware multi-task weighting (Eq. 12-18 of the paper).

The paper's scheme, verbatim:

    Loss = (1 - lambda) * Loss1 + lambda * Loss2                        (12)
    cosine_similarity = <grad Loss1, grad Loss2> / (|grad L1| |grad L2|) (13)
    clipped_similarity = min(cosine_similarity, 1.0)                     (14)
    adjustment_factor  = 0.01 * clipped_similarity                       (15)
    if gradients are similar:  weights unchanged                         (16)
    if gradients conflict:     w1 -= adj, w2 += adj                      (17, 18)

Two things are worth flagging, because they matter for anyone reproducing this:

1. Equations 14-15 make the *magnitude* of the adjustment proportional to the
   cosine similarity, but the adjustment is only applied when the gradients
   conflict, i.e. when the cosine similarity is negative.  Taken literally, a
   conflicting pair (cos < 0) yields a NEGATIVE adjustment_factor, which flips
   the signs of Eq. 17-18 and would *decrease* the weight of the drug-efficacy
   task -- the opposite of the stated intent ("we prioritize the drug efficacy
   prediction task by increasing its gradient weight").  We therefore use
   ``|adjustment_factor|``, matching the prose rather than the algebra, and
   expose ``literal_sign=True`` to reproduce the equations as written.

2. Eq. 14 clips from above at 1.0, which is a no-op for a cosine similarity
   (it is already bounded by 1).  It is implemented anyway for fidelity.

The weights are clamped to [0, 1] ("we ensure that the weights remain within a
reasonable range") and re-normalised so they continue to sum to 1.
"""

from __future__ import annotations

import torch

__all__ = ["GradientConflictWeighter", "flat_grad"]


def flat_grad(loss: torch.Tensor, params, retain_graph: bool = True) -> torch.Tensor:
    """Flattened gradient of ``loss`` w.r.t. ``params`` (zeros where unused)."""
    params = [p for p in params if p.requires_grad]
    grads = torch.autograd.grad(
        loss, params, retain_graph=retain_graph, allow_unused=True
    )
    return torch.cat(
        [
            (g if g is not None else torch.zeros_like(p)).reshape(-1)
            for g, p in zip(grads, params)
        ]
    )


class GradientConflictWeighter:
    """Adaptive task weights driven by gradient cosine similarity."""

    def __init__(
        self,
        lam: float = 0.8,
        step: float = 0.01,
        w_min: float = 0.0,
        w_max: float = 1.0,
        literal_sign: bool = False,
    ):
        # Eq. 12 initialisation: w1 = 1 - lambda, w2 = lambda, lambda = 0.8.
        self.w1 = 1.0 - lam
        self.w2 = lam
        self.step = step
        self.w_min = w_min
        self.w_max = w_max
        self.literal_sign = literal_sign
        self.history: list[dict[str, float]] = []

    def update(self, loss1: torch.Tensor, loss2: torch.Tensor, params) -> float:
        """Update the task weights in place; returns the cosine similarity."""
        g1 = flat_grad(loss1, params)
        g2 = flat_grad(loss2, params)
        denom = g1.norm() * g2.norm()
        if denom.item() == 0.0:
            cos = 0.0
        else:
            cos = float(torch.dot(g1, g2) / denom)

        clipped = min(cos, 1.0)               # Eq. 14
        adj = self.step * clipped             # Eq. 15
        if not self.literal_sign:
            adj = abs(adj)

        if cos < 0.0:                          # gradients conflict -> Eq. 17-18
            self.w1 -= adj
            self.w2 += adj
        # else: Eq. 16, weights unchanged.

        self.w1 = min(max(self.w1, self.w_min), self.w_max)
        self.w2 = min(max(self.w2, self.w_min), self.w_max)
        total = self.w1 + self.w2
        if total > 0:
            self.w1 /= total
            self.w2 /= total

        self.history.append({"cos": cos, "w1": self.w1, "w2": self.w2})
        return cos

    def combine(self, loss1: torch.Tensor, loss2: torch.Tensor) -> torch.Tensor:
        """Eq. 12 with the current (possibly adapted) weights."""
        return self.w1 * loss1 + self.w2 * loss2
