"""Local Composition Perplexity (LCP) sequence-design restraint.

Implements the penalty exactly as defined in Figure 1 of the prompt release
bundle (``prompts/Figure 1.jpg``, section I.2.3, shared with permission from
Richard Shuai):

.. math::

    C_3 = \\frac{L}{L-w+1} \\sum_{i=1}^{L-w+1}
          \\left(e^{\\hat S} - e^{S_i}\\right)^2 \\Delta\\!\\left(S_i < \\hat S\\right)

where

* ``S_i`` is the Shannon entropy (in **nats**) of the amino-acid composition of
  the length-``w`` window starting at position ``i``,
* ``e^{S_i}`` is that window's *perplexity*,
* ``\\hat S`` is a threshold entropy, and ``\\Delta(S_i < \\hat S)`` is an
  indicator that switches the quadratic penalty on only for windows whose
  entropy fell *below* the threshold.

Figure 1's chosen constants are ``w = 30`` and ``\\hat S = 2.32`` nats (the 5th
percentile of 30-residue local-window entropies in PDB sequences).

The prompt calls LCP a "mandatory sequence-design restraint against homopolymer
stretches" and requires ``lcp_score`` to be recorded for every sequence
designed.  ``lcp_score`` here is :math:`C_3`; lower is better, ``0.0`` means no
window fell below the threshold.
"""

from __future__ import annotations

import math
from collections import Counter
from dataclasses import dataclass
from typing import Iterable, Sequence

__all__ = [
    "DEFAULT_WINDOW",
    "DEFAULT_THRESHOLD_ENTROPY",
    "window_entropies",
    "lcp_score",
    "lcp_per_position_penalty",
    "LCPReport",
    "lcp_report",
]

#: Window length ``w`` used in Figure 1.
DEFAULT_WINDOW = 30

#: Threshold entropy ``\hat S`` in nats: the 5th percentile of 30-residue local
#: window entropies in PDB sequences, per Figure 1.
DEFAULT_THRESHOLD_ENTROPY = 2.32


def _shannon_nats(window: Sequence[str]) -> float:
    """Shannon entropy, in nats, of the composition of ``window``."""
    n = len(window)
    if n == 0:
        return 0.0
    counts = Counter(window)
    h = 0.0
    for c in counts.values():
        p = c / n
        h -= p * math.log(p)
    return h


def window_entropies(
    sequence: str, window: int = DEFAULT_WINDOW
) -> list[float]:
    """Entropies (nats) of every length-``window`` sliding window of ``sequence``.

    Returns ``L - w + 1`` values.  For sequences shorter than ``window`` the
    single whole-sequence window is used, so short peptides still get a score
    rather than silently scoring 0.
    """
    seq = sequence.strip().upper()
    if window <= 0:
        raise ValueError("window must be positive")
    if len(seq) < window:
        return [_shannon_nats(seq)] if seq else []
    return [
        _shannon_nats(seq[i : i + window])
        for i in range(len(seq) - window + 1)
    ]


def lcp_score(
    sequence: str,
    window: int = DEFAULT_WINDOW,
    threshold_entropy: float = DEFAULT_THRESHOLD_ENTROPY,
) -> float:
    """``C_3`` of Figure 1 for ``sequence``.  Lower is better; 0.0 = unpenalised."""
    seq = sequence.strip().upper()
    if not seq:
        return 0.0
    entropies = window_entropies(seq, window)
    if not entropies:
        return 0.0

    L = len(seq)
    n_windows = len(entropies)
    perp_hat = math.exp(threshold_entropy)

    total = 0.0
    for s_i in entropies:
        if s_i < threshold_entropy:  # the indicator Delta(S_i < S_hat)
            total += (perp_hat - math.exp(s_i)) ** 2

    # Prefactor L / (L - w + 1).  n_windows is exactly L - w + 1 for L >= w; for
    # the short-sequence fallback it is 1, which makes the prefactor L and keeps
    # the same "sum scaled to sequence length" meaning.
    return (L / n_windows) * total


def lcp_per_position_penalty(
    sequence: str,
    window: int = DEFAULT_WINDOW,
    threshold_entropy: float = DEFAULT_THRESHOLD_ENTROPY,
) -> list[float]:
    """Attribute ``C_3`` to residues, for use as a per-position design restraint.

    Figure 1 defines the penalty per *window*.  A sequence-design model needs a
    per-*position* gradient/bias, so each window's penalty is spread evenly over
    the residues it covers.  Summing the returned vector reproduces
    :func:`lcp_score` up to floating point.
    """
    seq = sequence.strip().upper()
    if not seq:
        return []
    entropies = window_entropies(seq, window)
    if not entropies:
        return [0.0] * len(seq)

    L = len(seq)
    n_windows = len(entropies)
    w = min(window, L)
    prefactor = L / n_windows
    perp_hat = math.exp(threshold_entropy)

    per_pos = [0.0] * L
    for i, s_i in enumerate(entropies):
        if s_i >= threshold_entropy:
            continue
        pen = prefactor * (perp_hat - math.exp(s_i)) ** 2
        share = pen / w
        for j in range(i, i + w):
            per_pos[j] += share
    return per_pos


@dataclass(frozen=True)
class LCPReport:
    """Diagnostics behind one ``lcp_score``."""

    sequence_length: int
    window: int
    threshold_entropy: float
    n_windows: int
    n_windows_penalised: int
    min_window_entropy: float
    lcp_score: float

    @property
    def fraction_penalised(self) -> float:
        return self.n_windows_penalised / self.n_windows if self.n_windows else 0.0


def lcp_report(
    sequence: str,
    window: int = DEFAULT_WINDOW,
    threshold_entropy: float = DEFAULT_THRESHOLD_ENTROPY,
) -> LCPReport:
    """:func:`lcp_score` plus the window statistics that produced it."""
    entropies = window_entropies(sequence, window)
    return LCPReport(
        sequence_length=len(sequence.strip()),
        window=window,
        threshold_entropy=threshold_entropy,
        n_windows=len(entropies),
        n_windows_penalised=sum(1 for s in entropies if s < threshold_entropy),
        min_window_entropy=min(entropies) if entropies else 0.0,
        lcp_score=lcp_score(sequence, window, threshold_entropy),
    )


def rank_by_lcp(sequences: Iterable[str], **kw) -> list[tuple[str, float]]:
    """Sequences paired with their LCP, best (lowest) first."""
    scored = [(s, lcp_score(s, **kw)) for s in sequences]
    scored.sort(key=lambda t: t[1])
    return scored
