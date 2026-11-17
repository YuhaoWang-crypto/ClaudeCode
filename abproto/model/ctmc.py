r"""Discrete sequence generation by a Continuous-Time Markov Chain (CTMC) with
EXACT motif hard constraints, as in Multiflow's discrete head (Campbell et al.,
"Generative Flows on Discrete State-Spaces", ICML 2024), extended with a
constrained rate matrix that *guarantees* zero developability liabilities.

Masking construction
---------------------
Alphabet = 20 amino acids + a MASK state M. The conditional path from source
(all-MASK) to data x_1 unmasks each position independently with schedule
kappa(t):

    p_t(x_t^i = x_1^i) = kappa(t),   p_t(x_t^i = M) = 1 - kappa(t),
    kappa(0) = 0,  kappa(1) = 1   (linear: kappa(t) = t).

A denoiser p_theta(x_1 | x_t, c) predicts the clean tokens. The generative
process only ever turns M into a concrete token (unmasking). For a masked
position i the transition rate M -> a is

    R_t(M -> a | i) = ( kappa'(t) / (1 - kappa(t)) ) * p_theta(x_1^i = a | x_t, c).

Euler step: in dt, position i unmasks with prob p_unmask = kappa'/(1-kappa) * dt,
and if it unmasks, token a ~ p_theta(. | x_t, c) restricted to allowed tokens.

Hard motif constraints (the key addition)
-----------------------------------------
A liability motif (N-glycosylation N-X-[S/T], deamidation N[GSTN], isomerization
D[GSTD], or any CDR cysteine) is only *created* once all positions of its window
are concrete. So we zero the rate of any transition M -> a at position i that
would COMPLETE a forbidden motif given the already-decided neighbours:

    R_t(M -> a | i) := 0   if assigning a at i completes a forbidden window.

This is an exact edit of the rate matrix: the generator can never reach a state
containing the motif. Contrast with rejection sampling, which only lowers its
probability. Feasibility is guaranteed for these sparse motifs (at the [S/T] or
N slot many tokens always remain), with a documented fallback otherwise.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Any, Callable, Optional

import numpy as np

AA = "ACDEFGHIKLMNPQRSTVWY"
MASK = len(AA)                 # index of the MASK token
AA_INDEX = {a: i for i, a in enumerate(AA)}

# A denoiser: (token_ids[L], t, cond) -> posterior over clean tokens, shape (L, 20).
Denoiser = Callable[[np.ndarray, float, Any], np.ndarray]


# --------------------------------------------------------------------------
# Motif constraints
# --------------------------------------------------------------------------
@dataclass
class MotifRule:
    """A forbidden motif as a window of per-position allowed-residue sets.

    ``window`` is a list of sets of amino acids; a length-L window is forbidden
    iff residue j is in window[j] for all j. Example N-glyc: [{'N'}, all\\{'P'},
    {'S','T'}].
    """

    name: str
    window: list[set[str]]

    def length(self) -> int:
        return len(self.window)

    def matches(self, residues: list[str]) -> bool:
        return all(r in allowed for r, allowed in zip(residues, self.window))


def _all_but(excluded: str) -> set[str]:
    return set(AA) - set(excluded)


DEFAULT_LIABILITIES: list[MotifRule] = [
    # N-linked glycosylation sequon N-X-[S/T], X != P
    MotifRule("nglyc", [{"N"}, _all_but("P"), {"S", "T"}]),
    # Asn deamidation hotspots N-[G/S/T/N]
    MotifRule("deamidation", [{"N"}, {"G", "S", "T", "N"}]),
    # Asp isomerization hotspots D-[G/S/T/D]
    MotifRule("isomerization", [{"D"}, {"G", "S", "T", "D"}]),
    # Free cysteine anywhere in a CDR (length-1 window)
    MotifRule("cys", [{"C"}]),
]


class MotifConstraint:
    """Checks, for a candidate assignment, whether it completes any forbidden
    motif given the currently decided tokens. Operates on a region (the CDR
    positions); positions outside ``region_mask`` are never constrained."""

    def __init__(self, rules: Optional[list[MotifRule]] = None):
        # Note: distinguish None (use defaults) from [] (explicitly no rules).
        self.rules = DEFAULT_LIABILITIES if rules is None else rules

    def allowed_tokens(self, tokens: np.ndarray, decided: np.ndarray, pos: int,
                       region_mask: np.ndarray) -> np.ndarray:
        """Return a boolean (20,) mask of amino acids permitted at ``pos``.

        A token is forbidden iff, combined with already-decided neighbours, it
        completes a forbidden window that lies entirely in the constrained region
        and whose other positions are all decided.
        """
        allowed = np.ones(len(AA), dtype=bool)
        if not region_mask[pos]:
            return allowed                      # framework position: unconstrained
        L = len(tokens)
        for rule in self.rules:
            w = rule.length()
            # Every window of this rule that contains ``pos``.
            for start in range(max(0, pos - w + 1), min(pos, L - w) + 1):
                idxs = list(range(start, start + w))
                if any(not region_mask[j] for j in idxs):
                    continue                    # window leaves the design region
                others = [j for j in idxs if j != pos]
                if any(not decided[j] for j in others):
                    continue                    # not all other slots decided yet
                offset = pos - start
                for cand in range(len(AA)):
                    if not allowed[cand]:
                        continue
                    residues = []
                    for j in idxs:
                        residues.append(AA[cand] if j == pos else AA[int(tokens[j])])
                    if rule.matches(residues):
                        allowed[cand] = False
        return allowed


# --------------------------------------------------------------------------
# Schedule
# --------------------------------------------------------------------------
def kappa(t: float) -> float:
    return float(np.clip(t, 0.0, 1.0))


def kappa_dot(t: float) -> float:
    return 1.0


# --------------------------------------------------------------------------
# Constrained CTMC sampler
# --------------------------------------------------------------------------
@dataclass
class ConstrainedCTMC:
    denoiser: Denoiser
    constraint: MotifConstraint = field(default_factory=MotifConstraint)
    steps: int = 100
    eps: float = 1e-6

    def sample(self, length: int, region_mask: np.ndarray, cond: Any = None,
               seed: Optional[int] = None,
               init_tokens: Optional[np.ndarray] = None) -> np.ndarray:
        r"""Generate a sequence of ``length`` tokens.

        ``region_mask`` (bool, length L) marks positions that are being designed
        (CDRs); positions outside are kept from ``init_tokens`` (the germline
        framework) and never unmasked/constrained. Returns token ids in [0,20).
        """
        rng = np.random.default_rng(seed)
        tokens = (np.full(length, MASK, dtype=int) if init_tokens is None
                  else np.array(init_tokens, dtype=int))
        decided = np.zeros(length, dtype=bool)
        # Design-region positions are ALWAYS (re)masked and regenerated; every
        # other position is frozen at its init value (the germline framework).
        for i in range(length):
            if region_mask[i]:
                tokens[i] = MASK
                decided[i] = False
            else:
                decided[i] = True

        dt = 1.0 / self.steps
        for k in range(self.steps):
            t = k * dt
            p_unmask = kappa_dot(t) / max(1.0 - kappa(t), self.eps) * dt
            post = self.denoiser(tokens, t, cond)      # (L, 20)
            for i in range(length):
                if decided[i] or not region_mask[i]:
                    continue
                if rng.random() >= min(p_unmask, 1.0):
                    continue
                tokens[i] = self._draw(tokens, decided, i, region_mask, post[i], rng)
                decided[i] = True

        # Final pass: force-decide any still-masked design position (t->1).
        if region_mask.any():
            post = self.denoiser(tokens, 1.0, cond)
            for i in range(length):
                if region_mask[i] and not decided[i]:
                    tokens[i] = self._draw(tokens, decided, i, region_mask, post[i], rng)
                    decided[i] = True
        return tokens

    def _draw(self, tokens, decided, pos, region_mask, posterior, rng) -> int:
        """Sample a token from the posterior restricted to motif-allowed tokens.

        This is where the rate matrix is edited: forbidden transitions get rate 0
        (posterior mass set to 0) and the rest is renormalised."""
        allowed = self.constraint.allowed_tokens(tokens, decided, pos, region_mask)
        p = np.array(posterior, dtype=float)
        p[~allowed] = 0.0
        s = p.sum()
        if s <= self.eps:
            # Degenerate: every token blocked (does not occur for the default
            # sparse liabilities). Fall back to any allowed token, else any.
            pool = np.where(allowed)[0]
            if pool.size == 0:
                pool = np.arange(len(AA))
            return int(rng.choice(pool))
        p /= s
        return int(rng.choice(len(AA), p=p))


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------
def tokens_to_seq(tokens: np.ndarray) -> str:
    return "".join(AA[int(t)] if int(t) != MASK else "-" for t in tokens)


def seq_to_tokens(seq: str) -> np.ndarray:
    return np.array([AA_INDEX.get(a, MASK) for a in seq], dtype=int)


def uniform_denoiser(tokens: np.ndarray, t: float, cond: Any = None) -> np.ndarray:
    """A trivial denoiser (uniform posterior) so the sampler runs with no trained
    model. Replace with the network's p_theta(x_1 | x_t, c)."""
    L = len(tokens)
    return np.full((L, len(AA)), 1.0 / len(AA))
