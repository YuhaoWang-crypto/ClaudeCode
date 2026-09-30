"""
bioif.registry -- routing, and choosing between competing models of one hop.

The payoff of typing the seams is that composition becomes a search problem
rather than an integration project. You register an adapter once; every chain
that needs that hop finds it.

Two adapters can consume and produce the same pair of quantities. That is the
normal case in practice -- several groups model the same step and disagree --
and it is where a registry earns its keep, because the choice can be made on
declared properties instead of on who wrote the pipeline:

  1. a candidate that REFUSES this particular claim is not a candidate;
     one that extrapolates is worse than one in domain;
  2. among those, the narrower calibrated interval wins -- an edge with a
     checked interval beats one with none, and algebra (width 0) beats both;
  3. ties go to the less assumption-bearing kind.

Rule 2 is only meaningful because `calibrated_width` is a declared property
of an adapter rather than a number in a paper. Note what it deliberately does
NOT do: it does not choose the model that gives the answer you like, and it
cannot see which model is *right* -- only which is narrower among models that
are all conformally valid. The residual disagreement between them is measured
separately (see real/transfer_models.disagreement_vs_width) and is not
something routing can resolve.

Path search prefers short paths, then paths with fewer assumption-bearing
bridges.
"""
from __future__ import annotations

import math
from collections import deque

from .adapter import Adapter, BRIDGE, COERCION, EMPIRICAL
from .core import Claim, EXTRAPOLATE, IN_DOMAIN, REFUSE, Quantity
from .quantities import LOSSLESS

_KIND_RANK = {COERCION: 0, EMPIRICAL: 1, BRIDGE: 2}
_STATUS_RANK = {IN_DOMAIN: 0, EXTRAPOLATE: 1, REFUSE: 2}


class Registry:
    def __init__(self, include_lossless: bool = True):
        self._by_input: dict[str, list[Adapter]] = {}
        self.adapters: list[Adapter] = []
        if include_lossless:
            for c in LOSSLESS:
                self.register(c)

    def register(self, adapter: Adapter) -> Adapter:
        self._by_input.setdefault(adapter.consumes.key, []).append(adapter)
        self.adapters.append(adapter)
        return adapter

    # -- competing models of one hop --------------------------------------
    def candidates(self, src: Quantity, dst: Quantity) -> list[Adapter]:
        """Every registered adapter for exactly this hop."""
        return [a for a in self._by_input.get(src.key, [])
                if a.produces.key == dst.key]

    def rank(self, cands: list[Adapter], claim: Claim | None = None,
             alpha: float = 0.2) -> list[tuple]:
        """Score competing adapters; lowest sorts first. Returns (score, adapter)."""
        rows = []
        for a in cands:
            status = a.domain(claim).status if claim is not None else IN_DOMAIN
            rows.append(((_STATUS_RANK[status], a.calibrated_width(alpha),
                          _KIND_RANK.get(a.kind, 9), a.name), a))
        rows.sort(key=lambda r: r[0])
        return rows

    def select(self, src: Quantity, dst: Quantity, claim: Claim | None = None,
               alpha: float = 0.2, allow_refused: bool = False) -> Adapter | None:
        """The registry's pick for one hop, given the claim that will cross it."""
        rows = self.rank(self.candidates(src, dst), claim, alpha)
        for score, a in rows:
            if allow_refused or score[0] < _STATUS_RANK[REFUSE]:
                return a
        return None

    # -- path search -------------------------------------------------------
    def route(self, src: Quantity, dst: Quantity, max_len: int = 8,
              claim: Claim | None = None, alpha: float = 0.2):
        """
        Breadth-first shortest path, tie-broken toward fewer bridges.

        With `claim`, parallel edges are collapsed by `select` using that
        claim, so a route is chosen with the actual data in hand: an edge
        that refuses this claim is not used, and among those that accept it
        the narrowest calibrated interval wins. Without a claim the same
        ordering applies minus the domain test.
        """
        if src.key == dst.key:
            return []
        best: list[Adapter] | None = None
        q: deque[tuple[str, list[Adapter]]] = deque([(src.key, [])])
        seen_cost: dict[str, int] = {src.key: 0}
        while q:
            node, path = q.popleft()
            if len(path) >= max_len:
                continue
            # group parallel edges by the quantity they produce, then pick one
            groups: dict[str, list[Adapter]] = {}
            for a in self._by_input.get(node, []):
                if any(a is p for p in path):
                    continue
                groups.setdefault(a.produces.key, []).append(a)
            for out_key, cands in groups.items():
                chosen = self._pick(cands, path, claim, alpha)
                if chosen is None:
                    continue
                npath = path + [chosen]
                if out_key == dst.key:
                    if best is None or _cost(npath) < _cost(best):
                        best = npath
                    continue
                cost = len(npath)
                if seen_cost.get(out_key, 10 ** 9) < cost:
                    continue
                seen_cost[out_key] = cost
                q.append((out_key, npath))
        return best

    def _pick(self, cands: list[Adapter], path: list[Adapter],
              claim: Claim | None, alpha: float) -> Adapter | None:
        if len(cands) == 1 and claim is None:
            return cands[0]
        probe = self._probe(claim, path) if claim is not None else None
        rows = self.rank(cands, probe, alpha)
        for score, a in rows:
            if score[0] < _STATUS_RANK[REFUSE]:
                return a
        return None

    @staticmethod
    def _probe(claim: Claim, path: list[Adapter]) -> Claim:
        """
        Push a cheap copy of the claim along the path so far, so each hop's
        domain test sees the context it would actually receive.
        """
        import random
        from .core import Estimate
        c = claim.derive(estimate=Estimate.point(claim.estimate.mean, 32))
        rng = random.Random(0)
        for a in path:
            # a probe is best-effort: if it no longer typechecks, or the hop
            # refuses, stop and let the caller rank on what is known rather
            # than raising out of a routing query
            if a.consumes.key != c.quantity.key or a.domain(c).status == REFUSE:
                return c
            c = a.apply(c, rng, noise=False)
        return c

    def describe(self, alpha: float = 0.2) -> str:
        lines = ["registered edges:"]
        seen: dict[tuple[str, str], int] = {}
        for a in self.adapters:
            seen[(a.consumes.key, a.produces.key)] = \
                seen.get((a.consumes.key, a.produces.key), 0) + 1
        for a in self.adapters:
            n = seen[(a.consumes.key, a.produces.key)]
            w = a.calibrated_width(alpha)
            wtxt = "-" if not math.isfinite(w) else f"{w:.2f}"
            comp = f"  [1 of {n}]" if n > 1 else ""
            lines.append(f"  [{a.kind:9s}] {a.name:34s} "
                         f"{a.consumes.name} -> {a.produces.name}"
                         f"   width@{alpha}={wtxt}{comp}")
        return "\n".join(lines)


def _cost(path: list[Adapter]) -> tuple[int, int]:
    return (len(path), sum(1 for a in path if a.kind == BRIDGE))
