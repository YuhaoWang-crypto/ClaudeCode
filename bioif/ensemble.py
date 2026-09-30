"""
bioif.ensemble -- running every route, not just the chosen one.

The registry picks one model per hop, on declared properties (§ registry).
That is the right thing to do and it is not sufficient, because the choice
itself carries uncertainty that none of the chosen model's diagnostics can
see: a conformal interval covers the residuals of the model you picked, and
says nothing about whether picking it was right.

So: run the alternatives too. The spread between the endpoints they produce
is a model-choice uncertainty estimate that costs no experiment -- only
compute -- and it can be compared directly against the width the selected
model advertises. When the spread is the larger of the two, the advertised
interval is understating what is actually unknown, and `RouteComparison`
says so rather than leaving it to be noticed.
"""
from __future__ import annotations

import statistics
from dataclasses import dataclass, field

from .adapter import Adapter
from .chain import ChainResult, run_chain
from .core import Claim, Quantity, REFUSE


@dataclass
class RouteComparison:
    """Every in-domain route from one claim to one endpoint quantity."""
    results: dict[str, ChainResult] = field(default_factory=dict)
    paths: dict[str, list[Adapter]] = field(default_factory=dict)
    selected: str = ""
    refused: dict[str, str] = field(default_factory=dict)

    @property
    def ok(self) -> bool:
        return bool(self.results)

    def endpoints(self) -> dict[str, float]:
        return {k: r.final.estimate.mean for k, r in self.results.items()
                if r.ok}

    def point_spread(self) -> float:
        """max - min over the routes' point answers: the cost of the choice."""
        v = list(self.endpoints().values())
        return max(v) - min(v) if len(v) > 1 else 0.0

    def selected_width(self, level: float = 0.80) -> float:
        r = self.results.get(self.selected)
        if not r or not r.ok:
            return float("inf")
        lo, hi = r.final.estimate.ci(level)
        return hi - lo

    def envelope(self, level: float = 0.80) -> tuple[float, float]:
        """Union of every route's interval -- the honest span across models."""
        los, his = [], []
        for r in self.results.values():
            if r.ok:
                lo, hi = r.final.estimate.ci(level)
                los.append(lo)
                his.append(hi)
        return (min(los), max(his)) if los else (float("nan"), float("nan"))

    def choice_exceeds_interval(self, level: float = 0.80) -> bool:
        """Is switching model a bigger move than the chosen model's own interval?"""
        return self.point_spread() > self.selected_width(level)

    def flag(self, level: float = 0.80) -> str | None:
        if not self.choice_exceeds_interval(level):
            return None
        return (f"[model-choice] {len(self.results)} in-domain routes disagree "
                f"by {self.point_spread():.3f}, more than the selected route's "
                f"{level:.0%} interval of {self.selected_width(level):.3f}; the "
                f"calibrated interval does not cover the choice of model")


def compare_routes(reg, src: Quantity, dst: Quantity, claim: Claim,
                   alpha: float = 0.2, seed: int = 7) -> RouteComparison:
    """
    Enumerate the competing first hops, route each to the endpoint, run all.

    Only the first hop is varied because that is where the competition is in
    this repo; the rest of each path is routed normally with the claim in
    hand, so a route that later hits a refusal drops out honestly.
    """
    out = RouteComparison()
    cands = [a for a in reg._by_input.get(src.key, [])]
    ranked = reg.rank(cands, claim, alpha)
    for score, a in ranked:
        label = a.name
        if score[0] >= 2:                     # REFUSE
            out.refused[label] = a.domain(claim).reason
            continue
        # route the remainder with the claim as it will look AFTER this hop,
        # so downstream domain tests see the context the transfer produces
        import random as _r
        after = a.apply(claim, _r.Random(0), noise=False)
        rest = reg.route(a.produces, dst, claim=after, alpha=alpha)
        if rest is None:
            out.refused[label] = f"no route from {a.produces.name} to {dst.name}"
            continue
        path = [a] + rest
        res = run_chain(path, claim, seed=seed)
        if not res.ok:
            out.refused[label] = res.refusal_reason
            continue
        out.paths[label] = path
        out.results[label] = res
        if not out.selected:
            out.selected = label               # ranked, so the first is the pick
    return out
