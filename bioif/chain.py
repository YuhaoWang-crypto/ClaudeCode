"""
bioif.chain -- execution, uncertainty propagation, and the weakest-link report.

A long causal chain of point predictions is not a long causal prediction; it
is a product of unvalidated maps that returns a number with no error bar. The
executor here does three things a naive `f(g(h(x)))` pipeline does not:

  1. It propagates a full sample distribution, so nonlinear hops (Hill
     functions, thresholds, log/linear changes) widen honestly instead of
     being linearised at the mean.

  2. It lets any link REFUSE. A missing dose, or a cell system outside the
     calibration set, stops the chain with a reason instead of producing a
     confident number nobody can trace.

  3. It attributes the variance of the final answer back to individual links,
     by re-running with each link's own noise switched off. That converts
     "the chain is uncertain" into "link 3 is 61% of the uncertainty", which
     is the only form of the statement you can act on.

(3) is what makes the whole exercise worth doing. The realistic goal for a
long chain is not an accurate endpoint prediction -- it is a defensible
ranking of which single measurement would most tighten the endpoint.
"""
from __future__ import annotations

import random
from dataclasses import dataclass, field

from .adapter import Adapter
from .core import Claim, Estimate, REFUSE, EXTRAPOLATE, IN_DOMAIN, Verdict


@dataclass
class LinkTrace:
    adapter: Adapter
    verdict: Verdict
    out: Claim | None
    var_share: float = 0.0          # fraction of final variance owned by link


@dataclass
class ChainResult:
    source: Claim
    traces: list[LinkTrace] = field(default_factory=list)
    final: Claim | None = None
    refused_at: str | None = None
    refusal_reason: str = ""
    source_var_share: float = 0.0

    @property
    def ok(self) -> bool:
        return self.final is not None

    def report(self) -> str:
        lines = []
        lines.append("SOURCE")
        lines.append("  " + self.source.describe().replace("\n", "\n  "))
        for t in self.traces:
            mark = {IN_DOMAIN: "ok", EXTRAPOLATE: "EXTRAPOLATED",
                    REFUSE: "REFUSED"}[t.verdict.status]
            lines.append(f"\nLINK  {t.adapter.name}  [{t.adapter.kind}]  ({mark})")
            if t.verdict.reason:
                lines.append(f"  ! {t.verdict.reason}")
            if t.out is not None:
                lines.append("  " + t.out.describe().replace("\n", "\n  "))
                lines.append(f"  variance share of final answer: "
                             f"{t.var_share * 100:4.1f}%")
        if self.refused_at:
            lines.append(f"\nCHAIN REFUSED at {self.refused_at}: "
                         f"{self.refusal_reason}")
            return "\n".join(lines)
        lines.append(f"\nsource measurement itself: "
                     f"{self.source_var_share * 100:4.1f}% of final variance")
        lines.append("  (shares need not sum to 100%: the hops are nonlinear, "
                     "so they interact)")
        lines.append("\nFINAL")
        lines.append("  " + self.final.describe().replace("\n", "\n  "))
        if self.final.flags:
            lines.append("  assumptions carried to this number:")
            for f in self.final.flags:
                lines.append(f"    - {f}")
        lines.append("  provenance: " +
                     " -> ".join(str(p) for p in self.final.provenance))
        return "\n".join(lines)


def _run_once(path: list[Adapter], source: Claim, seed: int,
              mute: int | None = None,
              mute_source: bool = False) -> tuple[Claim | None, list[Verdict]]:
    """
    Execute the path. `mute` silences link i's own noise; `mute_source`
    replaces the source distribution by its mean. Both are used for variance
    attribution, not for reporting.
    """
    rng = random.Random(seed)
    claim = source
    if mute_source:
        claim = claim.derive(
            estimate=Estimate.point(source.estimate.mean,
                                    len(source.estimate.samples)))
    verdicts: list[Verdict] = []
    for i, a in enumerate(path):
        v = a.domain(claim)
        verdicts.append(v)
        if v.status == REFUSE:
            return None, verdicts
        noise = (mute != i)
        inflate = v.inflate if (v.status == EXTRAPOLATE and noise) else 1.0
        claim = a.apply(claim, rng, noise=noise, inflate=inflate)
    return claim, verdicts


def run_chain(path: list[Adapter], source: Claim, seed: int = 0) -> ChainResult:
    """Execute a chain and attribute its output variance to its links."""
    res = ChainResult(source=source)

    # Full run, keeping every intermediate for the trace.
    rng = random.Random(seed)
    claim = source
    for a in path:
        v = a.domain(claim)
        if v.status == REFUSE:
            res.traces.append(LinkTrace(a, v, None))
            res.refused_at = a.name
            res.refusal_reason = v.reason
            return res
        inflate = v.inflate if v.status == EXTRAPOLATE else 1.0
        claim = a.apply(claim, rng, noise=True, inflate=inflate)
        res.traces.append(LinkTrace(a, v, claim))
    res.final = claim

    # Attribution: re-run with one source of stochasticity muted at a time.
    # share_i = 1 - Var(final | link i noiseless) / Var(final)
    total = claim.estimate.var
    if total > 0:
        for i, t in enumerate(res.traces):
            muted, _ = _run_once(path, source, seed, mute=i)
            if muted is not None:
                t.var_share = max(0.0, 1.0 - muted.estimate.var / total)
        muted, _ = _run_once(path, source, seed, mute_source=True)
        if muted is not None:
            res.source_var_share = max(0.0, 1.0 - muted.estimate.var / total)
    return res


# --------------------------------------------------------------------------
# Value of information: turn the attribution into an experiment ranking
# --------------------------------------------------------------------------

@dataclass
class Experiment:
    """A measurement that would collapse one link's uncertainty."""
    link: str
    description: str
    cost: float                     # arbitrary consistent units (e.g. k$)


def rank_experiments(res: ChainResult,
                     catalogue: dict[str, Experiment]) -> list[tuple]:
    """
    Rank candidate measurements by expected variance removed per unit cost.

    This is the deliverable of a long chain. Not "the predicted endpoint is
    X", but "of the six things you could measure, this one buys the most
    resolution per unit of effort, and these three buy almost nothing".
    """
    if not res.ok:
        return []
    total = res.final.estimate.var
    rows = []
    for t in res.traces:
        exp = catalogue.get(t.adapter.name)
        if exp is None or t.var_share <= 0:
            continue
        rows.append((t.var_share * total / exp.cost, t.var_share, exp))
    rows.sort(key=lambda r: -r[0])
    return rows
