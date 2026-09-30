"""
bioif.real.transfer_models -- competing models for one hop.

Until now every hop in this package had exactly one implementation, so
"which model" was never a question the interface had to answer. It is the
question that matters most in practice: several groups have a model for the
same step, they disagree, and a pipeline has to pick one.

Three predictors are offered for the assay-transfer hop, all fitted on the
same real ChEMBL pairs, all wrapped in the same split-conformal procedure so
their intervals are comparable:

  identity  y = x              -- reuse the number as measured
  shift     y = x + b          -- a single median offset between the assays
  linear    y = a + b*x        -- ordinary least squares

The point is not that one wins. All three are conformally *valid* by
construction; what differs is efficiency (interval width), and that is what
the router can legitimately choose on. What the router cannot choose away is
the disagreement between them -- see `disagreement_vs_width`, which measures
how far apart the models' point predictions sit relative to the interval each
one advertises. That gap is model-choice uncertainty, and conformal does not
cover it: it covers the residuals of the model you already picked.
"""
from __future__ import annotations

import math
import random
import statistics
from dataclasses import dataclass, field

from . import conformal as C


# --------------------------------------------------------------------------
# Predictors
# --------------------------------------------------------------------------

class Predictor:
    name = "predictor"

    def fit(self, train: list[C.Point]) -> "Predictor":
        return self

    def predict(self, x: float) -> float:
        raise NotImplementedError

    def describe(self) -> str:
        return self.name


class IdentityPredictor(Predictor):
    """Reuse the number as measured. The null model this repo keeps warning about."""
    name = "identity"

    def predict(self, x: float) -> float:
        return x


class ShiftPredictor(Predictor):
    """A single median offset between the two assays."""
    name = "shift"

    def __init__(self):
        self.b = 0.0

    def fit(self, train):
        self.b = statistics.median([p.y - p.x for p in train]) if train else 0.0
        return self

    def predict(self, x):
        return x + self.b

    def describe(self):
        return f"shift  (y = x {self.b:+.3f})"


class LinearPredictor(Predictor):
    """Ordinary least squares, which can also rescale, not just offset."""
    name = "linear"

    def __init__(self):
        self.a, self.b = 0.0, 1.0

    def fit(self, train):
        n = len(train)
        if n < 3:
            self.a, self.b = 0.0, 1.0
            return self
        mx = statistics.fmean(p.x for p in train)
        my = statistics.fmean(p.y for p in train)
        sxx = sum((p.x - mx) ** 2 for p in train)
        sxy = sum((p.x - mx) * (p.y - my) for p in train)
        self.b = sxy / sxx if sxx > 1e-12 else 1.0
        self.a = my - self.b * mx
        return self

    def predict(self, x):
        return self.a + self.b * x

    def describe(self):
        return f"linear (y = {self.a:+.3f} {self.b:+.3f}x)"


PREDICTORS = {"identity": IdentityPredictor, "shift": ShiftPredictor,
              "linear": LinearPredictor}


# --------------------------------------------------------------------------
# One conformalised model
# --------------------------------------------------------------------------

def signed_quantiles(residuals: list[float], alpha: float):
    """
    Two-sided split-conformal offsets from SIGNED calibration residuals.

    `conformal.py` uses the symmetric |residual| variant, which is the
    textbook form and what §6 of INTEROP.md reports. This module needs the
    two-sided signed form instead, for a reason that only showed up once
    three models were competing:

    an adapter that emits `predict(x) + resampled residual` has already
    absorbed any bias in `predict`, because the residuals carry that bias
    with the opposite sign. Its emitted distribution is therefore NOT
    centred on `predict(x)`, and a symmetric interval around `predict(x)`
    does not describe what it emits. Scoring such an adapter on the
    symmetric width punishes it for a bias its own output no longer has.

    The signed form fixes that: the interval is exactly the central
    (1-alpha) region of what the adapter emits. A consequence worth stating
    -- see demo_routing R2 -- is that a biased predictor and its
    bias-corrected twin become the SAME conformalised model.
    """
    n = len(residuals)
    if n == 0:
        return -math.inf, math.inf
    s = sorted(residuals)
    k_lo = math.floor((n + 1) * (alpha / 2.0))
    k_hi = math.ceil((n + 1) * (1.0 - alpha / 2.0))
    if k_lo < 1 or k_hi > n:
        return -math.inf, math.inf
    return s[k_lo - 1], s[k_hi - 1]


@dataclass
class ConformalModel:
    predictor: Predictor
    offsets: dict[float, tuple[float, float]]     # alpha -> (lo, hi), signed
    residuals: list[float]
    x_range: tuple[float, float]
    n_train: int
    n_cal: int
    pair: str = ""

    def predict(self, x: float) -> float:
        """The predictor's raw output, before conformalisation."""
        return self.predictor.predict(x)

    def center(self, x: float) -> float:
        """
        The centre of what this model actually emits.

        Any bias in `predict` is carried in the residuals, so the emitted
        distribution sits at predict(x) + median(residual). This, not
        `predict`, is what two models should be compared on.
        """
        med = statistics.median(self.residuals) if self.residuals else 0.0
        return self.predict(x) + med

    def interval(self, x: float, alpha: float) -> tuple[float, float]:
        lo, hi = self.offsets[alpha]
        c = self.predict(x)
        return c + lo, c + hi

    def width(self, alpha: float) -> float:
        lo, hi = self.offsets.get(alpha, (-math.inf, math.inf))
        return hi - lo

    def attainable(self, alpha: float) -> bool:
        return math.isfinite(self.width(alpha))


def fit_conformal(pts: list[C.Point], pair: str, predictor: str,
                  alphas=(0.2, 0.1), seed: int = 11,
                  fracs=(0.35, 0.65, 0.0)) -> ConformalModel:
    """Split conformal: fit on one half of the compounds, calibrate on the other."""
    sub = [p for p in pts if p.pair == pair]
    tr, ca, _ = C.split_by_compound(sub, random.Random(seed), fracs=fracs)
    pred = PREDICTORS[predictor]().fit(tr)
    res = [p.y - pred.predict(p.x) for p in ca]
    xs = [p.x for p in tr + ca] or [0.0]
    return ConformalModel(
        pred, {a: signed_quantiles(res, a) for a in alphas}, res,
        (min(xs), max(xs)), len(tr), len(ca), pair)


def coverage(model: ConformalModel, test: list[C.Point], alpha: float):
    if not test or not model.attainable(alpha):
        return float("nan"), float("inf")
    hit = 0
    for p in test:
        lo, hi = model.interval(p.x, alpha)
        hit += (lo <= p.y <= hi)
    return hit / len(test), model.width(alpha)


# --------------------------------------------------------------------------
# Comparing the competitors
# --------------------------------------------------------------------------

def compare_predictors(pts: list[C.Point], min_n: int = 16, alpha: float = 0.2,
                       repeats: int = 300, seed: int = 0,
                       fracs=(0.25, 0.60, 0.15)) -> dict:
    """
    Per assay pair, the coverage and width of each competing model.

    Same splits for all three, so the comparison is paired: any difference is
    the model, not the split.
    """
    names = [n for n in sorted({p.pair for p in pts})
             if len([p for p in pts if p.pair == n]) >= min_n]
    out: dict = {}
    for name in names:
        sub = [p for p in pts if p.pair == name]
        rec = {k: {"cov": [], "w": []} for k in PREDICTORS}
        for r in range(repeats):
            rng = random.Random(seed + r)
            tr, ca, te = C.split_by_compound(sub, rng, fracs=fracs)
            if not (tr and ca and te):
                continue
            for k, cls in PREDICTORS.items():
                pred = cls().fit(tr)
                res = [p.y - pred.predict(p.x) for p in ca]
                off = signed_quantiles(res, alpha)
                if not math.isfinite(off[1] - off[0]):
                    continue
                m = ConformalModel(pred, {alpha: off}, res, (0, 0),
                                   len(tr), len(ca))
                cov, w = coverage(m, te, alpha)
                rec[k]["cov"].append(cov)
                rec[k]["w"].append(w)
        out[name] = {"n": len(sub)}
        for k in PREDICTORS:
            d = rec[k]
            out[name][k] = {
                "cov": statistics.fmean(d["cov"]) if d["cov"] else float("nan"),
                "w": statistics.fmean(d["w"]) if d["w"] else float("inf"),
            }
    return out


def disagreement_vs_width(pts: list[C.Point], min_n: int = 16,
                          alpha: float = 0.2, seed: int = 11,
                          fracs=(0.35, 0.65, 0.0)) -> dict:
    """
    How far apart do the models' point predictions sit, compared with the
    interval each of them advertises?

    For every pair, fit all three models on the same split, then over the
    pair's own x values measure:
      * spread   = max - min over the three models' emitted CENTRES
                   (predict(x) + median residual -- see ConformalModel.center;
                   comparing raw predict(x) instead overstates the gap, because
                   a biased predictor's emitted distribution is not centred
                   there)
      * width    = the conformal interval width of the model the router picks
      * covered  = does the selected model's interval contain the others'
                   point predictions?

    Where `covered` is false, model choice is moving the answer further than
    the selected model's own interval admits. That is uncertainty the
    conformal guarantee does not address, because the guarantee is about the
    residuals of one chosen model, not about having chosen it.
    """
    names = [n for n in sorted({p.pair for p in pts})
             if len([p for p in pts if p.pair == n]) >= min_n]
    out: dict = {}
    for name in names:
        sub = [p for p in pts if p.pair == name]
        models = {k: fit_conformal(pts, name, k, (alpha,), seed, fracs)
                  for k in PREDICTORS}
        usable = {k: m for k, m in models.items() if m.attainable(alpha)}
        if len(usable) < 2:
            continue
        best = min(usable, key=lambda k: usable[k].width(alpha))
        spreads, contained = [], []
        for p in sub:
            preds = [m.center(p.x) for m in usable.values()]
            spreads.append(max(preds) - min(preds))
            lo, hi = usable[best].interval(p.x, alpha)
            contained.append(all(lo <= v <= hi for v in preds))
        out[name] = {
            "n": len(sub),
            "selected": best,
            "width": usable[best].width(alpha),
            "median_spread": statistics.median(spreads),
            "max_spread": max(spreads),
            "frac_contained": statistics.fmean([float(c) for c in contained]),
            "models": {k: m.predictor.describe() for k, m in usable.items()},
        }
    return out
