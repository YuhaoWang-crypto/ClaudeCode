"""
bioif.real.conformal -- a calibrated interval for the assay-transfer hop,
and an honest check of whether it covers.

Everywhere else in bioif, an adapter's uncertainty is asserted: a stub says
sd=0.55 and the chain believes it. Split conformal prediction replaces the
assertion with a finite-sample guarantee -- under exchangeability, an
interval built from held-out residuals covers the truth with probability at
least 1-alpha, with no assumption about the error's shape and none about the
predictor being any good.

The hop chosen is the one this package keeps running into: **can a potency
measured in assay A be reused where assay B is needed?** That is the §1A
axis-3 bottleneck in its most concrete form, it has real labels (compounds
measured in both assays), and a chain must answer it before reusing any
number at all.

Three things came out of running it on real data, in increasing order of
how much they matter:

  1. Conformal hits nominal coverage within a single assay pair, where the
     usual Gaussian +/- z*sd does not.
  2. Roughly half the candidate assay pairs on this target are not
     continuous measurements at all -- they are patent potency BINS at
     one-log spacing wearing a float column. Fitting a regression to them is
     a category error, so `classify_pair` refuses them.
  3. A conformal interval calibrated on one assay pair does NOT cover on a
     different assay pair. The guarantee is conditional on exchangeability,
     and two assay pairs are not exchangeable. This is measured, not argued.

Data: the committed ChEMBL snapshot, IC50 records only, so the readout-type
confound from INTEROP.md §5 defect A cannot leak in.
"""
from __future__ import annotations

import collections
import itertools
import math
import random
import statistics
from dataclasses import dataclass, field

from .chembl import load_snapshot

# --------------------------------------------------------------------------
# The transfer task
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class Point:
    compound: str
    pair: str
    x: float          # pChEMBL in the source assay
    y: float          # pChEMBL in the target assay


@dataclass(frozen=True)
class PairInfo:
    """A candidate transfer, with the data-quality verdict attached."""
    key: str
    source: str
    target: str
    n: int
    tie_x: float          # fraction of source values equal to the modal value
    n_distinct_x: int
    tie_y: float
    n_distinct_y: int
    bias: float
    sd: float
    status: str           # 'continuous' | 'binned'
    reason: str = ""


#: A column of floats is not evidence of a continuous measurement. Patent
#: potency bands ("IC50 < 0.1 uM") arrive as a float at the band midpoint, so
#: a handful of values repeat hundreds of times. Regression and conformal
#: intervals are both meaningless on that, and it is invisible unless checked.
TIE_LIMIT = 0.25
DISTINCT_FRAC = 0.40


def classify(vals_x: list[float], vals_y: list[float]) -> tuple[str, str]:
    n = len(vals_x)
    cx, cy = collections.Counter(round(v, 2) for v in vals_x), \
        collections.Counter(round(v, 2) for v in vals_y)
    tie_x = cx.most_common(1)[0][1] / n
    tie_y = cy.most_common(1)[0][1] / n
    if tie_x > TIE_LIMIT or tie_y > TIE_LIMIT or len(cx) < n * DISTINCT_FRAC:
        return "binned", (
            f"{tie_x:.0%} of source values share one number over "
            f"{len(cx)} distinct levels (n={n}); this looks like binned "
            f"potency bands, not a continuous measurement")
    return "continuous", ""


def discover_pairs(rows=None, min_n: int = 10) -> dict[str, PairInfo]:
    """Every assay pair on the target with >= min_n compounds measured in both."""
    rows = rows if rows is not None else load_snapshot()
    rows = [r for r in rows if r["standard_type"] == "IC50" and r["pchembl_value"]]
    by = collections.defaultdict(lambda: collections.defaultdict(list))
    for r in rows:
        by[r["molecule_chembl_id"]][r["assay_chembl_id"]].append(
            float(r["pchembl_value"]))
    shared = collections.defaultdict(list)
    for cmp_id, av in by.items():
        for a, b in itertools.combinations(sorted(av), 2):
            shared[(a, b)].append((cmp_id, statistics.fmean(av[a]),
                                   statistics.fmean(av[b])))
    out: dict[str, PairInfo] = {}
    for (a, b), v in shared.items():
        if len(v) < min_n:
            continue
        xs, ys = [t[1] for t in v], [t[2] for t in v]
        d = [y - x for x, y in zip(xs, ys)]
        status, reason = classify(xs, ys)
        cx = collections.Counter(round(t, 2) for t in xs)
        cy = collections.Counter(round(t, 2) for t in ys)
        key = f"{a[-6:]}->{b[-6:]}"
        out[key] = PairInfo(key, a, b, len(v),
                            cx.most_common(1)[0][1] / len(v), len(cx),
                            cy.most_common(1)[0][1] / len(v), len(cy),
                            statistics.fmean(d), statistics.pstdev(d),
                            status, reason)
    return dict(sorted(out.items(), key=lambda kv: -kv[1].n))


def build_points(pairs: dict[str, PairInfo] | None = None,
                 rows=None, continuous_only: bool = True) -> list[Point]:
    rows = rows if rows is not None else load_snapshot()
    pairs = pairs if pairs is not None else discover_pairs(rows)
    rows = [r for r in rows if r["standard_type"] == "IC50" and r["pchembl_value"]]
    by = collections.defaultdict(lambda: collections.defaultdict(list))
    for r in rows:
        by[r["molecule_chembl_id"]][r["assay_chembl_id"]].append(
            float(r["pchembl_value"]))
    out: list[Point] = []
    for key, info in pairs.items():
        if continuous_only and info.status != "continuous":
            continue
        for cmp_id, av in by.items():
            if info.source in av and info.target in av:
                out.append(Point(cmp_id, key,
                                 statistics.fmean(av[info.source]),
                                 statistics.fmean(av[info.target])))
    return out


# --------------------------------------------------------------------------
# Predictor: deliberately trivial, because conformal does not need it good
# --------------------------------------------------------------------------

def fit_shift(train: list[Point]) -> float:
    """The whole 'model': the median shift between the two assays."""
    return statistics.median([p.y - p.x for p in train]) if train else 0.0


def predict(x: float, shift: float) -> float:
    return x + shift


# --------------------------------------------------------------------------
# Split conformal
# --------------------------------------------------------------------------

def conformal_quantile(scores: list[float], alpha: float) -> float:
    """
    The finite-sample split-conformal quantile: the ceil((n+1)(1-alpha))-th
    smallest nonconformity score.

    When that index exceeds n, the guarantee simply cannot be met at this
    alpha with this much calibration data, and the correct return value is
    infinity. A method that quietly returned the largest observed score
    instead would be claiming a guarantee it does not have -- which is the
    same failure mode as an adapter that cannot refuse.
    """
    n = len(scores)
    if n == 0:
        return math.inf
    k = math.ceil((n + 1) * (1.0 - alpha))
    return math.inf if k > n else sorted(scores)[k - 1]


@dataclass
class Calibration:
    shift: float
    q: dict[float, float]
    n_train: int
    n_cal: int
    x_range: tuple[float, float]
    pair: str = ""

    def interval(self, x: float, alpha: float) -> tuple[float, float]:
        c = predict(x, self.shift)
        return c - self.q[alpha], c + self.q[alpha]

    def attainable(self, alpha: float) -> bool:
        return math.isfinite(self.q.get(alpha, math.inf))


def calibrate(train: list[Point], cal: list[Point],
              alphas=(0.2, 0.1)) -> Calibration:
    shift = fit_shift(train)
    scores = [abs(p.y - predict(p.x, shift)) for p in cal]
    xs = [p.x for p in train + cal] or [0.0]
    return Calibration(shift, {a: conformal_quantile(scores, a) for a in alphas},
                       len(train), len(cal), (min(xs), max(xs)))


def evaluate(test: list[Point], c: Calibration, alpha: float):
    """(empirical coverage, mean interval width) on the test set."""
    if not test or not c.attainable(alpha):
        return float("nan"), float("inf")
    hit = sum(c.interval(p.x, alpha)[0] <= p.y <= c.interval(p.x, alpha)[1]
              for p in test)
    return hit / len(test), 2.0 * c.q[alpha]


def gaussian_calibration(train: list[Point], alphas=(0.2, 0.1)) -> Calibration:
    """The usual thing people do instead: +/- z * sd of training residuals."""
    z = {0.2: 1.2816, 0.1: 1.6449, 0.05: 1.9600}
    shift = fit_shift(train)
    res = [p.y - predict(p.x, shift) for p in train]
    sd = statistics.pstdev(res) if len(res) > 1 else 0.0
    xs = [p.x for p in train] or [0.0]
    return Calibration(shift, {a: z[a] * sd for a in alphas},
                       len(train), 0, (min(xs), max(xs)))


# --------------------------------------------------------------------------
# Splits -- on COMPOUND, never on row (§1A axis 10)
# --------------------------------------------------------------------------

def split_by_compound(pts: list[Point], rng: random.Random,
                      fracs=(1 / 3, 1 / 3, 1 / 3)):
    cmps = sorted({p.compound for p in pts})
    rng.shuffle(cmps)
    n = len(cmps)
    i, j = int(fracs[0] * n), int((fracs[0] + fracs[1]) * n)
    a, b = set(cmps[:i]), set(cmps[i:j])
    return ([p for p in pts if p.compound in a],
            [p for p in pts if p.compound in b],
            [p for p in pts if p.compound not in a and p.compound not in b])


# --------------------------------------------------------------------------
# Experiments
# --------------------------------------------------------------------------

def experiment_within_pair(pts: list[Point], min_n: int = 12,
                           alphas=(0.2, 0.1), repeats: int = 300,
                           seed: int = 0) -> dict:
    """
    E1 -- calibrate and test inside one assay pair, so exchangeability holds.
    Conformal should hit nominal; the Gaussian baseline need not.
    """
    out: dict = {}
    names = sorted({p.pair for p in pts})
    for name in names:
        sub = [p for p in pts if p.pair == name]
        if len(sub) < min_n:
            continue
        rec = {a: {"cov": [], "w": [], "gcov": [], "gw": [], "unattainable": 0}
               for a in alphas}
        for r in range(repeats):
            rng = random.Random(seed + r)
            tr, ca, te = split_by_compound(sub, rng)
            if not (tr and ca and te):
                continue
            c = calibrate(tr, ca, alphas)
            g = gaussian_calibration(tr, alphas)
            for a in alphas:
                if not c.attainable(a):
                    rec[a]["unattainable"] += 1
                    continue
                cov, w = evaluate(te, c, a)
                rec[a]["cov"].append(cov)
                rec[a]["w"].append(w)
                gcov, gw = evaluate(te, g, a)
                rec[a]["gcov"].append(gcov)
                rec[a]["gw"].append(gw)
        out[name] = {"n": len(sub)}
        for a in alphas:
            d = rec[a]
            out[name][a] = {
                "cov": statistics.fmean(d["cov"]) if d["cov"] else float("nan"),
                "w": statistics.fmean(d["w"]) if d["w"] else float("inf"),
                "gcov": statistics.fmean(d["gcov"]) if d["gcov"] else float("nan"),
                "gw": statistics.fmean(d["gw"]) if d["gw"] else float("inf"),
                "unattainable": d["unattainable"] / max(repeats, 1),
            }
    return out


def experiment_cross_pair(pts: list[Point], alphas=(0.2, 0.1),
                          repeats: int = 100, seed: int = 0) -> dict:
    """
    E2 -- leave-one-assay-pair-out.

    Calibrate on every OTHER assay pair pooled, then test on the held-out
    pair. Conformal's guarantee is conditional on calibration and test data
    being exchangeable; two assay pairs have different systematic shifts, so
    they are not. The size of the coverage shortfall IS the size of the
    applicability-domain problem, and it is why the adapter below refuses an
    uncalibrated pair rather than widening for it.
    """
    out: dict = {}
    names = sorted({p.pair for p in pts})
    for held in names:
        sub = [p for p in pts if p.pair == held]
        rest = [p for p in pts if p.pair != held]
        if len(sub) < 10 or len(rest) < 30:
            continue
        rec = {a: [] for a in alphas}
        for r in range(repeats):
            rng = random.Random(seed + r)
            tr, ca, _ = split_by_compound(rest, rng, fracs=(0.5, 0.5, 0.0))
            c = calibrate(tr, ca, alphas)
            for a in alphas:
                if c.attainable(a):
                    rec[a].append(evaluate(sub, c, a)[0])
        out[held] = {"n": len(sub)}
        for a in alphas:
            out[held][a] = statistics.fmean(rec[a]) if rec[a] else float("nan")
    return out


def experiment_mondrian(pts: list[Point], min_n: int = 12, alpha: float = 0.2,
                        repeats: int = 300, seed: int = 0) -> dict:
    """
    E3 -- pooled (marginal) calibration vs per-pair (Mondrian) calibration,
    both scored per pair.

    Marginal conformal promises coverage averaged over the population. When
    the population is a mixture of assay pairs with different shifts, that
    average can be met while individual pairs are badly under-covered. The
    class-conditional fix uses the assay pair as the class -- i.e. exactly
    the context field the contract already carries.
    """
    names = [n for n in sorted({p.pair for p in pts})
             if len([p for p in pts if p.pair == n]) >= min_n]
    rec = {n: {"pooled": [], "mondrian": []} for n in names}
    marginal = []
    for r in range(repeats):
        rng = random.Random(seed + r)
        sp = {n: split_by_compound([p for p in pts if p.pair == n], rng)
              for n in names}
        ptr = [p for n in names for p in sp[n][0]]
        pca = [p for n in names for p in sp[n][1]]
        pooled = calibrate(ptr, pca, (alpha,))
        if not pooled.attainable(alpha):
            continue
        alltest = []
        for n in names:
            tr, ca, te = sp[n]
            if not (tr and ca and te):
                continue
            alltest += te
            rec[n]["pooled"].append(evaluate(te, pooled, alpha)[0])
            m = calibrate(tr, ca, (alpha,))
            if m.attainable(alpha):
                rec[n]["mondrian"].append(evaluate(te, m, alpha)[0])
        marginal.append(evaluate(alltest, pooled, alpha)[0])
    res = {n: {k: (statistics.fmean(v) if v else float("nan"))
               for k, v in d.items()} for n, d in rec.items()}
    res["_pooled_marginal"] = statistics.fmean(marginal) if marginal else float("nan")
    return res


def fit_deployable(pts: list[Point], pair: str, alphas=(0.2, 0.1),
                   seed: int = 11) -> Calibration:
    """The calibration the adapter ships with: half to fit, half to calibrate."""
    sub = [p for p in pts if p.pair == pair]
    tr, ca, _ = split_by_compound(sub, random.Random(seed), fracs=(0.5, 0.5, 0.0))
    c = calibrate(tr, ca, alphas)
    c.pair = pair
    return c
