#!/usr/bin/env python
"""Exploratory sweep over the modelling choices the report recommends.

The primary analysis in ``run_real_experiment.py`` is pre-specified and reported
once.  This script is the opposite: it tries many alternatives, and is therefore
**exploratory and hypothesis-generating, not confirmatory**.

Which means every variant here has to be reported, not just the best one --
otherwise this script would be committing precisely the model-multiplicity error
the project measures elsewhere.  Three defences are built in:

1. **Every variant is listed** in the output, in the order it was run.
2. **Each gets its own permutation null**, so "AUC 0.65" is read against what
   that same pipeline scores on shuffled labels rather than against 0.5.
3. **A multiplicity-aware summary**: with K variants tried, the expected maximum
   under the null is reported alongside the observed maximum. If the best
   variant does not clear that, the honest conclusion is that nothing was found.

The variants come from the assessment report's own recommendations: keep the
label continuous rather than dichotomised, model arrangement separately from
amount, harmonise across scanners, and check whether any single feature family
carries what signal exists.
"""

from __future__ import annotations

import argparse
import json
import sys
import warnings
from dataclasses import dataclass, field
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
warnings.filterwarnings("ignore")

import numpy as np  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402
from sklearn.base import clone  # noqa: E402
from sklearn.linear_model import Ridge  # noqa: E402
from sklearn.metrics import roc_auc_score  # noqa: E402
from sklearn.model_selection import KFold, StratifiedKFold  # noqa: E402
from sklearn.pipeline import Pipeline  # noqa: E402
from sklearn.preprocessing import StandardScaler  # noqa: E402

from mri_cd8_til.experiments import real_cohort  # noqa: E402
from mri_cd8_til.modeling.combat import ComBat, measure_batch_effect  # noqa: E402
from mri_cd8_til.modeling.pipeline import (  # noqa: E402
    UnivariateAUCSelector,
    make_lasso_logistic,
)

SEED = 20260810


@dataclass
class Variant:
    name: str
    rationale: str
    metric: str  # "auc" or "spearman"
    observed: float
    null_mean: float
    null_p95: float
    p_value: float
    n_effective: int
    detail: str = ""
    null: np.ndarray = field(repr=False, default_factory=lambda: np.array([]))

    @property
    def clears_null(self) -> bool:
        return self.observed > self.null_p95

    def describe(self) -> str:
        flag = "clears null" if self.clears_null else "within null"
        return (
            f"{self.name:<34} {self.metric:>8} {self.observed:+.3f}  "
            f"null {self.null_mean:+.3f} (p95 {self.null_p95:+.3f})  "
            f"p={self.p_value:.3f}  [{flag}]"
        )


def classifier(seed: int, k: int | None = 10) -> Pipeline:
    steps = [("scale", StandardScaler())]
    if k:
        steps.append(("select", UnivariateAUCSelector(k=k)))
    steps.append(("clf", make_lasso_logistic(C=1.0, random_state=seed)))
    return Pipeline(steps)


def regressor(seed: int) -> Pipeline:
    return Pipeline([("scale", StandardScaler()), ("reg", Ridge(alpha=10.0, random_state=seed))])


def cv_auc(model, X, y, n_repeats: int, seed: int) -> float:
    scores = []
    for repeat in range(n_repeats):
        cv = StratifiedKFold(5, shuffle=True, random_state=seed + repeat)
        for train, test in cv.split(X, y):
            if len(np.unique(y[test])) < 2:
                continue
            fitted = clone(model).fit(X[train], y[train])
            scores.append(roc_auc_score(y[test], fitted.decision_function(X[test])))
    return float(np.mean(scores)) if scores else float("nan")


def cv_spearman(model, X, values, n_repeats: int, seed: int) -> float:
    """Out-of-fold Spearman between predicted and observed continuous label.

    Pooling out-of-fold predictions before correlating (rather than averaging
    per-fold correlations) keeps the statistic interpretable at n=91, where a
    single fold holds ~18 patients and its correlation is close to noise.
    """
    scores = []
    for repeat in range(n_repeats):
        cv = KFold(5, shuffle=True, random_state=seed + repeat)
        predicted = np.empty_like(values, dtype=float)
        for train, test in cv.split(X):
            fitted = clone(model).fit(X[train], values[train])
            predicted[test] = fitted.predict(X[test])
        scores.append(spearmanr(predicted, values).statistic)
    return float(np.mean(scores))


def with_null(
    name: str,
    rationale: str,
    metric: str,
    compute,
    labels: np.ndarray,
    n_perm: int,
    seed: int,
    detail: str = "",
) -> Variant:
    observed = compute(labels)
    rng = np.random.default_rng(seed)
    null = np.array([compute(rng.permutation(labels)) for _ in range(n_perm)])
    null = null[np.isfinite(null)]
    p = float((np.sum(null >= observed) + 1) / (null.size + 1))
    return Variant(
        name=name,
        rationale=rationale,
        metric=metric,
        observed=observed,
        null_mean=float(null.mean()),
        null_p95=float(np.quantile(null, 0.95)),
        p_value=p,
        n_effective=int(labels.size),
        detail=detail,
        null=null,
    )


def expected_max_under_null(variants: list[Variant], n_draws: int = 20_000, seed: int = 0) -> float:
    """95th percentile of the maximum across K independent nulls.

    The bar the best variant has to clear once K variants have been tried. Built
    by resampling each variant's own null, so it reflects the actual spread of
    each pipeline rather than assuming a common distribution.
    """
    rng = np.random.default_rng(seed)
    usable = [v.null for v in variants if v.null.size > 5]
    if not usable:
        return float("nan")
    draws = np.array([[rng.choice(null) for null in usable] for _ in range(n_draws)])
    return float(np.quantile(draws.max(axis=1), 0.95))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--features", default="results/tcga_radiogenomics_features.xls")
    parser.add_argument("--labels", default="results/til_labels.json")
    parser.add_argument("--out", default="results/model_optimization.json")
    args = parser.parse_args()

    repeats = 2 if args.quick else 8
    n_perm = 40 if args.quick else 200

    print("=" * 80)
    print("EXPLORATORY OPTIMISATION SWEEP — every variant reported, each against its own null")
    print("=" * 80)

    cohort = real_cohort.build(args.features, args.labels)
    scanner_map = real_cohort.fetch_scanner_map(set(cohort.patient_ids))
    cohort = real_cohort.build(args.features, args.labels, scanner_map=scanner_map)
    X = cohort.X
    print(f"\nn={cohort.n}, {len(cohort.feature_names)} features, {len(set(cohort.scanner))} scanners")

    batch_before = measure_batch_effect(X, cohort.scanner)
    harmonised = ComBat().fit_transform(X, batch=cohort.scanner)
    batch_after = measure_batch_effect(harmonised, cohort.scanner)
    print(
        f"variance explained by scanner (median eta^2): "
        f"{batch_before.median_eta_squared:.3f} -> {batch_after.median_eta_squared:.3f} after ComBat"
    )

    log_til = np.log10(cohort.til_fraction + 1e-4)
    moran = cohort.til_moran
    variants: list[Variant] = []

    def add(name, rationale, metric, compute, labels, detail=""):
        v = with_null(name, rationale, metric, compute, labels, n_perm, SEED, detail)
        variants.append(v)
        print("  " + v.describe())

    print("\n[A] the label definition")
    add(
        "tertile split (primary)",
        "the pre-specified analysis: TIL-high vs the rest",
        "auc",
        lambda lab: cv_auc(classifier(SEED), X, lab, repeats, SEED),
        (cohort.til_fraction >= np.quantile(cohort.til_fraction, 2 / 3)).astype(int),
    )
    add(
        "median split",
        "balanced classes carry more power than a 1:2 tertile at this n",
        "auc",
        lambda lab: cv_auc(classifier(SEED), X, lab, repeats, SEED),
        (cohort.til_fraction >= np.median(cohort.til_fraction)).astype(int),
    )
    add(
        "continuous log TIL (regression)",
        "the report's recommendation: do not dichotomise away information",
        "spearman",
        lambda lab: cv_spearman(regressor(SEED), X, lab, repeats, SEED),
        log_til,
    )
    add(
        "extreme tertiles only",
        "drops the ambiguous middle third, the easiest version of the task",
        "auc",
        lambda lab: cv_auc(classifier(SEED), X[_extremes(cohort.til_fraction)], lab, repeats, SEED),
        (
            cohort.til_fraction[_extremes(cohort.til_fraction)]
            >= np.median(cohort.til_fraction)
        ).astype(int),
        detail=f"n={_extremes(cohort.til_fraction).sum()}",
    )

    print("\n[B] arrangement instead of amount")
    add(
        "Moran's I (regression)",
        "MRI texture may track how clustered the infiltrate is, not how much",
        "spearman",
        lambda lab: cv_spearman(regressor(SEED), X, lab, repeats, SEED),
        moran,
    )
    add(
        "Moran's I high vs low",
        "the same, dichotomised at the median",
        "auc",
        lambda lab: cv_auc(classifier(SEED), X, lab, repeats, SEED),
        (moran >= np.median(moran)).astype(int),
    )

    print("\n[C] harmonisation and feature families")
    add(
        "ComBat by scanner",
        "removes the scanner effect measured above",
        "auc",
        lambda lab: cv_auc(classifier(SEED), harmonised, lab, repeats, SEED),
        (cohort.til_fraction >= np.quantile(cohort.til_fraction, 2 / 3)).astype(int),
    )
    for family in ("kinetic", "texture", "morphology"):
        idx = cohort.family_indices(family)
        if idx.size < 2:
            continue
        add(
            f"{family} features only",
            "restricting to one family reduces the search space",
            "auc",
            lambda lab, idx=idx: cv_auc(classifier(SEED, k=None), X[:, idx], lab, repeats, SEED),
            (cohort.til_fraction >= np.quantile(cohort.til_fraction, 2 / 3)).astype(int),
            detail=f"{idx.size} features",
        )

    print("\n" + "=" * 80)
    print("MULTIPLICITY")
    print("=" * 80)
    auc_variants = [v for v in variants if v.metric == "auc"]
    bar = expected_max_under_null(auc_variants, seed=SEED)
    best = max(auc_variants, key=lambda v: v.observed)
    print(f"  {len(auc_variants)} AUC variants tried.")
    print(f"  95th percentile of the MAXIMUM under the null across them: {bar:.3f}")
    print(f"  best observed: {best.observed:.3f}  ({best.name})")
    verdict = (
        "the best variant clears the multiplicity-adjusted bar"
        if best.observed > bar
        else "NOTHING CLEARS IT — the best result is what trying this many variants buys on noise"
    )
    print(f"  -> {verdict}")

    print("\n  Every variant, as run:")
    for v in variants:
        print("    " + v.describe() + (f"   ({v.detail})" if v.detail else ""))

    payload = {
        "n_patients": cohort.n,
        "n_features": len(cohort.feature_names),
        "batch_eta2_before": batch_before.median_eta_squared,
        "batch_eta2_after": batch_after.median_eta_squared,
        "variants": [
            {
                "name": v.name,
                "rationale": v.rationale,
                "metric": v.metric,
                "observed": v.observed,
                "null_mean": v.null_mean,
                "null_p95": v.null_p95,
                "p_value": v.p_value,
                "clears_null": v.clears_null,
                "detail": v.detail,
            }
            for v in variants
        ],
        "multiplicity": {
            "n_auc_variants": len(auc_variants),
            "null_max_p95": bar,
            "best_observed": best.observed,
            "best_name": best.name,
            "clears_adjusted_bar": bool(best.observed > bar),
        },
    }
    path = Path(args.out)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2))
    print(f"\nwrote {path}")
    return 0


def _extremes(values: np.ndarray) -> np.ndarray:
    low, high = np.quantile(values, [1 / 3, 2 / 3])
    return (values <= low) | (values >= high)


if __name__ == "__main__":
    raise SystemExit(main())
