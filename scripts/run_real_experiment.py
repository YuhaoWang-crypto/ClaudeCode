#!/usr/bin/env python
"""Train and test on real data: MRI radiomics -> pathology TIL label, n=91.

Every earlier AUC in this project came from simulation. This one does not.
91 TCGA-BRCA patients with radiologist-supervised lesion segmentations, 36
published radiomic features, and TIL labels computed from their own diagnostic
slides.

The design is deliberately unflattering to itself:

* a **locked test set** split once, before anything is fitted, and touched once;
* every fitted thing -- standardisation, the label cut-point, feature selection,
  the model -- refitted inside each training fold;
* a **permutation null**, because with p=36 and n=64 training patients a
  respectable-looking AUC can appear from nothing;
* **leave-one-scanner-out**, the cross-centre estimate this cohort can support;
* the **published leaky protocol run on the same real data**, so the gap
  measured in simulation can be checked against reality.

Usage::

    python scripts/run_real_experiment.py
    python scripts/run_real_experiment.py --quick
"""

from __future__ import annotations

import argparse
import json
import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
warnings.filterwarnings("ignore")

import numpy as np  # noqa: E402
from sklearn.base import clone  # noqa: E402
from sklearn.ensemble import RandomForestClassifier  # noqa: E402
from sklearn.linear_model import LogisticRegression  # noqa: E402
from sklearn.metrics import roc_auc_score  # noqa: E402
from sklearn.model_selection import StratifiedKFold, train_test_split  # noqa: E402
from sklearn.pipeline import Pipeline  # noqa: E402
from sklearn.preprocessing import StandardScaler  # noqa: E402

from mri_cd8_til.experiments import real_cohort  # noqa: E402
from mri_cd8_til.modeling.metrics import bootstrap_auc_ci, calibration  # noqa: E402
from mri_cd8_til.modeling.pipeline import (  # noqa: E402
    UnivariateAUCSelector,
    columnwise_auc,
    make_lasso_logistic,
)

FEATURES = "results/tcga_radiogenomics_features.xls"
TIL_LABELS = "results/til_labels.json"
LABEL_QUANTILE = 2 / 3
TEST_FRACTION = 0.30
SEED = 20260810


def candidate_models(seed: int) -> dict[str, Pipeline]:
    """Five model families, spanning the bias/variance range sensibly at n=64.

    Nothing exotic: with 36 features and 64 training patients, a deep model
    would be fitting noise, and saying so is more useful than demonstrating it.
    """
    return {
        "lasso": Pipeline(
            [("scale", StandardScaler()), ("clf", make_lasso_logistic(C=1.0, random_state=seed))]
        ),
        "lasso_strong": Pipeline(
            [("scale", StandardScaler()), ("clf", make_lasso_logistic(C=0.1, random_state=seed))]
        ),
        "ridge": Pipeline(
            [
                ("scale", StandardScaler()),
                ("clf", LogisticRegression(C=0.5, max_iter=10_000, random_state=seed)),
            ]
        ),
        "lasso_top10": Pipeline(
            [
                ("scale", StandardScaler()),
                ("select", UnivariateAUCSelector(k=10)),
                ("clf", make_lasso_logistic(C=1.0, random_state=seed)),
            ]
        ),
        "random_forest": Pipeline(
            [
                ("scale", StandardScaler()),
                (
                    "clf",
                    RandomForestClassifier(
                        n_estimators=400, min_samples_leaf=4, max_features="sqrt",
                        random_state=seed, n_jobs=1,
                    ),
                ),
            ]
        ),
    }


def _scores(model, X: np.ndarray) -> np.ndarray:
    if hasattr(model, "decision_function"):
        return model.decision_function(X)
    return model.predict_proba(X)[:, 1]


def cv_auc(model, X, y, n_splits=5, n_repeats=20, seed=0) -> np.ndarray:
    """Repeated stratified CV AUC. Returns per-fold values, not a point."""
    out = []
    for repeat in range(n_repeats):
        cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed + repeat)
        for train, test in cv.split(X, y):
            if len(np.unique(y[test])) < 2:
                continue
            fitted = clone(model).fit(X[train], y[train])
            out.append(roc_auc_score(y[test], _scores(fitted, X[test])))
    return np.array(out)


def fold_safe_label(values: np.ndarray, train_idx: np.ndarray, quantile: float):
    """Cut-point fitted on training values only, applied to everyone.

    Fitting the tertile cut on the full cohort would let the test patients'
    labels influence which class the training patients land in -- a small leak,
    but exactly the kind this project exists to avoid.
    """
    cut = float(np.quantile(values[train_idx], quantile))
    return (values >= cut).astype(int), cut


def leaky_protocol(X, y, train_idx, test_idx, seed) -> float:
    """The reference study's protocol, run on this real data.

    Rank features by AUC over *everyone*, keep the top 10, fit on training,
    report on test.
    """
    scaler = StandardScaler().fit(X[train_idx])
    Xs = scaler.transform(X)
    selector = UnivariateAUCSelector(k=10).fit(Xs, y)  # <- sees the test labels
    Xsel = selector.transform(Xs)
    model = make_lasso_logistic(C=1.0, random_state=seed).fit(Xsel[train_idx], y[train_idx])
    return float(roc_auc_score(y[test_idx], model.decision_function(Xsel[test_idx])))


def leave_one_scanner_out(model, X, y, scanner, min_test=8) -> tuple[np.ndarray, list[str]]:
    scores, labels = [], []
    for held in np.unique(scanner):
        test = scanner == held
        train = ~test
        if test.sum() < min_test or len(np.unique(y[test])) < 2 or len(np.unique(y[train])) < 2:
            continue
        fitted = clone(model).fit(X[train], y[train])
        scores.append(roc_auc_score(y[test], _scores(fitted, X[test])))
        labels.append(f"{held} (n={int(test.sum())})")
    return np.array(scores), labels


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--features", default=FEATURES)
    parser.add_argument("--labels", default=TIL_LABELS)
    parser.add_argument("--out", default="results/real_experiment.json")
    args = parser.parse_args()

    repeats = 4 if args.quick else 20
    n_perm = 60 if args.quick else 400

    print("=" * 78)
    print("REAL TRAIN/TEST — TCGA-BRCA MRI radiomics -> pathology TIL label")
    print("=" * 78)

    print("\n[1] cohort")
    cohort = real_cohort.build(args.features, args.labels)
    scanner_map = real_cohort.fetch_scanner_map(set(cohort.patient_ids))
    cohort = real_cohort.build(args.features, args.labels, scanner_map=scanner_map)

    print(f"  n = {cohort.n} patients, {len(cohort.feature_names)} features")
    families: dict[str, list[str]] = {}
    for name in cohort.feature_names:
        families.setdefault(cohort.family_of(name), []).append(name)
    for family, names in sorted(families.items()):
        print(f"    {family:<22} {len(names)}")
    scanners, counts = np.unique(cohort.scanner, return_counts=True)
    print(f"  scanners: {len(scanners)}")
    for s, c in sorted(zip(scanners, counts), key=lambda x: -x[1])[:6]:
        print(f"    {c:3d}  {s}")

    y_all = cohort.binary_label(LABEL_QUANTILE)
    print(
        f"\n  label: TIL-high vs TIL-low at the upper tertile — "
        f"{int(y_all.sum())} high / {int((1 - y_all).sum())} low"
    )
    print(
        f"  TIL fraction: median {np.median(cohort.til_fraction):.4f}, "
        f"range [{cohort.til_fraction.min():.4f}, {cohort.til_fraction.max():.4f}]"
    )

    print("\n[2] locked test split (made once, before anything is fitted)")
    train_idx, test_idx = train_test_split(
        np.arange(cohort.n), test_size=TEST_FRACTION, stratify=y_all, random_state=SEED
    )
    y, cut = fold_safe_label(cohort.til_fraction, train_idx, LABEL_QUANTILE)
    print(f"  train {len(train_idx)} / test {len(test_idx)}; cut fitted on train only = {cut:.4f}")
    print(f"  train {int(y[train_idx].sum())} high; test {int(y[test_idx].sum())} high")

    X = cohort.X
    Xtr, ytr = X[train_idx], y[train_idx]

    print(f"\n[3] model selection inside the training set only ({repeats}x5-fold CV)")
    results = {}
    for name, model in candidate_models(SEED).items():
        folds = cv_auc(model, Xtr, ytr, n_repeats=repeats, seed=SEED)
        results[name] = folds
        print(
            f"    {name:<16} CV AUC {folds.mean():.3f} +/- {folds.std(ddof=1):.3f}  "
            f"[{np.quantile(folds, .05):.3f}, {np.quantile(folds, .95):.3f}]"
        )
    best_name = max(results, key=lambda k: results[k].mean())
    best = candidate_models(SEED)[best_name]
    print(f"  selected: {best_name} (chosen on training-set CV, never on the test set)")

    print("\n[4] permutation null — the same pipeline on shuffled labels")
    rng = np.random.default_rng(SEED)
    null = np.array(
        [cv_auc(best, Xtr, rng.permutation(ytr), n_repeats=1, seed=SEED + i).mean()
         for i in range(n_perm)]
    )
    observed = float(results[best_name].mean())
    p_value = float((np.sum(null >= observed) + 1) / (null.size + 1))
    print(
        f"    null AUC {null.mean():.3f} +/- {null.std(ddof=1):.3f} "
        f"(95th pct {np.quantile(null, .95):.3f})"
    )
    print(f"    observed {observed:.3f}  ->  permutation p = {p_value:.4f}")

    print("\n[5] the locked test set — touched once")
    fitted = clone(best).fit(Xtr, ytr)
    test_scores = _scores(fitted, X[test_idx])
    interval = bootstrap_auc_ci(y[test_idx], test_scores, n_boot=4000, seed=SEED)
    print(f"    {interval.describe()}")
    probabilities = (
        fitted.predict_proba(X[test_idx])[:, 1]
        if hasattr(fitted, "predict_proba")
        else 1 / (1 + np.exp(-test_scores))
    )
    calib = calibration(y[test_idx], probabilities)
    print(f"    {calib.describe()}")

    print("\n[6] the published leaky protocol, on this same real data")
    leaky = leaky_protocol(X, y, train_idx, test_idx, SEED)
    print(f"    select-on-everyone then split: test AUC {leaky:.3f}")
    print(f"    leakage-free result:           test AUC {interval.auc:.3f}")
    print(f"    difference: {leaky - interval.auc:+.3f}")

    print("\n[7] leave-one-scanner-out — the cross-scanner estimate")
    loso, loso_labels = leave_one_scanner_out(best, X, y, cohort.scanner)
    if loso.size:
        spread = loso.std(ddof=1) if loso.size > 1 else 0.0
        print(f"    mean AUC {loso.mean():.3f} +/- {spread:.3f} over {loso.size} scanner(s)")
        for label, score in zip(loso_labels, loso):
            print(f"      {label:<45} {score:.3f}")
    else:
        print("    not estimable: no scanner has enough patients of both classes")

    print("\n[8] which feature family carries the signal")
    family_scores = {}
    for family in sorted(families):
        idx = cohort.family_indices(family)
        if idx.size < 2:
            continue
        folds = cv_auc(best, Xtr[:, idx], ytr, n_repeats=repeats, seed=SEED)
        family_scores[family] = float(folds.mean())
        print(f"    {family:<22} CV AUC {folds.mean():.3f}  ({idx.size} features)")

    univariate = columnwise_auc(Xtr, ytr)
    order = np.argsort(np.abs(univariate - 0.5))[::-1][:6]
    print("\n    strongest single features (training set only):")
    for i in order:
        print(f"      {univariate[i]:.3f}  {cohort.feature_names[i]}")

    payload = {
        "n_patients": cohort.n,
        "n_features": len(cohort.feature_names),
        "feature_names": list(cohort.feature_names),
        "label": {
            "definition": "TIL-positive share of tissue patches, upper tertile",
            "cut_fitted_on_train": cut,
            "n_high": int(y.sum()),
            "n_low": int((1 - y).sum()),
        },
        "split": {"n_train": len(train_idx), "n_test": len(test_idx)},
        "cv_by_model": {
            k: {"mean": float(v.mean()), "std": float(v.std(ddof=1))} for k, v in results.items()
        },
        "selected_model": best_name,
        "permutation_null": {
            "mean": float(null.mean()),
            "p95": float(np.quantile(null, 0.95)),
            "observed": observed,
            "p_value": p_value,
        },
        "locked_test": {
            "auc": interval.auc,
            "ci_lower": interval.lower,
            "ci_upper": interval.upper,
            "n_positive": interval.n_positive,
            "n_negative": interval.n_negative,
            "brier": calib.brier,
            "calibration_slope": calib.slope,
            "calibration_intercept": calib.intercept,
        },
        "leaky_protocol_test_auc": leaky,
        "leave_one_scanner_out": {
            "mean": float(loso.mean()) if loso.size else None,
            "folds": dict(zip(loso_labels, loso.tolist())) if loso.size else {},
        },
        "family_cv_auc": family_scores,
        "top_univariate": {cohort.feature_names[i]: float(univariate[i]) for i in order},
    }
    path = Path(args.out)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2))
    print(f"\nwrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
