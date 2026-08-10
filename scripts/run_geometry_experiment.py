#!/usr/bin/env python
"""Can MRI predict the geometry-derived immunophenotype? Both of the paper's tasks.

The previous real experiment could only attempt the reference study's
``RM-whole`` task (inflamed vs rest), because the published TIL maps carry no
tumour boundary and ``RM-peri`` (desert vs excluded) is defined against one.
The geometry descriptors supply a boundary-free stand-in for that distinction,
so **both tasks are testable here for the first time**.

Two pre-specified primary analyses, each reported once:

``RM-whole``  inflamed vs (desert + excluded)
``RM-peri``   desert vs excluded, with inflamed cases dropped rather than
              silently coded as negatives

The imaging side is still the 91 patients with radiologist-supervised
segmentations -- geometry scales the *label* to 136, not the features. Dropping
inflamed cases makes ``RM-peri`` smaller again, and that shrinkage is reported
rather than hidden, because it is the binding constraint on what this can show.

Everything else matches ``run_real_experiment.py``: locked test set, everything
refitted inside folds, permutation nulls, and a multiplicity-adjusted bar over
however many targets end up being tried.
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
from sklearn.metrics import roc_auc_score  # noqa: E402
from sklearn.model_selection import StratifiedKFold, train_test_split  # noqa: E402
from sklearn.pipeline import Pipeline  # noqa: E402
from sklearn.preprocessing import StandardScaler  # noqa: E402

from mri_cd8_til.experiments import real_cohort  # noqa: E402
from mri_cd8_til.labels.schema import Immunophenotype, LabelTask, to_binary  # noqa: E402
from mri_cd8_til.modeling.metrics import bootstrap_auc_ci  # noqa: E402
from mri_cd8_til.modeling.pipeline import UnivariateAUCSelector, make_lasso_logistic  # noqa: E402

SEED = 20260810
TEST_FRACTION = 0.30


def models(seed: int) -> dict[str, Pipeline]:
    return {
        "lasso": Pipeline(
            [("scale", StandardScaler()), ("clf", make_lasso_logistic(C=1.0, random_state=seed))]
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


def _scores(model, X):
    return model.decision_function(X) if hasattr(model, "decision_function") else model.predict_proba(X)[:, 1]


def cv_auc(model, X, y, n_repeats: int, seed: int) -> np.ndarray:
    out = []
    for repeat in range(n_repeats):
        cv = StratifiedKFold(5, shuffle=True, random_state=seed + repeat)
        for train, test in cv.split(X, y):
            if len(np.unique(y[test])) < 2 or len(np.unique(y[train])) < 2:
                continue
            fitted = clone(model).fit(X[train], y[train])
            out.append(roc_auc_score(y[test], _scores(fitted, X[test])))
    return np.array(out)


def run_task(name: str, X: np.ndarray, y: np.ndarray, repeats: int, n_perm: int) -> dict:
    print(f"\n{'=' * 78}\n{name}\n{'=' * 78}")
    print(f"  n = {len(y)}  ({int(y.sum())} positive / {int((1 - y).sum())} negative)")
    if len(y) < 30 or min(int(y.sum()), int((1 - y).sum())) < 8:
        print("  too small to split and test honestly — reporting CV only")
        locked = None
        train_idx = np.arange(len(y))
        test_idx = np.array([], dtype=int)
    else:
        train_idx, test_idx = train_test_split(
            np.arange(len(y)), test_size=TEST_FRACTION, stratify=y, random_state=SEED
        )
        print(f"  locked split: train {len(train_idx)} / test {len(test_idx)}")
        locked = True

    Xtr, ytr = X[train_idx], y[train_idx]

    print("  model selection on the training set only:")
    cv_results = {}
    for model_name, model in models(SEED).items():
        folds = cv_auc(model, Xtr, ytr, repeats, SEED)
        cv_results[model_name] = folds
        print(f"    {model_name:<16} CV AUC {folds.mean():.3f} +/- {folds.std(ddof=1):.3f}")
    best_name = max(cv_results, key=lambda k: cv_results[k].mean())
    best = models(SEED)[best_name]
    observed = float(cv_results[best_name].mean())
    print(f"  selected: {best_name}")

    rng = np.random.default_rng(SEED)
    null = np.array(
        [cv_auc(best, Xtr, rng.permutation(ytr), 1, SEED + i).mean() for i in range(n_perm)]
    )
    null = null[np.isfinite(null)]
    p_value = float((np.sum(null >= observed) + 1) / (null.size + 1))
    print(
        f"  permutation null {null.mean():.3f} (95th pct {np.quantile(null, .95):.3f})  "
        f"observed {observed:.3f}  ->  p = {p_value:.4f}"
    )

    result = {
        "task": name,
        "n": int(len(y)),
        "n_positive": int(y.sum()),
        "cv_by_model": {k: float(v.mean()) for k, v in cv_results.items()},
        "selected_model": best_name,
        "cv_observed": observed,
        "null_mean": float(null.mean()),
        "null_p95": float(np.quantile(null, 0.95)),
        "p_value": p_value,
        "null": null.tolist(),
    }

    if locked and test_idx.size:
        fitted = clone(best).fit(Xtr, ytr)
        interval = bootstrap_auc_ci(y[test_idx], _scores(fitted, X[test_idx]), n_boot=4000, seed=SEED)
        print(f"  locked test set: {interval.describe()}")
        result["locked_test"] = {
            "auc": interval.auc,
            "ci_lower": interval.lower,
            "ci_upper": interval.upper,
            "n_positive": interval.n_positive,
            "n_negative": interval.n_negative,
        }
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--features", default="results/tcga_radiogenomics_features.xls")
    parser.add_argument("--labels", default="results/til_labels.json")
    parser.add_argument("--geometry", default="results/til_geometry.json")
    parser.add_argument("--out", default="results/geometry_experiment.json")
    args = parser.parse_args()

    repeats = 4 if args.quick else 20
    n_perm = 60 if args.quick else 300

    cohort = real_cohort.build(args.features, args.labels)
    geometry = json.loads(Path(args.geometry).read_text())
    by_patient = {r["case_id"]: r for r in geometry["records"]}

    keep = [i for i, pid in enumerate(cohort.patient_ids) if pid in by_patient]
    X = cohort.X[keep]
    patients = [cohort.patient_ids[i] for i in keep]
    phenotypes = [Immunophenotype(by_patient[p]["phenotype"]) for p in patients]

    print("=" * 78)
    print("MRI -> GEOMETRY-DERIVED IMMUNOPHENOTYPE")
    print("=" * 78)
    print(f"  imaging features available for {len(patients)} of "
          f"{geometry['n_patients']} phenotyped patients")
    print("  (geometry scales the LABEL to 136; the imaging features exist for 91)")
    for phenotype in Immunophenotype:
        count = sum(1 for p in phenotypes if p is phenotype)
        print(f"    {phenotype.value:<18} {count:3d}")

    tasks = []
    for task in (LabelTask.WHOLE, LabelTask.PERI):
        coded = [to_binary(p, task) for p in phenotypes]
        mask = np.array([c is not None for c in coded])
        y = np.array([c for c in coded if c is not None])
        label = (
            "RM-whole  ·  inflamed vs rest"
            if task is LabelTask.WHOLE
            else "RM-peri  ·  desert vs excluded  (inflamed cases dropped)"
        )
        tasks.append(run_task(label, X[mask], y, repeats, n_perm))

    print(f"\n{'=' * 78}\nMULTIPLICITY OVER THE TWO TASKS\n{'=' * 78}")
    rng = np.random.default_rng(SEED)
    nulls = [np.array(t["null"]) for t in tasks if len(t["null"]) > 5]
    if nulls:
        draws = np.array([[rng.choice(n) for n in nulls] for _ in range(20_000)])
        bar = float(np.quantile(draws.max(axis=1), 0.95))
        best = max(t["cv_observed"] for t in tasks)
        print(f"  bar across {len(nulls)} tasks: {bar:.3f}   best observed: {best:.3f}")
        print("  ->", "clears" if best > bar else "NOTHING CLEARS IT")
    else:
        bar, best = float("nan"), float("nan")

    payload = {
        "n_imaging": len(patients),
        "n_phenotyped": geometry["n_patients"],
        "phenotype_counts": {
            p.value: sum(1 for x in phenotypes if x is p) for p in Immunophenotype
        },
        "tasks": [{k: v for k, v in t.items() if k != "null"} for t in tasks],
        "multiplicity": {"bar": bar, "best_observed": best},
    }
    path = Path(args.out)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2))
    print(f"\nwrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
