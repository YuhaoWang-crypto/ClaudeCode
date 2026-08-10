#!/usr/bin/env python
"""What multi-centre data actually buys -- measured, not assumed.

Everything below holds the *true* per-feature effect size fixed, so any
difference between regimes is caused by the study design rather than by the
biology.  Four regimes, each evaluated with a leakage-free pipeline so that the
selection-leak measured in ``run_optimism_study.py`` cannot contaminate the
comparison:

``single_small``     one centre, n=182 -- the reference design's dimensions
``single_large``     one centre, n=900 -- isolates the effect of sample size alone
``pooled_random_cv`` many centres, random cross-validation
``loco_raw``         leave-one-centre-out, no harmonisation
``loco_combat``      leave-one-centre-out, with ComBat

The comparison ``single_large`` vs ``pooled_random_cv`` separates the two things
"multi-centre" bundles together: more patients, and more heterogeneity.  Only
the first is unambiguously good.

Two scenarios are run, because "centre effect" means two different things:

**nuisance**    scanners shift features, but phenotype prevalence is the same
                everywhere.  Centre is noise.
**confounded**  prevalence also varies by centre, so centre carries information
                about the outcome.  Here a pooled random split lets the model
                learn "which centre is this", and cross-centre performance
                collapses relative to the pooled estimate.

Real consortia are always somewhere between the two, and which one a given
cohort resembles is an empirical question about that cohort -- not something to
be assumed in either direction.
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
from sklearn.metrics import roc_auc_score  # noqa: E402

from mri_cd8_til.experiments.simulate import (  # noqa: E402
    SimulationConfig,
    effect_size_for_auc,
    simulate_cohort,
)
from mri_cd8_til.modeling.combat import ComBat, measure_batch_effect  # noqa: E402
from mri_cd8_til.modeling.pipeline import PipelineConfig, build_pipeline  # noqa: E402
from mri_cd8_til.modeling.validation import (  # noqa: E402
    leave_one_center_out,
    repeated_cv,
)

TRUE_FEATURE_AUC = 0.62

SCENARIOS = {
    "nuisance": {"shift": 0.6, "scale": 0.25, "prevalence_spread": 0.0},
    "confounded": {"shift": 0.6, "scale": 0.25, "prevalence_spread": 0.22},
}

REGIME_DESCRIPTIONS = {
    "single_small": "1 centre, n=182 (the reference design's size)",
    "single_large": "1 centre, n=N     (sample size alone)",
    "pooled_random_cv": "C centres pooled, random CV",
    "loco_raw": "leave-one-centre-out, raw features",
    "loco_combat": "leave-one-centre-out, ComBat harmonised",
}


def _config(n_patients, n_centers, shift, scale, spread, seed):
    return SimulationConfig(
        n_patients=n_patients,
        n_features=3332,
        n_informative=12,
        effect_size=effect_size_for_auc(TRUE_FEATURE_AUC),
        block_size=40,
        within_block_correlation=0.75,
        positive_fraction=0.35,
        n_centers=n_centers,
        center_shift_sd=shift,
        center_scale_sd=scale,
        center_prevalence_spread=spread,
        random_state=seed,
    )


def _loco_with_combat(X, y, center, pipeline_config):
    """Leave-one-centre-out, ComBat refitted inside every fold.

    ComBat sees only the training centres.  The held-out centre is standardised
    from its own features under the ``reference`` policy -- unsupervised, so no
    label leaks, but note that this still assumes the new centre arrives as a
    *batch*, not as one scan. Single-scan deployment is strictly harder.
    """
    scores = []
    for held_out in np.unique(center):
        test = center == held_out
        train = ~test
        if test.sum() < 10 or len(np.unique(y[test])) < 2 or len(np.unique(y[train])) < 2:
            continue
        harmoniser = ComBat(unseen_batch="reference").fit(X[train], batch=center[train])
        model = build_pipeline(pipeline_config).fit(
            harmoniser.transform(X[train], batch=center[train]), y[train]
        )
        decision = model.decision_function(harmoniser.transform(X[test], batch=center[test]))
        scores.append(float(roc_auc_score(y[test], decision)))
    return np.array(scores)


def run_scenario(name, params, args, pipeline_config):
    collected: dict[str, list[float]] = {k: [] for k in REGIME_DESCRIPTIONS}
    eta_before, eta_after = [], []

    for replicate in range(args.replicates):
        small = simulate_cohort(_config(182, 1, 0.0, 0.0, 0.0, replicate))
        collected["single_small"].append(
            repeated_cv(
                build_pipeline(pipeline_config), small.X, small.y, n_splits=5, n_repeats=2,
                random_state=replicate,
            ).mean
        )

        large = simulate_cohort(_config(args.n_pooled, 1, 0.0, 0.0, 0.0, 500 + replicate))
        collected["single_large"].append(
            repeated_cv(
                build_pipeline(pipeline_config), large.X, large.y, n_splits=5, n_repeats=2,
                random_state=replicate,
            ).mean
        )

        pooled = simulate_cohort(
            _config(
                args.n_pooled,
                args.n_centers,
                params["shift"],
                params["scale"],
                params["prevalence_spread"],
                1000 + replicate,
            )
        )
        collected["pooled_random_cv"].append(
            repeated_cv(
                build_pipeline(pipeline_config), pooled.X, pooled.y, n_splits=5, n_repeats=2,
                random_state=replicate,
            ).mean
        )
        collected["loco_raw"].append(
            leave_one_center_out(
                build_pipeline(pipeline_config), pooled.X, pooled.y, pooled.center
            ).mean
        )
        collected["loco_combat"].append(
            float(np.nanmean(_loco_with_combat(pooled.X, pooled.y, pooled.center, pipeline_config)))
        )

        eta_before.append(measure_batch_effect(pooled.X, pooled.center).median_eta_squared)
        harmonised = ComBat().fit_transform(pooled.X, batch=pooled.center)
        eta_after.append(measure_batch_effect(harmonised, pooled.center).median_eta_squared)

        if (replicate + 1) % 5 == 0:
            print(f"     ... {replicate + 1}/{args.replicates}")

    summary = {
        key: {
            "mean": float(np.nanmean(values)),
            "std": float(np.nanstd(values, ddof=1)) if len(values) > 1 else 0.0,
        }
        for key, values in collected.items()
    }
    summary["_batch_eta2_before"] = {"mean": float(np.mean(eta_before)), "std": 0.0}
    summary["_batch_eta2_after"] = {"mean": float(np.mean(eta_after)), "std": 0.0}
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates", type=int, default=20)
    parser.add_argument("--n-centers", type=int, default=6)
    parser.add_argument("--n-pooled", type=int, default=900)
    parser.add_argument("--out", default="results/multicenter_demo.json")
    args = parser.parse_args()
    pipeline_config = PipelineConfig(n_select=50)

    print("=" * 78)
    print("WHAT MULTI-CENTRE DATA BUYS")
    print("=" * 78)
    print(
        f"   True per-feature AUC held at {TRUE_FEATURE_AUC:.2f} in every regime.\n"
        f"   Pooled cohorts: {args.n_pooled} patients across {args.n_centers} unequal centres.\n"
        f"   Every regime uses a leakage-free pipeline, so nothing here is inflated by\n"
        f"   the selection leak measured in run_optimism_study.py.\n"
    )

    results = {}
    for name, params in SCENARIOS.items():
        print(f"\n   --- scenario: {name} (prevalence spread {params['prevalence_spread']:.2f}) ---")
        results[name] = run_scenario(name, params, args, pipeline_config)

    print("\n" + "=" * 78)
    print("RESULTS (mean AUC over replicates)")
    print("=" * 78)
    header = f"   {'regime':<20} {'nuisance':>16} {'confounded':>16}   description"
    print(header)
    print("   " + "-" * (len(header) + 12))
    for key, description in REGIME_DESCRIPTIONS.items():
        a = results["nuisance"][key]
        b = results["confounded"][key]
        print(
            f"   {key:<20} {a['mean']:>8.3f}+/-{a['std']:.3f} {b['mean']:>8.3f}+/-{b['std']:.3f}"
            f"   {description}"
        )

    print("\n   Reading the table:\n")
    for name in SCENARIOS:
        r = results[name]
        size_gain = r["single_large"]["mean"] - r["single_small"]["mean"]
        loco_gap = r["pooled_random_cv"]["mean"] - r["loco_raw"]["mean"]
        combat_gain = r["loco_combat"]["mean"] - r["loco_raw"]["mean"]
        best = r["loco_combat"]["mean"] - r["single_small"]["mean"]
        print(f"   [{name}]")
        print(f"     sample size alone (182 -> {args.n_pooled}):        {size_gain:+.3f} AUC")
        print(f"     pooled random CV over leave-one-centre-out: {loco_gap:+.3f} AUC of optimism")
        print(f"     ComBat, at the held-out centre:             {combat_gain:+.3f} AUC")
        print(
            f"     net: honest cross-centre AUC vs the reference design's internal "
            f"estimate: {best:+.3f}"
        )
        print(
            f"     variance explained by centre: "
            f"{r['_batch_eta2_before']['mean']:.3f} -> {r['_batch_eta2_after']['mean']:.3f}\n"
        )

    path = Path(args.out)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(
            {
                "true_feature_auc": TRUE_FEATURE_AUC,
                "n_replicates": args.replicates,
                "n_centers": args.n_centers,
                "n_pooled": args.n_pooled,
                "scenarios": results,
            },
            indent=2,
        )
    )
    print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
