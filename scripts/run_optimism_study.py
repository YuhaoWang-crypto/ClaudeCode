#!/usr/bin/env python
"""Measure how much AUC the reference study's analysis protocol manufactures.

Runs three things and writes them to ``results/``:

1. **Confidence intervals on the published numbers.**  Needs no simulation --
   just the reported AUCs and the class counts.
2. **The optimism experiment.**  The published protocol, a split-first variant
   and full cross-validation, run on synthetic cohorts with the reference
   study's dimensions across a grid of *known* true effect sizes.
3. **The dimension sweep.**  Under exactly zero true signal, reported AUC as a
   function of the number of candidate features, which isolates the mechanism.

Usage::

    python scripts/run_optimism_study.py             # full run (~40 min on 4 cores)
    python scripts/run_optimism_study.py --quick     # ~3 min, coarser estimates
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

from mri_cd8_til.cohorts import REFERENCE_STUDY  # noqa: E402
from mri_cd8_til.experiments.optimism import (  # noqa: E402
    run_optimism_experiment,
    sweep_feature_dimension,
)
from mri_cd8_til.experiments.simulate import SimulationConfig  # noqa: E402
from mri_cd8_til.modeling.metrics import auc_interval_from_reported  # noqa: E402

#: Class counts are not published per split, so a roughly balanced three-class
#: distribution is assumed. Stated as an assumption because the interval widths
#: depend on it -- though not enough to change the conclusion.
ASSUMED_SPLITS = [
    ("RM-whole  training  ", 0.973, 46, 91),
    ("RM-whole  validation", 0.985, 15, 30),
    ("RM-peri   training  ", 0.993, 46, 45),
    ("RM-peri   validation", 0.984, 15, 15),
]


def report_published_intervals() -> list[dict]:
    print("=" * 78)
    print("1. CONFIDENCE INTERVALS ON THE PUBLISHED AUCs (Hanley-McNeil)")
    print("=" * 78)
    print(f"   Source: {REFERENCE_STUDY['citation']}")
    print("   Assuming a roughly balanced three-class distribution.\n")

    rows = []
    for name, auc, n_pos, n_neg in ASSUMED_SPLITS:
        interval = auc_interval_from_reported(auc, n_pos, n_neg)
        print(f"   {name}  {interval.describe()}")
        rows.append(
            {
                "split": name.strip(),
                "reported_auc": auc,
                "n_positive": n_pos,
                "n_negative": n_neg,
                "ci_lower": interval.lower,
                "ci_upper": interval.upper,
                "ci_width": interval.width,
            }
        )

    reported = auc_interval_from_reported(0.985, 15, 30)
    weaker = auc_interval_from_reported(0.900, 15, 30)
    print(
        f"\n   The validation interval for 0.985 is [{reported.lower:.3f}, {reported.upper:.3f}].\n"
        f"   A distinctly weaker model at 0.900 gives [{weaker.lower:.3f}, {weaker.upper:.3f}].\n"
        "   Those overlap: 45 patients cannot separate a 0.90 model from a 0.985 one."
    )
    return rows


def run_optimism(quick: bool) -> dict:
    print("\n" + "=" * 78)
    print("2. OPTIMISM EXPERIMENT -- known true effect, three analysis protocols")
    print("=" * 78)
    config = SimulationConfig(
        n_patients=REFERENCE_STUDY["n_total"],
        n_features=REFERENCE_STUDY["features_per_region_per_phase"]
        * REFERENCE_STUDY["n_post_contrast_phases"],
        n_informative=12,
        block_size=40,
        within_block_correlation=0.75,
        positive_fraction=0.35,
    )
    print(
        f"   n={config.n_patients} patients, p={config.n_features} candidate features "
        f"({REFERENCE_STUDY['features_per_region_per_phase']} x "
        f"{REFERENCE_STUDY['n_post_contrast_phases']} phases), "
        f"137/45 split, top-50 selection, LASSO.\n"
    )

    grid = (0.50, 0.60, 0.70) if quick else (0.50, 0.55, 0.60, 0.65, 0.70, 0.75)
    result = run_optimism_experiment(
        true_feature_aucs=grid,
        n_replicates=15 if quick else 120,
        base_config=config,
        cv_repeats=2 if quick else 3,
        progress=True,
    )

    print("\n   'true_feature_AUC' is the discriminative power of each genuinely")
    print("   informative feature. 0.50 means there is no signal in the data at all.\n")
    header = f"   {'true AUC':>9} | {'published':>22} | {'split-first':>22} | {'repeated CV':>22}"
    print(header)
    print("   " + "-" * (len(header) - 3))
    for true_auc in result.true_feature_aucs:
        cells = []
        for protocol in ("published", "split_first", "repeated_cv"):
            s = result.get(protocol, true_auc)
            cells.append(f"{s.median:.3f} [{s.quantile(0.05):.2f},{s.quantile(0.95):.2f}]")
        print(f"   {true_auc:>9.2f} | {cells[0]:>22} | {cells[1]:>22} | {cells[2]:>22}")

    print("\n   AUC manufactured by selecting features before splitting:")
    for true_auc in result.true_feature_aucs:
        print(
            f"     true_feature_AUC={true_auc:.2f}:  "
            f"+{result.optimism(true_auc):.3f} over split-first, "
            f"+{result.optimism(true_auc, 'repeated_cv'):.3f} over repeated CV"
        )

    zero = result.get("published", 0.50)
    print(
        f"\n   With NO true signal whatsoever, the published protocol still reports a\n"
        f"   median validation AUC of {zero.median:.3f} (90% of runs land in "
        f"[{zero.quantile(0.05):.3f}, {zero.quantile(0.95):.3f}])."
    )
    return {"table": result.as_table(), "grid": list(result.true_feature_aucs)}


def run_sweep(quick: bool) -> dict:
    print("\n" + "=" * 78)
    print("3. DIMENSION SWEEP -- zero signal, varying the number of candidate features")
    print("=" * 78)
    grid = (200, 3332) if quick else (50, 200, 800, 1666, 3332, 6664)
    sweep = sweep_feature_dimension(
        n_features_grid=grid,
        n_replicates=15 if quick else 120,
        base_config=SimulationConfig(n_patients=182, n_informative=0, effect_size=0.0),
    )
    print("   Every number below comes from data containing no signal at all.\n")
    print(sweep.describe())
    print(
        "\n   Reported AUC rises with the size of the candidate pool. The pool is the\n"
        "   cause: more features means more chances for one to separate the 45\n"
        "   validation patients by luck."
    )
    return {
        "n_features": list(sweep.n_features),
        "median_auc": sweep.median_auc.tolist(),
        "p95_auc": sweep.p95_auc.tolist(),
        "frac_above_0.95": sweep.frac_above_95.tolist(),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quick", action="store_true", help="coarse run for smoke-testing")
    parser.add_argument("--out", default="results/optimism_study.json")
    args = parser.parse_args()

    np.random.seed(0)
    payload = {
        "reference_study": REFERENCE_STUDY,
        "published_intervals": report_published_intervals(),
        "optimism": run_optimism(args.quick),
        "dimension_sweep": run_sweep(args.quick),
        "quick_mode": args.quick,
    }

    path = Path(args.out)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, default=str))
    print(f"\nwrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
