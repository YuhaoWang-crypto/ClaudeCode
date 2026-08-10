#!/usr/bin/env python
"""How many targets has this project tried, and what bar does that set?

Each experiment corrects for the variants *it* tried.  Nobody corrects for the
variants tried across the whole project, and that is where a result quietly
becomes unreliable: a bar computed over 8 variants is the right bar only if 8 is
all that was ever attempted.

This script collects every binary target evaluated anywhere in the project and
recomputes the bar over all of them together.  It exists because the alternative
-- reporting the best target from a long search under its own local correction
-- is precisely the failure mode the project set out to measure.

Method: each target contributes its permutation null, summarised by mean and
95th percentile.  Where the full null was not retained, a normal approximation
recovers the spread as ``sd = (p95 - mean) / 1.645``.  That approximation is
stated rather than hidden; permutation nulls of a CV mean are close to normal by
the central limit theorem, and the bar it produces is within a few thousandths
of the empirical one where both are available.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402

RESULTS = Path("results")
Z95 = 1.6448536269514722


def collect() -> list[dict]:
    """Every binary AUC target the project has evaluated, with its null."""
    targets: list[dict] = []

    optimization = RESULTS / "model_optimization.json"
    if optimization.exists():
        payload = json.loads(optimization.read_text())
        for variant in payload["variants"]:
            if variant["metric"] != "auc":
                continue
            targets.append(
                {
                    "source": "model_optimization",
                    "name": variant["name"],
                    "observed": variant["observed"],
                    "null_mean": variant["null_mean"],
                    "null_p95": variant["null_p95"],
                }
            )

    geometry = RESULTS / "geometry_experiment.json"
    if geometry.exists():
        payload = json.loads(geometry.read_text())
        for task in payload["tasks"]:
            targets.append(
                {
                    "source": "geometry_experiment",
                    "name": task["task"].split("·")[0].strip(),
                    "observed": task["cv_observed"],
                    "null_mean": task["null_mean"],
                    "null_p95": task["null_p95"],
                    "locked_test_auc": task.get("locked_test", {}).get("auc"),
                    "locked_ci": (
                        task.get("locked_test", {}).get("ci_lower"),
                        task.get("locked_test", {}).get("ci_upper"),
                    ),
                }
            )

    real = RESULTS / "real_experiment.json"
    if real.exists():
        payload = json.loads(real.read_text())
        # The tertile target also appears in model_optimization; keep it once.
        if not any(t["name"].startswith("tertile") for t in targets):
            targets.append(
                {
                    "source": "real_experiment",
                    "name": "tertile TIL (primary)",
                    "observed": payload["permutation_null"]["observed"],
                    "null_mean": payload["permutation_null"]["mean"],
                    "null_p95": payload["permutation_null"]["p95"],
                    "locked_test_auc": payload["locked_test"]["auc"],
                    "locked_ci": (
                        payload["locked_test"]["ci_lower"],
                        payload["locked_test"]["ci_upper"],
                    ),
                }
            )
    return targets


def project_bar(targets: list[dict], n_draws: int = 200_000, seed: int = 0) -> float:
    """95th percentile of the maximum across every target's null."""
    rng = np.random.default_rng(seed)
    means = np.array([t["null_mean"] for t in targets])
    sds = np.maximum((np.array([t["null_p95"] for t in targets]) - means) / Z95, 1e-6)
    draws = rng.normal(means, sds, size=(n_draws, len(targets)))
    return float(np.quantile(draws.max(axis=1), 0.95))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default="results/project_multiplicity.json")
    args = parser.parse_args()

    targets = collect()
    if not targets:
        print("no results found — run the experiments first")
        return 1

    print("=" * 82)
    print("PROJECT-WIDE MULTIPLICITY ACCOUNTING")
    print("=" * 82)
    print(f"\nEvery binary target evaluated anywhere in this project ({len(targets)}):\n")
    print(f"  {'target':<34}{'source':<22}{'CV':>8}{'null p95':>10}")
    print("  " + "-" * 72)
    for target in sorted(targets, key=lambda t: -t["observed"]):
        print(
            f"  {target['name'][:33]:<34}{target['source']:<22}"
            f"{target['observed']:>8.3f}{target['null_p95']:>10.3f}"
        )

    bar = project_bar(targets)
    best = max(targets, key=lambda t: t["observed"])
    print(f"\n  bar over all {len(targets)} targets (95th pct of the max under the null): {bar:.3f}")
    print(f"  best cross-validated observation: {best['observed']:.3f}  ({best['name']})")
    clears = best["observed"] > bar
    print(f"  -> {'CLEARS the project-wide bar' if clears else 'does NOT clear the project-wide bar'}")

    with_test = [t for t in targets if t.get("locked_test_auc") is not None]
    if with_test:
        print("\n  Targets that also reached a locked test set:\n")
        for target in sorted(with_test, key=lambda t: -(t["locked_test_auc"] or 0)):
            low, high = target["locked_ci"]
            excludes = low is not None and low > 0.5
            print(
                f"    {target['name'][:36]:<38} test AUC {target['locked_test_auc']:.3f}  "
                f"95% CI [{low:.3f}, {high:.3f}]  "
                f"{'CI excludes chance' if excludes else 'CI includes chance'}"
            )
        print(
            "\n  A locked test set is not protected by the CV multiplicity bar -- it is a\n"
            "  separate, single evaluation. But it is also not immune to the search that\n"
            "  preceded it: the target it evaluates was chosen after seeing the others."
        )

    payload = {
        "n_targets": len(targets),
        "project_bar": bar,
        "best_cv": best["observed"],
        "best_name": best["name"],
        "clears_project_bar": bool(clears),
        "targets": targets,
    }
    path = Path(args.out)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2))
    print(f"\nwrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
