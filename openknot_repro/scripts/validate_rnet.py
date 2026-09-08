"""Check the RibonanzaNet reproduction against the release, and test the Fig. S4 claim.

Three questions, all answerable from `results/rnet_predictions.csv`:

1. Does this RNet reproduce the secondary structures in the release?  The
   checkpoints are bit-identical to the official ones, so any disagreement is
   in how pair probabilities are decoded into a structure.

2. Does the RNet_F1 filter the paper applies before its solid bars behave the
   same way here?  That is what actually matters downstream: the filter keeps
   designs at F1 >= 0.8, so what counts is agreement on that decision, not on
   the structure string.

3. Fig. S4: does an OpenKnot score computed from RNet-predicted reactivity
   track the experimental one, "especially for poorly performing designs"?
   This is the claim that licenses prospective filtering of designs without an
   experiment, so it is worth testing rather than assuming.

Usage:  python scripts/validate_rnet.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

RESULTS = Path(__file__).resolve().parent.parent / "results"
PREDICTIONS = RESULTS / "rnet_predictions.csv"
SUCCESS_CUTOFF = 90.0
RNET_F1_CUTOFF = 0.8


def spearman(a: np.ndarray, b: np.ndarray) -> float:
    ranks = lambda x: pd.Series(x).rank().to_numpy()  # noqa: E731
    ra, rb = ranks(a), ranks(b)
    return float(np.corrcoef(ra, rb)[0, 1])


def main() -> int:
    if not PREDICTIONS.exists():
        print("run scripts/rnet_predict.py first", file=sys.stderr)
        return 1
    data = pd.read_csv(PREDICTIONS)
    print(f"{len(data)} designs with RNet predictions "
          f"({data.groupby('round').size().to_dict()} per round)\n")

    # 1. structure agreement with the released RNet_structure column
    agreement = data["agreement_with_release"]
    identical = float((agreement > 0.999).mean())
    print("1. secondary structure vs the released RNet_structure column")
    print(f"   mean base-pair F1   {agreement.mean():.4f}")
    print(f"   median              {agreement.median():.4f}")
    print(f"   identical pair sets {100 * identical:.1f}%")

    # 2. the filter decision
    mine = data["rnet_F1"]
    released = data["released_RNet_F1"]
    keep_mine = mine >= RNET_F1_CUTOFF
    keep_released = released >= RNET_F1_CUTOFF
    same_decision = float((keep_mine == keep_released).mean())
    print("\n2. RNet_F1 (design vs target) compared with the released column")
    print(f"   Spearman            {spearman(mine, released):.4f}")
    print(f"   mean |difference|   {np.abs(mine - released).mean():.4f}")
    print(f"   exact agreement     {100 * float(np.isclose(mine, released).mean()):.1f}%")
    print(f"   same filter call at F1 >= {RNET_F1_CUTOFF}: {100 * same_decision:.1f}% "
          f"({int((keep_mine & keep_released).sum())} kept by both, "
          f"{int((~keep_mine & ~keep_released).sum())} dropped by both)")

    # 3. Fig. S4: simulated vs experimental OpenKnot score
    measured = data["target_openknot_score"]
    simulated = data["simulated_openknot_score"]
    usable = data[data["SN_filter"] == 1]
    print("\n3. simulated OpenKnot score (from predicted reactivity) vs experimental")
    print(f"   Pearson             {np.corrcoef(simulated, measured)[0, 1]:.4f}")
    print(f"   Spearman            {spearman(simulated, measured):.4f}")
    print(f"   Pearson, SN_filter=1 only ({len(usable)} designs): "
          f"{np.corrcoef(usable['simulated_openknot_score'], usable['target_openknot_score'])[0, 1]:.4f}")

    bins = [(0, 70), (70, 80), (80, 90), (90, 101)]
    print("\n   agreement by experimental score band (the Fig. S4 claim is that the")
    print("   simulation is most faithful for poor designs):")
    rows = []
    for low, high in bins:
        band = data[(measured >= low) & (measured < high)]
        if band.empty:
            continue
        signed = band["simulated_openknot_score"] - band["target_openknot_score"]
        rho = (
            spearman(band["simulated_openknot_score"], band["target_openknot_score"])
            if len(band) > 5
            else float("nan")
        )
        rows.append({
            "band": f"{low}-{high}",
            "n": len(band),
            "mean_abs_error": signed.abs().mean(),
            "mean_signed_error": signed.mean(),
            "spearman_within_band": rho,
        })
        print(f"     experimental {low:3d}-{high:3d}: n={len(band):4d}  "
              f"mean |error| {signed.abs().mean():6.2f}  "
              f"mean signed error {signed.mean():+6.2f}  "
              f"within-band Spearman {rho:+.3f}")

    # prospective filtering: does a low simulated score predict a low measured one?
    predicted_fail = simulated <= SUCCESS_CUTOFF
    actually_failed = measured <= SUCCESS_CUTOFF
    true_positive = int((predicted_fail & actually_failed).sum())
    precision = true_positive / max(int(predicted_fail.sum()), 1)
    recall = true_positive / max(int(actually_failed.sum()), 1)
    print(f"\n   using simulated score <= {SUCCESS_CUTOFF:.0f} to predict experimental failure:")
    print(f"     precision {precision:.3f}   recall {recall:.3f}   "
          f"(base rate {float(actually_failed.mean()):.3f})")

    summary = pd.DataFrame(rows)
    summary.to_csv(RESULTS / "rnet_validation_bands.csv", index=False)
    pd.DataFrame([{
        "n_designs": len(data),
        "structure_mean_pair_f1_vs_release": agreement.mean(),
        "structure_identical_fraction": identical,
        "rnet_f1_spearman_vs_release": spearman(mine, released),
        "rnet_f1_filter_agreement": same_decision,
        "simulated_vs_measured_pearson": float(np.corrcoef(simulated, measured)[0, 1]),
        "simulated_vs_measured_spearman": spearman(simulated, measured),
        "failure_prediction_precision": precision,
        "failure_prediction_recall": recall,
    }]).to_csv(RESULTS / "rnet_validation.csv", index=False)
    print(f"\nwrote {RESULTS}/rnet_validation.csv and rnet_validation_bands.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
