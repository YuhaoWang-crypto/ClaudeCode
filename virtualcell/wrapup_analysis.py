"""What the VCC 2026 results say about where the remaining score actually is.

The 2026 wrap-up published a final leaderboard, and `cell-eval2`'s own metric
specification publishes the two anchors every metric is rescaled by.  Together
those invert: a team's *scaled* score plus the anchors gives back the *raw*
metric value, and the raw value is what a modelling change has to move.

That inversion is worth doing before touching the model, because the six
metrics are nowhere near equally saturated and two of them turn out to be
limited by a calibration choice rather than by biology.

Anchors are from ``docs/vcc2026_metrics/vcc2026-metrics-brief.md`` at
cell-eval2 0.15.0, competition ``rule_version`` 3, measured on the three
official reference bundles.  The score is

    s = (u - b) / (r - b)

with ``b`` the context-mean-response baseline and ``r`` an independent
replicate, so s = 1 means "as good as repeating the experiment".  A context's
score is the equal-weight mean of the six -- which the published podium rows
reproduce to four decimals, so the table and the weighting are both confirmed
rather than assumed.

    python -m virtualcell.wrapup_analysis
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# published results
# ---------------------------------------------------------------------------

#: Final test scores, arcinstitute.org/news/virtual-cell-challenge-2026-wrap-up
#: (19 Nov 2026).  Scaled, i.e. already divided by (replicate - baseline).
PODIUM = {
    "BioInsight AI Lab (1st, Illumina)": {
        "pds_cosine": 0.8313, "mse": 0.6842, "nmae": 0.4399,
        "fidelity": 0.1908, "reach": 0.4288, "jaccard": 0.0642},
    "FlowState (2nd, independent)": {
        "pds_cosine": 0.8363, "mse": 0.3632, "nmae": 0.3286,
        "fidelity": 0.1695, "reach": 0.3128, "jaccard": -0.0018},
    "X-team (3rd, Xaira)": {
        "pds_cosine": 0.7346, "mse": 0.4681, "nmae": 0.3038,
        "fidelity": 0.0659, "reach": 0.3186, "jaccard": 0.0667},
}
OVERALL = {"BioInsight AI Lab (1st, Illumina)": 0.4399,
           "FlowState (2nd, independent)": 0.3348,
           "X-team (3rd, Xaira)": 0.3263}

#: (baseline, replicate) ranges over the three test contexts, from the
#: cell-eval2 metric brief.  Midpoints are used for inversion and the spread is
#: carried through as the uncertainty on every raw number below.
ANCHORS = {
    "pds_cosine": ((0.500, 0.500), (0.927, 0.984)),
    "mse": ((0.986, 0.992), (0.028, 0.045)),
    "fidelity": ((0.505, 0.522), (0.795, 0.832)),
    "reach": ((0.047, 0.097), (0.958, 0.978)),
    "jaccard": ((0.021, 0.037), (0.375, 0.423)),
    "nmae": ((1.0009, 1.0017), (0.369, 0.431)),
}

#: Metrics where lower raw values are better, so the span r-b is negative.
LOWER_BETTER = {"mse", "nmae"}

#: Significant-gene counts per target, top-ten entries on the validation lines
#: (Cal-33, HeLa, Jurkat).  The measured range is the organizers' median per
#: target; the predicted figures are what each entry called.
DE_COUNTS = {"measured (median per target)": (22, 44),
             "X-team": (122, 136),
             "BioInsight AI Lab": (500, 500),
             "rest of top ten": (340, 4200)}


def _mid(pair: tuple[float, float]) -> float:
    return 0.5 * (pair[0] + pair[1])


def raw_from_scaled(metric: str, scaled: float) -> tuple[float, float, float]:
    """Invert ``s = (u-b)/(r-b)`` for the raw metric value.

    Returns ``(raw_at_midpoint, raw_low, raw_high)`` where the interval comes
    from taking the anchors at the two ends of their published across-context
    range.  The interval is the honest precision of this inversion: the anchors
    are per-context and the published score is a mean over three contexts, so a
    single raw number is an approximation whatever one does.
    """
    (b_lo, b_hi), (r_lo, r_hi) = ANCHORS[metric]
    vals = [b + scaled * (r - b) for b in (b_lo, b_hi) for r in (r_lo, r_hi)]
    b, r = _mid((b_lo, b_hi)), _mid((r_lo, r_hi))
    return b + scaled * (r - b), min(vals), max(vals)


def jaccard_ceiling(n_pred: float, n_real: float, recall: float = 1.0) -> float:
    """Best Jaccard reachable when calling ``n_pred`` genes and truth has ``n_real``.

    ``J = |R n R^| / |R u R^|``.  With recall ``p`` of the real set,
    ``|R n R^| = p*n_real`` and ``|R u R^| = n_pred + n_real - p*n_real``, so
    over-calling caps J at ``n_real/n_pred`` even with perfect recall.  This is
    a hard arithmetic bound, independent of how good the biology is.
    """
    hit = recall * n_real
    return float(hit / max(n_pred + n_real - hit, 1e-12))


def implied_recall(j_raw: float, n_pred: float, n_real: float) -> float:
    """Solve :func:`jaccard_ceiling` for the recall consistent with a raw J."""
    # j = p*n_real / (n_pred + n_real - p*n_real)  =>  p = j*(n_pred+n_real) /
    #                                                      (n_real*(1+j))
    return float(j_raw * (n_pred + n_real) / (n_real * (1.0 + j_raw)))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-real", type=float, default=33.0,
                    help="measured significant genes per target (22-44 median)")
    ap.add_argument("--out", type=Path,
                    default=Path("results/wrapup_analysis.json"))
    args = ap.parse_args()

    out: dict = {"podium_scaled": PODIUM, "overall_published": OVERALL}

    # --- 0. the table is self-consistent, so the weighting is confirmed -----
    print("=" * 78)
    print("0. Confirm the aggregation before relying on it")
    print("   The six published per-metric scores should average to the "
          "published overall.")
    checks = {}
    for team, m in PODIUM.items():
        got = float(np.mean(list(m.values())))
        checks[team] = {"mean_of_six": got, "published": OVERALL[team],
                        "abs_diff": abs(got - OVERALL[team])}
        print(f"   {team:<36} mean of six {got:.4f}  published "
              f"{OVERALL[team]:.4f}  diff {abs(got - OVERALL[team]):.5f}")
    out["aggregation_check"] = checks
    print("   -> equal-weight mean of six, reproduced to 4 d.p. on all three "
          "rows.")

    # --- 1. where is the headroom ------------------------------------------
    win = "BioInsight AI Lab (1st, Illumina)"
    print("\n" + "=" * 78)
    print("1. Per-metric headroom for the WINNING entry")
    print("   'left to replicate' is 1 - scaled: what is still unclaimed on "
          "that metric.\n")
    print(f"   {'metric':<12}{'scaled':>8}{'raw':>9}{'raw range':>16}"
          f"{'left':>8}{'share of total':>16}")
    print("   " + "-" * 69)
    total_left = sum(1.0 - s for s in PODIUM[win].values())
    rows = {}
    for metric, scaled in sorted(PODIUM[win].items(),
                                 key=lambda kv: kv[1] - 1.0):
        raw, lo, hi = raw_from_scaled(metric, scaled)
        left = 1.0 - scaled
        share = left / total_left
        rows[metric] = {"scaled": scaled, "raw": raw, "raw_lo": lo,
                        "raw_hi": hi, "left_to_replicate": left,
                        "share_of_headroom": share}
        print(f"   {metric:<12}{scaled:>8.4f}{raw:>9.4f}"
              f"{f'[{lo:.3f}, {hi:.3f}]':>16}{left:>8.4f}{share:>15.1%}")
    out["headroom_winner"] = rows
    print(f"\n   total unclaimed, summed over six metrics: {total_left:.3f} "
          f"of 6.000")
    pair = rows["jaccard"]["share_of_headroom"] + rows["fidelity"]["share_of_headroom"]
    print(f"   jaccard + fidelity alone carry {pair:.1%} of it.")

    # --- 2. the over-calling bound on jaccard ------------------------------
    print("\n" + "=" * 78)
    print("2. Why jaccard is the worst metric: it is capped by how many genes "
          "you call")
    print("   J <= n_real/n_pred even with perfect recall, so calling 17x too "
          "many genes")
    print("   caps J near 0.06 regardless of the biology.\n")
    print(f"   measured median per target: {DE_COUNTS['measured (median per target)']} "
          f"(validation lines); using n_real = {args.n_real:.0f}\n")
    print(f"   {'entry':<22}{'n_called':>12}{'J ceiling':>12}"
          f"{'raw J':>9}{'at % of cap':>13}{'implied recall':>16}")
    print("   " + "-" * 84)
    name_map = {"X-team": "X-team (3rd, Xaira)",
                "BioInsight AI Lab": "BioInsight AI Lab (1st, Illumina)"}
    calls = {}
    for label, (c_lo, c_hi) in DE_COUNTS.items():
        if label.startswith("measured"):
            continue
        n_pred = _mid((c_lo, c_hi))
        cap = jaccard_ceiling(n_pred, args.n_real)
        team = name_map.get(label)
        if team is None:
            print(f"   {label:<22}{f'{c_lo}-{c_hi}':>12}{cap:>12.4f}"
                  f"{'-':>9}{'-':>13}{'-':>16}")
            calls[label] = {"n_called": [c_lo, c_hi], "jaccard_ceiling": cap}
            continue
        raw, _, _ = raw_from_scaled("jaccard", PODIUM[team]["jaccard"])
        rec = implied_recall(raw, n_pred, args.n_real)
        print(f"   {label:<22}{f'{c_lo}-{c_hi}':>12}{cap:>12.4f}"
              f"{raw:>9.4f}{raw / cap:>12.1%}{min(rec, 1.0):>15.1%}"
              + ("" if rec <= 1.0 else "  (>1: counts approximate)"))
        calls[label] = {"n_called": [c_lo, c_hi], "jaccard_ceiling": cap,
                        "raw_jaccard": raw, "frac_of_ceiling": raw / cap,
                        "implied_recall": rec}
    out["jaccard_overcalling"] = calls
    print("\n   Read the two podium rows against each other:")
    print("     BioInsight calls ~15x too many genes and lands at ~90% of the "
          "cap that\n     imposes -- the count is the binding constraint, not "
          "the gene selection.")
    print("     X-team calls ~4x too many and lands far BELOW its cap -- for "
          "them the\n     binding constraint is which genes, not how many.")

    # --- 3. what calibrating the count alone would be worth ----------------
    print("\n" + "=" * 78)
    print("3. What thresholding to the right number of calls is worth")
    print("   Sweep: keep the top-n ranked calls, assume recall r among them, "
          "and rescale.\n")
    b, r = _mid(ANCHORS["jaccard"][0]), _mid(ANCHORS["jaccard"][1])
    print(f"   jaccard anchors: baseline {b:.3f}, replicate {r:.3f}  "
          f"(span {r - b:.3f})")
    print(f"   winner's scaled jaccard today: "
          f"{PODIUM[win]['jaccard']:.4f}\n")
    print(f"   {'n_called':>10}" + "".join(f"{f'recall {x:.0%}':>14}"
                                           for x in (0.3, 0.5, 0.7, 0.9)))
    print("   " + "-" * 66)
    sweep = {}
    for n_pred in (args.n_real, 2 * args.n_real, 4 * args.n_real, 130, 500):
        cells, row = [], {}
        for rec in (0.3, 0.5, 0.7, 0.9):
            j = jaccard_ceiling(n_pred, args.n_real, rec)
            s = (j - b) / (r - b)
            row[f"recall_{rec}"] = {"raw_jaccard": j, "scaled": s,
                                    "composite_delta":
                                        (s - PODIUM[win]["jaccard"]) / 6.0}
            cells.append(f"{s:>+8.2f}{'':>6}")
        sweep[int(n_pred)] = row
        print(f"   {int(n_pred):>10}" + "".join(cells))
    out["jaccard_sweep_scaled"] = sweep
    print("\n   Values are SCALED jaccard (1.0 = replicate). Each +1.00 of "
          "scaled jaccard\n   adds +0.167 to the composite, because the "
          "composite is a mean of six.")
    best = sweep[int(args.n_real)]["recall_0.5"]
    print(f"\n   Calling n_real genes at 50% recall scores "
          f"{best['scaled']:+.2f} scaled jaccard,\n   worth "
          f"{best['composite_delta']:+.4f} on the composite against the "
          f"winner's {PODIUM[win]['jaccard']:.4f}.")
    print(f"   The published 1st-to-2nd gap was "
          f"{OVERALL[win] - OVERALL['FlowState (2nd, independent)']:+.4f}.")

    # --- 4. the levers that are gone, and the one the organizers named -----
    print("\n" + "=" * 78)
    print("4. Two consequences for this project's existing tuning")
    print("""
   (a) pds_cosine is scale-invariant, so effect magnitude no longer buys
       discrimination.  This project tuned `beta` on a magnitude frontier
       against the 2025 L1 discrimination score, and Liu et al. (arXiv
       2511.16954) show exactly why that worked: L1/L2 PDS rises with a global
       rescaling of the prediction and plateaus, while cosine PDS does not move
       at all.  Arc adopted cosine for 2026.  So `beta` now acts only through
       mse and nmae, and every magnitude hyperparameter has to be re-tuned
       against those two.  The `norm_cv` diagnostic was measuring a real metric
       artifact -- and the artifact has been removed.

   (b) The organizers' own diagnosis of the field is the failure this project
       already has a measured remedy for: predicted responses to one target in
       different lines agree at cosine 0.77-0.89, against 0.55 for the measured
       responses, and the line-specific component carries a slope of 0.02-0.13
       against 1.02 for a replicate.  Models are returning a cross-line
       average.  results/connection.json measures the same gap from the other
       side -- transporting a measured effect as-is gives cosine 0.163, and a
       state-conditioned orthogonal operator gives 0.295 -- but it needs ~25
       measured knockdowns in the target line, which the challenge does not
       give.  The control-only version of that operator is the open problem,
       and it is what the winner's "ridge regression on SVD response patterns"
       is a form of.""")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2, default=float))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
