"""Stage 5 - panel-level readout: hits per target, and selectivity across targets.

Screening one library against a whole panel buys something a single-target run
cannot give: for each compound, a predicted activity *profile*. That is the in
silico counterpart of what a profiling service sells, and it is where the
interesting failure modes live.

Two panel-level quantities are computed, and both are constructed so they cannot
flatter the panel:

**Selectivity** is reported as the gap between a compound's best target and its
second best, in log units, and only over targets where that compound is inside
the applicability domain. A compound predicted active everywhere is either a
genuine polypharmacology case or a promiscuous-prediction artifact, and the
`n_targets_high_ad` column lets a reader tell how much evidence the profile
rests on.

**Panel-wide frequent hitters** are flagged explicitly. A compound in the top
decile of many unrelated targets is more likely to be exploiting a fingerprint
shortcut than to be a real multi-target inhibitor, and the report names them
rather than letting them sit at the top of every list.

Hits are only called hits inside the applicability domain. Everything else is
written out but labelled.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"

TOP_DECILE = 0.90
HIT_MIN_PAFFINITY = 6.0  # ~1 uM predicted; below this a repurposing hit is not worth a plate
FREQUENT_HITTER_MIN_TARGETS = 4


def main() -> None:
    long = pd.read_csv(RESULTS / "panel_predictions_long.csv.gz")
    bench = {b["name"]: b for b in json.loads((RESULTS / "benchmark.json").read_text())}
    prov = {t["name"]: t for t in json.loads((RESULTS / "panel_provenance.json").read_text())["targets"]}
    controls = json.loads((RESULTS / "control_recall.json").read_text())

    targets = sorted(long["target"].unique())
    print(f"panel: {len(targets)} targets x {long['chembl_id'].nunique()} compounds\n")

    # --- per-target hit lists ----------------------------------------------
    hit_rows = []
    per_target = []
    for t in targets:
        sub = long[long["target"] == t]
        in_ad = sub[sub["ad_tier"].eq("high")]
        # Novel-to-this-target hits: exclude compounds already in the training
        # set, whose "prediction" is partly memorisation.
        cand = in_ad[(~in_ad["in_training"]) & (in_ad["pred_pAffinity"] >= HIT_MIN_PAFFINITY)]
        cand = cand.sort_values("pred_pAffinity", ascending=False)
        b = bench[t]
        m = b["models"][b["best_model"]]
        per_target.append(
            {
                "target": t,
                "bps_family": prov[t]["bps_family"],
                "model": b["best_model"],
                "scaffold_spearman": m["scaffold_spearman_mean"],
                "scaffold_rmse": m["scaffold_rmse_mean"],
                "noise_floor_log": b["noise_floor_log"],
                "n_train": b["n_compounds"],
                "n_in_ad_high": int(len(in_ad)),
                "n_hits": int(len(cand)),
                "n_hits_pains_flagged": int((cand["pains_alert"].fillna("") != "").sum()),
                "top_hit": (
                    f"{cand.iloc[0]['pref_name']} ({cand.iloc[0]['pred_pAffinity']:.2f})"
                    if len(cand)
                    else None
                ),
            }
        )
        cand = cand.assign(rank_in_target=np.arange(1, len(cand) + 1))
        hit_rows.append(cand)
        print(
            f"  {t:<10s} rho {m['scaffold_spearman_mean']:+.3f}  "
            f"RMSE {m['scaffold_rmse_mean']:.2f} vs noise {b['noise_floor_log']}  "
            f"AD-high {len(in_ad):>5d}  hits {len(cand):>4d}"
            + (f"  top: {cand.iloc[0]['pref_name']}" if len(cand) else "")
        )

    hits = pd.concat(hit_rows, ignore_index=True) if hit_rows else pd.DataFrame()
    hits.to_csv(RESULTS / "panel_hits.csv", index=False)
    pd.DataFrame(per_target).to_csv(RESULTS / "panel_target_summary.csv", index=False)

    # --- per-compound activity profile and selectivity ---------------------
    wide = long.pivot_table(index="chembl_id", columns="target", values="pred_pAffinity")
    tier = long.pivot_table(index="chembl_id", columns="target", values="ad_tier", aggfunc="first")
    meta = long.groupby("chembl_id").agg(
        pref_name=("pref_name", "first"),
        max_phase=("max_phase", "first"),
        pains_alert=("pains_alert", "first"),
    )

    # Only in-domain predictions may contribute to a selectivity claim.
    masked = wide.where(tier.eq("high"))
    n_high = tier.eq("high").sum(axis=1)

    # A compound can be in-domain on zero targets; idxmax raises on an all-NA
    # row, so those are held out of the profile columns and reported as such
    # rather than dropped silently.
    has_any = masked.notna().any(axis=1)
    best_t = pd.Series(index=masked.index, dtype=object)
    best_t.loc[has_any] = masked.loc[has_any].idxmax(axis=1)
    best_v = masked.max(axis=1)
    second = masked.apply(
        lambda r: r.dropna().nlargest(2).iloc[1] if r.notna().sum() >= 2 else np.nan, axis=1
    )
    print(
        f"  {int(has_any.sum())} of {len(masked)} compounds are in-domain on at "
        f"least one target; the rest have no usable profile"
    )

    # Frequent hitters: top-decile on many targets, computed per target so a
    # target with a generally higher predicted range cannot dominate.
    thresh = wide.quantile(TOP_DECILE)
    top_decile_count = (wide >= thresh).sum(axis=1)

    profile = meta.join(
        pd.DataFrame(
            {
                "n_targets_high_ad": n_high,
                "best_target": best_t,
                "best_pred_pAffinity": best_v.round(3),
                "second_best_pred": second.round(3),
                "selectivity_gap_log": (best_v - second).round(3),
                "n_targets_top_decile": top_decile_count,
                "frequent_hitter": top_decile_count >= FREQUENT_HITTER_MIN_TARGETS,
            }
        )
    ).join(wide.round(3), rsuffix="_pred")
    profile = profile.sort_values("best_pred_pAffinity", ascending=False)
    profile.to_csv(RESULTS / "compound_profiles.csv")

    fh_n = int(profile["frequent_hitter"].sum())
    print(
        f"\ncompound profiles: {len(profile)} compounds; "
        f"{fh_n} flagged frequent hitters (top decile on >= {FREQUENT_HITTER_MIN_TARGETS} targets)"
    )

    sel = profile[(profile["n_targets_high_ad"] >= 3) & (~profile["frequent_hitter"])]
    sel = sel.dropna(subset=["selectivity_gap_log"]).nlargest(10, "selectivity_gap_log")
    print("\nmost selective in-domain profiles (>=3 in-domain targets, not frequent hitters):")
    for cid, r in sel.iterrows():
        print(
            f"  {str(r['pref_name'])[:30]:<30s} best {r['best_target']:<9s} "
            f"{r['best_pred_pAffinity']:.2f}  gap {r['selectivity_gap_log']:.2f} log  "
            f"(in-domain on {int(r['n_targets_high_ad'])} targets)"
        )

    # --- control recall, panel-wide ----------------------------------------
    print("\nreference-inhibitor recall (did the virtual assay rank the real assay's controls high?)")
    recall_rows = []
    for t, recs in controls.items():
        found = [r for r in recs if r.get("in_library")]
        if not found:
            print(f"  {t:<10s} none of its controls are in the library")
            continue
        med = float(np.median([r["percentile_in_library"] for r in found]))
        n_top10 = sum(1 for r in found if r["percentile_in_library"] >= 90)
        n_train = sum(1 for r in found if r["in_training"])
        recall_rows.append(
            {"target": t, "n_controls_found": len(found), "median_percentile": round(med, 1),
             "n_in_top_decile": n_top10, "n_already_in_training": n_train}
        )
        print(
            f"  {t:<10s} {len(found)} controls, median percentile {med:5.1f}, "
            f"{n_top10} in top decile, {n_train} were in training"
        )
    pd.DataFrame(recall_rows).to_csv(RESULTS / "control_recall_summary.csv", index=False)

    print(f"\nwrote panel_hits.csv, panel_target_summary.csv, compound_profiles.csv, "
          f"control_recall_summary.csv")


if __name__ == "__main__":
    main()
