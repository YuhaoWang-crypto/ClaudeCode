"""Which scoring upgrade actually buys enrichment?

Measures, on the 1,320 experimentally labelled designs of the public release, how much
each proposed upgrade changes the ability to separate binders from non-binders:

  1. metric:    ipTM  ->  ipSAE_min
  2. ensemble:  one predictor -> z-scored ensemble of several
  3. self-consistency: adding sc-DockQ at quarter weight
  4. generator: hit rate by structure-generation method (descriptive, not controlled)

Reported as within-target average precision (chance = that target's hit rate) and as
the hit rate among the top-10 ranked designs of each target, which is what a wet-lab
budget actually feels.

Usage:
    python score_ablation.py --release-dir <protein_binder_design_data_release>
"""

import argparse
import itertools
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


def average_precision(scores, labels):
    order = np.argsort(-np.asarray(scores, dtype=float), kind="mergesort")
    y = np.asarray(labels, dtype=int)[order]
    if y.sum() == 0:
        return np.nan
    cum = np.cumsum(y)
    prec = cum / np.arange(1, len(y) + 1)
    return float((prec * y).sum() / y.sum())


def topk_hit_rate(scores, labels, k=10):
    order = np.argsort(-np.asarray(scores, dtype=float), kind="mergesort")
    y = np.asarray(labels, dtype=int)[order][:k]
    return float(y.mean())


def evaluate(df, score_col, targets, k=10):
    """Within-target AP and top-k hit rate, averaged over targets."""
    aps, tops, chance = [], [], []
    for t in targets:
        sub = df[df["target"] == t]
        s = sub[score_col].to_numpy(dtype=float)
        y = sub["binder_final"].astype(int).to_numpy()
        if np.isnan(s).all():
            continue
        s = np.nan_to_num(s, nan=np.nanmin(s))
        aps.append(average_precision(s, y))
        tops.append(topk_hit_rate(s, y, k))
        chance.append(y.mean())
    return float(np.mean(aps)), float(np.mean(tops)), float(np.mean(chance))


def zscore_within_target(df, cols, targets):
    """Per-target z-score of each column, then the row mean across columns."""
    out = pd.Series(index=df.index, dtype=float)
    for t in targets:
        m = df["target"] == t
        sub = df.loc[m, cols]
        z = (sub - sub.mean()) / sub.std(ddof=0).replace(0, np.nan)
        out.loc[m] = z.mean(axis=1)
    return out.fillna(out.mean())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--release-dir", required=True, type=Path)
    ap.add_argument("--top-k", type=int, default=10)
    args = ap.parse_args()

    t = args.release_dir / "tables"
    design = pd.read_parquet(t / "design_summary.parquet")
    cofold = pd.read_parquet(t / "insilico" / "cofold_predictions.parquet")

    labelled = design[design["binder_final"].notna()][["uuid", "target", "binder_final", "generator"]]

    # per design x predictor: best seed by ipSAE_min, and the same design's best ipTM
    agg = (cofold.groupby(["uuid", "cofolding_model"])
                 .agg(ipsae=("ipsae_min", "max"), iptm=("iptm_pae", "max"),
                      scdockq=("sc_dockq", "max"))
                 .reset_index())
    wide = agg.pivot(index="uuid", columns="cofolding_model",
                     values=["ipsae", "iptm", "scdockq"])
    wide.columns = [f"{a}_{b}" for a, b in wide.columns]
    df = labelled.merge(wide, on="uuid", how="left")

    targets = [x for x in df["target"].unique()
               if df.loc[df["target"] == x, "binder_final"].sum() >= 3
               and (~df.loc[df["target"] == x, "binder_final"].astype(bool)).sum() >= 3]
    print(f"{len(df)} labelled designs, {len(targets)} evaluable targets "
          f"(>=3 binders and >=3 non-binders)\n")

    PREDICTORS = ["boltz2", "chai1", "ptxv2", "ef2full", "ef2fast", "of3", "rf3", "afm3", "odde", "af3of3"]
    rows = []

    # --- 1. single predictor: ipTM vs ipSAE_min
    for p in PREDICTORS:
        for metric in ("iptm", "ipsae"):
            col = f"{metric}_{p}"
            if col not in df:
                continue
            a, k, c = evaluate(df, col, targets, args.top_k)
            rows.append(dict(scheme=f"single {p}", metric=metric.upper(), ap=a,
                             top_k=k, chance=c))

    # --- 2. ensembles (z-scored within target)
    ENSEMBLES = {
        "campaign trio (ef2full+ef2fast+ptxv2)": ["ef2full", "ef2fast", "ptxv2"],
        "boltz2+chai1+ptxv2": ["boltz2", "chai1", "ptxv2"],
        "boltz2+chai1": ["boltz2", "chai1"],
        "all 10 predictors": PREDICTORS,
    }
    for name, members in ENSEMBLES.items():
        for metric in ("iptm", "ipsae"):
            cols = [f"{metric}_{m}" for m in members if f"{metric}_{m}" in df]
            df["_ens"] = zscore_within_target(df, cols, targets)
            a, k, c = evaluate(df, "_ens", targets, args.top_k)
            rows.append(dict(scheme=name, metric=metric.upper(), ap=a, top_k=k, chance=c))

    # --- 3. + sc-DockQ at quarter weight (the paper's ranking score)
    for name, members in ENSEMBLES.items():
        cols = [f"ipsae_{m}" for m in members if f"ipsae_{m}" in df]
        dq = [f"scdockq_{m}" for m in members if f"scdockq_{m}" in df]
        z_ipsae = zscore_within_target(df, cols, targets)
        z_dq = zscore_within_target(df, dq, targets)
        df["_ens"] = (z_ipsae + 0.25 * z_dq) / 1.25
        a, k, c = evaluate(df, "_ens", targets, args.top_k)
        rows.append(dict(scheme=name + " + sc-DockQ(1/4)", metric="IPSAE", ap=a,
                         top_k=k, chance=c))

    res = pd.DataFrame(rows).sort_values("ap", ascending=False)
    res["ap"] = res["ap"].round(3)
    res["top_k"] = res["top_k"].round(3)
    res["chance"] = res["chance"].round(3)
    print(res.to_string(index=False))

    # --- paired comparison: ipSAE_min vs ipTM on the same predictor, per target
    print("\nPer-target paired delta, ipSAE_min - ipTM (same predictor, same designs):")
    deltas = {}
    for p in PREDICTORS:
        d = []
        for tg in targets:
            sub = df[df["target"] == tg]
            y = sub["binder_final"].astype(int).to_numpy()
            a_i = average_precision(np.nan_to_num(sub[f"iptm_{p}"].to_numpy(float), nan=0), y)
            a_s = average_precision(np.nan_to_num(sub[f"ipsae_{p}"].to_numpy(float), nan=0), y)
            d.append(a_s - a_i)
        deltas[p] = d
        w = stats.wilcoxon(d).pvalue if len(set(np.round(d, 9))) > 1 else float("nan")
        print(f"  {p:<8} mean {np.mean(d):+.3f}   wins {sum(x > 0 for x in d)}/{len(d)}   Wilcoxon p={w:.3g}")
    allд = np.concatenate(list(deltas.values()))
    print(f"  pooled over predictors: mean {allд.mean():+.3f}, "
          f"{(allд > 0).sum()}/{len(allд)} target-predictor pairs improve")

    # --- 4. generator: hit rate (descriptive)
    print("\nHit rate by structure-generation method (paper's caveat: Claude chose the "
          "method per target, so this is not a controlled comparison):")
    m = design[design["binder_final"].notna() & (design["vendor_agreement"] != "not_tested_either")]
    g = m.groupby("generator")["binder_final"].agg(["sum", "count"])
    g = g[g["count"] >= 100]
    g["hit_rate"] = (100 * g["sum"] / g["count"]).round(1)
    print(g.sort_values("hit_rate", ascending=False).to_string())

    # what the same designs look like once the best scoring scheme reorders them
    best = res.iloc[0]
    print(f"\nBest scheme: {best['scheme']} ({best['metric']}) — "
          f"AP {best['ap']} vs chance {best['chance']}, "
          f"top-{args.top_k} hit rate {best['top_k']:.0%} vs {best['chance']:.0%} overall")


if __name__ == "__main__":
    main()
