#!/usr/bin/env python
"""Paired, per-cell-line comparison of ProteinTalks against a no-proteome control.

Supplementary Table S5C of the published paper reports per-cell-line
leave-one-cell-line-out results for every model it benchmarked, including
ProteinTalks itself (listed as ``ppODE (SWA)`` and ``ppODE (Non-SWA)``), over
repeated runs. That lets us do the comparison the paper does not: a *paired*
test, cell line by cell line, of ProteinTalks against the drug-mean predictor
computed from the paper's own label matrix (Table S1, sheet C_Efficacy).

The drug-mean predictor uses no proteomics, no drug chemistry and no training:
for a held-out cell line it predicts each drug's efficacy rate across the other
17 cell lines.
"""

from __future__ import annotations

import argparse
import json
import os

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import accuracy_score, average_precision_score, roc_auc_score

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def drug_mean_per_cellline(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for cl in sorted(df["Cell"].unique()):
        te = df[df["Cell"] == cl]
        tr = df[df["Cell"] != cl]
        if te["y"].nunique() < 2:
            continue
        rate = tr.groupby("Pert_ID")["y"].mean()
        p = te["Pert_ID"].map(rate).fillna(tr["y"].mean()).values
        rows.append({
            "Cell_line": cl,
            "AUROC": roc_auc_score(te["y"], p),
            "AUPRC": average_precision_score(te["y"], p),
            "Accuracy": accuracy_score(te["y"], (p >= 0.5).astype(int)),
            "n": len(te),
        })
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table-s1", required=True)
    ap.add_argument("--table-s5", required=True)
    ap.add_argument("--out", default=os.path.join(HERE, "results", "paired_cellline_test.json"))
    args = ap.parse_args()

    lab = pd.read_excel(args.table_s1, sheet_name="C_Efficacy")
    lab["y"] = (lab["Effective"].astype(str).str.strip().str.upper() == "Y").astype(int)
    ctrl = drug_mean_per_cellline(lab).set_index("Cell_line")

    s5 = pd.read_excel(args.table_s5, sheet_name="C")
    # Average the paper's repeated runs within each (cell line, model).
    agg = s5.groupby(["Cell_line", "Model"])[["Accuracy", "AUROC", "AUPRC"]].mean()

    results = {"per_model": {}, "control": ctrl.reset_index().to_dict("records")}
    print(f"{'model':22s} {'AUROC':>8s} {'ctrl':>8s} {'diff':>8s} {'p(paired)':>11s}  "
          f"{'wins/n':>8s}")

    for model in s5["Model"].unique():
        sub = agg.xs(model, level="Model")
        common = sorted(set(sub.index) & set(ctrl.index))
        if len(common) < 5:
            continue
        a = sub.loc[common, "AUROC"].values.astype(float)
        b = ctrl.loc[common, "AUROC"].values.astype(float)
        d = a - b
        t = stats.wilcoxon(a, b) if len(common) >= 6 else (np.nan, np.nan)
        pval = float(t[1]) if not isinstance(t, tuple) or not np.isnan(t[1]) else float("nan")
        try:
            pval = float(stats.wilcoxon(a, b).pvalue)
        except Exception:
            pval = float("nan")
        results["per_model"][model] = {
            "n_cell_lines": len(common),
            "model_auroc": float(a.mean()),
            "control_auroc": float(b.mean()),
            "mean_diff": float(d.mean()),
            "wilcoxon_p": pval,
            "model_wins": int((d > 0).sum()),
            "per_cell_line": {c: {"model": float(x), "control": float(y)}
                              for c, x, y in zip(common, a, b)},
        }
        print(f"{model:22s} {a.mean():8.3f} {b.mean():8.3f} {d.mean():+8.3f} "
              f"{pval:11.4g}  {int((d>0).sum()):4d}/{len(common):<4d}")

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n[saved] {args.out}")
    print("\nControl = drug-mean predictor: for a held-out cell line, predict each")
    print("drug's efficacy rate across the other 17 cell lines. No proteomics used.")


if __name__ == "__main__":
    main()
