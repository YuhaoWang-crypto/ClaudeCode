#!/usr/bin/env python
"""The control experiment the paper does not report, run on the paper's own labels.

Input
-----
Supplementary Table S1, sheet ``C_Efficacy`` of the published Nature version:
1,116 rows of (cell line, perturbation ID, Effective Y/N), 18 cell lines x 62
perturbations, 373 "Y" / 743 "N". This is the exact ground truth of the drug
efficacy task, and it is openly downloadable.

Question
--------
How much of the reported performance is obtainable from the *label matrix
alone*, with no proteomics, no drug chemistry and no model?

Three predictors, all of which ignore the proteome entirely:

  drug-mean       P(effective) = rate for this drug across training cell lines
  cell-line-mean  P(effective) = rate for this cell line across training drugs
  additive        logistic regression on one-hot(drug) + one-hot(cell line)

evaluated under each of the paper's three splitting protocols:

  Setting 1  random 0.7/0.2/0.1 over conditions
  Setting 2  leave-one-cell-line-out   -> paper reports AUROC 0.953 +- 0.002
  Setting 3  leave-one-drug-out        -> paper reports AUROC 0.638 +- 0.018

The comparison for each setting is against the paper's own reported number for
that same setting, taken from Supplementary Tables S5A and S6A.
"""

from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, average_precision_score, roc_auc_score

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)

# Paper's reported ProteinTalks numbers, from the published supplementary tables.
PAPER = {
    "setting1": {"auroc": 0.960, "auprc": 0.854, "accuracy": 0.910,
                 "source": "main text (preprint v1 Results)"},
    "setting2": {"auroc": 0.953, "auprc": 0.891, "accuracy": 0.902,
                 "source": "Supplementary Table S5A (mean +- SE over folds)"},
    "setting3": {"auroc": 0.638, "auprc": 0.573, "accuracy": 0.775,
                 "source": "Supplementary Table S6A (mean +- SE over folds)"},
}
# Best non-ProteinTalks comparator reported by the paper in each setting.
PAPER_BEST_BASELINE = {
    "setting2": {"model": "KNN", "auroc": 0.805, "auprc": 0.631, "accuracy": 0.817},
    "setting3": {"model": "Linear regression", "auroc": 0.669, "auprc": 0.566,
                 "accuracy": 0.542},
}


def metrics(y, p):
    y = np.asarray(y).astype(int)
    p = np.asarray(p, dtype=float)
    out = {"n": int(len(y)), "positive_rate": float(y.mean())}
    if len(np.unique(y)) > 1:
        out["auroc"] = float(roc_auc_score(y, p))
        out["auprc"] = float(average_precision_score(y, p))
    else:
        out["auroc"] = float("nan")
        out["auprc"] = float("nan")
    out["accuracy"] = float(accuracy_score(y, (p >= 0.5).astype(int)))
    return out


def group_rate_predictor(df_tr, df_te, key: str):
    """Training-set positive rate of ``key``; global rate for unseen levels."""
    g = df_tr.groupby(key)["y"].mean()
    glob = float(df_tr["y"].mean())
    return df_te[key].map(g).fillna(glob).values


def additive_predictor(df_tr, df_te):
    """Logistic regression on one-hot drug + one-hot cell line. No proteome."""
    cats = pd.concat([df_tr, df_te])
    X = pd.get_dummies(cats[["Cell", "Pert_ID"]], columns=["Cell", "Pert_ID"])
    Xtr = X.iloc[: len(df_tr)].values.astype(float)
    Xte = X.iloc[len(df_tr) :].values.astype(float)
    clf = LogisticRegression(max_iter=2000, C=1.0)
    clf.fit(Xtr, df_tr["y"].values)
    return clf.predict_proba(Xte)[:, 1]


def run_setting1(df, n_rep: int, seed: int):
    rows = []
    rng = np.random.default_rng(seed)
    acc = {k: [] for k in ["drug-mean", "cell-line-mean", "additive"]}
    for _ in range(n_rep):
        perm = rng.permutation(len(df))
        n_tr = int(0.9 * len(df))          # 0.7 train + 0.2 val, both visible
        tr, te = df.iloc[perm[:n_tr]], df.iloc[perm[n_tr:]]
        acc["drug-mean"].append(metrics(te["y"], group_rate_predictor(tr, te, "Pert_ID")))
        acc["cell-line-mean"].append(metrics(te["y"], group_rate_predictor(tr, te, "Cell")))
        acc["additive"].append(metrics(te["y"], additive_predictor(tr, te)))
    for name, ms in acc.items():
        rows.append(_agg(name, ms))
    return rows


def run_leave_one_out(df, key: str, other: str):
    """Hold out each level of ``key``; the informative predictor is on ``other``."""
    acc = {k: [] for k in ["drug-mean", "cell-line-mean", "additive"]}
    for level in sorted(df[key].unique()):
        te = df[df[key] == level]
        tr = df[df[key] != level]
        if te["y"].nunique() < 2:
            continue
        acc["drug-mean"].append(metrics(te["y"], group_rate_predictor(tr, te, "Pert_ID")))
        acc["cell-line-mean"].append(metrics(te["y"], group_rate_predictor(tr, te, "Cell")))
        acc["additive"].append(metrics(te["y"], additive_predictor(tr, te)))
    return [_agg(name, ms) for name, ms in acc.items()]


def _agg(name, ms):
    out = {"predictor": name, "n_folds": len(ms)}
    for m in ["auroc", "auprc", "accuracy"]:
        v = np.array([x[m] for x in ms], dtype=float)
        v = v[~np.isnan(v)]
        out[m] = float(v.mean()) if len(v) else float("nan")
        out[f"{m}_se"] = float(v.std(ddof=1) / np.sqrt(len(v))) if len(v) > 1 else 0.0
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table-s1", required=True, help="TableS1_drug_cell_*.xlsx")
    ap.add_argument("--n-rep", type=int, default=200)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default=os.path.join(HERE, "results", "real_label_controls.json"))
    args = ap.parse_args()

    df = pd.read_excel(args.table_s1, sheet_name="C_Efficacy")
    df = df.rename(columns=str.strip)
    df["y"] = (df["Effective"].astype(str).str.strip().str.upper() == "Y").astype(int)
    print(f"[data] {len(df)} conditions | {df['Cell'].nunique()} cell lines | "
          f"{df['Pert_ID'].nunique()} perturbations | positives {df['y'].sum()} "
          f"({df['y'].mean():.3f})")

    out = {"n_conditions": int(len(df)), "paper_reported": PAPER}

    print("\n=== Setting 1: random split (paper: AUROC 0.960, AUPRC 0.854, ACC 0.910) ===")
    out["setting1"] = run_setting1(df, args.n_rep, args.seed)
    _show(out["setting1"], PAPER["setting1"])

    print("\n=== Setting 2: leave-one-cell-line-out "
          "(paper: AUROC 0.953, AUPRC 0.891, ACC 0.902) ===")
    out["setting2"] = run_leave_one_out(df, key="Cell", other="Pert_ID")
    _show(out["setting2"], PAPER["setting2"], PAPER_BEST_BASELINE["setting2"])

    print("\n=== Setting 3: leave-one-drug-out "
          "(paper: AUROC 0.638, AUPRC 0.573, ACC 0.775) ===")
    out["setting3"] = run_leave_one_out(df, key="Pert_ID", other="Cell")
    _show(out["setting3"], PAPER["setting3"], PAPER_BEST_BASELINE["setting3"])

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\n[saved] {args.out}")


def _show(rows, paper, best_baseline=None):
    print(f"  {'predictor':22s} {'AUROC':>16s} {'AUPRC':>16s} {'accuracy':>16s}")
    for r in rows:
        print(f"  {r['predictor']:22s} "
              f"{r['auroc']:.3f}+-{r['auroc_se']:.3f}   "
              f"{r['auprc']:.3f}+-{r['auprc_se']:.3f}   "
              f"{r['accuracy']:.3f}+-{r['accuracy_se']:.3f}   (n={r['n_folds']})")
    print(f"  {'ProteinTalks (paper)':22s} {paper['auroc']:.3f}"
          f"{'':12s}{paper['auprc']:.3f}{'':12s}{paper['accuracy']:.3f}")
    if best_baseline:
        print(f"  {'best baseline (paper)':22s} {best_baseline['auroc']:.3f}"
              f"{'':12s}{best_baseline['auprc']:.3f}{'':12s}"
              f"{best_baseline['accuracy']:.3f}   [{best_baseline['model']}]")


if __name__ == "__main__":
    main()
