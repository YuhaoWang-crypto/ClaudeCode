"""Transferability tests: does a model trained for one target/species carry to another?

Two claims in the report are tested here:
  - "cross-target selectivity holds" (NA config vs PA config rank compounds differently)
  - 13 separate per-species MIC models are worth having (vs one pooled antibacterial model)
"""
import os, numpy as np, pandas as pd
from scipy.stats import spearmanr
import evalkit as E

HERE = os.path.dirname(os.path.abspath(__file__))


def load(name):
    d = pd.read_csv(f"{HERE}/curated/{name}.csv")
    X, ok = E.fps(d.smiles.tolist())
    return d[ok].reset_index(drop=True), X[ok]


def fit_full(X, y):
    m = E.make_model(len(y))
    m.fit(X, y)
    return m


def transfer(src, dst):
    ds, Xs = load(src); dd, Xd = load(dst)
    shared = set(ds.smiles) & set(dd.smiles)
    mask_new = ~dd.smiles.isin(shared).values          # dst compounds unseen by src model
    m = fit_full(Xs, ds.p.values)
    pred = m.predict(Xd)
    nn = E.tanimoto_block(Xd, Xs).max(1)

    out = dict(src=src, dst=dst, n_dst=len(dd), n_shared=len(shared),
               rho_all=float(spearmanr(dd.p.values, pred).statistic),
               rho_unseen=float(spearmanr(dd.p.values[mask_new], pred[mask_new]).statistic),
               n_unseen=int(mask_new.sum()),
               rho_unseen_indomain=np.nan, n_unseen_indomain=0)
    sel = mask_new & (nn >= 0.4)
    if sel.sum() > 30:
        out["rho_unseen_indomain"] = float(spearmanr(dd.p.values[sel], pred[sel]).statistic)
        out["n_unseen_indomain"] = int(sel.sum())

    # intrinsic agreement ceiling: for compounds measured against BOTH endpoints,
    # how well does one endpoint's TRUE value predict the other's?
    if shared:
        a = ds.set_index("smiles").p
        b = dd.set_index("smiles").p
        both = pd.concat([a[list(shared)], b[list(shared)]], axis=1, keys=["src", "dst"]).dropna()
        if len(both) > 20:
            out["rho_truth_vs_truth"] = float(spearmanr(both.src, both.dst).statistic)
            out["n_truth_pairs"] = len(both)
    return out


if __name__ == "__main__":
    pairs = [("NA_influenza", "PA_influenza"), ("PA_influenza", "NA_influenza"),
             ("MIC_H_influenzae", "MIC_E_cloacae"), ("MIC_E_cloacae", "MIC_H_influenzae"),
             ("hERG", "MIC_E_cloacae")]
    rows = []
    for s, d in pairs:
        if not (os.path.exists(f"{HERE}/curated/{s}.csv") and os.path.exists(f"{HERE}/curated/{d}.csv")):
            continue
        r = transfer(s, d); rows.append(r)
        print(f"{s:18s} -> {d:18s} rho(all)={r['rho_all']:+.3f} "
              f"rho(unseen n={r['n_unseen']})={r['rho_unseen']:+.3f} "
              f"rho(unseen,in-domain n={r['n_unseen_indomain']})={r['rho_unseen_indomain']:+.3f} "
              f"| truth-vs-truth={r.get('rho_truth_vs_truth', float('nan')):+.3f}"
              f" (n={r.get('n_truth_pairs', 0)})", flush=True)
    pd.DataFrame(rows).to_csv(f"{HERE}/results/transfer.csv", index=False)
