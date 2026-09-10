import os, sys, json, time
import numpy as np, pandas as pd
from scipy.stats import spearmanr
import evalkit as E

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = HERE + "/results"; os.makedirs(OUT, exist_ok=True)

DATASETS = ["NA_influenza", "PA_influenza", "MIC_H_influenzae", "MIC_E_cloacae", "hERG"]


def load(name):
    d = pd.read_csv(f"{HERE}/curated/{name}.csv")
    X, ok = E.fps(d.smiles.tolist())
    d = d[ok].reset_index(drop=True); X = X[ok]
    return d, X, d.p.values.astype(np.float64)


def evaluate(name):
    t0 = time.time()
    d, X, y = load(name)
    kind = "ridge" if len(y) < 400 else "auto"
    print(f"\n=== {name}: n={len(y)} scaffolds={d.scaffold.nunique()} model={kind}", flush=True)

    clus = E.sphere_exclusion(X, sim_cut=0.4)
    print(f"    sphere-exclusion clusters (Tanimoto<0.4 between leaders): {len(set(clus))}", flush=True)

    splits = {
        "random": E.random_folds(len(y)),
        "scaffold": E.group_folds(d.scaffold.values),
        "cluster0.4": E.group_folds(clus),
    }
    ts = E.temporal_split(d.year.values)
    if ts:
        splits["temporal"] = ts

    rows, oofs = [], {}
    for sname, folds in splits.items():
        for qr in (False, True):
            oof, nn, ok = E.run_cv(X, y, folds, use_qrasar=qr, kind=kind)
            m = E.metrics(y[ok], oof[ok])
            m.update(dataset=name, split=sname, qrasar=qr,
                     n_train_frac=float(np.mean([len(tr) for tr, _ in folds]) / len(y)))
            rows.append(m)
            oofs[(sname, qr)] = (oof, nn, ok)
            print(f"    {sname:11s} qRASAR={int(qr)}  rho={m['spearman']:.3f} "
                  f"RMSE={m['rmse']:.3f} R2={m['r2']:.3f} EF5%={m['ef5']:.1f} "
                  f"<=1log={m['within1log']:.2f}", flush=True)

    # applicability-domain curve on the strictest available split
    key = ("cluster0.4", True)
    oof, nn, ok = oofs[key]
    ad = []
    edges = [0, .2, .3, .4, .5, .6, .7, 1.01]
    for lo, hi in zip(edges[:-1], edges[1:]):
        s = ok & (nn >= lo) & (nn < hi)
        if s.sum() >= 40:
            r = spearmanr(y[s], oof[s]).statistic
            ad.append(dict(dataset=name, lo=lo, hi=hi, n=int(s.sum()), spearman=float(r),
                           rmse=float(np.sqrt(np.mean((y[s] - oof[s]) ** 2)))))
    for a in ad:
        print(f"    AD nn[{a['lo']:.1f},{a['hi']:.1f}) n={a['n']:5d} rho={a['spearman']:.3f} "
              f"RMSE={a['rmse']:.3f}", flush=True)

    np.savez(f"{OUT}/{name}_oof.npz", y=y, oof=oof, nn=nn, ok=ok,
             smiles=np.array(d.smiles), clus=clus)
    print(f"    [{time.time()-t0:.0f}s]", flush=True)
    return rows, ad


if __name__ == "__main__":
    targets = sys.argv[1:] or DATASETS
    R, A = [], []
    for n in targets:
        if not os.path.exists(f"{HERE}/curated/{n}.csv"):
            print(f"skip {n} (not curated)"); continue
        r, a = evaluate(n); R += r; A += a
        pd.DataFrame(R).to_csv(f"{OUT}/splits_{'_'.join(targets)[:40]}.csv", index=False)
        pd.DataFrame(A).to_csv(f"{OUT}/ad_{'_'.join(targets)[:40]}.csv", index=False)
    print("\nDONE")
