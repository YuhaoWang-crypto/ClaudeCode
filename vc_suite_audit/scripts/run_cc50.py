"""CC50 pan-cytotoxicity: the suite's largest q-RASAR claim (+0.169 Spearman).

Uses the deployed configuration (RF, 300 trees, min_samples_leaf=10) so the comparison
is like-for-like, and tests it under scaffold, cluster and temporal splits.
"""
import os, time, json
import numpy as np, pandas as pd
from scipy.stats import spearmanr
from sklearn.ensemble import RandomForestRegressor
import evalkit as E

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = HERE + "/results"


def model():
    return RandomForestRegressor(n_estimators=int(os.environ.get("N_TREES", 300)),
                                 min_samples_leaf=10, n_jobs=4, random_state=E.SEED)


def run(folds, X, y, qr):
    oof = np.full(len(y), np.nan); nn = np.full(len(y), np.nan)
    for tr, te in folds:
        T_te = E.tanimoto_block(X[te], X[tr]); nn[te] = T_te.max(1)
        if qr:
            T_tr = E.tanimoto_block(X[tr], X[tr])
            Xtr = np.hstack([X[tr], E.qrasar_feats(T_tr, y[tr], loo=True)])
            Xte = np.hstack([X[te], E.qrasar_feats(T_te, y[tr], loo=False)])
            del T_tr
        else:
            Xtr, Xte = X[tr], X[te]
        m = model(); m.fit(Xtr, y[tr]); oof[te] = m.predict(Xte)
        del T_te
    ok = ~np.isnan(oof)
    return oof, nn, ok


if __name__ == "__main__":
    d = pd.read_csv(f"{HERE}/curated/CC50_pan.csv")
    X, ok = E.fps(d.smiles.tolist())
    d = d[ok].reset_index(drop=True); X = X[ok]; y = d.p.values
    print(f"CC50_pan n={len(y)} scaffolds={d.scaffold.nunique()}", flush=True)

    t0 = time.time()
    clus = E.sphere_exclusion(X, 0.4)
    print(f"clusters={len(set(clus))} [{time.time()-t0:.0f}s]", flush=True)

    splits = {"scaffold": E.group_folds(d.scaffold.values),
              "cluster0.4": E.group_folds(clus)}
    ts = E.temporal_split(d.year.values)
    if ts:
        splits["temporal"] = ts

    rows = []
    done = set()
    prev = f"{OUT}/splits_CC50.csv"
    if os.path.exists(prev):
        p = pd.read_csv(prev)
        rows = json.loads(p.to_json(orient="records"))
        done = {(r["split"], bool(r["qrasar"])) for r in rows}
        print(f"resuming; already done: {sorted(done)}", flush=True)

    for sname, folds in splits.items():
        for qr in (False, True):
            if (sname, qr) in done:
                continue
            t = time.time()
            oof, nn, k = run(folds, X, y, qr)
            m = E.metrics(y[k], oof[k]); m.update(dataset="CC50_pan", split=sname, qrasar=qr)
            rows.append(m)
            print(f"  {sname:11s} qRASAR={int(qr)} rho={m['spearman']:.3f} RMSE={m['rmse']:.3f} "
                  f"R2={m['r2']:.3f} EF5%={m['ef5']:.1f} [{time.time()-t:.0f}s]", flush=True)
            pd.DataFrame(rows).to_csv(f"{OUT}/splits_CC50.csv", index=False)
            if sname == "cluster0.4" and qr:
                ad = []
                for lo, hi in zip([0, .2, .3, .4, .5, .6, .7], [.2, .3, .4, .5, .6, .7, 1.01]):
                    s = k & (nn >= lo) & (nn < hi)
                    if s.sum() >= 40:
                        ad.append(dict(dataset="CC50_pan", lo=lo, hi=hi, n=int(s.sum()),
                                       spearman=float(spearmanr(y[s], oof[s]).statistic),
                                       rmse=float(np.sqrt(np.mean((y[s] - oof[s]) ** 2)))))
                pd.DataFrame(ad).to_csv(f"{OUT}/ad_CC50.csv", index=False)
                for a in ad:
                    print(f"    AD nn[{a['lo']:.1f},{a['hi']:.1f}) n={a['n']:6d} "
                          f"rho={a['spearman']:.3f} RMSE={a['rmse']:.3f}", flush=True)
    print("DONE_CC50")
