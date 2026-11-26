"""Biophysical baseline: ViennaRNA structure + simple k-mer/uAUG features -> gradient boosting on MRL.
Trained on the random library (as the CNN was), evaluated on the random test set and on real human 5'UTRs."""
import os, numpy as np, pandas as pd, RNA
from sklearn.ensemble import HistGradientBoostingRegressor
N = os.environ["NB_DATA"]; OUT = "results"
def feats(seqs):
    rows = []
    for s in seqs:
        s = s.upper(); rna = s.replace("T", "U")
        mfe = RNA.fold(rna)[1]; mfe3 = RNA.fold(rna[-20:])[1]       # whole UTR and 3' 20nt (near start codon)
        gc = (s.count("G") + s.count("C")) / len(s)
        uaug = s.count("ATG"); kozak = 1.0 if s.endswith(("ACC", "GCC", "A", "G")) else 0.0
        comp = [s.count(b) / len(s) for b in "ACGT"]
        rows.append([mfe, mfe3, gc, uaug, len(s)] + comp)
    return np.array(rows, float)
FN = ["MFE_full", "MFE_3prime20", "GC", "uAUG_count", "len", "fA", "fC", "fG", "fT"]

rnd = pd.read_csv(f"{N}/GSM3130435_egfp_unmod_1.csv.gz").sort_values("total_reads", ascending=False).reset_index(drop=True).iloc[:280000]
rng = np.random.default_rng(0)
tr = rnd.iloc[20000:]; tr = tr.iloc[rng.choice(len(tr), 60000, replace=False)]      # 60k random-train for the baseline
Xtr = feats(tr.utr.values); ytr = tr.rl.values
gb = HistGradientBoostingRegressor(max_iter=400, learning_rate=0.08, max_depth=6).fit(Xtr, ytr)

for name, path in [("random_test", f"{OUT}/pred_random_test.csv"), ("designed", f"{OUT}/pred_designed.csv")]:
    d = pd.read_csv(path); X = feats(d.utr.values)
    d["pred_baseline"] = gb.predict(X)
    for i, f in enumerate(FN): d[f] = X[:, i]
    d.to_csv(path, index=False)
import json
json.dump({"baseline_feature_importance_note": "HistGBM", "features": FN,
           "n_baseline_train": int(len(tr))}, open(f"{OUT}/baseline_info.json", "w"))
print("baseline done; features", FN)
