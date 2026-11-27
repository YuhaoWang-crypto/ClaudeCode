"""B. Learning curve: can a model learn SAR from very few client data points?
   C. Cross-endpoint transfer: does a model for one endpoint predict a different endpoint? (don't co-train)"""
import os, json, numpy as np, pandas as pd, pickle
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from sklearn.ensemble import RandomForestRegressor
from scipy.stats import spearmanr
RDLogger.DisableLog("rdApp.*")
C = os.environ["CB_DATA"]; OUT = "results"; rng = np.random.default_rng(0)
def fp(smi):
    m = Chem.MolFromSmiles(smi); return np.array(AllChem.GetMorganFingerprintAsBitVect(m, 2, 2048), np.int8) if m else None
def load(name):
    sp = pickle.load(open(f"{C}/{name}.pkl", "rb")); df = pd.concat([sp["train"], sp["valid"], sp["test"]], ignore_index=True)
    f = [fp(s) for s in df.Drug]; ok = [i for i, x in enumerate(f) if x is not None]
    return np.array([f[i] for i in ok]), df.Y.values[ok], df.Drug.values[ok]

# B. learning curve on Lipophilicity and Solubility (fix a held-out test, grow the training set)
S = {"learning_curve": {}}
for name in ["Lipophilicity_AstraZeneca", "Solubility_AqSolDB"]:
    X, y, _ = load(name); idx = rng.permutation(len(y)); te = idx[:600]; pool = idx[600:]
    curve = {}
    for n in [1, 3, 10, 30, 100, 300, 1000, len(pool)]:
        if n > len(pool): continue
        reps = []
        for r in range(5):
            tr = rng.choice(pool, n, replace=False)
            if len(set(y[tr])) < 2: reps.append(0.0); continue
            m = RandomForestRegressor(200, n_jobs=4, random_state=r).fit(X[tr], y[tr])
            reps.append(spearmanr(y[te], m.predict(X[te]))[0])
        curve[n] = round(float(np.nanmean(reps)), 3)
    S["learning_curve"][name] = curve
    print(name, "spearman by train-N:", curve)

# C. cross-endpoint transfer: train on endpoint A, predict endpoint B (molecules in both)
eps = ["Lipophilicity_AstraZeneca", "Solubility_AqSolDB", "Caco2_Wang", "Clearance_Hepatocyte_AZ"]
data = {e: load(e) for e in eps}
smi2y = {e: dict(zip(d[2], d[1])) for e, d in data.items()}
T = {}
for a in eps:
    Xa, ya, sa = data[a]; m = RandomForestRegressor(300, n_jobs=4, random_state=0).fit(Xa, ya)
    for b in eps:
        if a == b: continue
        Xb, yb, sb = data[b]; p = m.predict(Xb)
        T[f"{a[:10]}->{b[:10]}"] = round(float(spearmanr(yb, p)[0]), 3)
S["cross_endpoint_transfer_spearman"] = T
# also: pairwise correlation of the endpoints themselves on shared molecules
corr = {}
for i, a in enumerate(eps):
    for b in eps[i+1:]:
        sh = set(smi2y[a]) & set(smi2y[b])
        if len(sh) > 30:
            va = [smi2y[a][s] for s in sh]; vb = [smi2y[b][s] for s in sh]
            corr[f"{a[:10]}~{b[:10]}"] = {"n_shared": len(sh), "spearman": round(float(spearmanr(va, vb)[0]), 3)}
S["endpoint_label_correlation_on_shared_molecules"] = corr
json.dump(S, open(f"{OUT}/curve_transfer.json", "w"), indent=1)
print("\ntransfer:", json.dumps(T)); print("label corr:", json.dumps(corr))
