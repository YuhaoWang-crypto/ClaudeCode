"""QSAR / ADMET known-answer + applicability-domain test on 8 TDC endpoints.
Morgan-fingerprint models (RF/GBM). Compare random split vs scaffold split (new chemotypes) -
the matrix-row boundaries: 'check applicability domain' and 'one client point can't learn new SAR'."""
import os, json, numpy as np, pandas as pd, pickle
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds.MurckoScaffold import MurckoScaffoldSmiles
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score
RDLogger.DisableLog("rdApp.*")
C = os.environ["CB_DATA"]; OUT = "results"; rng = np.random.default_rng(0)
ENDPOINTS = {"Caco2_Wang": "reg", "Lipophilicity_AstraZeneca": "reg", "Solubility_AqSolDB": "reg",
             "Half_Life_Obach": "reg", "Clearance_Hepatocyte_AZ": "reg",
             "hERG": "clf", "CYP2D6_Veith": "clf", "BBB_Martins": "clf"}
def fp(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None: return None
    return np.array(AllChem.GetMorganFingerprintAsBitVect(m, 2, 2048), np.int8)
def scaffold(smi):
    try: return MurckoScaffoldSmiles(smi)
    except Exception: return smi
def featurize(df):
    f = [fp(s) for s in df.Drug]; ok = [i for i, x in enumerate(f) if x is not None]
    return np.array([f[i] for i in ok]), df.Y.values[ok], df.Drug.values[ok]

def fit_eval(Xtr, ytr, Xte, yte, kind):
    if kind == "reg":
        m = RandomForestRegressor(300, n_jobs=4, random_state=0).fit(Xtr, ytr); p = m.predict(Xte)
        return {"spearman": round(float(spearmanr(yte, p)[0]), 3),
                "pearson_r2": round(float(np.corrcoef(yte, p)[0, 1] ** 2), 3)}
    else:
        m = RandomForestClassifier(300, n_jobs=4, random_state=0).fit(Xtr, ytr); p = m.predict_proba(Xte)[:, 1]
        return {"AUROC": round(float(roc_auc_score(yte, p)), 3) if len(set(yte)) == 2 else None}

S = {}
for name, kind in ENDPOINTS.items():
    sp = pickle.load(open(f"{C}/{name}.pkl", "rb"))
    full = pd.concat([sp["train"], sp["valid"], sp["test"]], ignore_index=True)
    X, y, smi = featurize(full)
    scaf = np.array([scaffold(s) for s in smi])
    # scaffold split: hold out whole scaffolds (unseen chemotypes)
    uniq = pd.Series(scaf).value_counts().index.tolist(); rng.shuffle(uniq)
    cut = set(uniq[:int(0.2 * len(uniq))]) if len(uniq) > 5 else set(uniq[-1:])
    te_sc = np.array([s in cut for s in scaf]); 
    # random split of the same size
    idx = rng.permutation(len(y)); te_rand = np.zeros(len(y), bool); te_rand[idx[:int(te_sc.mean()*len(y))]] = True
    res = {"n": int(len(y)), "n_scaffolds": int(len(set(scaf))), "test_frac": round(float(te_sc.mean()), 2)}
    res["random_split"] = fit_eval(X[~te_rand], y[~te_rand], X[te_rand], y[te_rand], kind)
    res["scaffold_split"] = fit_eval(X[~te_sc], y[~te_sc], X[te_sc], y[te_sc], kind)
    S[name] = res
    met = "AUROC" if kind == "clf" else "spearman"
    print(f"{name:28s} {kind} random {met}={res['random_split'][met]}  scaffold {met}={res['scaffold_split'][met]}")
json.dump(S, open(f"{OUT}/qsar_splits.json", "w"), indent=1)
