"""Fair head-to-head on the SAME TDC canonical scaffold test split, same training data:
   ADMET-AI (pretrained D-MPNN ensemble) vs a simple Morgan-fingerprint RF baseline."""
import os, json, numpy as np, pandas as pd, pickle
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score
RDLogger.DisableLog("rdApp.*")
C = os.environ["CB_DATA"]; OUT = "results"
def fp(s):
    m = Chem.MolFromSmiles(s); return np.array(AllChem.GetMorganFingerprintAsBitVect(m, 2, 2048), np.int8) if m else None
def feat(df):
    f = [fp(s) for s in df.Drug]; ok = [i for i, x in enumerate(f) if x is not None]
    return np.array([f[i] for i in ok]), df.Y.values[ok]
MAP = {"Lipophilicity_AstraZeneca": "reg", "Solubility_AqSolDB": "reg", "Caco2_Wang": "reg",
       "Clearance_Hepatocyte_AZ": "reg", "hERG": "clf", "BBB_Martins": "clf", "CYP2D6_Veith": "clf"}
ai = json.load(open(f"{OUT}/admet_ai.json"))
S = {}
for name, kind in MAP.items():
    sp = pickle.load(open(f"{C}/{name}.pkl", "rb"))
    Xtr, ytr = feat(pd.concat([sp["train"], sp["valid"]])); Xte, yte = feat(sp["test"])
    if kind == "reg":
        m = RandomForestRegressor(300, n_jobs=4, random_state=0).fit(Xtr, ytr); p = m.predict(Xte)
        rf = {"spearman": round(float(spearmanr(yte, p)[0]), 3)}; met = "spearman"
    else:
        m = RandomForestClassifier(300, n_jobs=4, random_state=0).fit(Xtr, ytr); p = m.predict_proba(Xte)[:, 1]
        rf = {"AUROC": round(float(roc_auc_score(yte, p)), 3)}; met = "AUROC"
    S[name] = {"metric": met, "RF_fingerprint": rf[met], "ADMET_AI": ai.get(name, {}).get(met), "n_test": int(len(yte))}
    print(f"{name:28s} {met:8s} RF={rf[met]}  ADMET-AI={S[name]['ADMET_AI']}")
json.dump(S, open(f"{OUT}/fair_compare.json", "w"), indent=1)
