"""Run the shipped pretrained ADMET-AI model on the scaffold-held-out molecules and compare to the
simple same-endpoint RF. ADMET-AI is pretrained on TDC (so in-distribution here) - this checks whether
the shipped multi-endpoint tool ranks correctly, for the 'ADMET risk radar' promise."""
import os, json, numpy as np, pandas as pd, pickle
from rdkit import Chem, RDLogger
from admet_ai import ADMETModel
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score
RDLogger.DisableLog("rdApp.*")
C = os.environ["CB_DATA"]; OUT = "results"
# map TDC endpoint -> ADMET-AI output column
MAP = {"Lipophilicity_AstraZeneca": ("Lipophilicity_AstraZeneca", "reg"), "Solubility_AqSolDB": ("Solubility_AqSolDB", "reg"),
       "Caco2_Wang": ("Caco2_Wang", "reg"), "hERG": ("hERG", "clf"), "BBB_Martins": ("BBB_Martins", "clf"),
       "CYP2D6_Veith": ("CYP2D6_Veith", "clf"), "Clearance_Hepatocyte_AZ": ("Clearance_Hepatocyte_AZ", "reg")}
model = ADMETModel()
S = {}
for name, (col, kind) in MAP.items():
    sp = pickle.load(open(f"{C}/{name}.pkl", "rb")); te = sp["test"]
    valid = [(s, y) for s, y in zip(te.Drug, te.Y) if Chem.MolFromSmiles(s) is not None]
    smis = [s for s, _ in valid]; y = np.array([v for _, v in valid])
    pred = model.predict(smiles=smis)
    if col not in pred.columns:
        S[name] = {"error": "column not found", "available_like": [c for c in pred.columns if col.split('_')[0] in c][:3]}; continue
    p = pred[col].values
    if kind == "reg":
        S[name] = {"n_test": int(len(y)), "spearman": round(float(spearmanr(y, p)[0]), 3), "pearson_r2": round(float(np.corrcoef(y, p)[0, 1] ** 2), 3)}
    else:
        S[name] = {"n_test": int(len(y)), "AUROC": round(float(roc_auc_score(y, p)), 3) if len(set(y)) == 2 else None}
    print(name, S[name])
json.dump(S, open(f"{OUT}/admet_ai.json", "w"), indent=1)
