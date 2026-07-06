#!/usr/bin/env python3
"""
Part 1b - DprE1 v2 Random-Forest activity model (open-source reproduction).

Paper (Table 1): DprE1 v2, Random Forest, 406 molecules, ROC(test) = 0.92.
Descriptors listed in the paper:
  ALogP, Molecular Weight, N & O count, #H-donors, #H-acceptors,
  #aromatic rings, #rotatable bonds, molecular surface area, molecular PSA,
  molecular polar SASA, molecular solubility, ECFP_4.

We map those to RDKit equivalents and add a 1024-bit ECFP4 (Morgan r=2)
fingerprint, then train a RandomForestClassifier and report ROC AUC averaged
over three iterations (as the paper does) using stratified cross-validation.
"""
import json
import numpy as np
from pathlib import Path
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, Lipinski, Crippen, rdMolDescriptors
from rdkit.Chem import rdFingerprintGenerator
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import roc_auc_score, roc_curve

RDLogger.DisableLog("rdApp.*")

BASE = Path(__file__).resolve().parent.parent
OUT = BASE / "results"
DATA = json.load(open(OUT / "dprE1_dataset.json"))

_MORGAN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)

DESCRIPTOR_NAMES = [
    "ALogP", "MolWt", "N_O_count", "HDonors", "HAcceptors",
    "AromaticRings", "RotatableBonds", "MolSurfaceArea", "MolPSA",
    "PolarSASA_approx", "MolSolubility_ESOL",
]


def esol_logS(mol):
    """Delaney ESOL estimate of aqueous solubility (log S)."""
    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    rb = Descriptors.NumRotatableBonds(mol)
    aromatic_atoms = sum(1 for a in mol.GetAtoms() if a.GetIsAromatic())
    ap = aromatic_atoms / mol.GetNumHeavyAtoms() if mol.GetNumHeavyAtoms() else 0
    return 0.16 - 0.63 * logp - 0.0062 * mw + 0.066 * rb - 0.74 * ap


def descriptors(mol):
    n_o = sum(1 for a in mol.GetAtoms() if a.GetSymbol() in ("N", "O"))
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    return [
        Crippen.MolLogP(mol),                 # ALogP ~ Crippen logP
        Descriptors.MolWt(mol),
        n_o,
        Lipinski.NumHDonors(mol),
        Lipinski.NumHAcceptors(mol),
        rdMolDescriptors.CalcNumAromaticRings(mol),
        Descriptors.NumRotatableBonds(mol),
        rdMolDescriptors.CalcLabuteASA(mol),  # molecular surface area
        tpsa,                                 # molecular PSA
        tpsa * Descriptors.FractionCSP3(mol), # crude polar-SASA proxy
        esol_logS(mol),                       # molecular solubility
    ]


def featurize(records):
    X, y = [], []
    for rec in records:
        mol = Chem.MolFromSmiles(rec["smiles"])
        if mol is None:
            continue
        desc = descriptors(mol)
        fp = _MORGAN.GetFingerprint(mol)
        fp_arr = np.zeros((1024,), dtype=np.int8)
        from rdkit.DataStructs import ConvertToNumpyArray
        ConvertToNumpyArray(fp, fp_arr)
        X.append(np.concatenate([desc, fp_arr]))
        y.append(rec["active"])
    return np.asarray(X, dtype=float), np.asarray(y)


def main():
    X, y = featurize(DATA)
    print(f"Feature matrix: {X.shape}  (11 descriptors + 1024-bit ECFP4)")
    print(f"Actives: {int(y.sum())} / {len(y)}")

    aucs = []
    oof_pred = None
    for i, seed in enumerate([1, 2, 3]):  # three iterations, as in the paper
        clf = RandomForestClassifier(
            n_estimators=500, max_features="sqrt",
            class_weight="balanced", random_state=seed, n_jobs=-1,
        )
        cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
        proba = cross_val_predict(clf, X, y, cv=cv, method="predict_proba",
                                  n_jobs=-1)[:, 1]
        auc = roc_auc_score(y, proba)
        aucs.append(auc)
        print(f"  iteration {i+1} (seed={seed}): 5-fold CV ROC AUC = {auc:.3f}")
        if oof_pred is None:
            oof_pred = proba

    print("-" * 52)
    print(f"Mean ROC AUC (3 iterations): {np.mean(aucs):.3f} +/- {np.std(aucs):.3f}")
    print(f"Paper DprE1 v2 reference:    0.92")

    # persist ROC curve points + a full-data model for downstream scoring
    fpr, tpr, _ = roc_curve(y, oof_pred)
    json.dump(
        {"mean_auc": float(np.mean(aucs)), "std_auc": float(np.std(aucs)),
         "aucs": [float(a) for a in aucs],
         "roc": {"fpr": fpr.tolist(), "tpr": tpr.tolist()},
         "n_molecules": int(len(y)), "n_active": int(y.sum()),
         "descriptor_names": DESCRIPTOR_NAMES},
        open(OUT / "activity_model_metrics.json", "w"), indent=2,
    )

    # train final model on all data and save via joblib for candidate scoring
    final = RandomForestClassifier(
        n_estimators=500, max_features="sqrt",
        class_weight="balanced", random_state=1, n_jobs=-1,
    ).fit(X, y)
    import joblib
    joblib.dump(final, OUT / "activity_rf_model.joblib")
    print(f"Saved metrics + final model to {OUT}")


if __name__ == "__main__":
    main()
