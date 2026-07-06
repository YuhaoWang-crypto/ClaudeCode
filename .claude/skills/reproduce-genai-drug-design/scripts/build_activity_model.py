# Template: dataset prep + RF activity model (Parts 1a+1b).
# Adapt paths/thresholds/descriptors to the target paper. Two stages below,
# split them into separate files in a real project (see dprE1_reproduction/src).
# ============================ STAGE 1: prepare dataset ============================
#!/usr/bin/env python3
"""
Part 1a - Dataset preparation for the DprE1 activity model.

Reproduces the training set described in Chikhale et al. (chemRxiv 2026,
"Generative AI and Structure-Based Workflow ... DprE1 Inhibitor Candidates").

The Zenodo release ships:
  * List_of_DprE1_Inhibitors_with_activity_data.xlsx  - IUPAC + activity
    (mixed "MIC/IC50/IC90" column, NO SMILES)
  * Structures_of_ligands.zip - the drawn structures as ChemDraw .cdx, one
    folder per DOI, converted here to SMILES with OpenBabel (see convert_cdx.py)

Paper (DprE1 v2 model): 406 molecules with reported IC50, converted to pIC50,
active if pIC50 >= 5.75 (192 actives).

Empirically the IC50 entries are exactly those reported in molar units (uM/nM);
MIC entries are reported in ug/mL. We therefore keep molar entries as the IC50
set and label active at pIC50 >= 5.75.

Structure resolution priority:
  1. OpenBabel SMILES from the matching .cdx file (the authors' own drawing)
  2. OPSIN SMILES from the IUPAC name (fallback)
"""
import re
import json
import math
import openpyxl
from pathlib import Path
from py2opsin import py2opsin
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

BASE = Path(__file__).resolve().parent.parent
DATA = BASE / "data"
OUT = BASE / "results"
OUT.mkdir(exist_ok=True)

XLSX = DATA / "List_of_DprE1_Inhibitors_with_activity_data.xlsx"
CDX_SMILES = json.load(open(DATA / "cdx_smiles.json"))
ACTIVE_PIC50 = 5.75

MOLAR_UNITS = {"µm": 1e-6, "um": 1e-6, "μm": 1e-6, "nm": 1e-9, "m": 1.0}


def parse_activity(raw):
    """Return (IC50_in_M, censored) for molar (IC50) entries, else None."""
    if raw is None:
        return None
    s = str(raw).strip()
    if not s:
        return None
    censored = ">" if ">" in s else ("<" if "<" in s else "")
    m = re.search(r"([-+]?\d*\.?\d+)", s.replace(",", ""))
    if not m:
        return None
    val = float(m.group(1))
    unit = "".join(re.findall(r"[a-zA-Zµμ/]+", s)).lower()
    factor = None
    for u, f in MOLAR_UNITS.items():
        if unit == u or unit.replace("mic", "") == u:
            factor = f
            break
    if factor is None or val <= 0:
        return None
    return val * factor, censored


# ---- CDX structure lookup by (doi-folder, label) -----------------------------
# group cdx keys by folder for flexible label matching
_folders = {}
for key, smi in CDX_SMILES.items():
    folder, stem = key.split("||", 1)
    _folders.setdefault(folder, {})[stem] = smi


def cdx_lookup(doi, ligand_code):
    if not doi:
        return None
    folder = str(doi).replace("/", "")
    stems = _folders.get(folder)
    if not stems:
        # some folders drop a leading registrant chunk; try suffix match
        cand = [f for f in _folders if folder.endswith(f) or f.endswith(folder)]
        if len(cand) == 1:
            stems = _folders[cand[0]]
        else:
            return None
    label = str(ligand_code).split("_")[-1]
    # exact stem == full ligand_code, or stem == label, or trailing-token match
    if ligand_code in stems:
        return stems[ligand_code]
    if label in stems:
        return stems[label]
    for stem, smi in stems.items():
        if stem.split("_")[-1] == label:
            return smi
    return None


def main():
    wb = openpyxl.load_workbook(XLSX)
    rows = list(wb.active.iter_rows(values_only=True))[1:]

    # forward-fill DOI (merged cells -> only first row of each group is filled)
    records = []
    cur_doi = None
    for r in rows:
        doi = r[1]
        if doi:
            cur_doi = doi
        ligand_code, iupac, act = r[2], r[5], r[6]
        parsed = parse_activity(act)
        if parsed is None:
            continue
        ic50_M, censored = parsed
        records.append(dict(
            ligand_code=ligand_code, doi=cur_doi, iupac=(str(iupac).strip() if iupac else ""),
            raw_activity=str(act), ic50_M=ic50_M, censored=censored,
            pIC50=round(-math.log10(ic50_M), 3),
        ))
    print(f"Molar (IC50) entries parsed:        {len(records)}")

    # structure resolution: CDX first, OPSIN fallback
    need_opsin = []
    for rec in records:
        smi = cdx_lookup(rec["doi"], rec["ligand_code"])
        if smi:
            rec["smiles"] = smi
            rec["source"] = "cdx"
        else:
            rec["smiles"] = None
            need_opsin.append(rec)

    n_cdx = sum(1 for r in records if r["smiles"])
    print(f"Resolved from CDX (OpenBabel):      {n_cdx}")

    opsin_names = [r["iupac"] for r in need_opsin]
    if opsin_names:
        smi_list = py2opsin(opsin_names, output_format="SMILES")
        if isinstance(smi_list, str):
            smi_list = smi_list.splitlines()
        for rec, smi in zip(need_opsin, smi_list):
            smi = (smi or "").strip()
            if smi and Chem.MolFromSmiles(smi):
                rec["smiles"] = Chem.MolToSmiles(Chem.MolFromSmiles(smi))
                rec["source"] = "opsin"
    n_opsin = sum(1 for r in records if r.get("source") == "opsin")
    print(f"Resolved from IUPAC (OPSIN):        {n_opsin}")

    # validate + canonicalize
    valid = []
    for rec in records:
        if not rec["smiles"]:
            continue
        mol = Chem.MolFromSmiles(rec["smiles"])
        if mol is None:
            continue
        rec["smiles"] = Chem.MolToSmiles(mol)
        rec["active"] = int(rec["pIC50"] >= ACTIVE_PIC50)
        valid.append(rec)

    # dedup on canonical SMILES -> keep most potent
    best = {}
    for rec in valid:
        k = rec["smiles"]
        if k not in best or rec["pIC50"] > best[k]["pIC50"]:
            best[k] = rec
    final = list(best.values())

    n_active = sum(r["active"] for r in final)
    print("-" * 52)
    print(f"Unique molecules (dedup):           {len(final)}")
    print(f"  actives (pIC50 >= 5.75):          {n_active}")
    print(f"  inactives:                        {len(final) - n_active}")
    print(f"Paper reference (DprE1 v2):         406 molecules, 192 actives")

    json.dump(final, open(OUT / "dprE1_dataset.json", "w"), indent=2)
    print(f"Wrote {OUT / 'dprE1_dataset.json'}")


if __name__ == "__main__":
    main()

# ============================ STAGE 2: build RF model ============================
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
