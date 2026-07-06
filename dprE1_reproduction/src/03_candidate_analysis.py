#!/usr/bin/env python3
"""
Part 1c - Candidate analysis: reproduce paper Table 2 for TCA1 + GTD_9.1..9.10.

SMILES are transcribed from the preprint Appendix-2 (ChemDraw 22.2).
We recompute MolWt, ALogP(~Crippen), MolPSA(~TPSA) and #rotatable bonds and
compare against the paper's Discovery-Studio values, then score each candidate
with the reproduced DprE1 v2 RF model (Part 1b).

Note: ALogP and MolPSA in the paper come from Discovery Studio, whose algorithms
differ from RDKit's Crippen logP / TPSA, so small offsets are expected. MolWt and
rotatable-bond counts are algorithm-independent and should match closely.
"""
import json
import numpy as np
from pathlib import Path
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdMolDescriptors, Crippen
import joblib
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module
mdl = import_module("02_build_activity_model")

RDLogger.DisableLog("rdApp.*")
BASE = Path(__file__).resolve().parent.parent
OUT = BASE / "results"

# --- SMILES from Appendix-2 ----------------------------------------------------
CANDIDATES = {
    "TCA1":    "O=C(OCC)NC(C1=C(NC(C2=NC3=CC=CC=C3S2)=O)SC=C1)=O",
    "GTD_9.1": "O=C(NCC(F)F)CCC1=C(C2=CC(CCNC(C)=O)=CC(N)=C2C)C3=C(O[C@@]4([H])[C@]3([H])CCCC4)C=N1",
    "GTD_9.2": "CN1[C@@H](C[C@@H]2OCCC[C@@H]12)CNC(C(C)(c3cc(c(cn4)cc5c4occ5)c(NCC(F)F)c(N)c3)C)=O",
    "GTD_9.3": "Cc1c(C2=C(N3CCC(OCCC[C@@H]4C[C@H]5CCC[C@@H]5N4)CC3)COC2=O)nc(C(NCC(F)F)=O)cc1N",
    "GTD_9.4": "CCC(Cc1cc(c2c3c(O[C@@H]4CC[C@@H](C[C@@H]43)C)cnc2CCC(NCC(F)F)=O)c(C)c(CN)c1)=O",
    "GTD_9.5": "Cc(c(C(NCC(F)F)=O)nc1)c2c1[nH]cc2CCO[C@@H]3CCN(C4=Cc5c(CC4)c6c(CCO6)cc5N)C3",
    "GTD_9.6": "CC1=CC(N2CC[C@H](NC[C@@H]3[C@@H](C(C)(c4cc(CC(F)F)c(CO)c(N)c4)C)Cn5c3ncc5)C2)=NC1",
    "GTD_9.7": "CC(CC1CCOCC1)(c2cc(C3=NC(C(N4CC[C@H](C4)CNCc5ccccc5)=C3)=O)c(F)c(N)c2)C",
    "GTD_9.8": "Cc1nc2c(CNC2)c(C[C@@](C(NCC(F)F)=O)(c3cc(CC[C@H]4Cc5c(O4)cccc5)c(C)c(N)c3)C)n1",
    "GTD_9.9": "CN1[C@@H](C[C@H]2CCCC[C@H]12)CNC(C(C)(c3cc(c4ncc5c(CCO5)c4)c(CCC(F)F)c(N)c3)C)=O",
    "GTD_9.10": "Cc1c(C(NCC(F)F)=O)cc(C2=C(N3CCC(OCCC[C@H]4Cc5c(O4)cccc5)CC3)COC2=O)nc1N",
}

# --- paper Table 2 reference values (Discovery Studio) -------------------------
PAPER = {  # name: (ALogP, MolPSA, MolWt, RotBonds, OverallDS, DprE1v2)
    "TCA1":    (2.88, 154, 375, 5, 0.114, 0.148),
    "GTD_9.1": (2.97, 106, 501, 9, 1.00, 0.637),
    "GTD_9.2": (2.83, 106, 528, 8, 0.99, 0.630),
    "GTD_9.3": (2.23, 119, 548, 10, 0.983, 0.661),
    "GTD_9.4": (4.06, 94.3, 514, 10, 0.977, 0.538),
    "GTD_9.5": (3.31, 106, 538, 8, 0.963, 0.670),
    "GTD_9.6": (1.76, 91.7, 513, 9, 0.943, 0.634),
    "GTD_9.7": (4.37, 80, 519, 9, 0.941, 0.627),
    "GTD_9.8": (3.83, 102, 536, 9, 0.934, 0.601),
    "GTD_9.9": (4.96, 80.5, 527, 8, 0.927, 0.623),
    "GTD_9.10": (3.13, 116, 557, 10, 0.924, 0.569),
}


def main():
    model = joblib.load(OUT / "activity_rf_model.joblib")
    rows = []
    print(f"{'Name':9} | {'MolWt (calc/paper)':22} | {'RotB (c/p)':10} | "
          f"{'ALogP (c/p)':14} | {'TPSA/MolPSA':14} | {'RF P(active)':12}")
    print("-" * 95)
    for name, smi in CANDIDATES.items():
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            print(f"{name:9} | INVALID SMILES")
            continue
        mw = Descriptors.MolWt(mol)
        rb = Descriptors.NumRotatableBonds(mol)
        alogp = Crippen.MolLogP(mol)
        tpsa = rdMolDescriptors.CalcTPSA(mol)

        # featurize exactly as the model expects and predict
        X, _ = mdl.featurize([{"smiles": Chem.MolToSmiles(mol), "active": 0}])
        p_active = float(model.predict_proba(X)[0, 1])

        p = PAPER[name]
        rows.append(dict(name=name, smiles=Chem.MolToSmiles(mol),
                         molwt_calc=round(mw, 1), molwt_paper=p[2],
                         rotb_calc=rb, rotb_paper=p[3],
                         alogp_calc=round(alogp, 2), alogp_paper=p[0],
                         tpsa_calc=round(tpsa, 1), molpsa_paper=p[1],
                         rf_p_active=round(p_active, 3),
                         paper_dprE1v2=p[5]))
        print(f"{name:9} | {mw:7.1f} / {p[2]:<12} | {rb:2d} / {p[3]:<5} | "
              f"{alogp:5.2f} / {p[0]:<6} | {tpsa:5.1f} / {p[1]:<6} | {p_active:.3f}")

    # agreement stats on algorithm-independent props
    mw_err = np.mean([abs(r["molwt_calc"] - r["molwt_paper"]) for r in rows])
    rb_match = np.mean([r["rotb_calc"] == r["rotb_paper"] for r in rows])
    print("-" * 95)
    print(f"Mean |MolWt error| vs paper: {mw_err:.1f} Da")
    print(f"Rotatable-bond exact-match rate: {rb_match*100:.0f}%")
    print(f"All candidates predicted active by RF: "
          f"{all(r['rf_p_active'] >= 0.5 for r in rows if r['name'].startswith('GTD'))}")

    json.dump(rows, open(OUT / "candidate_analysis.json", "w"), indent=2)
    print(f"Wrote {OUT / 'candidate_analysis.json'}")


if __name__ == "__main__":
    main()
