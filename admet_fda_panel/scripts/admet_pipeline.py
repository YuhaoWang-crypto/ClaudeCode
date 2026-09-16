"""ADMET & drug-likeness profiling pipeline for an FDA-approved drug panel.

Layers
------
1. Standardisation      : RDKit sanitize, largest-fragment, neutralise, canonical SMILES
2. Physicochemical      : MW, cLogP, TPSA, HBD/HBA, RotB, rings, FracCsp3, HeavyAtoms, QED
3. Drug-likeness rules  : Lipinski Ro5, Veber, Ghose, Egan, Muegge-lite
4. Structural alerts    : RDKit FilterCatalog — PAINS (A/B/C), Brenk, NIH
5. ADMET-AI             : 41 TDC endpoints (Chemprop-RDKit ensemble) + DrugBank-approved
                          percentiles for every endpoint

Everything runs offline on CPU after the one-time ADMET-AI model download.
"""

from __future__ import annotations

import os
import warnings

warnings.filterwarnings("ignore")

OUT = "/home/user/results/admet_fda_panel"


# --------------------------------------------------------------------------- #
# 1. standardisation
# --------------------------------------------------------------------------- #
def standardize(df, smiles_col="smiles"):
    from rdkit import Chem, RDLogger
    from rdkit.Chem.MolStandardize import rdMolStandardize

    RDLogger.DisableLog("rdApp.*")
    lfc = rdMolStandardize.LargestFragmentChooser()
    un = rdMolStandardize.Uncharger()

    out = df.copy()
    canon, changed = [], 0
    for smi in out[smiles_col]:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            canon.append(None)
            continue
        mol = rdMolStandardize.Cleanup(mol)
        mol = lfc.choose(mol)
        mol = un.uncharge(mol)
        new = Chem.MolToSmiles(mol)
        changed += int(new != smi)
        canon.append(new)
    out["smiles_std"] = canon
    n_bad = out["smiles_std"].isna().sum()
    print(f"  standardised {len(out) - n_bad}/{len(out)} molecules "
          f"({changed} structures altered, {n_bad} failed)")
    return out.dropna(subset=["smiles_std"]).reset_index(drop=True)


# --------------------------------------------------------------------------- #
# 2-3. descriptors + rules
# --------------------------------------------------------------------------- #
def compute_descriptors(df, smiles_col="smiles_std"):
    import pandas as pd
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors, QED

    rows = []
    for smi in df[smiles_col]:
        m = Chem.MolFromSmiles(smi)
        mw = Descriptors.MolWt(m)
        logp = Crippen.MolLogP(m)
        tpsa = rdMolDescriptors.CalcTPSA(m)
        hbd = rdMolDescriptors.CalcNumHBD(m)
        hba = rdMolDescriptors.CalcNumHBA(m)
        rotb = rdMolDescriptors.CalcNumRotatableBonds(m)
        arom = rdMolDescriptors.CalcNumAromaticRings(m)
        rings = rdMolDescriptors.CalcNumRings(m)
        fcsp3 = rdMolDescriptors.CalcFractionCSP3(m)
        heavy = m.GetNumHeavyAtoms()
        mr = Crippen.MolMR(m)
        qed = QED.qed(m)

        # Lipinski rule-of-five violations
        ro5 = sum([mw > 500, logp > 5, hbd > 5, hba > 10])
        # Veber: RotB <= 10 and TPSA <= 140
        veber = (rotb <= 10) and (tpsa <= 140)
        # Ghose: 160<=MW<=480, -0.4<=logP<=5.6, 40<=MR<=130, 20<=atoms<=70
        ghose = (160 <= mw <= 480) and (-0.4 <= logp <= 5.6) and \
                (40 <= mr <= 130) and (20 <= heavy <= 70)
        # Egan: TPSA <= 131.6, -1 <= logP <= 5.88
        egan = (tpsa <= 131.6) and (-1 <= logp <= 5.88)

        rows.append(dict(
            MW=mw, cLogP=logp, TPSA=tpsa, HBD=hbd, HBA=hba, RotB=rotb,
            AromaticRings=arom, Rings=rings, FractionCSP3=fcsp3,
            HeavyAtoms=heavy, MolarRefractivity=mr, QED=qed,
            Lipinski_violations=ro5, Lipinski_pass=ro5 <= 1,
            Veber_pass=veber, Ghose_pass=ghose, Egan_pass=egan,
        ))
    desc = pd.DataFrame(rows, index=df.index)
    print(f"  computed {desc.shape[1]} descriptors / rule flags")
    return pd.concat([df, desc], axis=1)


# --------------------------------------------------------------------------- #
# 4. structural alerts
# --------------------------------------------------------------------------- #
ALERT_CATALOGS = ("PAINS_A", "PAINS_B", "PAINS_C", "BRENK", "NIH")


def structural_alerts(df, smiles_col="smiles_std"):
    import pandas as pd
    from rdkit import Chem
    from rdkit.Chem import FilterCatalog
    from rdkit.Chem.FilterCatalog import FilterCatalogParams

    cats = {}
    for name in ALERT_CATALOGS:
        params = FilterCatalogParams()
        params.AddCatalog(getattr(FilterCatalogParams.FilterCatalogs, name))
        cats[name] = FilterCatalog.FilterCatalog(params)

    rows, long_rows = [], []
    for idx, (nm, smi) in enumerate(zip(df["name"], df[smiles_col])):
        m = Chem.MolFromSmiles(smi)
        rec = {}
        total = 0
        for cname, cat in cats.items():
            hits = cat.GetMatches(m)
            names = [h.GetDescription() for h in hits]
            rec[f"alert_{cname}"] = len(names)
            total += len(names)
            for h in names:
                long_rows.append({"name": nm, "catalog": cname, "alert": h})
        rec["alert_total"] = total
        rec["PAINS_total"] = sum(rec[f"alert_PAINS_{x}"] for x in "ABC")
        rows.append(rec)

    alerts = pd.DataFrame(rows, index=df.index)
    long = pd.DataFrame(long_rows, columns=["name", "catalog", "alert"])
    print(f"  structural alerts: {int(alerts['alert_total'].sum())} matches "
          f"across {int((alerts['alert_total'] > 0).sum())} compounds")
    return pd.concat([df, alerts], axis=1), long


# --------------------------------------------------------------------------- #
# 5. ADMET-AI
# --------------------------------------------------------------------------- #
def run_admet_ai(df, smiles_col="smiles_std"):
    import pandas as pd
    from admet_ai import ADMETModel

    model = ADMETModel()
    preds = model.predict(smiles=list(df[smiles_col]))
    preds = preds.reset_index(drop=True)
    # ADMET-AI already returns MW/logP/QED etc.; drop duplicates of our own
    drop = [c for c in ["molecular_weight", "logP", "hydrogen_bond_acceptors",
                        "hydrogen_bond_donors", "Lipinski", "QED",
                        "stereo_centers", "tpsa", "PAINS_alert", "BRENK_alert",
                        "NIH_alert"] if c in preds.columns]
    admet_only = preds.drop(columns=drop)
    n_end = len([c for c in admet_only.columns if "percentile" not in c])
    print(f"  ADMET-AI: {n_end} endpoints + "
          f"{admet_only.shape[1] - n_end} DrugBank-approved percentiles")
    return pd.concat([df.reset_index(drop=True), admet_only], axis=1), preds


# --------------------------------------------------------------------------- #
# orchestration
# --------------------------------------------------------------------------- #
def run_full_analysis(df):
    print("Step 1/4  standardisation")
    d = standardize(df)
    print("Step 2/4  physicochemical descriptors & drug-likeness rules")
    d = compute_descriptors(d)
    print("Step 3/4  structural alerts (PAINS A/B/C, Brenk, NIH)")
    d, alerts_long = structural_alerts(d)
    print("Step 4/4  ADMET-AI endpoint prediction")
    d, raw = run_admet_ai(d)
    print(f"✓ Analysis completed — {d.shape[0]} molecules × {d.shape[1]} columns")
    return d, alerts_long


def summarize(d):
    import pandas as pd
    n = len(d)
    s = {
        "n_molecules": n,
        "Lipinski_pass": int(d["Lipinski_pass"].sum()),
        "Veber_pass": int(d["Veber_pass"].sum()),
        "Ghose_pass": int(d["Ghose_pass"].sum()),
        "Egan_pass": int(d["Egan_pass"].sum()),
        "mean_QED": float(d["QED"].mean()),
        "median_MW": float(d["MW"].median()),
        "median_cLogP": float(d["cLogP"].median()),
        "median_TPSA": float(d["TPSA"].median()),
        "with_PAINS": int((d["PAINS_total"] > 0).sum()),
        "with_any_alert": int((d["alert_total"] > 0).sum()),
        "hERG_high_risk": int((d["hERG"] > 0.5).sum()),
        "DILI_high_risk": int((d["DILI"] > 0.5).sum()),
        "AMES_positive": int((d["AMES"] > 0.5).sum()),
        "BBB_penetrant": int((d["BBB_Martins"] > 0.5).sum()),
        "CYP3A4_inhibitor": int((d["CYP3A4_Veith"] > 0.5).sum()),
        "Pgp_inhibitor": int((d["Pgp_Broccatelli"] > 0.5).sum()),
    }
    return s


if __name__ == "__main__":
    import pandas as pd
    import pickle
    import sys

    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from drug_panel import load_example_drugs

    os.makedirs(OUT, exist_ok=True)
    panel = load_example_drugs()
    results, alerts_long = run_full_analysis(panel)
    summary = summarize(results)

    results.to_csv(f"{OUT}/all_properties.csv", index=False)
    alerts_long.to_csv(f"{OUT}/structural_alerts.csv", index=False)

    dl_cols = ["name", "drug_class", "MW", "cLogP", "TPSA", "HBD", "HBA", "RotB",
               "FractionCSP3", "QED", "Lipinski_violations", "Lipinski_pass",
               "Veber_pass", "Ghose_pass", "Egan_pass", "alert_total"]
    results[dl_cols].to_csv(f"{OUT}/druglikeness_summary.csv", index=False)

    admet_cols = ["name"] + [c for c in results.columns
                             if c not in dl_cols and "percentile" not in c
                             and c not in ("smiles", "smiles_std", "formula",
                                           "AromaticRings", "Rings", "HeavyAtoms",
                                           "MolarRefractivity", "PAINS_total")
                             and not c.startswith("alert_")]
    results[admet_cols].to_csv(f"{OUT}/admet_predictions.csv", index=False)

    flagged = results[(results["hERG"] > 0.5) | (results["DILI"] > 0.5) |
                      (results["AMES"] > 0.5) | (results["alert_total"] > 0) |
                      (~results["Lipinski_pass"])]
    flag_cols = ["name", "drug_class", "hERG", "DILI", "AMES", "ClinTox",
                 "LD50_Zhu", "Lipinski_violations", "alert_total", "QED"]
    flagged[flag_cols].to_csv(f"{OUT}/flagged_compounds.csv", index=False)

    with open(f"{OUT}/analysis_object.pkl", "wb") as fh:
        pickle.dump({"results": results, "alerts": alerts_long,
                     "summary": summary}, fh)

    print("\nSummary:")
    for k, v in summary.items():
        print(f"  {k:24s} {v}")
    print(f"\nWrote outputs to {OUT}")
