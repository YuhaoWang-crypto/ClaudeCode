"""Stage 1 - build the CDK2 pIC50 label set from ChEMBL.

Design decisions that matter for honesty of the downstream benchmark:

* Binding assays only (assay_type=B) and standard_relation '=' - censored
  values ('>', '<') carry no usable regression signal and would otherwise
  pile up at the assay detection limit.
* pChEMBL is used as the label. ChEMBL only computes it for '=' / nM / >0,
  so it is already the cleaned-up potency axis.
* Rows carrying a data_validity_comment (unit errors, out-of-range, suspected
  author error) are dropped rather than trusted.
* Duplicate measurements of the same compound are aggregated by median, and
  compounds whose replicate measurements span more than DUP_SPREAD_MAX log
  units are dropped - they are assay noise, not labels.
* Structures are desalted to the largest organic fragment and canonicalised
  before de-duplication, so salt forms collapse onto one label.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import pandas as pd
import requests
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")

API = "https://www.ebi.ac.uk/chembl/api/data/activity.json"
TARGET = "CHEMBL301"  # Cyclin-dependent kinase 2, Homo sapiens (UniProt P24941)
ACTIVITY_TYPES = ("IC50", "Ki")
PAGE = 1000

# Cleaning thresholds
DUP_SPREAD_MAX = 2.0  # log units between replicate measurements of one compound
MW_RANGE = (150.0, 700.0)
ALLOWED_ELEMENTS = set("H B C N O F Si P S Cl Se Br I".split())

OUT = Path(__file__).resolve().parents[1] / "data"


def fetch_activities(activity_type: str) -> list[dict]:
    """Page through the ChEMBL activity endpoint for one activity type."""
    params = {
        "target_chembl_id": TARGET,
        "standard_type": activity_type,
        "assay_type": "B",
        "standard_relation": "=",
        "pchembl_value__isnull": "false",
        "limit": PAGE,
        "offset": 0,
        "format": "json",
    }
    rows: list[dict] = []
    while True:
        for attempt in range(4):
            try:
                r = requests.get(API, params=params, timeout=120)
                r.raise_for_status()
                break
            except Exception as exc:  # network flake - back off and retry
                if attempt == 3:
                    raise
                print(f"    retry {attempt + 1} after {exc}")
                time.sleep(2**attempt)
        payload = r.json()
        batch = payload.get("activities", [])
        rows.extend(batch)
        total = payload.get("page_meta", {}).get("total_count", 0)
        print(f"  {activity_type}: {len(rows)}/{total}")
        if len(batch) < PAGE or len(rows) >= total:
            break
        params["offset"] += PAGE
    return rows


def standardise(smiles: str) -> str | None:
    """Desalt to the largest fragment, neutralise, return canonical SMILES."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = rdMolStandardize.Cleanup(mol)
    mol = rdMolStandardize.LargestFragmentChooser().choose(mol)
    if mol is None:
        return None
    mol = rdMolStandardize.Uncharger().uncharge(mol)
    if mol.GetNumHeavyAtoms() == 0:
        return None
    if any(a.GetSymbol() not in ALLOWED_ELEMENTS for a in mol.GetAtoms()):
        return None
    mw = Descriptors.MolWt(mol)
    if not MW_RANGE[0] <= mw <= MW_RANGE[1]:
        return None
    return Chem.MolToSmiles(mol)


def scaffold_of(smiles: str) -> str:
    """Bemis-Murcko scaffold, used later for the scaffold split."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ""
    try:
        core = MurckoScaffold.GetScaffoldForMol(mol)
        smi = Chem.MolToSmiles(core)
    except Exception:
        return ""
    # Acyclic molecules give an empty scaffold; bucket them by themselves so
    # they cannot leak across the split as one giant pseudo-scaffold.
    return smi if smi else f"ACYCLIC::{smiles}"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    raw: list[dict] = []
    for atype in ACTIVITY_TYPES:
        raw.extend(fetch_activities(atype))
    df = pd.DataFrame(raw)
    print(f"\nraw activity rows: {len(df)}")
    df.to_csv(OUT / "chembl_raw_activities.csv.gz", index=False)

    keep = [
        "molecule_chembl_id",
        "canonical_smiles",
        "pchembl_value",
        "standard_type",
        "standard_value",
        "standard_units",
        "assay_chembl_id",
        "document_chembl_id",
        "data_validity_comment",
        "target_chembl_id",
    ]
    df = df[[c for c in keep if c in df.columns]].copy()

    n0 = len(df)
    df = df[df["data_validity_comment"].isna()]
    print(f"drop flagged data_validity_comment: {n0} -> {len(df)}")

    df = df.dropna(subset=["canonical_smiles", "pchembl_value"])
    df["pchembl_value"] = df["pchembl_value"].astype(float)

    # ChEMBL caps pChEMBL sanity at ~1 mM .. 1 pM; anything outside is noise.
    df = df[(df["pchembl_value"] >= 3.0) & (df["pchembl_value"] <= 12.0)]
    print(f"after pChEMBL range filter: {len(df)}")

    print("standardising structures ...")
    df["smiles"] = [standardise(s) for s in df["canonical_smiles"]]
    df = df.dropna(subset=["smiles"])
    print(f"after structure standardisation: {len(df)}")

    # Aggregate replicate measurements per standardised structure.
    grp = df.groupby("smiles")["pchembl_value"]
    agg = pd.DataFrame(
        {
            "pchembl": grp.median(),
            "n_meas": grp.size(),
            "spread": grp.max() - grp.min(),
        }
    ).reset_index()

    n1 = len(agg)
    agg = agg[agg["spread"] <= DUP_SPREAD_MAX]
    print(f"drop compounds with replicate spread > {DUP_SPREAD_MAX}: {n1} -> {len(agg)}")

    # Carry one representative ChEMBL id per structure for traceability.
    rep = df.drop_duplicates("smiles").set_index("smiles")["molecule_chembl_id"]
    agg["molecule_chembl_id"] = agg["smiles"].map(rep)
    agg["scaffold"] = [scaffold_of(s) for s in agg["smiles"]]

    agg = agg.sort_values("pchembl", ascending=False).reset_index(drop=True)
    agg.to_csv(OUT / "cdk2_dataset.csv", index=False)

    summary = {
        "target": TARGET,
        "activity_types": list(ACTIVITY_TYPES),
        "n_raw_rows": int(n0),
        "n_compounds": int(len(agg)),
        "n_scaffolds": int(agg["scaffold"].nunique()),
        "pchembl_min": float(agg["pchembl"].min()),
        "pchembl_max": float(agg["pchembl"].max()),
        "pchembl_mean": float(agg["pchembl"].mean()),
        "pchembl_std": float(agg["pchembl"].std()),
        "largest_scaffold_frac": float(
            agg["scaffold"].value_counts().iloc[0] / len(agg)
        ),
    }
    (OUT / "dataset_summary.json").write_text(json.dumps(summary, indent=2))
    print("\n" + json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
