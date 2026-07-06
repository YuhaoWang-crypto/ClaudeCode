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
