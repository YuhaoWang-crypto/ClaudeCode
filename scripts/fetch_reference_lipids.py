"""Build data/reference_lipids.csv by pulling canonical SMILES from PubChem.

We deliberately do NOT hard-code SMILES for these complex lipids (easy to get
subtly wrong). Instead we resolve them by name against the PubChem PUG-REST API,
so every structure is authoritative and reproducible.

Run:  python scripts/fetch_reference_lipids.py
"""

from __future__ import annotations

import csv
import sys
import time
from pathlib import Path

import requests

# (name, role) for reference LNP components + known ionizable lipids.
REFERENCE = [
    ("DLin-MC3-DMA", "ionizable_lipid"),   # MC3, Onpattro
    ("SM-102", "ionizable_lipid"),         # Moderna COVID vaccine
    ("ALC-0315", "ionizable_lipid"),       # Pfizer/BioNTech COVID vaccine
    ("DLin-KC2-DMA", "ionizable_lipid"),
    ("C12-200", "ionizable_lipid"),
    ("cholesterol", "sterol"),
    ("DSPC", "helper_phospholipid"),
    ("DOPE", "helper_phospholipid"),
]

# NOTE: PubChem PUG-REST renamed its SMILES properties (2025): the isomeric
# structure is now 'SMILES' and the flat one is 'ConnectivitySMILES'. The old
# 'CanonicalSMILES'/'IsomericSMILES' names return empty.
PUG = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{}"
       "/property/SMILES,ConnectivitySMILES,MolecularFormula,MolecularWeight/JSON")


def resolve(name: str) -> dict | None:
    url = PUG.format(requests.utils.quote(name))
    try:
        r = requests.get(url, timeout=30)
    except requests.RequestException as e:
        print(f"  {name}: request error {e}", file=sys.stderr)
        return None
    if r.status_code != 200:
        print(f"  {name}: HTTP {r.status_code}", file=sys.stderr)
        return None
    props = r.json().get("PropertyTable", {}).get("Properties", [])
    return props[0] if props else None


def main() -> int:
    out_path = Path(__file__).resolve().parents[1] / "data" / "reference_lipids.csv"
    out_path.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    for name, role in REFERENCE:
        p = resolve(name)
        if not p:
            print(f"skip {name}")
            continue
        smiles = p.get("SMILES") or p.get("ConnectivitySMILES", "")
        rows.append({
            "name": name,
            "role": role,
            "smiles": smiles,
            "formula": p.get("MolecularFormula", ""),
            "mw": p.get("MolecularWeight", ""),
            "cid": p.get("CID", ""),
        })
        print(f"ok   {name:16s} {p.get('MolecularFormula','')}")
        time.sleep(0.25)  # be polite to PubChem

    with out_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["name", "role", "smiles", "formula", "mw", "cid"])
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {len(rows)} rows -> {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
