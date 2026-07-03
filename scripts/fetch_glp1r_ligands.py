"""Build a GLP1R (Problem-B) ligand dataset directly from the ChEMBL REST API.

GLP1R = CHEMBL1784 (UniProt P43220). We pull all potent activities
(pChEMBL >= threshold), attach canonical SMILES, and classify each ligand as a
small molecule vs. peptide/large by molecular weight — the two feed different
targeting-ligand tracks (small-molecule conjugate vs. peptide conjugate).

This uses the public ChEMBL API with pagination, so it is reproducible outside
this session and not subject to any MCP output-size cap.

Run:  python scripts/fetch_glp1r_ligands.py [--target CHEMBL1784] [--min-pchembl 6]
"""

from __future__ import annotations

import argparse
import csv
import sys
import time
from pathlib import Path

import requests
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors

RDLogger.DisableLog("rdApp.*")

BASE = "https://www.ebi.ac.uk/chembl/api/data/activity.json"
SM_MW_CUTOFF = 900.0  # below -> small-molecule conjugate candidate


def fetch(target: str, min_pchembl: float) -> list[dict]:
    rows, offset, limit = [], 0, 1000
    while True:
        params = {
            "target_chembl_id": target,
            "pchembl_value__gte": min_pchembl,
            "limit": limit,
            "offset": offset,
        }
        r = requests.get(BASE, params=params, timeout=60)
        r.raise_for_status()
        data = r.json()
        acts = data.get("activities", [])
        rows.extend(acts)
        total = data.get("page_meta", {}).get("total_count", len(rows))
        offset += limit
        print(f"  fetched {len(rows)}/{total}", file=sys.stderr)
        if offset >= total or not acts:
            break
        time.sleep(0.2)
    return rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", default="CHEMBL1784")
    ap.add_argument("--name", default="GLP1R")
    ap.add_argument("--min-pchembl", type=float, default=6.0)
    args = ap.parse_args()

    out_dir = Path(__file__).resolve().parents[1] / "data" / "targets" / args.name
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"querying ChEMBL for {args.name} ({args.target}), pChEMBL >= {args.min_pchembl}")
    acts = fetch(args.target, args.min_pchembl)

    # dedupe by molecule, keep max pchembl
    best: dict[str, dict] = {}
    for a in acts:
        smi = a.get("canonical_smiles")
        cid = a.get("molecule_chembl_id")
        pch = a.get("pchembl_value")
        if not smi or not cid or pch is None:
            continue
        pch = float(pch)
        if cid not in best or pch > best[cid]["pchembl"]:
            best[cid] = {"molecule_chembl_id": cid, "smiles": smi,
                         "standard_type": a.get("standard_type"), "pchembl": pch}

    rows = []
    for r in best.values():
        m = Chem.MolFromSmiles(r["smiles"])
        mw = Descriptors.MolWt(m) if m else None
        r["mw"] = round(mw, 1) if mw else ""
        r["class"] = "small_molecule" if (mw and mw < SM_MW_CUTOFF) else "peptide_or_large"
        rows.append(r)

    rows.sort(key=lambda x: x["pchembl"], reverse=True)
    out = out_dir / f"{args.name}_ligands.csv"
    with out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["molecule_chembl_id", "class", "pchembl",
                                          "standard_type", "mw", "smiles"])
        w.writeheader()
        w.writerows(rows)

    n_sm = sum(r["class"] == "small_molecule" for r in rows)
    print(f"\nwrote {len(rows)} unique ligands -> {out}")
    print(f"  small molecules: {n_sm}   peptide/large: {len(rows) - n_sm}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
