"""Mine DrugCLIP genome-wide virtual-screening predictions for a target.

DrugCLIP ("drug-the-whole-genome", Gao/Lan, Science 2025) publishes, for (almost)
every human protein, a set of predicted small-molecule binders per binding pocket
— each with a DrugCLIP retrieval score, an AutoDock-style docking score, and the
nearby pocket residues. The website is a Vue SPA backed by a FastAPI at
``https://dtwgapi.yanyanlan.com`` (endpoints reverse-engineered from the JS
bundle):

  GET  /protein/{uniprot_id}     -> metadata
  GET  /complexes/{uniprot_id}   -> [{pocket_id, ligand_id, drugclip_score,
                                      docking_score, nearby_residues, source, ...}]
  POST /get_smiles {complex_list} -> "<smiles>\\t<ligand_origin_id>\\n..." per row

This script pulls the predicted ligands for a UniProt ID, resolves their SMILES,
and writes a ranked CSV. These are ready-made Problem-B (targeting-ligand)
candidates for that receptor.

Run:  python scripts/fetch_drugclip.py --uniprot P43220 --name GLP1R
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import requests

API = "https://dtwgapi.yanyanlan.com"


def fetch_complexes(uniprot: str) -> list[dict]:
    r = requests.get(f"{API}/complexes/{uniprot}", timeout=60)
    r.raise_for_status()
    return r.json()


def resolve_smiles(ligand_ids: list[str]) -> dict[str, str]:
    """Map ligand_origin_id -> SMILES via the get_smiles endpoint."""
    r = requests.post(f"{API}/get_smiles", json={"complex_list": ligand_ids}, timeout=120)
    r.raise_for_status()
    blob = r.json()
    text = blob.get("data", blob)
    text = text.get("smiles") if isinstance(text, dict) else text
    mapping: dict[str, str] = {}
    for line in str(text).splitlines():
        parts = line.split("\t")
        if len(parts) >= 2:
            mapping[parts[1].strip()] = parts[0].strip()
    return mapping


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--uniprot", default="P43220")
    ap.add_argument("--name", default="GLP1R")
    args = ap.parse_args()

    out_dir = Path(__file__).resolve().parents[1] / "data" / "targets" / args.name
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"DrugCLIP: fetching predicted ligands for {args.name} ({args.uniprot})")
    comps = fetch_complexes(args.uniprot)
    print(f"  {len(comps)} predicted pocket-ligand complexes")

    smi = resolve_smiles([c["ligand_id"] for c in comps])

    rows = []
    for c in comps:
        rows.append({
            "ligand_id": c.get("ligand_id_origin"),
            "smiles": smi.get(c.get("ligand_id_origin"), ""),
            "drugclip_score": round(c.get("drugclip_score", float("nan")), 4),
            "docking_score": round(c.get("docking_score", float("nan")), 3),
            "source": c.get("source"),
            "pocket_id": c.get("pocket_id"),
            "pocket_group": c.get("pocket_group"),
            "nearby_residues": ";".join(c.get("nearby_residues", [])),
        })
    # rank by DrugCLIP score (retrieval confidence), highest first
    rows.sort(key=lambda x: x["drugclip_score"], reverse=True)

    out = out_dir / f"{args.name}_drugclip_predicted.csv"
    with out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["ligand_id", "smiles", "drugclip_score",
                                          "docking_score", "source", "pocket_id",
                                          "pocket_group", "nearby_residues"])
        w.writeheader()
        w.writerows(rows)

    n_smi = sum(bool(r["smiles"]) for r in rows)
    print(f"wrote {len(rows)} predicted ligands ({n_smi} with SMILES) -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
