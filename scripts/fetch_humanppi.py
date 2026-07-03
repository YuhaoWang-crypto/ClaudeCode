"""Mine predicted protein-protein interaction partners from the Grishin/Cong-lab
humanPPI database (Science 2025, doi:10.1126/science.adt1630).

The site is a Flask app at http://prodata.swmed.edu/humanPPI with a DataTables
JSON endpoint (reverse-engineered from the results page JS):

  GET /humanPPI/autocomplete?query=<gene>        -> resolve gene -> UniProt entry
  GET /humanPPI/data/<uniprot>?filter=<bool>&length=N&start=0
        -> {data: [{Interacting Partner, Gene Name, Protein Name, AF_Score,
                    RF_Score, DCA_Score, Subcellular Locality, confident, pdb...}],
            recordsTotal: N}

  filter=true  -> high-confidence subset; filter=false -> all predicted partners.

A full bulk dump of all predicted pairs is also available at
https://conglab.swmed.edu/humanPPI/downloads/final_predictions.tar.gz

For Problem B this gives, for a target receptor, its predicted interactors —
useful for peptide/PPI-derived targeting (which membrane proteins engage it) and
for picking surface handles.

Run:  python scripts/fetch_humanppi.py --uniprot P43220 --name GLP1R [--filter]
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import requests

BASE = "http://prodata.swmed.edu/humanPPI"


def resolve_gene(gene: str) -> str | None:
    r = requests.get(f"{BASE}/autocomplete", params={"query": gene}, timeout=30)
    r.raise_for_status()
    hits = r.json()
    return hits[0]["entry"] if hits else None


def fetch_partners(uniprot: str, high_confidence_only: bool) -> list[dict]:
    r = requests.get(f"{BASE}/data/{uniprot}",
                     params={"filter": str(high_confidence_only).lower(),
                             "length": 100000, "start": 0, "draw": 1},
                     timeout=60)
    r.raise_for_status()
    return r.json().get("data", [])


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--uniprot", default="P43220")
    ap.add_argument("--name", default="GLP1R")
    ap.add_argument("--filter", action="store_true",
                    help="high-confidence subset only (default: all predicted)")
    args = ap.parse_args()

    out_dir = Path(__file__).resolve().parents[1] / "data" / "targets" / args.name
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"humanPPI: predicted partners for {args.name} ({args.uniprot}), "
          f"{'high-confidence' if args.filter else 'all'}")
    rows = fetch_partners(args.uniprot, args.filter)
    print(f"  {len(rows)} predicted partners")

    # order columns sensibly; sort by AlphaFold interface score (best first)
    def score(r):
        try:
            return float(r.get("AF_Score", 0) or 0)
        except ValueError:
            return 0.0
    rows.sort(key=score, reverse=True)

    cols = ["Interacting Partner", "Gene Name", "Protein Name",
            "Subcellular Locality", "AF_Score", "RF_Score", "DCA_Score",
            "confident", "Function", "pdb"]
    out = out_dir / f"{args.name}_humanppi_partners.csv"
    with out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {len(rows)} partners -> {out}")

    # quick highlight: membrane/surface partners (relevant for targeting)
    surf = [r for r in rows if "membrane" in str(r.get("Subcellular Locality", "")).lower()]
    print(f"  of which membrane-localized: {len(surf)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
