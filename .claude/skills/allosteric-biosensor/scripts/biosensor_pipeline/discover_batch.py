"""
discover_batch.py -- mine many analyte -> small-receptor pairs at once.

Runs Tier-A PDB mining (discover.py) across a panel of analytes, resolving each
name to a SMILES via PubChem, then compiles a single catalog of candidate
receptor domains ranked short-first (the paper's minimal-binder principle).

    python3 -m biosensor_pipeline.discover_batch                 # default panel
    python3 -m biosensor_pipeline.discover_batch cortisol biotin # custom names

Writes biosensor_out/receptor_catalog.json + prints a size-ranked table.
Uses only public RCSB + PubChem REST (no MCP, no key).
"""

from __future__ import annotations
import json
import os
import sys
import urllib.parse
import urllib.request

from .discover import discover

OUT = os.path.join(os.getcwd(), "biosensor_out")

# analytes likely to have small, folded binding domains in the PDB
DEFAULT_PANEL = [
    # steroids / hormones
    "cortisol", "progesterone", "testosterone", "estradiol", "aldosterone",
    "dexamethasone", "17-hydroxyprogesterone", "cholic acid", "thyroxine",
    # vitamins / cofactors
    "biotin", "folic acid", "riboflavin", "retinoic acid", "vitamin D3",
    # neurotransmitters / metabolites
    "dopamine", "serotonin", "bilirubin", "prostaglandin E2",
    # drugs / xenobiotics
    "fentanyl", "caffeine", "theophylline", "colchicine", "warfarin",
    "methotrexate", "digoxigenin",
    # dyes / fluorogens
    "fluorescein",
]


def smiles_for(name: str) -> str | None:
    url = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
           + urllib.parse.quote(name)
           + "/property/CanonicalSMILES/JSON")
    try:
        with urllib.request.urlopen(url, timeout=25) as r:
            d = json.loads(r.read())
        props = d["PropertyTable"]["Properties"][0]
        return props.get("CanonicalSMILES") or props.get("ConnectivitySMILES")
    except Exception:
        return None


def main(argv):
    names = argv[1:] if len(argv) > 1 else DEFAULT_PANEL
    os.makedirs(OUT, exist_ok=True)
    catalog = []
    print(f"mining {len(names)} analytes (Tier-A PDB co-crystal)…\n")
    for nm in names:
        sm = smiles_for(nm)
        if not sm:
            print(f"  [skip] {nm}: no SMILES")
            continue
        try:
            res = discover(sm, name=nm, max_len=200, min_len=45, top=5)
        except Exception as e:
            print(f"  [err ] {nm}: {e}")
            continue
        n = res["n_candidates"]
        best = res["candidates"][0] if res["candidates"] else None
        tag = f"{best['entity_id']} {best['length']}aa {best['name'][:34]}" if best else "-"
        print(f"  {nm:24s} candidates={n:3d}  smallest: {tag}")
        for c in res["candidates"]:
            catalog.append({
                "analyte": nm, "ligand_ccd": c["ligand_ccd"],
                "entity_id": c["entity_id"], "pdb": c["pdb"],
                "length": c["length"], "receptor_name": c["name"], "seq": c["seq"],
            })

    catalog.sort(key=lambda r: r["length"])
    with open(os.path.join(OUT, "receptor_catalog.json"), "w") as fh:
        json.dump(catalog, fh, indent=1)

    print(f"\n{'='*78}\nGLOBAL CATALOG — {len(catalog)} analyte×receptor pairs, ranked short-first\n{'='*78}")
    print(f"  {'len':>4s}  {'entity':8s}  {'ligand':6s}  {'analyte':22s}  receptor")
    for r in catalog[:45]:
        print(f"  {r['length']:4d}  {r['entity_id']:8s}  {r['ligand_ccd']:6s}  "
              f"{r['analyte']:22s}  {r['receptor_name'][:38]}")
    print(f"\nwrote {OUT}/receptor_catalog.json ({len(catalog)} pairs)")
    print("Pick any row -> its 'seq' is ready; derive loop sites with "
          "structure.annotate(pdb,'A',seq); add Receptor+System; run_repro + Boltz.")


if __name__ == "__main__":
    main(sys.argv)
