"""
discover.py -- mine a receptor for an arbitrary analyte from the PDB.

Given an analyte (SMILES), find deposited structures whose protein chain is
co-crystallized with that ligand (a ready-made receptor), rank candidates by the
paper's "small, rigid, minimal binder" preference, and emit the sequence +
structure-derived loop sites so build_chimera can graft it into the reporter.

    python3 -m biosensor_pipeline.discover "<SMILES>" [--name NAME] [--max-len 200]

This is TIER A of receptor mining (PDB co-crystal). Tier B (known binding
protein via ChEMBL/literature) and Tier C (de-novo design via Boltz protein
design) are described in the skill's reference/adding-a-system.md.

Uses the public RCSB Search + Data APIs (no key). ✅ discovery is real data;
whether a mined receptor yields a working switch is ⚠️ until Boltz-validated
and, ultimately, assayed.
"""

from __future__ import annotations
import json
import sys
import urllib.request

SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
DATA = "https://data.rcsb.org/rest/v1/core"


def _post(url, payload):
    req = urllib.request.Request(url, data=json.dumps(payload).encode(),
                                 headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=30) as r:
        if r.status == 204:
            return {}
        return json.loads(r.read())


def _get(url):
    with urllib.request.urlopen(url, timeout=30) as r:
        return json.loads(r.read())


def ligand_ccds(smiles: str, max_ccds: int = 5) -> list[str]:
    """Find PDB chemical-component IDs matching the analyte (graph-relaxed)."""
    q = {
        "query": {"type": "terminal", "service": "chemical",
                  "parameters": {"value": smiles, "type": "descriptor",
                                 "descriptor_type": "SMILES",
                                 "match_type": "graph-relaxed-stereo"}},
        "return_type": "mol_definition",
        "request_options": {"paginate": {"start": 0, "rows": max_ccds}},
    }
    d = _post(SEARCH, q)
    return [r["identifier"] for r in d.get("result_set", [])]


def entries_with_ccd(ccd: str, rows: int = 25) -> list[str]:
    """PDB entries that contain a given ligand CCD."""
    q = {
        "query": {"type": "terminal", "service": "text_chem",
                  "parameters": {"attribute": "rcsb_chem_comp_container_identifiers.comp_id",
                                 "operator": "exact_match", "value": ccd}},
        "return_type": "entry",
        "request_options": {"paginate": {"start": 0, "rows": rows},
                            "results_content_type": ["experimental"]},
    }
    d = _post(SEARCH, q)
    return [r["identifier"] for r in d.get("result_set", [])]


def entities_with_ccd(ccd: str, rows: int = 25) -> list[str]:
    """Polymer-entity ids for the protein chains in entries containing the CCD."""
    ents = []
    for pdb in entries_with_ccd(ccd, rows=rows):
        try:
            d = _get(f"{DATA}/entry/{pdb}")
            ids = d.get("rcsb_entry_container_identifiers", {}).get("polymer_entity_ids", [])
            ents.extend(f"{pdb}_{e}" for e in ids)
        except Exception:
            continue
    return ents


def entity_info(entity_id: str) -> dict:
    """Sequence + length + name for a polymer entity id like '1DBB_1'."""
    pdb, ent = entity_id.split("_")
    d = _get(f"{DATA}/polymer_entity/{pdb}/{ent}")
    poly = d.get("entity_poly", {})
    names = d.get("rcsb_polymer_entity", {})
    return {
        "entity_id": entity_id, "pdb": pdb,
        "seq": poly.get("pdbx_seq_one_letter_code_can", "").replace("\n", ""),
        "length": d.get("rcsb_polymer_entity_container_identifiers", {}).get("entity_id") and
                  len(poly.get("pdbx_seq_one_letter_code_can", "").replace("\n", "")),
        "name": names.get("pdbx_description", ""),
        "type": poly.get("rcsb_entity_polymer_type", ""),
    }


def discover(smiles: str, name: str = "analyte", max_len: int = 220, top: int = 10) -> dict:
    """Return ranked receptor candidates for an analyte SMILES."""
    ccds = ligand_ccds(smiles)
    cand, seen = [], set()
    for ccd in ccds:
        for ent in entities_with_ccd(ccd):
            if ent in seen:
                continue
            seen.add(ent)
            try:
                info = entity_info(ent)
            except Exception:
                continue
            if info["type"] != "Protein" or not info["seq"]:
                continue
            info["ligand_ccd"] = ccd
            cand.append(info)
    # rank: prefer SHORT chains (paper's minimal-binder principle), drop huge ones
    cand = [c for c in cand if c["length"] and c["length"] <= max_len]
    cand.sort(key=lambda c: c["length"])
    return {"analyte": name, "smiles": smiles, "ligand_ccds": ccds,
            "n_candidates": len(cand), "candidates": cand[:top]}


def main(argv):
    if len(argv) < 2:
        print("usage: python3 -m biosensor_pipeline.discover '<SMILES>' [--name NAME] [--max-len N]")
        return
    smiles = argv[1]
    name = "analyte"; max_len = 220
    if "--name" in argv:
        name = argv[argv.index("--name") + 1]
    if "--max-len" in argv:
        max_len = int(argv[argv.index("--max-len") + 1])
    res = discover(smiles, name=name, max_len=max_len)
    print(f"analyte={name}  ligand CCDs={res['ligand_ccds']}  "
          f"receptor candidates (<= {max_len} aa): {res['n_candidates']}")
    for c in res["candidates"]:
        print(f"  {c['entity_id']:8s} len={c['length']:4d}  ligand={c['ligand_ccd']:4s}  {c['name'][:55]}")
    print("\nNext: pick a candidate, fetch its sequence (already shown), derive loop "
          "sites with structure.annotate(), add a Receptor+System to systems.py, then "
          "run_repro + Boltz. Small, single-domain chains make the best receptors.")


if __name__ == "__main__":
    main(sys.argv)
