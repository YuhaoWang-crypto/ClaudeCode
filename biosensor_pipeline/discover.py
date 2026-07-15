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
GRAPHQL = "https://data.rcsb.org/graphql"


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


def _gql(query: str) -> dict:
    req = urllib.request.Request(GRAPHQL, data=json.dumps({"query": query}).encode(),
                                 headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=45) as r:
        return json.loads(r.read())


def entities_for_entries(pdb_ids: list[str]) -> list[dict]:
    """Batch-fetch every protein chain of many entries (seq/length/name) in one call."""
    if not pdb_ids:
        return []
    ids = "[" + ",".join(f'"{p}"' for p in pdb_ids) + "]"
    q = ("{entries(entry_ids:%s){rcsb_id polymer_entities{rcsb_id "
         "entity_poly{pdbx_seq_one_letter_code_can rcsb_sample_sequence_length "
         "rcsb_entity_polymer_type} rcsb_polymer_entity{pdbx_description} "
         "polymer_entity_instances{rcsb_ligand_neighbors{ligand_comp_id}}}}}" % ids)
    d = _gql(q)
    out = []
    for e in (d.get("data", {}).get("entries") or []):
        for pe in (e.get("polymer_entities") or []):
            ep = pe.get("entity_poly") or {}
            nm = (pe.get("rcsb_polymer_entity") or {}).get("pdbx_description", "")
            contacted = set()
            for inst in (pe.get("polymer_entity_instances") or []):
                for ln in (inst.get("rcsb_ligand_neighbors") or []):
                    contacted.add(ln["ligand_comp_id"])
            out.append({
                "entity_id": pe["rcsb_id"], "pdb": e["rcsb_id"],
                "seq": (ep.get("pdbx_seq_one_letter_code_can") or "").replace("\n", ""),
                "length": ep.get("rcsb_sample_sequence_length"),
                "type": ep.get("rcsb_entity_polymer_type", ""),
                "name": nm,
                "contacts": contacted,
            })
    return out


def discover(smiles: str, name: str = "analyte", max_len: int = 220,
             min_len: int = 45, top: int = 10, rows_per_ccd: int = 25) -> dict:
    """Return ranked receptor candidates for an analyte SMILES.

    min_len excludes co-crystallized peptide fragments (coactivators etc.) that
    are too short to be a stand-alone receptor; a real minimal binder is a
    folded domain (~50-150 aa, per the paper). Candidates are deduped by protein
    name (keep the shortest chain per distinct protein) and ranked short-first.
    """
    ccds = set(ligand_ccds(smiles))
    by_name = {}
    for ccd in ccds:
        entries = entries_with_ccd(ccd, rows=rows_per_ccd)
        for info in entities_for_entries(entries):
            if info["type"] != "Protein" or not info["seq"]:
                continue
            L = info["length"]
            if not L or L < min_len or L > max_len:
                continue
            # CONTACT VERIFICATION: this chain must actually touch one of the
            # analyte's ligand codes -- removes G-proteins / ribosomal chains /
            # detergents that merely co-occur in the same entry.  ✅
            hit = info["contacts"] & ccds
            if not hit:
                continue
            info["ligand_ccd"] = sorted(hit)[0]
            key = info["name"].strip().lower() or info["entity_id"]
            if key not in by_name or L < by_name[key]["length"]:
                by_name[key] = info
    cand = sorted(by_name.values(), key=lambda c: c["length"])
    for c in cand:
        c["contacts"] = sorted(c.get("contacts", []))   # JSON-friendly
    return {"analyte": name, "smiles": smiles, "ligand_ccds": sorted(ccds),
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
