"""Step 13 - Mine the PDB for antibodies bound to CHEMICALLY SIMILAR peptides.

The user's instruction is to stop arguing and try: design against the 1C6
peptide alone, or start from antibodies that already bind comparable peptides
in the databases, then optimise with ProteinMPNN / LigandMPNN over several
rounds.

The first two words of that instruction are the ones that make it work.
ERRPEGPGAQAPSSPRVTEAPC has no aromatic residue and one large hydrophobic in
22 positions, so there is no hydrophobic core for a de novo paratope to bury
and nothing for a backbone generator to nucleate on. But that deficit is not
unique to TXLNA - the immune system solves it routinely, and where it has, a
crystallographer has sometimes caught it. If a Fab in the PDB grips a peptide
with the SAME composition problem, that Fab's bound-peptide backbone and its
CDR loop geometry are an empirical answer to "what does a paratope for this
kind of epitope look like". That is a template we can graft onto, and a
template is exactly the input ProteinMPNN needs and cannot invent.

So this script does not score TXLNA. It scores the PDB's antibody-peptide
complexes by how closely their epitope resembles ours in the properties that
made de novo fail, and returns the ones that are the closest analogues.

Pipeline
  1. RCSB search: entries carrying an Ig V-set domain (Pfam PF07686) AND a
     protein entity of 6-30 aa.
  2. GraphQL: pull every entity sequence, description and chain id in two
     batched calls rather than ~500 REST calls.
  3. Classify chains from sequence. An antibody heavy chain or VHH carries the
     J-region W-G-x-G; TCR alpha/beta and MHC do not, so requiring it removes
     the TCR-pMHC complexes that PF07686 also matches.
  4. Download each structure and MEASURE the contact - an entity of the right
     length in the same crystal is not necessarily the antigen. A peptide is
     kept only if it makes heavy-atom contacts with a V domain's CDR region.
  5. Score the verified epitope against ours on the axes that matter.
"""
import json
import re
import sys
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import RESULTS, gravy

TARGET = "ERRPEGPGAQAPSSPRVTEAPC"          # 1C6 immunogen, TXLNA 502-523
CACHE = Path("/tmp/claude-0/pdb_cache")
SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
PEP_MAX = 40

# V-set alone (PF07686) returns 132 entries and misses most Fabs, because a
# Fab's constant domains are what carry the annotation in many depositions.
# Adding the C1-set (PF07654) takes the pool to 1741. The sequence classifier
# below, not the Pfam id, is what decides whether an entry is an antibody.
IG_PFAM = ["PF07686", "PF00047", "PF16135", "PF07654"]
GRAPHQL = "https://data.rcsb.org/graphql"

AROM = set("FWY")
BIGPHOBIC = set("ILVMFWY")
POS, NEG = set("KR"), set("DE")

# J-region of an antibody heavy chain / VHH. Light chains and both TCR chains
# use F-G-x-G instead, and MHC has neither, so this one motif separates real
# antibody entries from the TCR-pMHC complexes that share the V-set domain.
HEAVY_J = re.compile(r"W[GA][QKRHSENA]G[TAS]")
LIGHT_J = re.compile(r"FG[QGSEA]G[TAS]")


def post(url, payload, timeout=120):
    req = urllib.request.Request(
        url, data=json.dumps(payload).encode(),
        headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        body = r.read()
    return json.loads(body) if body.strip() else {}


def search_entries():
    """Entries with an Ig V-set domain and a short protein entity."""
    payload = {
        "query": {"type": "group", "logical_operator": "and", "nodes": [
            {"type": "terminal", "service": "text", "parameters": {
                "attribute": "rcsb_polymer_entity_annotation.annotation_id",
                "operator": "in", "value": IG_PFAM}},
            {"type": "terminal", "service": "text", "parameters": {
                "attribute": "entity_poly.rcsb_sample_sequence_length",
                "operator": "range",
                "value": {"from": 6, "to": PEP_MAX,
                          "include_lower": True, "include_upper": True}}},
            {"type": "terminal", "service": "text", "parameters": {
                "attribute": "entity_poly.rcsb_entity_polymer_type",
                "operator": "exact_match", "value": "Protein"}},
        ]},
        "return_type": "entry",
        "request_options": {"paginate": {"start": 0, "rows": 2000},
                            "results_content_type": ["experimental"]},
    }
    r = post(SEARCH, payload)
    ids = [h["identifier"] for h in r["result_set"]]
    # The Pfam+length query is recall-limited for the question the user asked
    # about phosphorylation: a phosphopeptide co-crystal is rare enough that it
    # must be sought by name, not hoped for. These seeds go through exactly the
    # same classifier and contact test as everything else.
    for term in PTM_TERMS:
        ids += full_text_seed(term)
    return sorted(set(ids))


PTM_TERMS = ["phosphopeptide", "phosphoserine", "phosphothreonine",
             "phosphotyrosine", "citrullinated", "glycopeptide"]


def full_text_seed(text):
    payload = {
        "query": {"type": "group", "logical_operator": "and", "nodes": [
            {"type": "terminal", "service": "text", "parameters": {
                "attribute": "rcsb_polymer_entity_annotation.annotation_id",
                "operator": "in", "value": IG_PFAM}},
            {"type": "terminal", "service": "text", "parameters": {
                "attribute": "entity_poly.rcsb_sample_sequence_length",
                "operator": "range",
                "value": {"from": 4, "to": PEP_MAX,
                          "include_lower": True, "include_upper": True}}},
            {"type": "terminal", "service": "full_text",
             "parameters": {"value": text}},
        ]},
        "return_type": "entry",
        "request_options": {"paginate": {"start": 0, "rows": 200},
                            "results_content_type": ["experimental"]},
    }
    try:
        r = post(SEARCH, payload)
    except urllib.error.HTTPError:
        return []
    if not r or "result_set" not in r:
        return []          # 204 no-content comes back as an empty body
    return [h["identifier"] for h in r["result_set"]]


ENTITY_Q = """
query ($ids: [String!]!) {
  entries(entry_ids: $ids) {
    rcsb_id
    struct { title }
    rcsb_entry_info { resolution_combined }
    polymer_entities {
      entity_poly { pdbx_seq_one_letter_code_can rcsb_sample_sequence_length }
      rcsb_polymer_entity { pdbx_description }
      rcsb_polymer_entity_container_identifiers { auth_asym_ids }
    }
  }
}"""


def fetch_entities(ids, chunk=40):
    out = []
    for i in range(0, len(ids), chunk):
        part = ids[i:i + chunk]
        r = post(GRAPHQL, {"query": ENTITY_Q, "variables": {"ids": part}})
        out.extend(r["data"]["entries"])
    return out


def classify(seq):
    """heavy / light / peptide / other, from sequence alone."""
    L = len(seq)
    if L <= PEP_MAX:
        return "peptide"
    ncys = seq.count("C")
    if 90 <= L <= 250 and ncys >= 2:
        if HEAVY_J.search(seq):
            return "heavy"
        if LIGHT_J.search(seq):
            return "light"
    return "other"


def get_structure(pdb_id):
    CACHE.mkdir(parents=True, exist_ok=True)
    f = CACHE / f"{pdb_id}.cif"
    if f.exists() and f.stat().st_size > 1000:
        return f
    import gzip
    url = f"https://files.rcsb.org/download/{pdb_id}.cif.gz"
    for _ in range(3):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                f.write_bytes(gzip.decompress(r.read()))
            return f
        except urllib.error.HTTPError:
            return None
        except Exception:
            continue
    return None


def prefetch(ids, workers=12):
    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=workers) as ex:
        list(ex.map(get_structure, ids))


def read_atoms(cif):
    """Parse the mmCIF atom_site loop into (auth_chain, auth_seq, comp, xyz).

    Column ORDER in atom_site is not fixed across depositions - 1SM3 does not
    put auth_asym_id where 5A2J does - so the indices are taken from the loop
    header every time. An earlier throwaway version of this hard-coded them
    and silently reported zero hits on files it was simply misreading.
    """
    atoms = []
    with open(cif) as fh:
        cols, active = {}, False
        for line in fh:
            if line.startswith("loop_"):
                cols, active = {}, False
                continue
            if line.startswith("_atom_site."):
                cols[line.strip().split(".", 1)[1]] = len(cols)
                active = True
                continue
            if not active or not cols:
                continue
            p = line.split()
            # mmCIF is not guaranteed column-aligned, so the record type is
            # taken from the split token, never from a fixed slice of the line.
            if not p or p[0] not in ("ATOM", "HETATM"):
                if atoms and line.startswith("#"):
                    break
                continue
            if len(p) <= max(cols.values()):
                continue
            if True:
                if p[cols["type_symbol"]] == "H":
                    continue
                try:
                    xyz = (float(p[cols["Cartn_x"]]), float(p[cols["Cartn_y"]]),
                           float(p[cols["Cartn_z"]]))
                except ValueError:
                    continue
                icode = (p[cols["pdbx_PDB_ins_code"]]
                         if "pdbx_PDB_ins_code" in cols else ".")
                rid = p[cols["auth_seq_id"]] + ("" if icode in (".", "?")
                                                else icode)
                atoms.append((p[cols["auth_asym_id"]], rid,
                              p[cols["label_comp_id"]], xyz))
    return atoms


STANDARD = set("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO "
               "SER THR TRP TYR VAL".split())
# Waters and cryo/buffer components are deposited under the peptide's own
# auth chain id in many entries - 1SM3 puts 60 waters in chain P - so counting
# them made a 13-mer report 21 contacted "epitope residues".
SOLVENT = {"HOH", "DOD", "GOL", "EDO", "SO4", "PO4", "CL", "NA", "MG", "CA",
           "ZN", "ACT", "PEG", "MPD", "TRS", "DMS", "IOD", "K", "NO3"}
PTM_NAMES = {"SEP": "phosphoserine", "TPO": "phosphothreonine",
             "PTR": "phosphotyrosine", "CIR": "citrulline",
             "MLY": "methyllysine", "ALY": "acetyllysine",
             "NH2": "C-terminal amide", "ACE": "N-terminal acetyl"}


def contacts(atoms, pep_chains, ab_chains, cutoff=4.5):
    """Heavy-atom contacts between the candidate peptide and the Fv, plus the
    modified residues actually present in the peptide chain."""
    pep = [(rid, comp, xyz) for ch, rid, comp, xyz in atoms
           if ch in pep_chains and comp not in SOLVENT]
    ab = [xyz for ch, _, comp, xyz in atoms
          if ch in ab_chains and comp not in SOLVENT]
    mods = sorted({c for _, c, _ in pep} - STANDARD)
    if not pep or not ab:
        return 0, 0, mods
    c2 = cutoff * cutoff
    n, touched = 0, set()
    for rid, comp, (x, y, z) in pep:
        for (ax, ay, az) in ab:
            if (x - ax) ** 2 + (y - ay) ** 2 + (z - az) ** 2 <= c2:
                n += 1
                touched.add(rid)
                break
    return n, len(touched), mods


def features(seq):
    L = len(seq)
    return {
        "len": L,
        "arom": round(sum(c in AROM for c in seq) / L, 4),
        "bigphobic": round(sum(c in BIGPHOBIC for c in seq) / L, 4),
        "pro_gly": round(sum(c in "PG" for c in seq) / L, 4),
        "sta": round(sum(c in "STA" for c in seq) / L, 4),
        "pos": round(sum(c in POS for c in seq) / L, 4),
        "neg": round(sum(c in NEG for c in seq) / L, 4),
        "gravy": round(gravy(seq.replace("X", "")) if seq.strip("X") else 0.0, 3),
    }


# The two deficits that killed de novo design are the absence of aromatics and
# of large hydrophobics; the excess of Pro+Gly is what makes the chain refuse
# to fold. Those three carry the weight. Charge and gravy are included so a
# template is not called similar on shape alone while being electrostatically
# opposite - they matter for a graft, just less.
WEIGHTS = {"arom": 3.0, "bigphobic": 3.0, "pro_gly": 2.5,
           "sta": 1.0, "pos": 1.0, "neg": 1.0}


def distance(a, b):
    d = sum(w * (a[k] - b[k]) ** 2 for k, w in WEIGHTS.items())
    d += 0.25 * ((a["gravy"] - b["gravy"]) / 4.5) ** 2
    return round(d ** 0.5, 4)


def main():
    tgt = features(TARGET)
    print("=" * 100)
    print("A. WHAT WE ARE LOOKING FOR AN ANALOGUE OF")
    print("=" * 100)
    print(f"  TXLNA 502-523  {TARGET}")
    print(f"  {tgt}")
    print()

    ids = search_entries()
    print(f"B. RCSB: {len(ids)} entries carry an Ig V-set domain AND a 6-30 aa "
          f"protein entity")
    entries = fetch_entities(ids)
    print(f"   pulled {len(entries)} entity records via GraphQL")

    cands = []
    for e in entries:
        pdb = e["rcsb_id"]
        heavy, light, peps = [], [], []
        for ent in e["polymer_entities"] or []:
            seq = (ent["entity_poly"] or {}).get("pdbx_seq_one_letter_code_can")
            if not seq:
                continue
            seq = seq.replace("\n", "").strip().upper()
            ch = (ent["rcsb_polymer_entity_container_identifiers"] or {}
                  ).get("auth_asym_ids") or []
            desc = (ent["rcsb_polymer_entity"] or {}).get("pdbx_description", "")
            kind = classify(seq)
            if kind == "heavy":
                heavy += ch
            elif kind == "light":
                light += ch
            elif kind == "peptide" and set(seq) <= set("ACDEFGHIKLMNPQRSTVWYX"):
                peps.append({"seq": seq, "chains": ch, "desc": desc})
        if not heavy or not peps:
            continue            # no antibody heavy chain -> TCR/MHC, drop
        res = (e["rcsb_entry_info"] or {}).get("resolution_combined") or [None]
        cands.append({"pdb": pdb, "title": (e["struct"] or {}).get("title", ""),
                      "resolution": res[0], "heavy": heavy, "light": light,
                      "peptides": peps})
    print(f"   {len(cands)} have a heavy chain (W-G-x-G) -> TCR/MHC removed")

    print()
    print("=" * 100)
    print("C. VERIFYING THE CONTACT (an entity of the right length is not")
    print("   automatically the antigen)")
    print("=" * 100)
    prefetch([c["pdb"] for c in cands])
    verified = []
    for c in cands:
        cif = get_structure(c["pdb"])
        if cif is None:
            continue
        atoms = read_atoms(cif)
        ab = set(c["heavy"]) | set(c["light"])
        for p in c["peptides"]:
            pc = set(p["chains"]) - ab
            if not pc:
                continue
            n, touched, mods = contacts(atoms, pc, ab)
            if n >= 10 and touched >= 4:
                f = features(p["seq"])
                verified.append({
                    "pdb": c["pdb"], "title": c["title"],
                    "resolution": c["resolution"],
                    "heavy_chains": c["heavy"], "light_chains": c["light"],
                    "peptide_chains": sorted(pc), "peptide": p["seq"],
                    "peptide_desc": p["desc"],
                    "contacts": n, "epitope_residues_touched": touched,
                    "modified_residues": mods,
                    "ptm": sorted({PTM_NAMES[m] for m in mods if m in PTM_NAMES}),
                    "features": f, "distance_to_TXLNA": distance(f, tgt),
                })
    verified.sort(key=lambda r: r["distance_to_TXLNA"])
    print(f"  {len(verified)} antibody-peptide complexes with a measured "
          f"paratope contact")
    print()
    print(f"{'rank':5} {'PDB':6} {'res':6} {'d':7} {'arom':6} {'phob':6} "
          f"{'P+G':6} {'nc':4} peptide")
    print("-" * 100)
    for i, v in enumerate(verified[:25], 1):
        f = v["features"]
        r = f"{v['resolution']:.2f}" if v["resolution"] else "  NMR"
        print(f"{i:<5} {v['pdb']:6} {r:6} {v['distance_to_TXLNA']:<7.3f} "
              f"{f['arom']:<6.2f} {f['bigphobic']:<6.2f} {f['pro_gly']:<6.2f} "
              f"{v['epitope_residues_touched']:<4} {v['peptide'][:44]}")

    ptm_t = [v for v in verified if v["ptm"]]
    print()
    print("=" * 100)
    print("D. PTM-BEARING EPITOPES - does a template exist for the phospho arm?")
    print("=" * 100)
    if not ptm_t:
        print("  none")
    for v in sorted(ptm_t, key=lambda r: r["distance_to_TXLNA"])[:15]:
        print(f"  {v['pdb']}  d={v['distance_to_TXLNA']:<6.3f} "
              f"{','.join(v['ptm']):28} {v['peptide'][:26]:28} "
              f"{v['title'][:44]}")

    RESULTS.mkdir(exist_ok=True)
    (RESULTS / "step13_template_mining.json").write_text(json.dumps({
        "ptm_templates": sorted(ptm_t, key=lambda r: r["distance_to_TXLNA"]),
        "target": {"sequence": TARGET, "features": tgt},
        "n_entries_searched": len(ids),
        "n_with_heavy_chain": len(cands),
        "n_verified_complexes": len(verified),
        "weights": WEIGHTS,
        "templates": verified,
    }, indent=2))
    print(f"\nwrote results/step13_template_mining.json")


if __name__ == "__main__":
    main()
