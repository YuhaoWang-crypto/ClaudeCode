#!/usr/bin/env python3
"""
From PDB template -> binding-site (interface) residues.

This is the structural layer that turns a predicted PPI (with a PDB template,
from the humanPPI dataset) into the actual interface: which residues of each
protein form the binding site. That is the starting point for a PPI inhibitor /
interface-targeting design.

Method (standard interface definition):
  * take the two chains of the complex,
  * for every heavy atom of chain A, find heavy atoms of chain B within a cutoff
    (default 5.0 A; a tighter 3.5 A pass flags likely H-bond / salt-bridge core),
  * a residue is an INTERFACE residue if it has any such cross-chain contact,
  * rank interface residues by number of cross-chain contacts = interface
    HOTSPOTS (the residues most worth targeting / mutating).

Structure source: RCSB / PDBe. Those hosts are NOT in this sandbox's allowlist
(HTTP 403), so pass a local file with --pdb, or add the host to egress settings.
The candidate pairs and their template PDB IDs are in
data/structure_templates_immune.tsv (run with --list to see them).

Examples:
  python3 interface_residues.py --pdb 2beg.pdb --chains A,B        # local file
  python3 interface_residues.py --pdb-id 4ZQK --pair PDCD1 CD274   # PD-1/PD-L1
  python3 interface_residues.py --list                            # candidates

Deps: biopython
"""

import argparse
import gzip
import os
import urllib.request
from collections import Counter, defaultdict

from Bio.PDB import PDBParser, MMCIFParser, NeighborSearch
from Bio.PDB.Polypeptide import is_aa

HERE = os.path.dirname(os.path.abspath(__file__))
TEMPLATES = os.path.join(HERE, "data", "structure_templates_immune.tsv")
RESULTS = os.path.join(HERE, "results")


# --------------------------------------------------------------------------
# Structure loading: local file, or gated download from RCSB / PDBe
# --------------------------------------------------------------------------
def fetch_structure_file(pdb_id, dest_dir="/tmp"):
    pid = pdb_id.lower()
    urls = [
        f"https://files.rcsb.org/download/{pdb_id.upper()}.cif",
        f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb",
        f"https://www.ebi.ac.uk/pdbe/entry-files/download/pdb{pid}.ent",
    ]
    for url in urls:
        dest = os.path.join(dest_dir, os.path.basename(url))
        try:
            urllib.request.urlretrieve(url, dest)
            if os.path.getsize(dest) > 200:
                print(f"      downloaded {url}")
                return dest
        except Exception as e:
            print(f"      {url} -> {e.__class__.__name__} "
                  f"({getattr(e, 'code', '')})")
    return None


def load_model(pdb_id=None, local=None):
    if local:
        path = local
    else:
        path = fetch_structure_file(pdb_id)
        if not path:
            raise SystemExit(
                f"\n!! Could not fetch {pdb_id} (structure hosts are blocked here).\n"
                f"   Download {pdb_id}.pdb where RCSB is reachable and rerun with\n"
                f"   --pdb {pdb_id}.pdb")
    opener = gzip.open if path.endswith(".gz") else open
    parser = MMCIFParser(QUIET=True) if path.endswith((".cif", ".cif.gz")) \
        else PDBParser(QUIET=True)
    with opener(path, "rt") as fh:
        structure = parser.get_structure("x", fh)
    return structure[0], path


# --------------------------------------------------------------------------
# Interface computation
# --------------------------------------------------------------------------
def heavy_atoms(chain):
    atoms = []
    for res in chain:
        if not is_aa(res, standard=True):
            continue
        for atom in res:
            if atom.element != "H":
                atoms.append(atom)
    return atoms


def auto_chains(model):
    prot = [c.id for c in model if any(is_aa(r, standard=True) for r in c)]
    return prot


def compute_interface(model, c1, c2, cutoff=5.0, core=3.5):
    a1 = heavy_atoms(model[c1])
    a2 = heavy_atoms(model[c2])
    if not a1 or not a2:
        raise SystemExit(f"!! chain {c1} or {c2} has no standard residues")
    ns = NeighborSearch(a2)
    contacts = []                      # (res1_id, res2_id, dist)
    cnt1, cnt2 = Counter(), Counter()  # cross-chain contact counts per residue
    core1, core2 = set(), set()
    resname = {}
    for atom in a1:
        for nb in ns.search(atom.coord, cutoff):
            d = atom - nb
            r1, r2 = atom.get_parent(), nb.get_parent()
            id1 = (r1.id[1], r1.resname)
            id2 = (r2.id[1], r2.resname)
            resname[(c1, id1[0])] = id1[1]
            resname[(c2, id2[0])] = id2[1]
            contacts.append((id1[0], id2[0], d))
            cnt1[id1[0]] += 1
            cnt2[id2[0]] += 1
            if d <= core:
                core1.add(id1[0]); core2.add(id2[0])
    return contacts, cnt1, cnt2, core1, core2, resname


def fmt_residues(cnt, core, resname, chain, top_hot=8):
    items = sorted(cnt.items())          # by residue number
    hot = {r for r, _ in cnt.most_common(top_hot)}
    out = []
    for resnum, n in items:
        rn = resname.get((chain, resnum), "?")
        flags = ("*" if resnum in hot else " ") + ("c" if resnum in core else " ")
        out.append(f"{rn}{resnum}{flags}({n})")
    return out


# --------------------------------------------------------------------------
# Dataset lookup
# --------------------------------------------------------------------------
def load_templates():
    pairs = {}
    if not os.path.exists(TEMPLATES):
        return pairs
    with open(TEMPLATES) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("name1") or not line.strip():
                continue
            c = line.rstrip("\n").split("\t")
            pairs[frozenset((c[0], c[1]))] = {
                "rfprob": c[2], "afprob": c[3], "pdbtemp": c[4],
                "exact": c[5], "ortho": c[6], "homo": c[7]}
    return pairs


def list_candidates():
    pairs = load_templates()
    rows = sorted(pairs.items(), key=lambda kv: -float(kv[1]["rfprob"]))
    print(f"{'pair':22}{'RFprob':>8}{'AFprob':>8}  {'tmpl':9} exact/ortho/homo templates")
    for k, v in rows:
        a, b = tuple(k)
        tmpl = v["exact"] if v["exact"] != "none" else (
            v["ortho"] if v["ortho"] != "none" else v["homo"])
        print(f"{a+'–'+b:22}{v['rfprob']:>8}{v['afprob']:>8}  {v['pdbtemp']:9} {tmpl}")


# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", help="local structure file (.pdb/.cif[.gz])")
    ap.add_argument("--pdb-id", help="PDB ID to download (needs RCSB reachable)")
    ap.add_argument("--chains", help="two chain IDs, e.g. A,B (else auto-detect)")
    ap.add_argument("--pair", nargs=2, metavar=("P1", "P2"),
                    help="gene-name pair (for labelling / template lookup)")
    ap.add_argument("--cutoff", type=float, default=5.0)
    ap.add_argument("--core", type=float, default=3.5)
    ap.add_argument("--list", action="store_true", help="list structure-ready candidates")
    args = ap.parse_args()

    if args.list:
        list_candidates()
        return

    # resolve PDB id from the dataset if only a pair is given
    if args.pair and not args.pdb and not args.pdb_id:
        pairs = load_templates()
        rec = pairs.get(frozenset(args.pair))
        if rec and rec["exact"] != "none":
            args.pdb_id = rec["exact"].split(",")[0]
            print(f"[lookup] {args.pair[0]}–{args.pair[1]} exact template -> {args.pdb_id}")
        else:
            raise SystemExit(f"!! no exact template for {args.pair} in dataset")

    print("=" * 74)
    label = "–".join(args.pair) if args.pair else (args.pdb_id or args.pdb)
    print(f"INTERFACE / BINDING-SITE RESIDUES | {label}")
    print("=" * 74)
    model, path = load_model(args.pdb_id, args.pdb)

    if args.chains:
        c1, c2 = args.chains.split(",")
    else:
        ch = auto_chains(model)
        if len(ch) < 2:
            raise SystemExit("!! need >=2 protein chains; specify --chains")
        c1, c2 = ch[0], ch[1]
        print(f"[auto] using chains {c1}, {c2} (of {ch})")

    contacts, cnt1, cnt2, core1, core2, resname = compute_interface(
        model, c1, c2, args.cutoff, args.core)

    print(f"\nstructure: {os.path.basename(path)}  chains {c1}/{c2}  "
          f"cutoff={args.cutoff}A  core(H-bond)={args.core}A")
    print(f"cross-chain contacts: {len(contacts)}   "
          f"interface residues: {c1}={len(cnt1)}, {c2}={len(cnt2)}")
    print(f"\n[{c1}] binding-site residues  (*=hotspot  c=core<{args.core}A  (n)=#contacts):")
    print("   " + "  ".join(fmt_residues(cnt1, core1, resname, c1)))
    print(f"\n[{c2}] binding-site residues:")
    print("   " + "  ".join(fmt_residues(cnt2, core2, resname, c2)))

    # closest residue-residue contacts = the core interface interactions
    best = defaultdict(lambda: 1e9)
    for r1, r2, d in contacts:
        best[(r1, r2)] = min(best[(r1, r2)], d)
    print(f"\nclosest cross-chain residue pairs (binding core):")
    for (r1, r2), d in sorted(best.items(), key=lambda kv: kv[1])[:10]:
        n1 = resname.get((c1, r1), "?"); n2 = resname.get((c2, r2), "?")
        print(f"   {n1}{r1} ({c1})  ––  {n2}{r2} ({c2})   {d:.2f} A")

    os.makedirs(RESULTS, exist_ok=True)
    slug = "".join(ch if ch.isalnum() else "_" for ch in label).strip("_")
    out = os.path.join(RESULTS, f"interface_{slug}.tsv")
    with open(out, "w") as fh:
        fh.write(f"# interface residues for {label} from {os.path.basename(path)} "
                 f"chains {c1}/{c2}, cutoff {args.cutoff}A\n")
        fh.write("chain\tresnum\tresname\tn_contacts\tcore\thotspot\n")
        for cnt, core, ch in [(cnt1, core1, c1), (cnt2, core2, c2)]:
            hot = {r for r, _ in cnt.most_common(8)}
            for resnum, n in sorted(cnt.items()):
                fh.write(f"{ch}\t{resnum}\t{resname.get((ch,resnum),'?')}\t{n}"
                         f"\t{'yes' if resnum in core else 'no'}"
                         f"\t{'yes' if resnum in hot else 'no'}\n")
    print(f"\nwrote -> results/{os.path.basename(out)}")


if __name__ == "__main__":
    main()
