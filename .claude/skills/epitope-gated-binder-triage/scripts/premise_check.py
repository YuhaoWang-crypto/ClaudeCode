#!/usr/bin/env python3
"""Species-divergence premise check for a cross-reactivity argument.

Answers, reproducibly and from primary sources, the question that justifies
building a species-specific reagent at all:

  "How different is my target's ectodomain between species, and -- separately --
   how different is the functional interface I intend to block?"

These two numbers usually point in opposite directions, and conflating them is
the single easiest way to write a premise a reviewer can dismantle. In the
TGFbR2 case the ectodomain was 81.8% identical while the ligand interface was
11/12 identical: the divergence that justifies a species-specific antibody is
*peripheral*, not at the functional site. The defensible claim is therefore
about where antibody epitopes sit, not about a species-divergent binding site.

Fetches sequences from UniProt and the complex from RCSB, so it re-derives
everything rather than trusting a number someone typed into a slide.

Usage
-----
    python premise_check.py \\
        --uniprot-a P37173 --uniprot-b Q62312 \\
        --domain 24-166 \\
        --pdb 1KTZ --pdb-chain B \\
        --numbering-ref b
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from collections import defaultdict

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def fetch(url: str) -> str:
    out = subprocess.run(["curl", "-sS", url], capture_output=True, text=True)
    if out.returncode != 0:
        sys.exit(f"fetch failed: {url}\n{out.stderr}")
    return out.stdout


def uniprot_seq(acc: str) -> str:
    fasta = fetch(f"https://rest.uniprot.org/uniprotkb/{acc}.fasta")
    return "".join(l.strip() for l in fasta.splitlines() if not l.startswith(">"))


def pdb_contacts(pdb_id: str, chain: str, cutoff: float = 4.5):
    """Residues of `chain` within cutoff of any other chain. -> {resnum: aa}"""
    text = fetch(f"https://files.rcsb.org/download/{pdb_id}.pdb")
    res: dict[tuple[str, int], list] = defaultdict(list)
    names: dict[tuple[str, int], str] = {}
    for l in text.splitlines():
        if not l.startswith("ATOM"):
            continue
        rn = l[17:20].strip()
        if rn not in THREE_TO_ONE:
            continue
        key = (l[21], int(l[22:26]))
        res[key].append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
        names[key] = THREE_TO_ONE[rn]
    c2 = cutoff * cutoff
    mine = [k for k in res if k[0] == chain]
    other = [k for k in res if k[0] != chain]
    hits = {}
    for a in mine:
        found = False
        for b in other:
            for p in res[a]:
                for q in res[b]:
                    if sum((p[i] - q[i]) ** 2 for i in range(3)) < c2:
                        found = True
                        break
                if found:
                    break
            if found:
                break
        if found:
            hits[a[1]] = names[a]
    return hits


def derive_offset(pdb_res: dict[int, str], ref: str) -> int:
    """Shift such that ref[pdbnum + offset - 1] == pdb residue, maximising matches."""
    best = (None, -1)
    for off in range(0, 60):
        ok = sum(
            1 for n, aa in pdb_res.items()
            if 0 < n + off <= len(ref) and ref[n + off - 1] == aa
        )
        if ok > best[1]:
            best = (off, ok)
    return best[0], best[1]


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--uniprot-a", required=True, help="e.g. human accession")
    p.add_argument("--uniprot-b", required=True, help="e.g. mouse accession")
    p.add_argument("--domain", required=True, help="START-END in UniProt numbering")
    p.add_argument("--pdb", help="experimental complex, e.g. 1KTZ")
    p.add_argument("--pdb-chain", default="B", help="the target's chain in that complex")
    p.add_argument("--numbering-ref", choices=["a", "b"], default="b",
                   help="which UniProt entry the reported epitope numbering follows")
    a = p.parse_args()

    sa, sb = uniprot_seq(a.uniprot_a), uniprot_seq(a.uniprot_b)
    start, end = (int(x) for x in a.domain.split("-"))

    da, db = sa[start - 1:end], sb[start - 1:end]
    if len(da) != len(db):
        print("NOTE: domain slices differ in length -- sequences are not ungapped-alignable "
              "over this region; align properly before treating positions as interchangeable.",
              file=sys.stderr)
    n = min(len(da), len(db))
    ident = sum(x == y for x, y in zip(da[:n], db[:n]))
    subs = [(i + start, x, y) for i, (x, y) in enumerate(zip(da[:n], db[:n])) if x != y]

    print(f"== ectodomain {a.domain} ({a.uniprot_a} vs {a.uniprot_b}) ==")
    print(f"length {n}, identical {ident}, identity {100*ident/n:.1f}%")
    print(f"substitutions ({len(subs)}): " + ", ".join(f"{p_}{x}>{y}" for p_, x, y in subs))

    if not a.pdb:
        return 0

    hits = pdb_contacts(a.pdb, a.pdb_chain)
    # Derive the offset against BOTH entries and keep the better fit. The PDB
    # belongs to one species; validating against the other would flag genuine
    # species differences as mapping errors, which is exactly backwards.
    off_a, m_a = derive_offset(hits, sa)
    off_b, m_b = derive_offset(hits, sb)
    if m_a >= m_b:
        off, matched, src = off_a, m_a, a.uniprot_a
    else:
        off, matched, src = off_b, m_b, a.uniprot_b
    print(f"\n== interface from {a.pdb} chain {a.pdb_chain} ==")
    print(f"derived numbering offset +{off}; best fit is {src} "
          f"({matched}/{len(hits)} residue identities reproduced)")
    print(f"  [{a.uniprot_a}: {m_a}/{len(hits)} at +{off_a}"
          f" | {a.uniprot_b}: {m_b}/{len(hits)} at +{off_b}]")
    if matched != len(hits):
        print(f"WARNING: offset +{off} does not reproduce every residue identity even "
              f"against its best-matching entry ({src}) -- the mapping is unsound, "
              f"do not proceed.", file=sys.stderr)
    if off_a != off_b:
        print(f"NOTE: the two entries need different offsets, so they are not "
              f"ungapped-alignable here; positions are not interchangeable.",
              file=sys.stderr)

    print(f"\n{'pdb':>6}  {'uniprot':>7}  {a.uniprot_a}/{a.uniprot_b}")
    diff = []
    for pn in sorted(hits):
        un = pn + off
        x, y = sa[un - 1], sb[un - 1]
        flag = "  <-- DIFFERS" if x != y else ""
        if x != y:
            diff.append(f"{un}{x}>{y}")
        print(f"{hits[pn]}{pn:>5}  {un:>7}  {x}/{y}{flag}")
    print(f"\ninterface residues: {len(hits)}, differing between species: {len(diff)} "
          f"({', '.join(diff) if diff else 'none'})")
    print(f"interface identity {100*(len(hits)-len(diff))/len(hits):.1f}% "
          f"vs ectodomain {100*ident/n:.1f}%")
    print("\nIf the interface is markedly more conserved than the ectodomain, the "
          "cross-reactivity argument must rest on peripheral epitope divergence, "
          "NOT on a species-divergent binding site.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
