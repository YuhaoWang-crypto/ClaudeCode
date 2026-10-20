#!/usr/bin/env python3
"""Stage 5: group the survivors into subfamilies and find each one's neighbours.

The precedent's highest-value output was not a score, it was a **clade**: all
three experimentally validated enzymes fell into one clade of five, which turns
"19 structural candidates" into "prioritise this branch, and use the sister
branches as controls". This script reconstructs that grouping.

    python3 cluster.py out/run.domains.fasta \
        --table out/run.annotated.tsv \
        --references O64642 A0A0P0VE25 A0A0P0XSC2 \
        --max-seqs 80 --out out/run

IMPORTANT -- this is **distance-based clustering, not maximum-likelihood
phylogeny.** The source paper built an ML tree (MEGA, JTT+gamma, bootstrap).
Average-linkage clustering on pairwise identity is a cheaper approximation
that recovers subfamily structure but carries no branch support and should
never be reported as a phylogeny. Where the clade call actually matters, rerun
with a real tree builder (IQ-TREE/RAxML on a MAFFT alignment) and say which
one you used.

Distances are computed on the **catalytic domain only** (the region that
aligned to the query, as emitted by `mine.py --emit-fasta`); full-length
sequences mix in unrelated domains and distort the grouping.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import sys

import numpy as np

try:
    from Bio import Align
    from Bio.Align import substitution_matrices
except ImportError:                                      # pragma: no cover
    sys.exit("needs biopython: pip install biopython")

from scipy.cluster.hierarchy import fcluster, linkage, to_tree   # noqa: E402
from scipy.spatial.distance import squareform                     # noqa: E402


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            elif line:
                buf.append(line)
    if name:
        seqs[name] = "".join(buf)
    return seqs


def make_aligner():
    a = Align.PairwiseAligner()
    a.substitution_matrix = substitution_matrices.load("BLOSUM62")
    a.open_gap_score = -11
    a.extend_gap_score = -1
    a.mode = "global"
    return a


def pairwise_identity(aligner, s1, s2):
    """Fraction identity over the shorter sequence (robust for remote pairs)."""
    try:
        aln = aligner.align(s1, s2)[0]
    except Exception:                                    # noqa: BLE001
        return 0.0
    a, b = aln[0], aln[1]
    matches = sum(1 for x, y in zip(a, b) if x == y and x != "-")
    return matches / max(1, min(len(s1), len(s2)))


def newick(node, labels):
    if node.is_leaf():
        return f"{labels[node.id]}:{node.dist:.4f}"
    left = newick(node.get_left(), labels)
    right = newick(node.get_right(), labels)
    return f"({left},{right}):{node.dist:.4f}"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fasta", help="catalytic domains from mine.py --emit-fasta")
    ap.add_argument("--table", help="annotated candidates.tsv, for labels")
    ap.add_argument("--references", nargs="*", default=[],
                    help="accessions with known/validated function -- the "
                         "clade containing these is the one to prioritise")
    ap.add_argument("--include", nargs="*", default=[],
                    help="force these accessions into the set without marking "
                         "them as references -- use it to keep the sample "
                         "taxonomically spread rather than top-score-heavy")
    ap.add_argument("--per-lineage", type=int, default=0,
                    help="if set, take at most this many per lineage_group "
                         "(needs --table); the cheapest way to avoid a tree "
                         "built almost entirely from the query's own kingdom")
    ap.add_argument("--max-seqs", type=int, default=80,
                    help="cap the set; all-vs-all is quadratic")
    ap.add_argument("--cut", type=float, default=0.55,
                    help="distance at which to cut the dendrogram into clades")
    ap.add_argument("--min-len", type=int, default=120)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    seqs = read_fasta(args.fasta)
    meta = {}
    if args.table:
        for r in csv.DictReader(open(args.table), delimiter="\t"):
            meta.setdefault(r["accession"], r)

    # Rank by catalytic-core quality, but always keep the references.
    def rank(acc):
        r = meta.get(acc)
        if not r:
            return 9e9
        try:
            return float(r.get("core_rmsd", 9e9))
        except (TypeError, ValueError):
            return 9e9

    pool = [a for a, s in seqs.items() if len(s) >= args.min_len]
    chosen = [a for a in args.references if a in pool]
    for acc in args.include:
        if acc in pool and acc not in chosen:
            chosen.append(acc)

    per_lineage = {}
    for acc in sorted(pool, key=rank):
        if len(chosen) >= args.max_seqs:
            break
        if acc in chosen:
            continue
        if args.per_lineage:
            grp = meta.get(acc, {}).get("lineage_group", "unknown")
            if per_lineage.get(grp, 0) >= args.per_lineage:
                continue
            per_lineage[grp] = per_lineage.get(grp, 0) + 1
        chosen.append(acc)
    missing = [a for a in args.references if a not in pool]
    if missing:
        print(f"references absent from the FASTA: {', '.join(missing)}")

    print(f"clustering {len(chosen)} catalytic domains "
          f"({len(chosen)*(len(chosen)-1)//2} pairwise alignments)")
    aligner = make_aligner()
    n = len(chosen)
    D = np.zeros((n, n))
    for k, (i, j) in enumerate(itertools.combinations(range(n), 2)):
        d = 1.0 - pairwise_identity(aligner, seqs[chosen[i]], seqs[chosen[j]])
        D[i, j] = D[j, i] = d
        if k % 250 == 0 and k:
            print(f"  {k} pairs", flush=True)

    Z = linkage(squareform(D, checks=False), method="average")
    labels = fcluster(Z, t=args.cut, criterion="distance")

    with open(f"{args.out}.tree.nwk", "w") as fh:
        fh.write(newick(to_tree(Z), chosen) + ";\n")
    np.savetxt(f"{args.out}.distance.tsv", D, delimiter="\t",
               header="\t".join(chosen), fmt="%.4f")

    groups = {}
    for acc, cl in zip(chosen, labels):
        groups.setdefault(int(cl), []).append(acc)

    ref_clades = {int(cl) for acc, cl in zip(chosen, labels)
                  if acc in args.references}

    rows = []
    print(f"\n{len(groups)} clades at cut distance {args.cut}")
    for cl in sorted(groups, key=lambda c: (-len(groups[c]), c)):
        members = groups[cl]
        mark = "  <== contains reference(s)" if cl in ref_clades else ""
        print(f"\nclade {cl}: {len(members)} members{mark}")
        for acc in sorted(members, key=rank):
            r = meta.get(acc, {})
            flag = "*" if acc in args.references else " "
            print(f"  {flag} {acc:12s} core={r.get('core_rmsd','?'):>6s} "
                  f"{r.get('lineage_group','?')[:13]:13s} "
                  f"{r.get('tax_name','?')[:24]:24s} "
                  f"{r.get('description','')[:40]}")
            rows.append({
                "accession": acc, "clade": cl,
                "with_reference": cl in ref_clades,
                "is_reference": acc in args.references,
                "core_rmsd": r.get("core_rmsd", ""),
                "lineage_group": r.get("lineage_group", ""),
                "tax_name": r.get("tax_name", ""),
                "known_class": r.get("known_class", ""),
                "description": r.get("description", ""),
            })

    with open(f"{args.out}.clades.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    n_with = sum(1 for r in rows if r["with_reference"] and not r["is_reference"])
    print(f"\nwrote {args.out}.clades.tsv, {args.out}.tree.nwk, "
          f"{args.out}.distance.tsv")
    print(f"{n_with} non-reference proteins share a clade with a reference "
          f"-- those are the priority candidates.")
    print("Distance-based clustering, NOT an ML phylogeny: no branch support. "
          "Label it accordingly.")


if __name__ == "__main__":
    main()
