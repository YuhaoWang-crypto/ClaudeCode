#!/usr/bin/env python3
"""
M1 - Sequence QC, standardisation and humanness.

QC (format, composition, completeness) is table stakes. The addition here is
*humanness*, which the epitope count alone cannot express: for an antibody-
derived ligand the dominant determinant of clinical immunogenicity is how far
the framework has drifted from the human germline that the recipient's
repertoire is already tolerised to. A VHH with 12 predicted DR binders in a
framework 92% identical to human IGHV3-23 is a different risk from the same
count in a bacterial scaffold with no human counterpart.

Reports per sequence:
  * composition QC (canonical residues, ambiguity codes, length, MW, pI)
  * % identity to human germline IGHV3-23 (global alignment, BLOSUM62)
  * VHH hallmark residues at Kabat 37/44/45/47 for the nanobody-type entries
"""
import csv
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import read_fasta, read_metadata, data_path, results_path  # noqa: E402

CANONICAL = set("ACDEFGHIKLMNPQRSTVWY")
AMBIGUOUS = set("BZXUO*")

# Average residue masses (Da), monomer minus water.
MW = {"A": 71.08, "R": 156.19, "N": 114.10, "D": 115.09, "C": 103.14,
      "E": 129.12, "Q": 128.13, "G": 57.05, "H": 137.14, "I": 113.16,
      "L": 113.16, "K": 128.17, "M": 131.19, "F": 147.18, "P": 97.12,
      "S": 87.08, "T": 101.10, "W": 186.21, "Y": 163.18, "V": 99.13}

# pKa set (EMBOSS) for a quick isoelectric point.
PKA = {"Nterm": 8.6, "Cterm": 3.6, "C": 8.5, "D": 3.9, "E": 4.1,
       "H": 6.5, "K": 10.8, "R": 12.5, "Y": 10.1}


def molecular_weight(seq):
    return sum(MW.get(a, 110.0) for a in seq) + 18.02


def isoelectric_point(seq):
    counts = {a: seq.count(a) for a in "CDEHKRY"}

    def charge(ph):
        q = 1.0 / (1.0 + 10 ** (ph - PKA["Nterm"]))
        q -= 1.0 / (1.0 + 10 ** (PKA["Cterm"] - ph))
        for a in "KRH":
            q += counts[a] / (1.0 + 10 ** (ph - PKA[a]))
        for a in "DECY":
            q -= counts[a] / (1.0 + 10 ** (PKA[a] - ph))
        return q

    lo, hi = 0.0, 14.0
    for _ in range(100):
        mid = (lo + hi) / 2
        if charge(mid) > 0:
            lo = mid
        else:
            hi = mid
    return round((lo + hi) / 2, 2)


def global_identity(a, b):
    """Needleman-Wunsch with BLOSUM62; returns (%identity over aligned cols)."""
    from Bio import Align
    from Bio.Align import substitution_matrices
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    aligner.mode = "global"
    aln = aligner.align(a, b)[0]
    s1, s2 = aln[0], aln[1]
    aligned = [(x, y) for x, y in zip(s1, s2) if x != "-" and y != "-"]
    if not aligned:
        return 0.0
    ident = sum(x == y for x, y in aligned)
    return 100.0 * ident / len(aligned)


def vhh_hallmarks(seq):
    """
    VHH hallmark tetrad in FR2 (Kabat 37/44/45/47). Camelid VHHs carry
    F/Y-E-R-G/L where a human VH carries V-G-L-W. Located by the conserved
    'W(F/V)RQ' / 'WVRQ' FR2 motif rather than by residue counting, which
    breaks on indel-containing CDR1s.
    """
    i = -1
    for motif in ("WFRQ", "WVRQ", "WYRQ", "WLRQ"):
        i = seq.find(motif)
        if i >= 0:
            break
    if i < 0:
        return None
    # motif W is Kabat 36; 37 is the next residue, 44/45/47 follow the
    # 'APGK...' stretch: use offsets from the motif start.
    try:
        h37 = seq[i + 1]
        h44, h45 = seq[i + 8], seq[i + 9]
        h47 = seq[i + 11]
    except IndexError:
        return None
    camelid = (h37 in "FYL", h44 == "E", h45 == "R", h47 in "GLF")
    return {"K37": h37, "K44": h44, "K45": h45, "K47": h47,
            "camelid_hallmarks": sum(camelid), "of": 4}


def main():
    seqs = read_fasta(data_path("sequences.fasta"))
    meta = read_metadata()
    germline = seqs["HumanVH3_23_germline"]

    rows = []
    for sid, seq in seqs.items():
        bad = sorted(set(seq) - CANONICAL)
        hm = vhh_hallmarks(seq) if len(seq) > 90 else None
        row = {
            "id": sid,
            "role": meta[sid]["role"],
            "length": len(seq),
            "source": meta[sid]["source"],
            "qc_canonical_only": not bad,
            "non_canonical": ",".join(bad) or "-",
            "ambiguous_codes": ",".join(sorted(set(seq) & AMBIGUOUS)) or "-",
            "mw_kda": round(molecular_weight(seq) / 1000.0, 2),
            "pI": isoelectric_point(seq),
            "pct_id_human_IGHV3_23": round(global_identity(seq, germline), 1),
            "vhh_hallmarks": (f"{hm['K37']}{hm['K44']}{hm['K45']}{hm['K47']} "
                              f"({hm['camelid_hallmarks']}/4 camelid)") if hm else "n/a",
        }
        rows.append(row)

    out = results_path("m1_sequence_qc.tsv")
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    with open(results_path("m1_sequence_qc.json"), "w") as f:
        json.dump(rows, f, indent=2)

    hdr = f"{'id':24s} {'role':22s} {'len':>4s} {'kDa':>6s} {'pI':>5s} {'%id VH3-23':>10s}  hallmarks"
    print(hdr)
    print("-" * len(hdr))
    for r in rows:
        print(f"{r['id']:24s} {r['role']:22s} {r['length']:4d} {r['mw_kda']:6.2f} "
              f"{r['pI']:5.2f} {r['pct_id_human_IGHV3_23']:9.1f}%  {r['vhh_hallmarks']}")
    fails = [r["id"] for r in rows if not r["qc_canonical_only"]]
    print(f"\nQC: {len(rows)-len(fails)}/{len(rows)} sequences pass composition checks"
          + (f"; failures: {fails}" if fails else ""))


if __name__ == "__main__":
    main()
