#!/usr/bin/env python3
"""
M4 - Self / pre-existing-tolerance filter.

The single largest source of false positives in a naive HLA-DR screen of an
antibody-derived ligand: most predicted binders sit in *framework* regions
whose 9-mer cores also occur in the human proteome - in human germline V
domains above all. Those cores are seen by a repertoire that has been
negatively selected against them, so counting them as immunogenic risk
inflates the score of every VHH, scFv and Fab-derived ligand equally and
destroys the ability to rank them.

Two screens, both against UniProt Swiss-Prot Homo sapiens:

  exact       the predicted 9-mer core occurs verbatim in a human protein
              -> treat as tolerised (weight 0 by default)
  tcr_face    the residues that face the TCR (P2, P3, P5, P7, P8 of the core)
              match a human 9-mer at the same positions -> a cross-reactive
              human-like TCR face, the JanusMatrix concept in simplified form
              -> heavily down-weighted (weight 0.35 by default)

This is a *screen*, not the published JanusMatrix algorithm: it does not
require the human counterpart to bind the same allele. It is therefore
conservative in the direction of calling more peptides tolerised, and every
flagged core is written out with its human hit so the call can be checked.

Output: results/m4_core_tolerance.tsv (one row per distinct predicted core).
"""
import csv
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import load_config, data_path, results_path  # noqa: E402


def iter_proteome(path):
    name, buf = None, []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if name:
                    yield name, "".join(buf)
                name, buf = line[1:].strip().split()[0], []
            else:
                buf.append(line.strip())
    if name:
        yield name, "".join(buf)


def main():
    cfg = load_config()
    faces = cfg["tolerance_filter"]["tcr_face_positions"]
    proteome = os.path.join(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))), cfg["tolerance_filter"]["proteome"])

    # ---- collect every distinct predicted core (any EL/WB-or-better call) ---
    cores = defaultdict(set)          # core -> {sequence ids it came from}
    with open(results_path("m3_binding_long.tsv")) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["call_el"] in ("SB", "WB"):
                cores[row["core"]].add(row["id"])
    cores = {c: v for c, v in cores.items() if len(c) == 9}
    print(f"{len(cores)} distinct 9-mer binding cores to screen")

    def face(k):
        return "".join(k[p - 1] for p in faces)

    face_index = defaultdict(set)
    for c in cores:
        face_index[face(c)].add(c)

    exact_hit = {}      # core -> human protein accession
    face_hit = {}       # core -> (human protein, human 9mer)
    n_prot = n_res = 0
    for acc, seq in iter_proteome(proteome):
        n_prot += 1
        n_res += len(seq)
        for i in range(len(seq) - 8):
            k = seq[i:i + 9]
            if k in cores and k not in exact_hit:
                exact_hit[k] = acc
            fk = face(k)
            if fk in face_index:
                for c in face_index[fk]:
                    if c not in face_hit:
                        face_hit[c] = (acc, k)
    print(f"screened {n_prot} human proteins / {n_res:,} residues")

    w_exact = cfg["tolerance_filter"]["discount_exact"]
    w_face = cfg["tolerance_filter"]["discount_tcrface"]

    out = results_path("m4_core_tolerance.tsv")
    n_e = n_f = 0
    with open(out, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["core", "from_sequences", "tolerance_class",
                    "human_hit_protein", "human_hit_9mer", "weight"])
        for c in sorted(cores):
            if c in exact_hit:
                cls, prot, hk, wt = "exact_human_9mer", exact_hit[c], c, w_exact
                n_e += 1
            elif c in face_hit:
                cls, (prot, hk), wt = "human_tcr_face", face_hit[c], w_face
                n_f += 1
            else:
                cls, prot, hk, wt = "foreign", "-", "-", 1.0
            w.writerow([c, ",".join(sorted(cores[c])), cls, prot, hk, wt])

    print(f"\n  exact human 9-mer core : {n_e:4d}  ({100*n_e/len(cores):.1f}%)  weight {w_exact}")
    print(f"  human TCR-face match   : {n_f:4d}  ({100*n_f/len(cores):.1f}%)  weight {w_face}")
    print(f"  foreign                : {len(cores)-n_e-n_f:4d}  "
          f"({100*(len(cores)-n_e-n_f)/len(cores):.1f}%)  weight 1.0")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
