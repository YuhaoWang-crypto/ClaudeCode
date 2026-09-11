#!/usr/bin/env python3
"""
M7 - B-cell / ADA layer.

The measured endpoint for a leached affinity ligand is an anti-drug *antibody*
assay. Antibodies are made by B cells, and a B cell only class-switches and
affinity-matures with help from a CD4 T cell that recognises a peptide from
the *same* protein. So the peptides that matter most are the ones where a
predicted HLA-DR epitope and a predicted B-cell epitope sit in the same region
of the fold: T-help and B-cell recognition co-localised.

A DR-only screen cannot see this. This module adds:
  * BepiPred-2.0 per-residue linear B-cell epitope propensity (IEDB REST),
  * contiguous B-cell epitope regions, and
  * T/B coincidence regions - a foreign HLA-DR cluster from M5 overlapping a
    B-cell region - which are the ones worth carrying into a wet-lab ADA or
    PBMC assay first.

Linear B-cell prediction is the weakest model in this pipeline (most real ADA
epitopes are conformational). Treat its output as a prioritisation aid, not a
risk number, and never as a standalone claim.
"""
import csv
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (load_config, read_fasta, data_path, results_path,  # noqa: E402
                    iedb_bcell)

THRESHOLD = 0.5      # BepiPred-2.0 default epitope threshold
MIN_LEN = 6          # minimum contiguous length to call a region


def regions(scores, threshold=THRESHOLD, min_len=MIN_LEN):
    out, start = [], None
    for i, s in enumerate(scores):
        if s >= threshold and start is None:
            start = i
        elif s < threshold and start is not None:
            if i - start >= min_len:
                out.append((start, i - 1))
            start = None
    if start is not None and len(scores) - start >= min_len:
        out.append((start, len(scores) - 1))
    return out


def main():
    cfg = load_config()
    seqs = read_fasta(data_path("sequences.fasta"))

    clusters = []
    with open(results_path("m5_clusters.tsv")) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["tolerance_class"] != "all_tolerised":
                clusters.append(r)

    per_residue, region_rows, coincidence = {}, [], []
    for sid, seq in seqs.items():
        rows = iedb_bcell(cfg["prediction"]["bcell_endpoint"],
                          cfg["prediction"]["bcell_method"], seq)
        scores = [s for _, _, s in rows]
        per_residue[sid] = scores
        print(f"  {sid:24s} {len(scores):4d} residues scored", flush=True)

        regs = regions(scores)
        for a, b in regs:
            region_rows.append({"id": sid, "start": a + 1, "end": b + 1,
                                "length": b - a + 1, "peptide": seq[a:b + 1],
                                "mean_score": round(sum(scores[a:b + 1]) / (b - a + 1), 3)})
        for c in clusters:
            if c["id"] != sid:
                continue
            cs, ce = int(c["start"]), int(c["end"])
            for a, b in regs:
                ov = min(ce, b + 1) - max(cs, a + 1) + 1
                if ov > 0:
                    coincidence.append({
                        "id": sid,
                        "t_cluster": f"{cs}-{ce}",
                        "t_peak_core": c["peak_core"],
                        "t_pop_presenting": c["pop_presenting"],
                        "b_region": f"{a+1}-{b+1}",
                        "overlap_aa": ov,
                        "region_peptide": seq[a:b + 1],
                    })

    with open(results_path("m7_bcell_per_residue.json"), "w") as f:
        json.dump(per_residue, f)
    for name, data in (("m7_bcell_regions.tsv", region_rows),
                       ("m7_tb_coincidence.tsv", coincidence)):
        with open(results_path(name), "w", newline="") as f:
            if not data:
                f.write("(none)\n")
                continue
            w = csv.DictWriter(f, fieldnames=list(data[0]), delimiter="\t")
            w.writeheader()
            w.writerows(data)

    print(f"\n{len(region_rows)} linear B-cell regions, "
          f"{len(coincidence)} T/B coincidence windows")
    for c in sorted(coincidence, key=lambda c: -float(c["t_pop_presenting"]))[:12]:
        print(f"  {c['id']:24s} T {c['t_cluster']:>9s} ({c['t_peak_core']}) "
              f"x B {c['b_region']:>9s}  overlap {c['overlap_aa']:2d} aa  "
              f"pop {float(c['t_pop_presenting'])*100:5.1f}%")


if __name__ == "__main__":
    main()
