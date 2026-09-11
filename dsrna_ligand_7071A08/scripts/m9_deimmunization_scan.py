#!/usr/bin/env python3
"""
M9 - Anchor-position deimmunisation scan of the dominant foreign cluster
     (optional module; only meaningful if the ligand can be re-engineered).

For the highest-risk foreign epitope in the test article, every one of the 19
substitutions is placed at each MHC-II anchor pocket position (P1, P4, P6, P9)
of the binding core, the mutated 15-mer is re-scored across the full panel, and
the variants are ranked by how much population-weighted presentation they
remove.

Two caveats are enforced in the output rather than buried in the discussion:
  * BLOSUM62 score of the substitution is reported, because a change that
    abolishes DR binding but breaks the fold is not a design.
  * Whether the substituted residue is the one human germline IGHV3-23 carries
    at the aligned position is reported, because germline-matching
    substitutions are the ones least likely to create a *new* epitope or
    destabilise the framework.

This does not model the effect on ligand-target affinity. Any candidate has to
go back through a binding/stability screen before it means anything.
"""
import csv
import json
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (load_config, read_fasta, data_path, results_path,  # noqa: E402
                    iedb_mhcii, CoverageModel, population_weights)

ANCHORS = [1, 4, 6, 9]     # 1-based positions in the 9-mer binding core
AA = "ACDEFGHIKLMNPQRSTVWY"


def main():
    cfg = load_config()
    pcfg = cfg["prediction"]
    seqs = read_fasta(data_path("sequences.fasta"))
    with open(data_path("drb1_allele_frequencies.json")) as f:
        tables = json.load(f)
    presenting_fraction = CoverageModel(tables, population_weights(cfg)).weighted
    with open(results_path("m2_panel_alleles.txt")) as f:
        panel = [l.strip() for l in f if l.strip()]

    # ---- pick the target: worst foreign cluster of the test article -------
    test_id = None
    with open(results_path("m5_ligand_summary.tsv")) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["role"] == "test_article":
                test_id = r["id"]
    target = None
    with open(results_path("m5_clusters.tsv")) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["id"] == test_id and r["tolerance_class"] != "all_tolerised":
                if target is None or float(r["pop_presenting"]) > float(target["pop_presenting"]):
                    target = r
    if target is None:
        sys.exit("no foreign cluster in the test article - nothing to deimmunise")

    seq = seqs[test_id]
    peptide = target["peak_peptide"]
    core = target["peak_core"]
    off = peptide.find(core)
    pep_start = seq.find(peptide)
    print(f"target: {test_id} cluster {target['start']}-{target['end']}  "
          f"peptide {peptide}  core {core}  "
          f"pop_presenting {float(target['pop_presenting'])*100:.1f}%")

    # ---- build variants ---------------------------------------------------
    variants = {"WT": peptide}
    meta = {"WT": ("-", 0, "-", core)}
    for p in ANCHORS:
        wt_aa = core[p - 1]
        for aa in AA:
            if aa == wt_aa:
                continue
            new_core = core[:p - 1] + aa + core[p:]
            new_pep = peptide[:off + p - 1] + aa + peptide[off + p:]
            name = f"P{p}_{wt_aa}{pep_start + off + p}{aa}"
            variants[name] = new_pep
            meta[name] = (f"{wt_aa}->{aa}", p, wt_aa, new_core)

    # BLOSUM62 for the conservativeness column
    from Bio.Align import substitution_matrices
    bl = substitution_matrices.load("BLOSUM62")

    # germline residue at the aligned position, for the "germline-matching" column
    from Bio import Align
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = bl
    aligner.open_gap_score, aligner.extend_gap_score, aligner.mode = -11, -1, "global"
    germ = seqs["HumanVH3_23_germline"]
    aln = aligner.align(seq, germ)[0]
    gmap, i, j = {}, 0, 0
    for a, b in zip(aln[0], aln[1]):
        if a != "-" and b != "-":
            gmap[i] = b
        if a != "-":
            i += 1
        if b != "-":
            j += 1

    # ---- score ------------------------------------------------------------
    names = list(variants)
    print(f"scoring {len(names)} peptides x {len(panel)} DR molecules")
    from concurrent.futures import ThreadPoolExecutor

    def chunk(xs, n):
        for i in range(0, len(xs), n):
            yield xs[i:i + n]

    jobs = [(blk, dict((n, variants[n]) for n in nblk))
            for blk in chunk(panel, pcfg["alleles_per_request"])
            for nblk in chunk(names, 40)]

    def run(job):
        blk, subset = job
        rows = iedb_mhcii(pcfg["endpoint"], "netmhciipan_el", subset, blk, 15)
        print(f"  ok {len(subset)} peptides x {len(blk)} alleles", flush=True)
        return rows

    per_variant = defaultdict(dict)
    with ThreadPoolExecutor(max_workers=4) as ex:
        for rows in ex.map(run, jobs):
            for r in rows:
                rank = float(r["rank"])
                prev = per_variant[r["id"]].get(r["allele"], 9e9)
                per_variant[r["id"]][r["allele"]] = min(prev, rank)

    sb = pcfg["sb_rank"]
    out = []
    for name in names:
        hits = [a for a, rk in per_variant[name].items() if rk < sb]
        pop = presenting_fraction(hits)
        sub, p, wt_aa, new_core = meta[name]
        abs_pos = pep_start + off + p - 1 if name != "WT" else None
        out.append({
            "variant": name, "substitution": sub, "core": new_core,
            "n_sb_alleles": len(hits),
            "pop_presenting": round(pop, 4),
            "blosum62": int(bl[wt_aa][name[-1]]) if name != "WT" else "-",
            "germline_residue": gmap.get(abs_pos, "-") if abs_pos is not None else "-",
            "germline_match": (name[-1] == gmap.get(abs_pos)) if abs_pos is not None else "-",
        })
    wt = next(r for r in out if r["variant"] == "WT")
    for r in out:
        r["delta_pop_presenting"] = round(r["pop_presenting"] - wt["pop_presenting"], 4)
        r["delta_sb"] = r["n_sb_alleles"] - wt["n_sb_alleles"]
    out.sort(key=lambda r: (r["pop_presenting"], -(r["blosum62"] if isinstance(r["blosum62"], int) else 0)))

    with open(results_path("m9_deimmunization_scan.tsv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0]), delimiter="\t")
        w.writeheader()
        w.writerows(out)

    print(f"\nWT: {wt['n_sb_alleles']} SB alleles, "
          f"{wt['pop_presenting']*100:.1f}% presenting\n")
    print(f"{'variant':16s} {'sub':8s} {'core':10s} {'SB':>3s} {'pop%':>7s} "
          f"{'dPop':>7s} {'BL62':>5s}  germline")
    print("-" * 74)
    for r in out[:15]:
        if r["variant"] == "WT":
            continue
        print(f"{r['variant']:16s} {r['substitution']:8s} {r['core']:10s} "
              f"{r['n_sb_alleles']:3d} {r['pop_presenting']*100:6.1f}% "
              f"{r['delta_pop_presenting']*100:+6.1f}% {r['blosum62']:5} "
              f" {r['germline_residue']} {'(match)' if r['germline_match'] is True else ''}")


if __name__ == "__main__":
    main()
