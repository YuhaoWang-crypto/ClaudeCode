#!/usr/bin/env python3
"""
M5 - Epitope consolidation, promiscuity, and population-weighted risk scoring.

Counting strong binders is the wrong unit. A 15-mer scan reports the same
epitope up to seven times (every frame that contains the core), and it treats
a hit on DRB1*07:01 - carried by ~23% of Europeans - as equal to a hit on
DRB1*13:03, carried by ~1%. Two ligands with identical SB counts can differ
several-fold in the fraction of patients who could actually present them.

This module therefore:
  1. collapses overlapping 15-mers into *epitopes*, keyed on (sequence,
     binding core), keeping the best rank and its position;
  2. computes, for each epitope, the fraction of the weighted US/EU population
     carrying at least one DR molecule predicted to present it
     (Hardy-Weinberg phenotypic frequency over the presenting DRB1 set);
  3. applies the M4 tolerance weight so framework cores shared with the human
     proteome stop dominating the score;
  4. rolls epitopes into positional clusters; and
  5. emits one headline number per ligand:

     pIRS  = population-weighted immunogenic risk score
           = 100/L * SUM over foreign epitopes of
                     (weighted US/EU presenting fraction) * (tolerance weight)

     i.e. "population-weighted presentable foreign epitope content per 100
     residues". It is a *relative* scale - interpretable only against the
     benchmark ligands and controls run in the same batch (M6).

Also reported, and easier to explain to a non-specialist:
     pop_at_risk = fraction of the weighted US/EU population carrying at least
     one DR molecule predicted to present at least one *foreign* epitope of
     this ligand.
"""
import csv
import json
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (load_config, read_fasta, read_metadata, data_path,  # noqa: E402
                    results_path, CoverageModel, population_weights)

CLUSTER_GAP = 8      # residues; epitope cores closer than this join one cluster


def main():
    cfg = load_config()
    seqs = read_fasta(data_path("sequences.fasta"))
    meta = read_metadata()
    with open(data_path("drb1_allele_frequencies.json")) as f:
        tables = json.load(f)
    model = CoverageModel(tables, population_weights(cfg))
    presenting_fraction = model.weighted

    tol = {}
    with open(results_path("m4_core_tolerance.tsv")) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            tol[r["core"]] = (r["tolerance_class"], float(r["weight"]))

    # ---- collapse 15-mer frames into epitopes -----------------------------
    # epitope key: (id, core). value: best rank per allele + best position.
    ep = defaultdict(lambda: {"alleles_sb": {}, "alleles_wb": {},
                              "best_rank": 9e9, "pos": None, "peptide": None})
    with open(results_path("m3_binding_long.tsv")) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["call_consensus"] == "-":
                continue
            e = ep[(r["id"], r["core"])]
            rank = float(r["el_rank"])
            bucket = "alleles_sb" if r["call_consensus"] == "SB" else "alleles_wb"
            prev = e[bucket].get(r["allele"], 9e9)
            e[bucket][r["allele"]] = min(prev, rank)
            if rank < e["best_rank"]:
                e["best_rank"] = rank
                e["pos"] = int(r["start"]) + r["peptide"].find(r["core"])
                e["peptide"] = r["peptide"]

    # ---- per-epitope table ------------------------------------------------
    rows = []
    for (sid, core), e in ep.items():
        if not e["alleles_sb"]:
            continue                       # epitope = >=1 consensus strong binder
        cls, wt = tol.get(core, ("foreign", 1.0))
        sb = sorted(e["alleles_sb"])
        drb345 = [a for a in sb if not a.startswith("HLA-DRB1")]
        rows.append({
            "id": sid, "core": core, "pos": e["pos"], "peptide": e["peptide"],
            "best_el_rank": round(e["best_rank"], 3),
            "n_sb_alleles": len(sb),
            "n_sb_drb1": len(sb) - len(drb345),
            "sb_alleles": ";".join(a.replace("HLA-", "") for a in sb),
            "drb345_sb": ";".join(a.replace("HLA-", "") for a in drb345) or "-",
            "n_wb_alleles": len(e["alleles_wb"]),
            "tolerance_class": cls, "tolerance_weight": wt,
            "pop_presenting": round(presenting_fraction(sb), 4),
        })
    rows.sort(key=lambda r: (r["id"], r["pos"]))

    with open(results_path("m5_epitopes.tsv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    # ---- clusters ---------------------------------------------------------
    clusters = []
    by_seq = defaultdict(list)
    for r in rows:
        by_seq[r["id"]].append(r)
    for sid, rs in by_seq.items():
        rs = sorted(rs, key=lambda r: r["pos"])
        cur = [rs[0]]
        for r in rs[1:]:
            if r["pos"] - cur[-1]["pos"] <= CLUSTER_GAP:
                cur.append(r)
            else:
                clusters.append((sid, cur))
                cur = [r]
        clusters.append((sid, cur))

    crows = []
    for sid, cs in clusters:
        alleles = sorted({a for c in cs for a in c["sb_alleles"].split(";")})
        alleles_full = ["HLA-" + a for a in alleles]
        peak = min(cs, key=lambda c: c["best_el_rank"])
        foreign = [c for c in cs if c["tolerance_class"] == "foreign"]
        crows.append({
            "id": sid,
            "start": min(c["pos"] for c in cs),
            "end": max(c["pos"] for c in cs) + 8,
            "n_epitopes": len(cs),
            "n_foreign_epitopes": len(foreign),
            "peak_core": peak["core"], "peak_peptide": peak["peptide"],
            "peak_el_rank": peak["best_el_rank"],
            "max_allele_breadth": max(c["n_sb_alleles"] for c in cs),
            "union_sb_alleles": len(alleles),
            "pop_presenting": round(presenting_fraction(alleles_full), 4),
            "tolerance_class": ("all_tolerised" if not foreign else
                                "mixed" if len(foreign) < len(cs) else "foreign"),
        })
    crows.sort(key=lambda r: (r["id"], r["start"]))
    with open(results_path("m5_clusters.tsv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(crows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(crows)

    # ---- per-ligand summary ----------------------------------------------
    summary = []
    for sid, seq in seqs.items():
        rs = by_seq.get(sid, [])
        L = len(seq)
        foreign = [r for r in rs if r["tolerance_class"] == "foreign"]
        pirs = 100.0 / L * sum(r["pop_presenting"] * r["tolerance_weight"] for r in rs)
        pirs_raw = 100.0 / L * sum(r["pop_presenting"] for r in rs)
        union_foreign = sorted({("HLA-" + a) for r in foreign
                                for a in r["sb_alleles"].split(";")})
        summary.append({
            "id": sid, "role": meta[sid]["role"], "length": L,
            "n_epitopes": len(rs),
            "n_foreign_epitopes": len(foreign),
            "n_tolerised_epitopes": len(rs) - len(foreign),
            "max_promiscuity": max((r["n_sb_alleles"] for r in rs), default=0),
            # kept at 4 dp: downstream fold-changes divide these, and rounding
            # to 2 dp first moves the ratio in the second decimal
            "pIRS": round(pirs, 4),
            "pIRS_no_tolerance_filter": round(pirs_raw, 4),
            "pop_at_risk": round(presenting_fraction(union_foreign), 4),
            "n_clusters": sum(1 for c in crows if c["id"] == sid),
            "n_foreign_clusters": sum(1 for c in crows
                                      if c["id"] == sid and c["tolerance_class"] != "all_tolerised"),
        })
    summary.sort(key=lambda r: -r["pIRS"])
    with open(results_path("m5_ligand_summary.tsv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(summary[0]), delimiter="\t")
        w.writeheader()
        w.writerows(summary)

    print(f"{len(rows)} epitopes, {len(crows)} clusters\n")
    hdr = (f"{'id':24s} {'role':22s} {'ep':>3s} {'fgn':>4s} {'tol':>4s} "
           f"{'maxProm':>7s} {'pIRS':>7s} {'raw':>7s} {'pop@risk':>9s}")
    print(hdr); print("-" * len(hdr))
    for s in summary:
        print(f"{s['id']:24s} {s['role']:22s} {s['n_epitopes']:3d} "
              f"{s['n_foreign_epitopes']:4d} {s['n_tolerised_epitopes']:4d} "
              f"{s['max_promiscuity']:7d} {s['pIRS']:7.2f} "
              f"{s['pIRS_no_tolerance_filter']:7.2f} {s['pop_at_risk']*100:8.1f}%")


if __name__ == "__main__":
    main()
