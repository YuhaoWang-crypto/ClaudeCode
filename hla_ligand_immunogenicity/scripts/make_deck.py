#!/usr/bin/env python3
"""
Collect everything the slide deck needs into results/deck_data.json, then hand
off to scripts/make_deck.js (pptxgenjs) to render report.pptx.

Keeping the data assembly in Python means the deck and the HTML report read the
same numbers out of the same result files.
"""
import csv
import json
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (load_config, read_fasta, read_metadata, data_path,  # noqa: E402
                    results_path, figures_path, ROOT)


def tsv(name):
    p = results_path(name)
    if not os.path.exists(p):
        return []
    with open(p) as f:
        if f.readline().startswith("(none)"):
            return []
        f.seek(0)
        return list(csv.DictReader(f, delimiter="\t"))


def jsn(name):
    p = results_path(name)
    return json.load(open(p)) if os.path.exists(p) else {}


def main():
    cfg = load_config()
    seqs = read_fasta(data_path("sequences.fasta"))
    meta = read_metadata()
    panel = jsn("m2_panel.json")
    rank = tsv("m6_calibrated_ranking.tsv")
    suit = jsn("m6_system_suitability.json")
    val = jsn("m4_filter_validation.json")
    clusters = tsv("m5_clusters.tsv")
    epitopes = tsv("m5_epitopes.tsv")
    summary = {r["id"]: r for r in tsv("m5_ligand_summary.tsv")}
    coincide = tsv("m7_tb_coincidence.tsv")
    grid = tsv("m8_exposure_grid.tsv")
    deimm = tsv("m9_deimmunization_scan.tsv")
    binding = tsv("m3_binding_long.tsv")
    qc = tsv("m1_sequence_qc.tsv")

    test_id = next(r["id"] for r in summary.values() if r["role"] == "test_article")
    test_clusters = [c for c in clusters if c["id"] == test_id]
    foreign = [c for c in test_clusters if c["tolerance_class"] != "all_tolerised"]
    top = max(foreign, key=lambda c: float(c["pop_presenting"])) if foreign else None

    data = {
        "test_id": test_id,
        "test_length": len(seqs[test_id]),
        "test_source": meta[test_id]["source"],
        "test_seq": seqs[test_id],
        "panel": panel,
        "qc": qc,
        "rank": rank,
        "suitability": suit,
        "filter_validation": val,
        "summary": summary[test_id],
        "test_rank": next(r for r in rank if r["id"] == test_id),
        "anchor_rank": next(r for r in rank if r["id"] == cfg["benchmarks"]["anchor_low"]),
        "clusters": test_clusters,
        "top_cluster": top,
        "epitopes": sorted([e for e in epitopes if e["id"] == test_id],
                           key=lambda e: -int(e["n_sb_alleles"]))[:8],
        "coincidence": sorted([c for c in coincide if c["id"] == test_id],
                              key=lambda c: -float(c["t_pop_presenting"]))[:5],
        "exposure": grid,
        "deimm": deimm[:10],
        "deimm_wt": next((r for r in deimm if r["variant"] == "WT"), None),
        "n_el_sb": sum(1 for r in binding if r["call_el"] == "SB"),
        "n_cons_sb": sum(1 for r in binding if r["call_consensus"] == "SB"),
        "ba_confirm_rank": cfg["prediction"]["ba_confirm_rank"],
        "figures": {k: figures_path(v) for k, v in {
            "panel": "fig1_panel_coverage.png",
            "landscape": "fig2_binding_landscape.png",
            "ranking": "fig3_calibrated_ranking.png",
            "tb": "fig4_tb_coincidence.png",
            "deimm": "fig5_deimmunization.png"}.items()
            if os.path.exists(figures_path(v))},
        "out": os.path.join(ROOT, "report.pptx"),
    }
    p = results_path("deck_data.json")
    with open(p, "w") as f:
        json.dump(data, f, indent=1)
    print(f"wrote {p}")

    js = os.path.join(ROOT, "scripts", "make_deck.js")
    subprocess.run(["node", js, p], check=True, cwd=ROOT)


if __name__ == "__main__":
    main()
