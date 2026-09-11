#!/usr/bin/env python3
"""
M8 - Exposure context: turn a sequence-level score into a dose-level risk call.

A leached affinity ligand is a process-related impurity, and impurity risk is
the product of *intrinsic* immunogenic potential (M1-M7) and *how much of it a
patient actually receives*. A ligand with a high pIRS delivered at 20 ng per
dose is a different regulatory conversation from the same ligand at 100 ug.

This module builds the exposure grid the risk call is made on:

    ug ligand per dose = leachate (ng ligand / mg product) x dose (mg) / 1000

and bands it against the exposure regime in which rProtein A leachate has an
established clinical safety record. The bands are a documented internal
convention for triage, not a regulatory threshold - no agency publishes a
numeric leachate immunogenicity limit, and the ICH Q6B / EMA expectation is
that leachate is controlled to a justified, consistently achieved level.
"""
import csv
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import load_config, results_path  # noqa: E402

# Triage bands on ug ligand delivered per dose.
BANDS = [
    (0.1,   "negligible", "Below the exposure at which protein impurities have "
                          "shown measurable ADA induction in the clinic"),
    (1.0,   "low",        "Comparable to well-controlled rProtein A leachate "
                          "exposure in approved mAb products"),
    (10.0,  "moderate",   "Above typical qualified leachate exposure; justify "
                          "with ligand-specific data"),
    (1e18,  "elevated",   "Immunologically meaningful protein dose; treat the "
                          "ligand as a co-administered antigen"),
]


def band(ug):
    for limit, name, note in BANDS:
        if ug < limit:
            return name, note
    return BANDS[-1][1], BANDS[-1][2]


def main():
    cfg = load_config()
    ex = cfg["exposure"]

    grid = []
    for ppm in ex["leachate_ppm"]:
        for dose in ex["dose_mg"]:
            ug = ppm * dose / 1000.0
            nmol = ug / (ex["ligand_mw_kda"] * 1000.0) * 1e3   # nmol
            name, note = band(ug)
            grid.append({"leachate_ng_per_mg": ppm, "dose_mg": dose,
                         "ug_ligand_per_dose": round(ug, 4),
                         "nmol_per_dose": round(nmol, 5),
                         "exposure_band": name, "basis": note})

    with open(results_path("m8_exposure_grid.tsv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(grid[0]), delimiter="\t")
        w.writeheader()
        w.writerows(grid)

    # Combine with the intrinsic score for the test article.
    summary = {}
    with open(results_path("m5_ligand_summary.tsv")) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            summary[r["id"]] = r

    test = [r for r in summary.values() if r["role"] == "test_article"]
    anchor = summary.get(cfg["benchmarks"]["anchor_low"])
    out = {"exposure_grid": grid}
    if test and anchor:
        t = test[0]
        ratio = float(t["pIRS"]) / max(float(anchor["pIRS"]), 1e-9)
        out["intrinsic"] = {
            "test_article": t["id"],
            "pIRS": float(t["pIRS"]),
            "benchmark": anchor["id"],
            "benchmark_pIRS": float(anchor["pIRS"]),
            "fold_vs_benchmark": round(ratio, 2),
        }
    with open(results_path("m8_exposure_context.json"), "w") as f:
        json.dump(out, f, indent=2)

    print(f"{'ng/mg':>6s} {'dose mg':>8s} {'ug/dose':>10s} {'nmol/dose':>10s}  band")
    print("-" * 56)
    for g in grid:
        print(f"{g['leachate_ng_per_mg']:6d} {g['dose_mg']:8d} "
              f"{g['ug_ligand_per_dose']:10.4f} {g['nmol_per_dose']:10.5f}  "
              f"{g['exposure_band']}")
    if "intrinsic" in out:
        i = out["intrinsic"]
        print(f"\nintrinsic: {i['test_article']} pIRS {i['pIRS']} = "
              f"{i['fold_vs_benchmark']}x {i['benchmark']} ({i['benchmark_pIRS']})")


if __name__ == "__main__":
    main()
