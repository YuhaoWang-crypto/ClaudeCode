#!/usr/bin/env python3
"""
M6 - Benchmark calibration and system-suitability check.

An absolute pIRS has no meaning on its own; nobody can say whether 3.1 is a
lot. The previous report ran the ligand alone and then reasoned about the raw
counts, which is how "13 strong binders" ends up sounding alarming without any
way to know it is unremarkable for a 126-residue non-human protein.

So the ligand is run *in a batch with controls*, exactly as a bioassay would
be, and the batch is only reportable if the controls behave:

  system suitability
    positive controls  (tetanus toxin p2 / p30 universal epitope regions)
        must show high promiscuity - the panel and thresholds can detect a
        known promiscuous DR epitope.
    self negative controls (human germline VH3-23, human serum albumin)
        must fall to a low score once the M4 tolerance filter is applied -
        the filter is doing its job and not merely deleting everything.
    tolerance-filter discrimination
        the self controls must drop *more* than the foreign ligands do.

  calibration
    every ligand is expressed as a fold-change against ProteinA_Z, the ligand
    unit of alkali-stable rProtein A resins, whose leachate has decades of
    clinical exposure at controlled levels. That is the only defensible
    reference point available without proprietary clinical data.

Risk bands are assigned on the calibrated fold-change, and are an internal
triage convention - they are not a regulatory classification.
"""
import csv
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import load_config, results_path  # noqa: E402

BANDS = [(1.0, "comparable-to-benchmark"),
         (2.0, "modestly-elevated"),
         (4.0, "elevated"),
         (1e18, "high")]


def band(fold):
    for lim, name in BANDS:
        if fold <= lim:
            return name
    return BANDS[-1][1]


def main():
    cfg = load_config()
    rows = {}
    with open(results_path("m5_ligand_summary.tsv")) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            for k in ("pIRS", "pIRS_no_tolerance_filter", "pop_at_risk"):
                r[k] = float(r[k])
            for k in ("length", "n_epitopes", "n_foreign_epitopes",
                      "n_tolerised_epitopes", "max_promiscuity"):
                r[k] = int(r[k])
            rows[r["id"]] = r

    b = cfg["benchmarks"]
    anchor = rows[b["anchor_low"]]
    pos = [r for r in rows.values() if r["role"] == "positive_control"]
    selfs = [r for r in rows.values() if r["role"] == "negative_control_self"]
    ligands = [r for r in rows.values()
               if r["role"] in ("test_article", "benchmark_ligand",
                                "clinical_anchor", "class_comparator")]

    def drop(r):
        raw = r["pIRS_no_tolerance_filter"]
        return 0.0 if raw <= 0 else 100.0 * (raw - r["pIRS"]) / raw

    checks = []
    checks.append({
        "check": "positive controls promiscuous",
        "detail": "; ".join(f"{r['id']} max breadth {r['max_promiscuity']}/24 DR" for r in pos),
        "pass": all(r["max_promiscuity"] >= 6 for r in pos),
    })
    checks.append({
        "check": "self controls suppressed by tolerance filter",
        "detail": "; ".join(f"{r['id']} pIRS {r['pIRS']:.2f} (raw {r['pIRS_no_tolerance_filter']:.2f})"
                            for r in selfs),
        "pass": all(r["pIRS"] <= 0.35 * anchor["pIRS"] for r in selfs),
    })
    mean_self_drop = sum(drop(r) for r in selfs) / max(len(selfs), 1)
    mean_lig_drop = sum(drop(r) for r in ligands) / max(len(ligands), 1)
    checks.append({
        "check": "tolerance filter discriminates self from foreign",
        "detail": f"mean score drop: self {mean_self_drop:.1f}% vs ligands {mean_lig_drop:.1f}%",
        "pass": mean_self_drop > mean_lig_drop + 10,
    })
    checks.append({
        "check": "benchmark anchor scored",
        "detail": f"{anchor['id']} pIRS {anchor['pIRS']:.2f}",
        "pass": anchor["pIRS"] > 0,
    })

    out_rows = []
    for r in sorted(rows.values(), key=lambda r: -r["pIRS"]):
        fold = r["pIRS"] / max(anchor["pIRS"], 1e-9)
        out_rows.append({
            "id": r["id"], "role": r["role"], "length": r["length"],
            "pIRS": round(r["pIRS"], 2),
            "pIRS_raw": round(r["pIRS_no_tolerance_filter"], 2),
            "tolerance_drop_pct": round(drop(r), 1),
            "fold_vs_ProteinA_Z": round(fold, 2),
            "risk_band": band(fold) if r["role"] in
                         ("test_article", "benchmark_ligand", "clinical_anchor",
                          "class_comparator") else "n/a (control)",
            "n_foreign_epitopes": r["n_foreign_epitopes"],
            "max_promiscuity": r["max_promiscuity"],
            "pop_at_risk_pct": round(r["pop_at_risk"] * 100, 1),
        })

    with open(results_path("m6_calibrated_ranking.tsv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out_rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)
    with open(results_path("m6_system_suitability.json"), "w") as f:
        json.dump({"checks": checks, "batch_valid": all(c["pass"] for c in checks)},
                  f, indent=2)

    print("system suitability")
    for c in checks:
        print(f"  [{'PASS' if c['pass'] else 'FAIL'}] {c['check']}")
        print(f"         {c['detail']}")
    print(f"\n  batch reportable: {all(c['pass'] for c in checks)}\n")

    hdr = (f"{'id':24s} {'role':22s} {'pIRS':>7s} {'raw':>7s} {'drop%':>6s} "
           f"{'xPrA-Z':>7s} {'pop@risk':>9s}  band")
    print(hdr); print("-" * (len(hdr) + 10))
    for r in out_rows:
        print(f"{r['id']:24s} {r['role']:22s} {r['pIRS']:7.2f} {r['pIRS_raw']:7.2f} "
              f"{r['tolerance_drop_pct']:6.1f} {r['fold_vs_ProteinA_Z']:7.2f} "
              f"{r['pop_at_risk_pct']:8.1f}%  {r['risk_band']}")


if __name__ == "__main__":
    main()
