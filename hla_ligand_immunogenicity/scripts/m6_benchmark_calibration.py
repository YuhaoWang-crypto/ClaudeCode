#!/usr/bin/env python3
"""
M6 - Benchmark calibration and system suitability.

An absolute pIRS has no meaning on its own; nobody can say whether 3.1 is a
lot. A ligand run alone produces exactly that kind of number, and "13 strong
binders" then sounds alarming without any way to know it is unremarkable for a
126-residue non-human protein.

So the ligand is run *in a batch with controls*, as a bioassay would be, and
the controls do two different jobs:

  CALIBRATION (what the numbers mean)
    ProteinA_Z - the ligand unit of alkali-stable rProtein A resins, whose
        leachate has decades of controlled clinical exposure - sets the
        intrinsic-risk reference. Every ligand is expressed as a fold-change
        against it. That is the only defensible anchor available without
        proprietary clinical data.
    The tetanus toxin universal T-helper epitopes p2 (830-844) and p30
        (947-967) - peptides that drive CD4 responses in most donors in vitro -
        set the *promiscuity* reference. Whatever DR breadth the predictor
        assigns them is the practical ceiling of the method, and a test-article
        epitope can then be reported as more or less promiscuous than a
        textbook universal epitope rather than as a bare count.

  SUITABILITY (whether the batch is reportable)
    the universal epitopes must be detected with multi-allele breadth at the
        weak-binder tier;
    the human self controls must collapse once the M4 tolerance filter is
        applied; and
    the filter must drop the self controls further than it drops the foreign
        ligands, i.e. it discriminates rather than deletes.

Risk bands are assigned on the calibrated fold-change and are an internal
triage convention, not a regulatory classification.
"""
import csv
import json
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import load_config, read_fasta, data_path, results_path  # noqa: E402

BANDS = [(1.0, "comparable-to-benchmark"),
         (2.0, "modestly-elevated"),
         (4.0, "elevated"),
         (1e18, "high")]


def band(fold):
    for lim, name in BANDS:
        if fold <= lim:
            return name
    return BANDS[-1][1]


def epitope_breadth(binding, sid, window, sb, wb):
    """Best EL %Rank per DR molecule over the 15-mers inside `window`."""
    frames = {window[i:i + 15] for i in range(max(len(window) - 14, 1))}
    best = {}
    for r in binding:
        if r["id"] != sid or r["peptide"] not in frames or r["el_rank"] in ("", "None"):
            continue
        v = float(r["el_rank"])
        if v < best.get(r["allele"], 9e9):
            best[r["allele"]] = v
    return {
        "n_alleles": len(best),
        "n_sb": sum(1 for v in best.values() if v < sb),
        "n_wb": sum(1 for v in best.values() if v < wb),
        "best_rank": round(min(best.values()), 3) if best else None,
    }


def main():
    cfg = load_config()
    b = cfg["benchmarks"]
    sb_t, wb_t = cfg["prediction"]["sb_rank"], cfg["prediction"]["wb_rank"]

    rows = {}
    with open(results_path("m5_ligand_summary.tsv")) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            for k in ("pIRS", "pIRS_no_tolerance_filter", "pop_at_risk"):
                r[k] = float(r[k])
            for k in ("length", "n_epitopes", "n_foreign_epitopes",
                      "n_tolerised_epitopes", "max_promiscuity"):
                r[k] = int(r[k])
            rows[r["id"]] = r
    with open(results_path("m3_binding_long.tsv")) as f:
        binding = list(csv.DictReader(f, delimiter="\t"))

    anchor = rows[b["anchor_low"]]
    selfs = [r for r in rows.values() if r["role"] == "negative_control_self"]
    ligands = [r for r in rows.values()
               if r["role"] in ("test_article", "benchmark_ligand",
                                "clinical_anchor", "class_comparator")]
    test = next(r for r in rows.values() if r["role"] == "test_article")

    # ---- promiscuity calibration against the universal epitopes -----------
    universal = {}
    for sid, window in b["positive_control_epitopes"].items():
        universal[sid] = epitope_breadth(binding, sid, window, sb_t, wb_t)
    ceiling_sb = max(u["n_sb"] for u in universal.values())
    ceiling_wb = max(u["n_wb"] for u in universal.values())
    n_panel = universal[list(universal)[0]]["n_alleles"]

    calibration = {
        "universal_epitopes": universal,
        "panel_size": n_panel,
        "sb_threshold": sb_t, "wb_threshold": wb_t,
        "universal_epitope_ceiling_sb": ceiling_sb,
        "universal_epitope_ceiling_wb": ceiling_wb,
        "test_article": test["id"],
        "test_peak_promiscuity_sb": test["max_promiscuity"],
        "test_peak_vs_universal_sb": (round(test["max_promiscuity"] / ceiling_sb, 2)
                                      if ceiling_sb else None),
        "interpretation": (
            f"Experimentally universal T-helper epitopes reach {ceiling_sb}/{n_panel} DR "
            f"molecules at EL %Rank < {sb_t:g} and {ceiling_wb}/{n_panel} at < {wb_t:g} on this "
            f"panel. The strong-binder tier is therefore a high-specificity, low-sensitivity "
            f"criterion: it does not reproduce the textbook universal epitopes, so a peptide that "
            f"does clear it across many molecules is more promiscuous than they are, not merely "
            f"'a strong binder'."),
    }
    with open(results_path("m6_promiscuity_calibration.json"), "w") as f:
        json.dump(calibration, f, indent=2)

    # ---- system suitability ----------------------------------------------
    def drop(r):
        raw = r["pIRS_no_tolerance_filter"]
        return 0.0 if raw <= 0 else 100.0 * (raw - r["pIRS"]) / raw

    min_wb = b["positive_control_min_breadth_wb"]
    checks = [{
        "check": "universal epitopes detected with multi-allele breadth",
        "detail": "; ".join(f"{k} {v['n_wb']}/{v['n_alleles']} DR at %Rank<{wb_t:g} "
                            f"({v['n_sb']} at <{sb_t:g}, best {v['best_rank']})"
                            for k, v in universal.items()),
        "pass": all(v["n_wb"] >= min_wb for v in universal.values()),
    }, {
        "check": "self controls suppressed by tolerance filter",
        "detail": "; ".join(f"{r['id']} pIRS {r['pIRS']:.2f} "
                            f"(raw {r['pIRS_no_tolerance_filter']:.2f})" for r in selfs),
        "pass": all(r["pIRS"] <= 0.35 * anchor["pIRS"] for r in selfs),
    }]
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

    # ---- calibrated ranking ----------------------------------------------
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

    print("promiscuity calibration")
    for k, v in universal.items():
        print(f"  {k:16s} {v['n_sb']:2d}/{v['n_alleles']} DR at %Rank<{sb_t:g}   "
              f"{v['n_wb']:2d}/{v['n_alleles']} at <{wb_t:g}   best {v['best_rank']}")
    print(f"  test article peak promiscuity: {test['max_promiscuity']}/{n_panel} at "
          f"%Rank<{sb_t:g}  ->  {calibration['test_peak_vs_universal_sb']}x the universal-epitope "
          f"ceiling\n")

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
