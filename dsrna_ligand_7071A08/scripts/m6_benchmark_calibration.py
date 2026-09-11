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


MIN_OVERLAP = 9      # a 15-mer must cover a whole binding core of the epitope


def epitope_breadth(binding, sid, window, seq, sb, wb, panel_size):
    """
    Best EL %Rank per DR molecule over the scanned 15-mers that cover `window`.

    Matched by position overlap, not by string equality: a defined epitope is
    rarely exactly 15 residues (HA306-318 is 13), so requiring the scanned
    15-mer to equal the epitope silently matches nothing and reports a breadth
    of zero for the best-evidenced epitope in the set.
    """
    w0 = seq.find(window)
    if w0 < 0:
        return {"n_alleles": 0, "n_sb": 0, "n_wb": 0, "best_rank": None,
                "note": "epitope window not found in the sequence"}
    w0, w1 = w0 + 1, w0 + len(window)          # 1-based inclusive
    best = {}
    for r in binding:
        if r["id"] != sid or r["el_rank"] in ("", "None"):
            continue
        s, e = int(r["start"]), int(r["end"])
        if min(e, w1) - max(s, w0) + 1 < MIN_OVERLAP:
            continue
        v = float(r["el_rank"])
        if v < best.get(r["allele"], 9e9):
            best[r["allele"]] = v
    return {
        "n_alleles": panel_size,
        "n_scored": len(best),
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
    seqs = read_fasta(data_path("sequences.fasta"))
    with open(results_path("m2_panel_alleles.txt")) as f:
        n_panel = sum(1 for l in f if l.strip())

    anchor = rows[b["anchor_low"]]
    selfs = [r for r in rows.values() if r["role"] == "negative_control_self"]
    ligands = [r for r in rows.values()
               if r["role"] in ("test_article", "benchmark_ligand",
                                "clinical_anchor", "class_comparator")]
    test = next(r for r in rows.values() if r["role"] == "test_article")

    # ---- promiscuity calibration against the universal epitopes -----------
    universal = {}
    for sid, window in b["positive_control_epitopes"].items():
        universal[sid] = epitope_breadth(binding, sid, window, seqs[sid],
                                         sb_t, wb_t, n_panel)
    ceiling_sb = max(u["n_sb"] for u in universal.values())
    ceiling_wb = max(u["n_wb"] for u in universal.values())

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
        "universal_epitope_range_sb": [min(u["n_sb"] for u in universal.values()),
                                       ceiling_sb],
        "universal_epitope_range_wb": [min(u["n_wb"] for u in universal.values()),
                                       ceiling_wb],
        "interpretation": (
            f"Experimentally universal T-helper epitopes reach "
            f"{min(u['n_sb'] for u in universal.values())}-{ceiling_sb} of {n_panel} DR molecules "
            f"at EL %Rank < {sb_t:g}, and "
            f"{min(u['n_wb'] for u in universal.values())}-{ceiling_wb} at < {wb_t:g}. The test "
            f"article's dominant core reaches {test['max_promiscuity']}/{n_panel} at the "
            f"strong-binder tier - "
            + ("above" if test["max_promiscuity"] > ceiling_sb else
               "at the top of" if test["max_promiscuity"] == ceiling_sb else
               "within") +
            f" that range. Read the scale in both directions: the strong-binder tier recovers only "
            f"{100*ceiling_sb/n_panel:.0f}% of the DR molecules a universal epitope is known to be "
            f"presented by, so peptides below it are unflagged rather than cleared."),
    }
    with open(results_path("m6_promiscuity_calibration.json"), "w") as f:
        json.dump(calibration, f, indent=2)

    # ---- boundary controls: where the tolerance filter's assumption breaks --
    # These do not pass or fail the batch. They measure two things the filter
    # cannot do, so the limitation is a number in the output rather than a
    # sentence in the discussion.
    tol = {}
    with open(results_path("m4_core_tolerance.tsv")) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            tol[r["core"]] = r
    boundary = {}
    for sid, window in b.get("boundary_control_epitopes", {}).items():
        if sid not in rows:
            continue
        br = epitope_breadth(binding, sid, window, seqs[sid], sb_t, wb_t, n_panel)
        # which predicted cores inside this window survived the filter
        cores_in = [c for c, r in tol.items()
                    if sid in r["from_sequences"].split(",") and c in
                    {window[i:i + 9] for i in range(max(len(window) - 8, 1))}]
        classes = [tol[c]["tolerance_class"] for c in cores_in]
        boundary[sid] = {
            "role": rows[sid]["role"],
            "epitope_window": window,
            "dr_breadth_sb": br["n_sb"], "dr_breadth_wb": br["n_wb"],
            "best_rank": br["best_rank"],
            "predicted_cores_in_window": len(cores_in),
            "cores_called_tolerised": sum(1 for c in classes if c != "foreign"),
            "cores_called_foreign": sum(1 for c in classes if c == "foreign"),
            "ligand_pIRS": rows[sid]["pIRS"],
            "ligand_pIRS_unfiltered": rows[sid]["pIRS_no_tolerance_filter"],
        }
    calibration["boundary_controls"] = boundary
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

    if boundary:
        print("boundary controls (measured limitations, not pass/fail)")
        for sid, v in boundary.items():
            print(f"  {sid:26s} {v['role']:26s} DR breadth {v['dr_breadth_sb']}/{n_panel} SB, "
                  f"{v['dr_breadth_wb']}/{n_panel} WB")
            print(f"  {'':26s} {v['cores_called_tolerised']}/{v['predicted_cores_in_window']} "
                  f"predicted cores in the epitope window were called tolerised; "
                  f"pIRS {v['ligand_pIRS']:.2f} (unfiltered {v['ligand_pIRS_unfiltered']:.2f})")
        print()

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
