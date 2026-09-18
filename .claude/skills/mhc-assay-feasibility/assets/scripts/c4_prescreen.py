#!/usr/bin/env python3
"""
C4 - The tetramer go/no-go pre-screen (RFP Step 2).

The RFP asks a CRO to custom-synthesise two HLA alleles and "8 uniblock
complexes" from 5 candidate 9-mers across HLA-A*32:01 and HLA-B*57:01. A
peptide-MHC monomer that does not fold does not become a tetramer, and the
failure is discovered after the money and the weeks are spent. This module
ranks the 10 possible (peptide, allele) complexes before any of that.

The decision is deliberately asymmetric. A complex wrongly built is money; a
complex wrongly skipped is a lost epitope nobody will go back for. So the rule
is: build anything the predictor does not confidently reject, and let the
confident rejections carry the calibrated NPV from C3 as their warrant.

Candidates come from data/candidate_peptides.tsv when the sponsor supplies one.
Absent that, a stand-in panel is drawn from benchmark peptides HELD OUT of the
C3 calibration, so the demo makes real calls against known answers instead of
scoring peptides it was tuned on.

Output: results/c4_prescreen.tsv, results/c4_prescreen.json
"""
import csv
import json
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import data_path, results_path, iedb_mhci  # noqa: E402

ENDPOINT = "https://tools-cluster-interface.iedb.org/tools_api/mhci/"
ALLELES = ["HLA-A*32:01", "HLA-B*57:01"]
N_CANDIDATES = 5          # the RFP's 5 uniblock peptides
COMPLEXES_REQUESTED = 8   # the RFP's "8 uniblock complexes"
SEED = 20260918
# Classic affinity standard used by the field and by the comparator report.
IC50_STANDARD = 500.0


def tsv(name, path=None):
    p = path or results_path(name)
    if not os.path.exists(p):
        return []
    with open(p) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def load_candidates():
    """Sponsor sequences if present, else a held-out stand-in panel."""
    supplied = data_path("candidate_peptides.tsv")
    if os.path.exists(supplied):
        rows = tsv(None, supplied)
        return [(r["id"], r["peptide"].strip().upper()) for r in rows], "sponsor-supplied"

    cal = json.load(open(results_path("c3_calibration.json")))
    heldout = set(tuple(x) for x in cal.get("held_out", []))
    if not heldout:
        raise SystemExit("no candidates and no held-out set in c3_calibration.json")
    bench = {(r["allele"], r["peptide"]): int(r["label"]) for r in tsv("c1_benchmark.tsv")}
    rng = random.Random(SEED)
    # one peptide per candidate slot, balanced so the demo shows both calls
    peps = sorted({p for _, p in heldout})
    pos = [p for p in peps if any(bench.get((a, p)) == 1 for a in ALLELES)]
    neg = [p for p in peps if p not in pos]
    rng.shuffle(pos)
    rng.shuffle(neg)
    chosen = (pos[:3] + neg[:2])[:N_CANDIDATES]
    rng.shuffle(chosen)
    return [(f"uniblock_{i+1}", p) for i, p in enumerate(chosen)], "held-out benchmark stand-in"


def score(peptides):
    out = {}
    for method, key in (("netmhcpan_el", "el"), ("netmhcpan_ba", "ba")):
        rows = iedb_mhci(ENDPOINT, method, peptides, ALLELES)
        for r in rows:
            d = out.setdefault((r["allele"], r["peptide"]), {})
            d[f"{key}_rank"] = float(r["percentile_rank"])
            if "ic50" in r:
                d["ic50_nM"] = float(r["ic50"])
    return out


def main():
    cands, provenance = load_candidates()
    cal = json.load(open(results_path("c3_calibration.json")))
    cuts = {a: cal["per_allele"][a]["recommended_skip_cut"] for a in ALLELES}
    npvs = {a: cal["per_allele"][a]["recommended_skip_npv"] for a in ALLELES}

    peps = [p for _, p in cands]
    scored = score(peps)
    bench = {(r["allele"], r["peptide"]): int(r["label"]) for r in tsv("c1_benchmark.tsv")}

    rows = []
    for cid, pep in cands:
        for allele in ALLELES:
            s = scored.get((allele, pep), {})
            el = s.get("el_rank")
            cut = cuts[allele]
            call = "BUILD" if el is not None and el < cut else "SKIP"
            rows.append({
                "candidate": cid, "peptide": pep, "allele": allele,
                "el_rank": el, "ba_rank": s.get("ba_rank"),
                "ic50_nM": s.get("ic50_nM"),
                "passes_ic50_500": (s.get("ic50_nM") is not None
                                    and s["ic50_nM"] < IC50_STANDARD),
                "skip_cut": cut, "call": call,
                "warrant": (f"NPV {npvs[allele]:.2f} at %Rank>={cut}"
                            if call == "SKIP" else ""),
                "known_label": bench.get((allele, pep)),
            })

    rows.sort(key=lambda r: (r["el_rank"] if r["el_rank"] is not None else 99))
    build = [r for r in rows if r["call"] == "BUILD"]

    # Where the requested budget lands
    budget = rows[:COMPLEXES_REQUESTED]
    checked = [r for r in rows if r["known_label"] is not None]
    agree = sum(1 for r in checked
                if (r["call"] == "BUILD") == (r["known_label"] == 1))

    out = {
        "alleles": ALLELES, "candidate_provenance": provenance,
        "n_candidates": len(cands), "complexes_possible": len(rows),
        "complexes_requested_by_rfp": COMPLEXES_REQUESTED,
        "skip_cuts": cuts, "skip_npv": npvs,
        "n_build": len(build), "n_skip": len(rows) - len(build),
        "ic50_standard_nM": IC50_STANDARD,
        "n_passing_ic50_standard": sum(1 for r in rows if r["passes_ic50_500"]),
        "agreement_with_known_labels": (
            {"n_checkable": len(checked), "n_agree": agree,
             "accuracy": round(agree / len(checked), 3) if checked else None}),
        "rows": rows,
    }
    with open(results_path("c4_prescreen.json"), "w") as f:
        json.dump(out, f, indent=1)
    cols = ["candidate", "peptide", "allele", "el_rank", "ba_rank", "ic50_nM",
            "passes_ic50_500", "skip_cut", "call", "warrant", "known_label"]
    with open(results_path("c4_prescreen.tsv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    print(f"candidates: {len(cands)} ({provenance})")
    print(f"possible complexes: {len(rows)}   RFP asks for {COMPLEXES_REQUESTED}\n")
    print(f"{'candidate':12s} {'peptide':11s} {'allele':13s} {'EL%':>6s} {'BA%':>6s} "
          f"{'IC50':>9s} {'<500':>5s} {'call':>6s} {'truth':>6s}")
    print("-" * 84)
    for r in rows:
        t = {1: "bind", 0: "non", None: "-"}[r["known_label"]]
        print(f"{r['candidate']:12s} {r['peptide']:11s} {r['allele']:13s} "
              f"{r['el_rank']:6.2f} {r['ba_rank']:6.2f} {r['ic50_nM']:9.0f} "
              f"{'Y' if r['passes_ic50_500'] else 'n':>5s} {r['call']:>6s} {t:>6s}")
    print(f"\nBUILD {len(build)} / SKIP {len(rows)-len(build)}")
    if checked:
        print(f"agreement with IEDB labels on the {len(checked)} checkable complexes: "
              f"{agree}/{len(checked)}")
    print(f"passing the classic IC50<{IC50_STANDARD:.0f} nM standard: "
          f"{out['n_passing_ic50_standard']}/{len(rows)}")
    print(f"\nwrote {results_path('c4_prescreen.tsv')}")


if __name__ == "__main__":
    main()
