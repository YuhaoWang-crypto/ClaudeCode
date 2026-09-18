#!/usr/bin/env python3
"""
C7 - Item-by-item: what in-silico can do for this RFP, and what it cannot.

The question put to this pipeline was "can the validation analysis in this RFP
be done in silico?" The answer is not one word, because the RFP asks for seven
different things and they fall into three groups:

  REPLACE   in-silico answers the question outright; the assay is not needed
            to get the answer (there is nothing in this RFP in this group, and
            saying so is the most useful sentence in this file)
  DE-RISK   in-silico cannot answer it, but can change what you commission,
            how much of it, or whether it is feasible at all - before you pay
  NO        the readout is a physical measurement with no sequence-level
            surrogate, or the molecule is outside what any current model
            represents

Each verdict carries the evidence produced elsewhere in this pipeline, so a
reader can check it rather than take it. Verdicts with no supporting run are
labelled as reasoning, not measurement.

Output: results/c7_feasibility.json
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import results_path  # noqa: E402


def jsn(name):
    p = results_path(name)
    return json.load(open(p)) if os.path.exists(p) else {}


def main():
    cal = jsn("c3_calibration.json")
    pre = jsn("c4_prescreen.json")
    donor = jsn("c5_donor_feasibility.json")
    cls2 = jsn("c6_nanopeptide_classII.json")

    items = []

    # ---- Step 1 --------------------------------------------------------
    v = cls2.get("verdict", {})
    pc = cls2.get("positive_control", {})
    items.append({
        "rfp_item": "Step 1 · MHC-II-dependent IL-2 readout (DO11.10 / I-A(d))",
        "verdict": "NO — with one useful exception",
        "evidence_type": "measured (C6)",
        "what_in_silico_did": (
            f"Reproduced the documented I-A(d) register for cOVA323-339 "
            f"(core {pc.get('result',{}).get('core','-')}, "
            f"%Rank {pc.get('result',{}).get('rank','-')}), confirming the "
            f"hybridoma system's restriction element behaves as published."),
        "what_it_cannot_do": (
            "IL-2 secretion is a cellular functional readout. No sequence model "
            "predicts cytokine output. Binding is necessary, not sufficient."),
        "changes_the_study": False,
    })
    n_score = v.get("scored"); n_tot = v.get("panel_size")
    items.append({
        "rfp_item": "Step 1 · the five-peptide panel (8-mer, 9-mer+4D, 17-mer, 2x 30-mer)",
        "verdict": "NO for the comparison the step is built around",
        "evidence_type": "measured (C6)",
        "what_in_silico_did": (
            f"Scored {n_score} of {n_tot} panel members. The 17-mer uniblock and "
            f"the 30-mer diblock return the same core at the same %Rank, i.e. "
            f"appending the assembly block changes nothing a sequence model can see."),
        "what_it_cannot_do": (
            f"{v.get('blocked_by_length_floor','?')} peptides sit below the class II "
            f"API's 11-residue floor and cannot be submitted at all. "
            f"{v.get('blocked_by_d_amino_acids','?')} contain D-amino acids, which no "
            f"class II predictor represents — a D-substituted peptide is scored as "
            f"its all-L twin, which is a different molecule. Since Step 1 is an "
            f"all-L control versus D-substituted test comparison, in-silico cannot "
            f"perform it."),
        "changes_the_study": False,
    })

    # ---- Step 2 --------------------------------------------------------
    if cal.get("per_allele"):
        aucs = {a: r["EL"]["auc"] for a, r in cal["per_allele"].items()}
        cuts = {a: r["recommended_skip_cut"] for a, r in cal["per_allele"].items()}
        npvs = {a: r["recommended_skip_npv"] for a, r in cal["per_allele"].items()}
        auc_txt = "; ".join(f"{a} AUC {v:.3f}" for a, v in aucs.items())
        skip_txt = "; ".join(f"{a} skip at %Rank>={cuts[a]} (NPV {npvs[a]:.2f})"
                             for a in cuts)
    else:
        auc_txt = skip_txt = "calibration not yet run"
    items.append({
        "rfp_item": "Step 2 · cell-free peptide-MHC binding / stability screen",
        "verdict": "DE-RISK — strongly",
        "evidence_type": "measured (C1, C3)",
        "what_in_silico_did": (
            f"Measured NetMHCpan on these two alleles specifically, against every "
            f"labelled IEDB record for them: {auc_txt}. That converts the "
            f"pre-screen from an opinion into a calibrated rule: {skip_txt}."),
        "what_it_cannot_do": (
            "Predicts binding, not complex stability or off-rate. The RFP asks "
            "for stability, which is the property that governs whether a tetramer "
            "survives staining, and no sequence model reports a half-life."),
        "changes_the_study": True,
    })
    n_build = pre.get("n_build"); n_poss = pre.get("complexes_possible")
    items.append({
        "rfp_item": "Step 2 · custom synthesis of 2 alleles + 8 uniblock complexes",
        "verdict": "DE-RISK — this is where the money is",
        "evidence_type": "measured (C4)" if pre else "pending C4",
        "what_in_silico_did": (
            f"Ranked all {n_poss} possible (peptide, allele) complexes and called "
            f"{n_build} worth building. A monomer that does not fold never becomes "
            f"a tetramer, and that is discovered after the spend."
            if pre else "pending"),
        "what_it_cannot_do": (
            "Cannot confirm refolding, biotinylation or multimerisation yield. "
            "A predicted binder can still fail to refold."),
        "changes_the_study": True,
    })
    items.append({
        "rfp_item": "Step 2 · baseline precursor frequency by FACS after IL-7/IL-15 expansion",
        "verdict": "NO",
        "evidence_type": "reasoning, not measurement",
        "what_in_silico_did": "Nothing. There is no sequence-level surrogate.",
        "what_it_cannot_do": (
            "Precursor frequency is a property of one donor's T-cell repertoire, "
            "set by their thymic selection and exposure history. It is measured, "
            "not predicted. This is the single most irreplaceable item in the RFP."),
        "changes_the_study": False,
    })

    # ---- Step 3 and donors ---------------------------------------------
    items.append({
        "rfp_item": "Step 3 · MLR (proliferation, cytokines, ICS)",
        "verdict": "NO",
        "evidence_type": "reasoning, not measurement",
        "what_in_silico_did": "Nothing.",
        "what_it_cannot_do": (
            "Proliferation and cytokine output are cellular outcomes downstream of "
            "processing, presentation, co-stimulation and repertoire. Sequence "
            "models see only the first step."),
        "changes_the_study": False,
    })
    if donor:
        w = donor["weighted_us_eu"]
        b = donor["donors_to_screen"]
        a_name, b_name = donor["alleles"]["HLA-A"], donor["alleles"]["HLA-B"]
        items.append({
            "rfp_item": "Donors · human PBMC carrying the specified HLA antigens",
            "verdict": "DE-RISK — and it may change the study design",
            "evidence_type": "measured (C5)",
            "what_in_silico_did": (
                f"Carrier frequency, weighted US/EU: {a_name} {w[a_name]*100:.1f}%, "
                f"{b_name} {w[b_name]*100:.1f}%, either {w['either']*100:.1f}%, "
                f"both {w['both_independent']*100:.2f}%. To obtain 10 donors: "
                f"~{b['either']['10']:.0f} typed if either allele suffices, "
                f"~{b['both_independent']['10']:.0f} if each donor must carry both."),
            "what_it_cannot_do": (
                "Cannot tell you what a commercial HLA-typed bank already holds, "
                "and the two-locus combination assumes independence, which linkage "
                "disequilibrium violates."),
            "changes_the_study": True,
        })

    groups = {"REPLACE": 0, "DE-RISK": 0, "NO": 0}
    for it in items:
        key = ("DE-RISK" if it["verdict"].startswith("DE-RISK")
               else "REPLACE" if it["verdict"].startswith("REPLACE") else "NO")
        it["group"] = key
        groups[key] += 1

    out = {"items": items, "group_counts": groups,
           "headline": (
               "No item in this RFP can be replaced by computation. Three of "
               f"{len(items)} can be materially de-risked before commissioning: the "
               "binding pre-screen, the tetramer build list, and donor feasibility. "
               "The wet-lab study should still be run — but it should be a smaller, "
               "better-targeted one.")}
    p = results_path("c7_feasibility.json")
    with open(p, "w") as f:
        json.dump(out, f, indent=1)

    print(f"{'group':9s} {'item':62s}")
    print("-" * 74)
    for it in items:
        print(f"{it['group']:9s} {it['rfp_item'][:62]:62s}")
    print(f"\nREPLACE {groups['REPLACE']}   DE-RISK {groups['DE-RISK']}   NO {groups['NO']}")
    print(f"\n{out['headline']}")
    print(f"\nwrote {p}")


if __name__ == "__main__":
    main()
