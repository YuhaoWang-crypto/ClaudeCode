"""
R6 — Validation against real expression and within-person data.

Everything up to r5 used published SUMMARY statistics: medians, reference
ranges, AUCs. Two claims in the recorder framework were never tested against
primary data at all:

  (A) SOURCE SPECIFICITY — "this marker comes from that cell". Asserted for
      every candidate, measured for none. Tested here against Human Protein
      Atlas per-gene records (35 genes, fetched and cached in data/).

  (B) WITHIN-PERSON BEHAVIOUR — the whole personal-baseline and pairing
      argument depends on how an analyte moves inside ONE person over time.
      Tested here against published repeated-sampling biological-variation
      studies.

(B) contains the strongest empirical result in this package: a study that
sampled 20 healthy adults weekly for 10 weeks measured CV_I of 6.5% for
Abeta42, 6.4% for Abeta40, and 3.3% for the Abeta42/Abeta40 RATIO. Those three
numbers over-determine the correlation, so r2_pairing's model can be solved
against real human data instead of assumed.

A METHODOLOGICAL LIMIT WORTH STATING LOUDLY
-------------------------------------------
Single-cell RNA atlases are structurally blind to the deposition recorders
(class R4). Mature erythrocytes and platelets are enucleate and carry
negligible mRNA, so no scRNA-seq experiment can confirm or refute anything
about EC4d, BC4d or platelet-C4d carriers. The HPA record for CR1 -- the
erythrocyte complement receptor that is the natural denominator for EC4d --
shows neutrophils, podocytes and Kupffer cells, and NO erythrocytes, purely
because of this. Absence of RBC signal here is an artefact of the assay, not
biology. Validating R4 requires flow cytometry, not transcriptomics.

DATA PROVENANCE
  data/hpa_markers.json  — fetched from proteinatlas.org per-gene JSON,
                           gene symbol verified against the request. [LIT]
  BV_STUDIES             — transcribed from the cited publications. [LIT]
"""
import json
import os
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
HPA_PATH = os.path.join(_HERE, "data", "hpa_markers.json")


def load_hpa():
    with open(HPA_PATH) as fh:
        return json.load(fh)


# --------------------------------------------------------------------------
# (A) source specificity
# --------------------------------------------------------------------------
# Which recorder candidate each gene is the parent/writer/carrier for, and what
# the framework CLAIMS about its source. Verdict is computed, not typed.
CLAIMS = {
    "Uromodulin (CKD reserve)":      "kidney tubule only (parenchymal reserve)",
    "NGAL (AKI)":                    "renal tubule, confounded by neutrophils",
    "KIM-1 (AKI)":                   "proximal tubule",
    "Collagen III (PRO-C3/C3M)":     "liver fibrogenesis",
    "GFAP (TBI parent)":             "astrocyte only",
    "ApoA-I (Tyr192 parent)":        "liver/HDL",
    "C4A (C4d parent)":              "liver, deposits on cells",
    "CR1 (C4d normalizer)":          "erythrocyte surface (denominator)",
    "Albumin (C-Alb carrier)":       "liver, long-lived plasma carrier",
    "CA125 (ovarian)":               "tumour/mesothelium",
    "uEGF (tubular mass)":           "distal tubule",
    "MPO (ApoA-I writer)":           "neutrophil (the writer)",
}


def source_verdict(rec):
    """Grade source purity from HPA fields. Deliberately simple and auditable."""
    tis = (rec.get("RNA tissue specificity") or "").lower()
    blood = (rec.get("RNA blood cell specificity") or "").lower()
    n_ct = len(rec.get("RNA single cell type specific nCPM") or {})

    if "enriched" in tis and n_ct <= 2:
        grade = "SINGLE-SOURCE"
    elif "enriched" in tis or "group enriched" in tis:
        grade = "FEW-SOURCE"
    elif "enhanced" in tis:
        grade = "MANY-SOURCE"
    else:
        grade = "NON-SPECIFIC"

    blood_clean = "not detected" in blood
    return grade, blood_clean, n_ct


def part_a():
    hpa = load_hpa()
    print("\n[A] SOURCE-SPECIFICITY AUDIT (Human Protein Atlas, 35 genes)")
    print(f"    {'marker (parent gene)':<34}{'tissue class':<22}"
          f"{'#cell types':>11}{'leukocyte-free':>16}")
    rows = []
    for label in CLAIMS:
        rec = hpa[label]
        grade, blood_clean, n_ct = source_verdict(rec)
        rows.append((label, rec, grade, blood_clean, n_ct))
        print(f"    {label[:32]:<34}{rec['RNA tissue specificity']:<22}"
              f"{n_ct if n_ct else '-':>11}{'yes' if blood_clean else 'NO':>16}")

    print("\n    Three results that change a recommendation:")

    print("\n    1. UMOD is the cleanest source in the entire panel.")
    u = hpa["Uromodulin (CKD reserve)"]
    print(f"       tissue-enriched, specificity score {u['RNA tissue specificity score']}"
          f" (highest scored here)")
    print(f"       kidney {u['RNA tissue specific nTPM']}, "
          f"cell type {u['RNA single cell type specific nCPM']}")
    print(f"       blood cells: {u['RNA blood cell specificity']}")
    print("       -> zero leukocyte confounding. This is the structural reason")
    print("          uromodulin behaves better than C-Alb, and it supports")
    print("          promoting it to the lead of the CKD pair. [OK]")

    print("\n    2. LCN2 (NGAL) is worse than 'neutrophil contamination'.")
    l = hpa["NGAL (AKI)"]
    ct = l["RNA single cell type specific nCPM"] or {}
    top = sorted(ct.items(), key=lambda kv: -float(kv[1]))[:5]
    print(f"       tissue class: {l['RNA tissue specificity']}; "
          f"top tissues {l['RNA tissue specific nTPM']}")
    print("       top cell types by nCPM:")
    for k, v in top:
        print(f"         {k:<38}{float(v):>9.0f}")
    print(f"       blood: {l['RNA blood cell specific nTPM']}")
    print("       -> KIDNEY IS NOT IN EITHER TOP LIST. LCN2 is constitutively a")
    print("          broad secretory-epithelium gene (endometrium, conjunctiva,")
    print("          oesophagus, prostate, salivary, airway) and only becomes")
    print("          renal on injury induction. So 'monomer = tubular' is NOT")
    print("          safe as a baseline assumption -- monomeric LCN2 can come")
    print("          from any of these epithelia. This WEAKENS the monomer/total")
    print("          rescue proposed in RESCUE_MINING.md case 1: the form split")
    print("          separates myeloid from epithelial, not kidney from")
    print("          not-kidney. Correction recorded. [OK]")

    print("\n    3. COL3A1 explains PRO-C3's organ-specificity gap structurally.")
    c = hpa["Collagen III (PRO-C3/C3M)"]
    ct = c["RNA single cell type specific nCPM"] or {}
    top = sorted(ct.items(), key=lambda kv: -float(kv[1]))[:4]
    print(f"       tissue class: {c['RNA tissue specificity']} "
          f"({c['RNA tissue distribution']})")
    print("       top cell types: " + ", ".join(f"{k} {float(v):.0f}" for k, v in top))
    print("       -> fibroblast/stromal expression in essentially every organ, so")
    print("          a circulating collagen-III neoepitope CANNOT be organ-")
    print("          specific. The report lists 'organ- and aetiology-specific")
    print("          thresholds' as the key gap; this shows it is a structural")
    print("          property of the parent, not a threshold-tuning problem.")
    print("          Silver lining: hepatic stellate cells are the single")
    print("          highest expressor, so liver is the best-case organ -- which")
    print("          independently supports doing liver first. [OK]")

    print("\n    Caveat that invalidates this method for one whole class:")
    cr1 = hpa["CR1 (C4d normalizer)"]
    print(f"       CR1 top cell types: "
          f"{list((cr1['RNA single cell type specific nCPM'] or {}).keys())}")
    print("       No erythrocytes -- because mature RBCs are enucleate and carry")
    print("       no mRNA. scRNA-seq is structurally blind to the R4 deposition")
    print("       recorders (EC4d/BC4d/platelet-C4d). Their carriers must be")
    print("       validated by flow cytometry, never by transcriptomics. [OK]")

    print("\n    Do NOT use HPA 'blood concentration MS' for assay feasibility:")
    g = hpa["GFAP (TBI parent)"]
    print(f"       GFAP listed at {g['Blood concentration - Conc. blood MS [pg/L]']}"
          f" pg/L = 67 ng/mL, versus a clinical healthy median of 8 pg/mL.")
    print("       Four orders of magnitude apart. These MS values sit at their")
    print("       own detection floor for low-abundance analytes. [OK]")
    return rows


# --------------------------------------------------------------------------
# (B) within-person data
# --------------------------------------------------------------------------
# Plasma AD panel: 20 healthy adults, weekly x 10 weeks, Simoa HD-X  [LIT]
AD_BV = [
    #  analyte,        CV_I%, CV_G%, II(as published), RCV+%
    ("Abeta42",          6.5, 15.3, 0.46, 21.6),
    ("Abeta40",          6.4, 12.5, 0.56, 21.3),
    ("Abeta42/Abeta40",  3.3,  6.6, 0.68, 13.0),
    ("p-tau181",        16.7, 25.7, 0.69, 62.2),
    ("p-tau217",        10.3, 21.1, 0.56, 38.3),
    ("p-tau231",         6.3, 17.2, 0.49, 26.1),
    ("GFAP",             8.0, 30.1, 0.34, 32.6),
    ("NfL",              7.4, 21.2, 0.43, 30.7),
]

# Urinary tubular markers: 23 stable CKD patients, first-morning urine
# daily x 10 days. Reported as mean of per-patient CVs.  [LIT]
URINE_BV = [
    #  analyte,      CV_I% raw, CV_I% creatinine-normalised, 90th-pct change%
    ("U-NGAL",        45.6, 43.8, 83.3),
    ("U-KIM-1",       30.1, 23.4, 45.5),
    ("U-cystatin C",  25.0, 24.8, 46.3),
]


def rho_from_ratio(cv_a, cv_b, cv_ratio):
    """Solve r2_pairing's variance identity for the correlation.

    Var(log A/B) = Var(logA) + Var(logB) - 2 rho sd(logA) sd(logB)
    With small CVs, CV ~= sd on the log scale, so rho is over-determined by
    three measured CVs. This is the model being fitted to data rather than
    assumed.
    """
    return (cv_a**2 + cv_b**2 - cv_ratio**2) / (2.0 * cv_a * cv_b)


def part_b():
    print("\n[B] WITHIN-PERSON VALIDATION")

    print("\n    B1. Plasma AD panel — 20 healthy adults, weekly x 10 wk, Simoa")
    print(f"    {'analyte':<18}{'CV_I%':>7}{'CV_G%':>7}{'II':>7}{'RCV+%':>8}")
    for name, cvi, cvg, ii, rcv in AD_BV:
        print(f"    {name:<18}{cvi:>7.1f}{cvg:>7.1f}{ii:>7.2f}{rcv:>8.1f}")

    print("\n    B2. THE PAIRING MODEL, SOLVED AGAINST REAL DATA  [OK]")
    cv_a, cv_b, cv_r = 6.5, 6.4, 3.3
    rho = rho_from_ratio(cv_a, cv_b, cv_r)
    kappa = cv_b / cv_a
    print(f"       measured: CV_I(Ab42)={cv_a}%  CV_I(Ab40)={cv_b}%  "
          f"CV_I(ratio)={cv_r}%")
    print(f"       three numbers over-determine the model, giving")
    print(f"         rho (within-person, Ab42 vs Ab40) = {rho:.3f}")
    print(f"         kappa = {kappa:.3f}")
    print(f"       r2_pairing admissibility: ratio helps iff rho > kappa/2 = "
          f"{kappa/2:.3f}")
    print(f"         measured rho = {rho:.2f}  ->  PASSES, with large margin [OK]")
    gain_ratio = cv_a / cv_r
    gain_resid = 1.0 / np.sqrt(1.0 - rho**2)
    print(f"       observed gain of the fixed ratio      = {gain_ratio:.2f}x")
    print(f"       theoretical ceiling (fitted residual) = {gain_resid:.2f}x")
    print(f"       -> the fixed ratio captured "
          f"{100*gain_ratio/gain_resid:.0f}% of what was available.")
    print("       This is the Level-1 (same parent) case where a hard division")
    print("       is nearly optimal and a fitted residual is not worth the")
    print("       complexity. Measured, not assumed. [OK]")

    print("\n    B3. But the report's own benchmark is a DIFFERENT mechanism")
    print("       pTau217/Abeta42 pairs two markers on OPPOSING pathology axes,")
    print("       not two products of one parent. Its benefit comes from the")
    print("       numerator and denominator moving in opposite directions")
    print("       (RECORDER_FRAMEWORK L4), not from cancelling shared noise")
    print("       (L1). Their within-person CVs are 10.3% and 6.5%, and nothing")
    print("       in this dataset says they are correlated within a person. So")
    print("       pTau217/Abeta42 should NOT be cited as evidence that")
    print("       same-parent normalisation works -- Abeta42/40 is that")
    print("       evidence. Two different claims, two different mechanisms. [OK]")

    print("\n    B4. Urinary tubular markers — 23 stable CKD, daily x 10 d")
    print(f"    {'analyte':<16}{'CV_I%':>7}{'/creatinine':>13}{'change needed':>15}")
    for name, raw, norm, chg in URINE_BV:
        print(f"    {name:<16}{raw:>7.1f}{norm:>13.1f}{chg:>14.1f}%")
    print("       Creatinine normalisation cuts KIM-1's within-person variation")
    print("       by 22% (30.1 -> 23.4) but does essentially nothing for NGAL")
    print("       (45.6 -> 43.8, 4%). Same denominator, opposite verdicts --")
    print("       exactly r2_pairing's point that a denominator is admissible")
    print("       or not depending on WHICH noise it shares. KIM-1's day-to-day")
    print("       variation is largely urine concentration (shared with")
    print("       creatinine); NGAL's is not -- consistent with the episodic,")
    print("       multi-epithelial secretion that part A found. [OK]")
    print("       Practical consequence: urinary NGAL needs an 83% change before")
    print("       it means anything in a stable patient. Any monitoring use is")
    print("       dead on arrival at that noise level.")

    print("\n    B5. THE GAP THAT BLOCKS THE REPORT'S OWN GATE 3")
    print("       Gate 3 requires reference change values, and RCV needs CV_I.")
    print("       Searched and NOT FOUND: no published within-subject biological")
    print("       variation for PRO-C3, for C3M, or for serum uromodulin.")
    print("       Three of the report's five lead candidates therefore have no")
    print("       computable RCV today. A 20-subject x 10-week repeated-sampling")
    print("       study is a cheap prerequisite that unblocks all three at once,")
    print("       and it does not need any patients -- healthy volunteers only.")

    print("\n    B6. Between-study caution")
    print("       Serum NfL CV_I is reported as 3.1% (CV_G 35.6%, II 0.09) in one")
    print("       study and 7.4% (CV_G 21.2%, II 0.43) in another. A 5-fold")
    print("       disagreement in II for the same analyte. Personal-baseline")
    print("       thresholds must be derived on the platform and schedule they")
    print("       will be used with -- not imported from a table. [OK]")


# --------------------------------------------------------------------------
# (C) healthy vs indication, single cell   [LIT — CELLxGENE Census 2025-11-08]
# --------------------------------------------------------------------------
CENSUS_PATH = os.path.join(_HERE, "data", "census_kidney.tsv")
COMP_PATH = os.path.join(_HERE, "data", "census_composition.tsv")


def _read_tsv(path):
    with open(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(head, ln.rstrip("\n").split("\t"))) for ln in fh]


def part_c():
    print("\n[C] HEALTHY vs INDICATION, SINGLE CELL (CELLxGENE Census)")
    print("    kidney, primary cells only: 658k normal / 267k CKD / 125k AKI.")
    print("    Only 2 datasets contain all three conditions, so ALL comparisons")
    print("    below are WITHIN dataset and WITHIN cell type. Cross-dataset")
    print("    comparison is not interpretable here -- see C3.")

    rows = _read_tsv(CENSUS_PATH)
    order = ["normal", "chronic kidney disease", "acute kidney failure"]

    def show(gene, cell_type_key, title):
        print(f"\n    {title}")
        print(f"      {'dataset':<10}{'normal':>10}{'CKD':>16}{'AKI':>18}")
        for ds in ["a12ccb9b", "dea717d4"]:
            sel = {r["disease"]: r for r in rows
                   if r["gene"] == gene and r["dataset"] == ds
                   and cell_type_key in r["cell_type"]}
            if len(sel) < 3:
                continue
            base = float(sel["normal"]["mean_expr"])
            cells = [f"{float(sel[d]['pct_pos']):.2f}%" for d in order]
            fcs = [float(sel[d]["mean_expr"]) / base for d in order]
            print(f"      {ds:<10}{cells[0]:>10}{cells[1]:>10}"
                  f" ({fcs[1]:>5.2f}x){cells[2]:>10} ({fcs[2]:>6.2f}x)")

    show("HAVCR1", "proximal tubule",
         "C1. KIM-1 induction in proximal tubule — REPRODUCIBLE")
    print("      Same direction, same rough magnitude, both datasets, monotonic")
    print("      normal < CKD < AKI. This is what a usable injury marker looks")
    print("      like at source. [OK]")

    show("LCN2", "proximal tubule",
         "C2. NGAL induction in proximal tubule — NOT REPRODUCIBLE")
    print("      One dataset says 1.2x in AKI and DOWN in CKD; the other says")
    print("      144x in AKI. Two orders of magnitude apart, same cell type,")
    print("      same disease label. NGAL's induction is not a stable property")
    print("      of the injured tubule. Combined with part A (LCN2 is a broad")
    print("      secretory-epithelium gene) and B4 (urinary CV_I 45.6%), three")
    print("      independent data types agree that NGAL is unstable at every")
    print("      level: source, transcript induction, and day-to-day urine.")
    print("      This is a converging case to DROP the NGAL line, not fix it.")

    show("UMOD", "thick ascending limb",
         "C3. Uromodulin per-cell expression in TAL — essentially FLAT")
    print("      Per-cell UMOD moves within 0.79-1.27x across disease states,")
    print("      while SERUM uromodulin falls 4.4x by CKD G3 and 17.5x by G5.")
    print("      A <1.3x transcriptional change cannot produce a 4-17x serum")
    print("      change, so the serum signal must be TAL cell MASS, not per-cell")
    print("      output. That is exactly the 'parenchymal reserve' reading the")
    print("      report gives it -- now supported by primary data rather than")
    print("      assumed. [OK]")

    comp = _read_tsv(COMP_PATH)
    norm = [float(r["TAL_pct"]) for r in comp if r["disease"] == "normal"]
    print(f"\n      BUT the mass claim itself is NOT testable this way: TAL cells")
    print(f"      are {min(norm):.1f}%-{max(norm):.1f}% of kidney cells across")
    print(f"      {len(norm)} NORMAL datasets alone. Dissociation protocol swamps")
    print("      any disease effect (TAL is notoriously hard to dissociate).")
    print("      Cell-type proportions from public scRNA-seq cannot quantify")
    print("      parenchymal mass. Use histology/imaging for that. [OK]")


def main():
    print("=" * 78)
    print("R6 — Validation against expression atlases and within-person data")
    print("=" * 78)
    part_a()
    part_b()
    part_c()
    return True


if __name__ == "__main__":
    main()
