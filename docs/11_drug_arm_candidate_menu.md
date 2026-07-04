# Drug-Arm Candidate Menu (ChEMBL-verified)

A ready-to-pick menu of approved kinase/pathway inhibitors for the chemogenomic
drug arm (`08_drug_arm_sop.md` §4), with **potency, target, and approval status
queried live from ChEMBL v34** (snapshot 2026-07-04). Pick the row whose pathway
matches your cell line's driver; the "expected screen readout" column says what
the kinome screen should surface around it.

> **Provenance:** all IC50, target, mechanism, and approval fields below are from
> **ChEMBL v34** (queried live). Values are representative potent measurements —
> use them to bracket your **own** dose-response (`08` Step 1), not as the screen
> concentration itself.

## Verified candidate table

| Pathway | Drug (ChEMBL ID) | Target gene(s) | Mechanism | Representative IC50 (ChEMBL) | Approved (1st-in-class) |
|---------|------------------|----------------|-----------|------------------------------|-------------------------|
| **BRAF-MAPK** | Vemurafenib (CHEMBL1229517) | **BRAF** (V600E) | BRAF^V600E inhibitor | biochem 3.2–31 nM; A375 pERK 33–150 nM | 2011 ✓ |
| **MEK-MAPK** | Trametinib (CHEMBL2103875) | **MAP2K1/MAP2K2** (MEK1/2) | allosteric MEK1/2 inhibitor | MEK2 1.6 nM, MEK1 3.4 nM; COLO205 GI50 1 nM | 2013 ✓ |
| **EGFR** | Osimertinib (CHEMBL3353410) | **EGFR** (T790M) | 3rd-gen covalent EGFR TKI | EGFR L858R/T790M 2.5–7 nM (H1975 cell 2.5 nM) | 2015 |
| **KRAS** | Sotorasib (CHEMBL4535757) | **KRAS** (G12C) | covalent KRAS^G12C inhibitor | KRAS 30–68 nM; MIA PaCa-2 GI50 5 nM | 2021 ✓ |
| **PI3K** | Alpelisib (CHEMBL2396661) | **PIK3CA** (p110α) | PI3Kα-selective inhibitor | PI3Kα 4.6–5 nM | 2019 |
| **AKT** | Capivasertib (CHEMBL2325741) | **AKT1/AKT2/AKT3** | pan-AKT inhibitor | AKT2 8 nM, AKT3 8 nM (also PKA 7 nM, S6K 6 nM) | 2023 ✓ |
| **CDK4/6** | Palbociclib (CHEMBL189963) | **CDK4/CDK6** | CDK4/6 inhibitor | CDK4 11 nM, CDK6 16 nM | 2015 ✓ |

## How to use each row (matched background + expected readout)

| Drug | Matched cell background | Is the drug target in the kinome library? | Expected screen readout |
|------|-------------------------|-------------------------------------------|-------------------------|
| Vemurafenib | BRAF^V600E melanoma (A375) | **Yes** (BRAF) — see `10_case_study_A375_BRAFi.md` | Sensitizers = MEK/ERK/RTK/PI3K/CDK4-6; resistance sparse (KO under-powered) |
| Trametinib | BRAF/KRAS-mutant (melanoma, CRC, NSCLC) | **Yes** (MAP2K1/2) | KO of RTK-feedback + PI3K/AKT nodes sensitize; ERK-reactivators (e.g. RAF1) modulate |
| Osimertinib | EGFR-mutant NSCLC (PC9, H1975) | **Yes** (EGFR) | Sensitizers = parallel RTKs (MET, ERBB family), MAPK/PI3K nodes; MET/AXL relate to resistance |
| Sotorasib | KRAS^G12C NSCLC/CRC (MIA PaCa-2, NCI-H358) | **No** (KRAS is a GTPase, not a kinase) — but its **pathway kinases are** | Sensitizers = RTK-feedback (EGFR, SHP2-linked), MEK/ERK, PI3K — adaptive-reactivation biology |
| Alpelisib | PIK3CA-mutant breast | **Only in the extended library** (PIK3CA is a lipid kinase — include the ~20 lipid kinases) | Sensitizers = AKT/mTOR axis, RTKs; feedback to MAPK |
| Capivasertib | PIK3CA/AKT1-mutant, PTEN-null | **Yes** (AKT1/2/3) | Sensitizers = PDPK1, mTOR, RTKs; MAPK cross-talk |
| Palbociclib | RB-intact, ER+ breast / CDK4/6-dependent | **Yes** (CDK4, CDK6) | Sensitizers = CDK2/cyclin-E bypass, mitotic kinases; RB-pathway modulators |

## Notes that affect design

- **Kinome-library coverage of the drug's own target is a built-in positive
  control.** If the drug inhibits a kinase in your library (BRAF, MEK, EGFR,
  AKT, CDK4/6), its knockout should interact with the drug in the expected
  direction — a sanity check that the drug arm is working.
- **Sotorasib is the exception:** KRAS is a **GTPase, not a kinase**, so it is
  *not* in a kinome library. Use it when your interest is the **downstream/parallel
  kinase network** around KRAS (RTK-SHP2 feedback, MEK/ERK, PI3K), not KRAS itself.
- **PIK3CA (alpelisib target) is a lipid kinase** — only present if you adopt the
  ~538-gene (protein + lipid) target set (`data/kinome_targets.md`).
- **Capivasertib is a pan-AKT inhibitor with measurable PKA/S6K activity** — expect
  broader pathway effects; interpret AKT-axis hits with that polypharmacology in mind.
- Every value here is a *literature/biochemical* anchor. The screen concentration
  is always the **IC20–IC30 of proliferation measured in your own line** (`08` §1).
