# Worked Case Study: A375 (BRAF^V600E melanoma) + vemurafenib drug arm

A **fully specified** instantiation of the generic workflow in `07/08/09`, with
potency and clinical numbers filled in from **ChEMBL v34** and
**ClinicalTrials.gov**. Swap the cell line/drug for your own system, keeping the
same structure.

> **Data provenance:** IC50/mechanism below are from **ChEMBL v34**
> (vemurafenib = CHEMBL1229517; target BRAF = CHEMBL5145, UniProt P15056 V600E);
> clinical dose and combination landscape from **ClinicalTrials.gov**. Both
> queried live in this session — not recalled from memory.

---

## 1. System spec

| Item | Setting |
|------|---------|
| Cell line | **A375** (human cutaneous melanoma), **BRAF^V600E**, doubling ~17–24 h |
| Engineering | Build/validate **A375-Cas9** stable line first (`04_library_construction.md`) |
| Library | Kinome ~3,508 guides (518 kinases × 6 + 400 controls; `05_screen_and_sequencing.md`) |
| Drug arm | **vemurafenib (PLX4032 / Zelboraf)**, BRAF^V600E inhibitor |
| Genetic-context alternative | For a genotype synthetic-lethal design, use A375 (BRAF^V600E) vs. a BRAF-WT melanoma line as an isogenic-style pair |

**Why this is a good exemplar:** A375 is **oncogene-addicted** to BRAF — BRAF
knockout is itself a strong dropout — and vemurafenib applies pathway (MAPK)
pressure, so kinases whose loss changes MAPK/drug dependence become visible.

---

## 2. Real potency data (ChEMBL v34) — sets the dose-response range

| Assay level | IC50 | Note |
|-------------|------|------|
| Biochemical BRAF^V600E | **3.2 – 31 nM** (pChEMBL 7.5–8.5, multiple sources) | intrinsic target potency |
| A375 cellular pERK inhibition | **33 nM** (WB, 90 min) / **150 nM** (ELISA, 72 h) / 190–260 nM (others) | **target-engagement** level |
| Anti-proliferative GI50 (A375) | typically **~0.5–1 µM** (not in this set; measure your own) | drives the screen dose |

> Key distinction: **target-engagement IC50 (~0.1 µM) ≠ anti-proliferative GI50
> (~0.5–1 µM)**. Screen pressure is set by the anti-proliferative curve, so the
> Step 3 dose-response must be measured in-house — do not transplant ChEMBL's
> pERK values directly.

---

## 3. Dose selection (concretizing `08` Step 1)

1. **Dose-response range:** 8–10 points, **10 nM to 10 µM** (3-fold dilutions,
   spanning the cellular IC50 ~0.1 µM by two logs each way).
2. **72 h CellTiter-Glo** → fit GI50; expect ~0.5–1 µM (confirm empirically).
3. **Take IC20–IC30** as the screen concentration `C_screen` — for this system
   likely **~0.2–0.5 µM**.
4. **Chronic validation:** sustain `C_screen` for 10–14 days so the drug arm
   grows at ~70–80% of the DMSO arm's rate; adjust accordingly.

---

## 4. Concrete samples and drugZ matrix

Sample naming (per `templates/sample_sheet.csv`, with `DRUG` → `VEM`):

```
plasmid | T0_1..3 | DMSO_1..3 | VEM_1..3        (10 samples)
```

### Line A · Essentiality (vehicle vs. T0, BAGEL2)

```bash
BAGEL.py fc -i kinome_screen.count.txt -o kinome_fc -c T0_1,T0_2,T0_3
BAGEL.py bf -i kinome_fc.foldchange -o results/essentiality_bagel.tsv \
  -e CEGv2.txt -n NEGv2.txt -c DMSO_1,DMSO_2,DMSO_3
```

### Line B · Synthetic lethality (vemurafenib vs. vehicle, drugZ)

```bash
python drugz.py -i kinome_screen.count.txt -o results/synlethal_drugz.tsv \
  -c DMSO_1,DMSO_2,DMSO_3 -x VEM_1,VEM_2,VEM_3
# normZ < 0: KO sensitizes A375 to vemurafenib (sensitizer / synthetic lethal)
# normZ > 0: KO confers resistance
```

---

## 5. Expected hits for this system (positive controls + interpretation prior)

### Positive controls (confirm the screen worked first)

- **Essentiality arm:** pan-essential kinases (PLK1, CDK1, AURKB) + **BRAF
  itself** strongly deplete — A375 is BRAF-addicted, so BRAF KO drops out even
  without drug. If these don't drop, do not trust new hits yet.
- **Drug arm:** expect signal in the MAPK-sensitization direction (below).

### Sensitizers / synthetic lethal (normZ < 0, deplete more with drug)

Direction is reliable; specific genes to be confirmed with your data + literature.
These categories are the biological basis of the clinical combinations:

| Category | Candidate kinases | Clinical parallel |
|----------|-------------------|-------------------|
| Vertical MAPK reinforcement | MAP2K1/MAP2K2 (MEK), MAPK1/MAPK3 (ERK) | **BRAFi + MEKi** (vemurafenib+cobimetinib approved) |
| Adaptive RTK feedback (ERK rebound) | EGFR, ERBB2/3, FGFR, IGF1R, PDGFRB, MET | bypass reactivation drives resistance |
| Parallel survival pathway | PIK3CA, AKT1/2/3, PDPK1, MTOR | BRAFi + PI3K/AKT synergy |
| Cell cycle | CDK4, CDK6 | **BRAFi + CDK4/6i** combination |
| Survival/adhesion | SRC family, PTK2 (FAK) | sustains survival signaling |

> Note on MEK/ERK: they are **baseline-essential in the vehicle arm too**; drugZ
> (drug vs. vehicle) subtracts baseline essentiality and keeps only the
> drug-specific **increment** — which is exactly why "drug vs. T0" is wrong here.

### Resistance / enrichment (normZ > 0, KO → survive) — expect sparse

**Methodological point:** clinical BRAFi resistance is mostly **gain-of-function**
(RTK upregulation, MAP3K8/COT overexpression, CRAF compensation) — invisible to a
knockout screen — while loss-of-function resistance events (NF1, PTEN, CIC, MED12)
are mostly **non-kinase**. So a **kinome-KO screen is under-powered for the
resistance direction** — a design limitation, not a failure. To map
activation-driven resistance systematically, use **overexpression/ORF or CRISPRa**
screens instead of knockout.

---

## 6. Clinical relevance (ClinicalTrials.gov, live query)

- Registered trials of vemurafenib + the MEK inhibitor cobimetinib: **≥36**
  (e.g., NCT01656642, NCT04722575). BRAFi+MEKi is now standard combination
  therapy — the clinical confirmation that this screen's "sensitizer direction =
  MAPK feedback-reactivation nodes" is the right axis.
- Practical use: mapping hits back onto **approved/investigational combinations**
  quickly triages which sensitizer hits are translationally interesting
  (validated MEK/CDK4-6 vs. novel nodes).

---

## 7. One-page checklist (this system)

1. Build A375-Cas9, validate activity → clone kinome library.
2. Vemurafenib dose-response (10 nM–10 µM) → set `C_screen` (expect ~0.2–0.5 µM)
   → chronic validation.
3. Transduce (MOI 0.3, ≥13.5M cells) → puro select → T0 (3 reps).
4. Split DMSO arm / vemurafenib@`C_screen` arm, 3 reps each, maintain ≥3.5M cells,
   dose to ~day 21.
5. Harvest → gDNA (~23 µg/sample) → 2-step PCR → sequence (~44M reads).
6. BAGEL2 (essentiality, with BRAF positive control) + drugZ (sensitizer
   direction: MEK/ERK/PI3K/CDK4-6 prior).
7. Map hits onto the BRAFi combination landscape for translational triage.
