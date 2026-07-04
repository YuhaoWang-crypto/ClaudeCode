# Drug-Arm SOP: Dose Pre-experiment, Dosing Scheme, and Hit Interpretation

This turns the "drug arm" of `07_two_screen_designs.md` into something directly
executable. The worked example uses **A375 (BRAF^V600E melanoma) + vemurafenib
(PLX4032)** — a template; swap cell line, drug, and dose for your own (placeholders
marked `<...>`). See `10_case_study_A375_BRAFi.md` for the fully data-filled version.

> The example cell line/drug are just anchors. What actually determines drug-arm
> success is the **dose-response you measure in your own line**, not literature IC
> values — so Step 1 is mandatory.

**Example drug check (queried live from ChEMBL v34):**

| Item | Value |
|------|-------|
| Drug | Vemurafenib (Zelboraf; research codes PLX4032 / RG7204) |
| ChEMBL ID | CHEMBL1229517 |
| Mechanism | **BRAF (serine/threonine kinase) inhibitor**, direct binding, for the **V600E** mutant |
| Target | BRAF (CHEMBL5145, UniProt P15056 V600E) |
| Status | Approved (2011 FDA, first-in-class; USAN stem `-rafenib` = RAF kinase inhibitor) |

---

## 1. Dose pre-experiment (set IC20–IC30) — mandatory

Goal: find a **sub-lethal concentration** that applies mild, sustained selective
pressure over the whole screen — neither collapsing the population nor leaving
sensitizer/resistance guides unmoved.

### 1a. Short-term dose-response (72 h)

| Item | Setting |
|------|---------|
| Plate | 96-well (or 384), ≥3 replicate wells/dose |
| Seeding | `<2000–5000>` cells/well, attach 24 h (tune so 72 h stays sub-confluent) |
| Dose series | 8–10 points, 2–3-fold serial dilution, spanning two logs around expected IC50 |
| Controls | equal-volume DMSO vehicle wells + no-cell blanks |
| Incubation | 72 h |
| Readout | CellTiter-Glo (ATP luminescence) or resazurin, normalized to vehicle |
| Fit | 4-parameter logistic → read GI50/IC50 and **IC20–IC30** |

### 1b. Chronic dose validation (critical, often skipped)

A 72 h IC differs from the equivalent inhibition under **sustained dosing over the
~3-week screen**. Validate with a small growth curve that the chosen dose lowers
the **proliferation rate** by only ~20–30% long-term:

- 6-well, `<chosen IC20–IC30>` vs. DMSO, 3 replicate wells each.
- Count/passage every 2–3 days for 10–14 days; plot cumulative population doublings.
- Target: drug arm doubling rate ~**70–80%** of the vehicle arm. Fine-tune the
  final concentration accordingly.

**Output:** a final concentration `<C_screen>` for the screen drug arm.

---

## 2. Dosing scheme

| Item | Setting |
|------|---------|
| Start | at T0 (post-puro), split half to vehicle, half to drug arm |
| Concentration | drug arm constant `<C_screen>`; vehicle arm equal-volume DMSO |
| Refresh | re-dose to `<C_screen>` at every media change/passage (every 2–3 days) |
| Duration | to ~day 21 (≈ 8–10 population doublings) |
| Coverage | both arms maintain **≥1000× cells/guide** (drug arm grows slower — increase starting cells or extend to accumulate doublings) |
| Replicates | ≥3 biological replicates per arm |

---

## 3. drugZ comparison matrix (synthetic-lethal/resistance readout)

drugZ compares **drug vs. vehicle at the same T_end** (not vs. T0 — that folds in
baseline essentiality):

| Comparison | control (-c) | treatment (-x) | Direction |
|------------|--------------|----------------|-----------|
| Synthetic lethal / sensitizer | T_end(DMSO) ×3 | T_end(Drug) ×3 | **more depleted with drug** → KO makes drug more lethal |
| Resistance | same | same | **enriched with drug** → KO confers resistance |

drugZ reports both directions in one run (normZ negative = sensitizer, positive =
resistance) with rank and FDR.

---

## 4. Expected hit categories (double as internal positive controls)

Using the BRAFi example — **categories are reliable; specific genes to be confirmed
with your data and the literature**:

- **Resistance / enrichment (KO → survive):** nodes that reactivate MAPK. The
  classic category is **RTK-mediated bypass reactivation** (EGFR, PDGFRB, IGF1R,
  FGFR family) and RAF/MEK-level compensation (RAF1/CRAF, MAP2K1, MAP3K8/COT).
- **Sensitizer / dropout (KO → more killing):** survival/feedback nodes
  **synthetic-lethal with MAPK inhibition** (cell-cycle kinases sustaining
  proliferation, components limiting pathway feedback reactivation).

> Interpretation: first check whether these **known categories** appear in the
> expected direction — if they do, dose and pressure window are reasonable and new
> hits are credible; if both directions are flat, the dose is likely too low or
> the duration too short.

### Candidate kinase-inhibitor menu (pick a drug arm by pathway)

| Pathway/context | Common drug-arm candidates | Typical cell background |
|-----------------|----------------------------|-------------------------|
| BRAF-MAPK | vemurafenib/dabrafenib (BRAFi), trametinib (MEKi) | BRAF^V600E melanoma (e.g., A375) |
| EGFR | osimertinib, gefitinib | EGFR-mutant NSCLC (e.g., PC9) |
| KRAS downstream | MEKi (selumetinib), SHP2i, KRAS^G12C inhibitors | KRAS-mutant lung/colorectal |
| PI3K-AKT-mTOR | alpelisib (PI3Kα), mTOR inhibitors | PIK3CA-mutant |
| CDK4/6 | palbociclib | ER+ breast, RB-intact |

> ChEMBL / clinical-trial sources can pull candidate inhibitors' mechanism,
> potency (IC50/Ki), and clinical stage for your target — see the fully worked
> `10_case_study_A375_BRAFi.md`, which does this for vemurafenib.

---

## 5. Sample sheet (see `templates/sample_sheet.csv`)

Single-drug example = **10 samples**: plasmid 1 + T0 ×3 + T_end(DMSO) ×3 +
T_end(Drug) ×3. Each added dose/agent = +3 samples. Indices and FASTQ paths are in
the template file.
