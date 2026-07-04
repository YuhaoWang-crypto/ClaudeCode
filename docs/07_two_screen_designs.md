# Two Screen Types in One Experiment: Essentiality + Chemogenomic/Synthetic-lethal

Bottom line: **you don't need two experiments.** Using one kinome library, one
Cas9 line, and one transduction, set up a **vehicle (DMSO) arm** and a **drug
arm** in parallel to get both results:

- **Essentiality / dropout screen** = vehicle arm **T_end vs. T0** (which kinases
  this line depends on).
- **Chemogenomic / synthetic-lethal screen** = **drug arm vs. vehicle arm** (same
  T_end) (which kinase knockouts change drug sensitivity).

The vehicle arm does double duty: it is both the essentiality experiment and the
control for the drug arm.

---

## 1. Arms and timepoints

```
                        ┌─ vehicle (DMSO) ──passage/vehicle──► T_end (DMSO)
plasmid ─►transduce─►puro─►T0 ─┤
                        └─ drug ─────────────passage/drug─────► T_end (Drug)
```

| Sample | Timing | Purpose |
|--------|--------|---------|
| plasmid | pre-virus | library input reference, QC |
| T0 | post-selection (~day 7) | shared reference for both analyses |
| T_end (DMSO) | ~day 21 | essentiality (vs. T0) + synthetic-lethal control (vs. drug) |
| T_end (Drug) | ~day 21 | synthetic-lethal (vs. DMSO) |

- **Replicates:** ≥3 biological replicates per arm.
- **Shared T0:** one T0 set (3 reps); both arms split from the same transduced
  cells before dosing.

### Sample count and sequencing (exact)

The calculator conservatively counts T0 per-arm (13 samples). With a shared T0:

```
plasmid 1 + T0 ×3 + T_end(DMSO) ×3 + T_end(Drug) ×3 = 10 samples
```

- ~4.4M reads/sample → **~44M reads** total, well under one NovaSeq lane.
- Each added drug dose/agent → +3 samples, +~13M reads.

Other scale numbers (≥3.5M cells/sample maintained, ~23 µg gDNA/sample, ~13.5M
cells to plate at MOI 0.3) are in `05_screen_and_sequencing.md`. Recompute with:

```bash
python3 scripts/coverage_calculator.py --genes 518 --guides-per-gene 6 \
  --controls 400 --coverage 1000 --arms 2 --timepoints 2 --replicates 3
```

---

## 2. The critical variable: setting up the drug arm

Essentiality needs no extra setup; the **synthetic-lethal arm lives or dies on
dose and pressure window**:

1. **Dose-response pre-experiment:** measure the dose-response curve in your line
   and pick a **sub-lethal dose — usually IC20–IC30** (~20–30% growth inhibition).
   Too high collapses the population and loses coverage; too low gives no
   selective pressure.
2. **Sustained dosing** across the whole passaging window (refresh at each media
   change), to ~8–10 doublings.
3. **Maintain coverage:** even as the drug arm grows slower, hold **≥1000×
   cells/guide** — increase starting cells or extend time to accumulate enough
   doublings.
4. **Matched DMSO:** vehicle arm gets the same volume of DMSO.

### Genetic-background synthetic lethality (optional alternative)

If your question is "which kinases become essential in a given genotype" (not a
drug), replace the drug arm with an **isogenic pair**: screen the same library in
**wild-type vs. mutant** lines (e.g., TP53⁺ᐟ⁺ vs. TP53⁻ᐟ⁻, or KRAS-mutant vs.
corrected) and compare dropout profiles. Same structure — "two arms" becomes "two
cell lines."

---

## 3. Running the two analyses

From one guide-count matrix (`06_analysis.md` §1), run two analysis lines:

### A. Essentiality / dropout (vehicle arm)

| Step | Method |
|------|--------|
| Comparison | T_end(DMSO) **vs.** T0 |
| Primary | **BAGEL2** (calibrated Bayes Factor via CEGv2/NEGv2 priors) |
| Secondary | **MAGeCK RRA** (gene-level LFC/rank) |
| Correction | **Chronos/CERES** copy-number bias |

### B. Chemogenomic / synthetic-lethal (drug vs. vehicle)

| Step | Method |
|------|--------|
| Comparison | T_end(Drug) **vs.** T_end(DMSO) (same timepoint, removes pure essentiality) |
| Primary | **drugZ** (tuned for drug-vs-control synthetic-lethal/resistance) |
| Secondary | **MAGeCK MLE** (multi-condition design matrix) |
| Interpretation | more depleted with drug → **sensitizer (synthetic lethal)**; enriched with drug → **resistance** |

> Key: synthetic lethality is **drug vs. vehicle**, not drug vs. T0 — the latter
> folds in kinases that are simply essential. The vehicle arm exists precisely to
> subtract that baseline essentiality.

---

## 4. Deliverables and naming (large files referenced by name)

| Product | Note |
|---------|------|
| `results/guide_counts.tsv` | guide×sample count matrix (small, linkable) |
| `results/essentiality_bagel.tsv` | essentiality BF/gene-level (small, linkable) |
| `results/synlethal_drugz.tsv` | synthetic-lethal drugZ results (small, linkable) |
| `results/qc_report.html` | QC report |
| raw FASTQ, full count BAMs | large — **referenced by name** in summaries, not linked |

---

## 5. One-page flow

1. Build/validate Cas9 line → clone kinome library (`04_library_construction.md`).
2. Dose-response pre-experiment → IC20–IC30.
3. Large-scale transduction (MOI 0.3, ≥13.5M cells) → puro select → collect T0 (3 reps).
4. Split: DMSO arm + drug arm, 3 reps each, maintain ≥3.5M cells/sample, dose to ~day 21.
5. Harvest both arms → gDNA (~23 µg/sample) → 2-step PCR → sequence (~44M reads).
6. Count matrix → Line A BAGEL2 (essentiality) + Line B drugZ (synthetic lethal)
   → QC → hit interpretation and validation (multi-guide concordance + cDNA rescue).
