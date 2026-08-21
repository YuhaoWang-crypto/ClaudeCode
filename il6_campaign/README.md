# De novo miniprotein binders against human IL-6 — an in-silico demo campaign

A scaled-down run of the workflow from *Autonomous de novo protein binder design with Claude*
(Anthropic, Aug 2026), applied to a target that was **not** among the paper's 16: human
interleukin-6 (UniProt P05231, mature 183-mer).

Everything here is a **prediction**. Nothing was expressed or measured. The paper's own
conclusion applies with full force: a confident co-fold is a useful filter, not evidence of
binding, and experimental screening is the only way to learn whether a design actually binds.

## What was run

| Stage | What | Tool | Scale |
|---|---|---|---|
| Target & epitope | IL-6 chain of the IL-6/IL-6Rα/gp130 hexamer (PDB 1P9M), epitopes derived from the experimental complex (5 Å contacts) | gemmi / Biopython | 2 epitopes |
| Backbone + sequence design | epitope-conditioned miniprotein design, 55–85 aa, no Cys, hydrophobic fraction ≤ 0.42, N-glyc motifs excluded | Boltz-2 protein design API | 4 runs × 120 = **480 designs** |
| Ranking | interface confidence (ipTM) + minimum interface PAE | design-time model | all 480 |
| Interface analysis | buried surface area, atom-pair counts, H-bonds, salt bridges, epitope recall | Biopython Shrake-Rupley | top 12 per run |
| Independent confirmation | re-co-folding of the delivered sequence with the target, **no design template, binder without an MSA**, 5 samples | Boltz-2 structure+binding | 6 designs |
| Specificity control | lead binder co-folded against an unrelated antigen (CLEC12A ectodomain), the paper's off-target control | Boltz-2 | 1 |

Total compute cost: **≈ USD 25** (design $0.05/sequence, co-fold $0.02/sample). Wall clock ≈ 50 min.

## Two epitopes, derived not guessed

IL-6 signals by assembling with IL-6Rα (site I) and then gp130 (site II/III). Both footprints
were computed directly from 1P9M rather than taken from memory:

- **site I** — IL-6 residues contacting IL-6Rα: R30, L33, K54, N61, K66, E69, C73, F74, Q75, F78, E172, Q175, L178, R179, A180, R182, Q183. Blocking it prevents the first receptor engagement (the mechanism of siltuximab/sirukumab).
- **site II** — IL-6 residues contacting gp130: L19, R24, K27, Q28, R30, Y31, D34, E110, Q111, R113, A114, M117, S118, V121, Q124, F125, K128. Blocking it leaves IL-6Rα binding intact but stops signal assembly.

Numbering is 1P9M chain-B numbering, which is the mature-sequence position **+ 1** (so the
canonical site I arginines R179/R182 are mature positions 178/181). Getting this off by one
silently shifts every epitope claim, so it is checked in code, not assumed.

## A real bug, and what it cost

The first two runs specified the epitope by 0-indexed position among the **resolved** residues of
1P9M chain B (163 of them). The API indexes into the **entity** sequence (186 residues, including
the ones the crystal does not resolve). The designs came back beautifully confident — and bound
a surface ~18 residues away from the intended one: site I recall 0.06.

Re-running with corrected indices (`index = mature_position + 2`) fixed it:

| Run | median ipTM | designs > 0.9 | epitope recall of top 12 (median) |
|---|---|---|---|
| site I, buggy indices | 0.677 | 6 | 0.06 |
| site I, corrected | 0.677 | 6 | **0.53** (best 0.71) |
| site II, buggy indices | 0.734 | 19 | 0.29 |
| site II, corrected | **0.838** | **31** | **0.82** (best 0.94) |

The lesson is the paper's in miniature: the confidence score did not notice. Both the buggy and
corrected site I runs score identically (median ipTM 0.677) while binding entirely different
surfaces. Only an explicit geometric check against the experimental complex caught it.

## Lead candidates

`data/candidates.csv` / `.fasta` hold the top 12 designs of each corrected run with full metrics.
The two leads:

| | **IL6-S2-01** (site II, gp130 blocker) | **IL6-S1-01** (site I, IL-6Rα blocker) |
|---|---|---|
| Length | 61 aa, all-α | 61 aa, all-α |
| Sequence | `SSNPLVRLAEELLRELEENPESEYAEVLLESAELILKDLEEKYPEEAARLRERLERLKKRL` | `KQKERREKAYNEEVEALQKELGLSKSLAKKILTAIETGDKSLLPKEQEELYEKALEIIEKL` |
| Design-time ipTM / min iPAE | 0.966 / 0.71 Å | 0.917 / 1.53 Å |
| Independent co-fold ipTM (5 samples) | **0.948 – 0.955** | 0.870 – 0.897 |
| Complex pLDDT | 0.90 | 0.92 |
| Boltz binding score | **0.765** | 0.211 |
| Buried surface (co-fold) | 2508 Å² | 1532 Å² |
| H-bonds / salt bridges | 16 / 9 | 15 / 5 |
| Epitope recall | 0.88 of the gp130 footprint | 0.65 of the IL-6Rα footprint |
| Designed pose reproduced by the independent co-fold | 24/32 interface residues (Jaccard 0.63) | 15/17 (Jaccard 0.75) |

For scale, the **natural** interfaces in 1P9M bury 1454 Å² (IL-6·IL-6Rα, K_D ≈ nM) and 1352 Å²
(IL-6·gp130) with 3 and 5 H-bonds. Both designs bury as much or more surface with more polar
contacts — necessary for a good binder, nowhere near sufficient as proof of one.

IL6-S1-01 engages the canonical site I hot spot directly: salt bridges from its Glu13 to IL-6
Arg179 and Arg182 (1P9M numbering), the two arginines that dominate IL-6Rα recognition.

## Specificity control

The paper runs every design against an unrelated antigen (CLEC12A). Same test here, on IL6-S1-01:

| | ipTM (5 samples) | Boltz binding score |
|---|---|---|
| vs IL-6 | 0.870 – 0.897 | 0.211 |
| vs CLEC12A ectodomain | 0.745 – 0.847 | 0.00067 |

The binding score separates target from off-target by ~300×; **ipTM barely separates them at
all** (0.89 vs 0.79). That is the single most important caveat in this report: interface
confidence on a designed complex is not a specificity assay. A prediction model asked to dock two
proteins will dock them.

## Re-scoring pass: orthogonal judge + ipSAE_min + a DockQ gate

The ranking above used design-time ipTM from the model that produced the designs. Three
upgrades were then applied to the 32-design shortlist (union of the top 12 by ipTM and
the top 12 by ipSAE_min in each run):

1. **ipSAE_min instead of ipTM**, computed with the reference implementation
   (`ipsae.py` v4, Dunbrack lab — it reads Boltz PAE natively), pae cutoff 10 Å, dist 15 Å.
2. **Chai-1 as an orthogonal judge**, run on a Modal A10G from sequence alone (no design
   template, no MSA), returning its own PAE so the same ipSAE_min applies. Ensemble score
   = mean of the per-run z-scores of the two predictors.
3. **DockQ as a gate, not a ranking term** — real DockQ (fnat + LRMSD + iRMSD) between the
   design model and the Chai-1 prediction, threshold 0.23 ("acceptable" pose).

### What it changed

| | site I run | site II run |
|---|---|---|
| ipSAE_min (Boltz), median over 120 designs | 0.03 | 0.25 |
| Shortlisted designs passing the DockQ gate | **5 / 18** | **11 / 14** |
| ipTM top-5 still in the gated top-5 | 2 / 5 | 3 / 5 |
| Designs Chai-1 scores at ipSAE_min = 0 that Boltz scored 0.46–0.83 | 8 | 2 |

Three findings worth keeping:

- **The designer over-scores its own work.** Ten shortlisted designs carry Boltz ipSAE_min
  of 0.46–0.83 and get exactly **0.00** from Chai-1, with DockQ 0.01–0.11. The single
  highest Boltz ipSAE_min in the whole site I run (0.830) is one of them. This is the
  athlete-and-referee problem made visible, and it is why the orthogonal judge is the
  upgrade that earns its cost.
- **The site II design ranked #1 by ipTM fails the gate.** Chai-1 gives it ipTM 0.383 and
  DockQ 0.108; it drops from rank 1 to rank 10 and is cut. Both of the leads delivered
  earlier survive: IL6-S1-01 is #1 in the gated site I list (DockQ 0.44) and IL6-S2-01
  is #2 in site II (DockQ 0.38).
- **ipSAE_min alone did not obviously pick better designs.** On site I it reshuffled half
  the top 12 (Spearman ρ with ipTM = 0.66, top-12 overlap 6/12), and its top-5 mean
  epitope recall was higher than ipTM's (0.46 vs 0.28) — but at top-12 the two are
  identical (0.24) and neither metric correlates with epitope recall (ρ = +0.09 and
  +0.02). On site II the two metrics agree (ρ = 0.87) and both pick perfectly
  on-epitope designs. The metric swap is free and directionally right; the evidence that
  it helps comes from the paper's labelled data, not from this run.
- **The gate is not a substitute for the epitope check.** Two site I designs pass DockQ
  comfortably while contacting almost none of the IL-6Rα footprint (recall 0.00 and 0.12):
  a reproducible pose on the wrong surface. Rank on the ensemble, gate on DockQ, and
  *separately* require epitope recall.

DockQ on the earlier Boltz confirmations of the two leads: IL6-S1-01 **0.936** ("high",
fnat 0.95, iRMSD 0.47 Å), IL6-S2-01 **0.641** ("medium", fnat 0.71, iRMSD 2.53 Å) — both
far more informative than the interface-residue Jaccard proxy used in the first pass.

Added cost: ~USD 2 of Modal A10G time for 32 Chai-1 folds; ipSAE and DockQ are CPU-only.

## Honest limits

- No wet-lab data. Hit rate here is unknown and unknowable from these numbers; the paper's
  measured rate for comparable campaigns was 27% overall, 49% for top-ranked designs.
- Boltz-2's affinity head is calibrated on **small molecules**, not protein–protein complexes, so
  the "binding score" is used here only as a relative ranking signal — no K_D is quoted, and none
  should be inferred.
- One design run per epitope; run-to-run variance is not characterised.
- No novelty screen against the PDB/UniRef was run (the paper's designs were 98% novel at 30%
  identity); these sequences are idealised helical bundles and should be screened before use.
- Every structure in `structures/` is a prediction, including the "confirmation" co-folds.

## Reproduce

```bash
pip install gemmi biopython numpy matplotlib
python analyze_interface.py <complex.cif> --target-chain A --binder-chain X --offset 1
python render_complex.py  <complex.cif> --target-chain A --binder-chain X --offset 1 --epitope siteII --out fig.png
python campaign_summary.py --designs designs.json --iface v2_iface.json --out figures/campaign_summary.png
```

Design and co-fold jobs were run through the Boltz-2 API (`boltz_start_protein_design`,
`boltz_start_structure_and_binding`); job IDs are recorded in `data/run_log.md`.

## Files

```
analyze_interface.py     interface metrics + epitope-overlap for any two-chain complex
batch_interface.py       the same over every design in a run
run_ipsae.py             drives the reference ipsae.py over a directory of Boltz outputs
chai_score_modal.py      Chai-1 scoring on a Modal GPU (orthogonal judge)
dockq.py                 DockQ between a design model and an independent prediction
rescore_ensemble.py      two-predictor z-score ensemble + DockQ gate
rescoring_figure.py      the re-scoring figure
render_complex.py        structure figure (backbone tube + contact map), no PyMOL needed
campaign_summary.py      three-panel campaign summary
data/candidates.csv      top 12 designs per corrected run, full metrics
data/candidates.fasta    the same sequences
data/run_log.md          Boltz job IDs and settings
data/ipsae_scores.csv    ipSAE_min, pDockQ, LIS for all 240 corrected-run designs
data/interface_all.csv   interface metrics + epitope recall for all 240
data/rescored_shortlist.csv  the 32-design shortlist with both predictors, DockQ and ranks
structures/              design models and independent co-folds of the two leads (mmCIF)
figures/                 lead complexes and campaign summary
```
