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
render_complex.py        structure figure (backbone tube + contact map), no PyMOL needed
campaign_summary.py      three-panel campaign summary
data/candidates.csv      top 12 designs per corrected run, full metrics
data/candidates.fasta    the same sequences
data/run_log.md          Boltz job IDs and settings
structures/              design models and independent co-folds of the two leads (mmCIF)
figures/                 lead complexes and campaign summary
```
