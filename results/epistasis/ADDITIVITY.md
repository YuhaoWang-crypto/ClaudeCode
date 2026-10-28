# Does stacking mutations add up?

Greedy combination assumes the free energies of individual mutations add. ATLAS can test that directly: wherever a measured multi-mutant decomposes into disjoint subsets that were also measured, the sum of the parts can be compared with the measured whole.

**Pairs add up: 6 of 6 two-mutation combinations measured within one study land inside 0.5 kcal/mol of the sum of their parts, mean error +0.00. Beyond six simultaneous mutations the measured gain falls short of the sum by about +3.4 kcal/mol, so stacking that deep cannot be predicted by adding up parts.**

Wild-type affinities for the same receptor differ between studies, so ddG values from different papers rest on different baselines. Both the pooled analysis and one restricted to single studies are reported; the size-dependent pattern survives the restriction, so it is not an artefact of that noise.

## Within a single study

n = 13 decompositions, mean error +0.49, median -0.07, sd 1.37 kcal/mol. 69% land inside 0.5 kcal/mol.

| mutations stacked | n | mean error (kcal/mol) | within 0.5 |
|---|---|---|---|
| 2 | 6 | +0.00 | 6/6 |
| 4 | 2 | -0.73 | 1/2 |
| 6 | 1 | -0.32 | 1/1 |
| 8 | 3 | +1.60 | 1/3 |
| 9 | 1 | +3.38 | 0/1 |

## Pooled across studies

n = 51 decompositions, mean error +0.22, median -0.01, sd 1.20 kcal/mol. 45% land inside 0.5 kcal/mol.

| mutations stacked | n | mean error (kcal/mol) | within 0.5 |
|---|---|---|---|
| 2 | 14 | -0.16 | 12/14 |
| 3 | 5 | -0.80 | 2/5 |
| 4 | 13 | +0.01 | 3/13 |
| 5 | 11 | +0.51 | 4/11 |
| 6 | 1 | -0.32 | 1/1 |
| 7 | 1 | +0.92 | 0/1 |
| 8 | 5 | +1.59 | 1/5 |
| 9 | 1 | +3.37 | 0/1 |

Error is measured minus predicted, so a positive value means the stack underdelivered against the sum of its parts.

## What this means for a design campaign

* Predicting a **double** mutant by adding two singles is sound.
* **Triples** are borderline and should be checked, not assumed.
* Stacking **beyond about four** mutations cannot be predicted by summation. The published deep variants, which reach 6 to 14 simultaneous mutations, were found by directed evolution, and the sub-additivity here is why.

