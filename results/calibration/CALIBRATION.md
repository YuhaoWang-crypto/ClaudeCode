# Calibrating the triage score against measured ddG

## Data

| source | entries |
|---|---|
| ATLAS (TCR-side single mutants with ddG) | 273 |
| SKEMPI 2.0 (`TCR/pMHC` subset, single mutants) | 189 |
| overlapping between the two | 94 |
| **pooled unique mutations** | **368** |
| successfully mapped onto structures and scored | 261 |
| structures covered | 25 |

Of the scored mutations, 46 tighten binding (ddG < -0.2 kcal/mol) and 121 weaken it (ddG > 0.5). ddG is positive for weaker binding throughout.

## 1. The hand-weighted score, as shipped

| metric | value |
|---|---|
| Spearman vs measured ddG | -0.082 |
| Pearson vs measured ddG | +0.006 |
| AUC, ranking stabilising mutations first | +0.542 |
| AUC, ranking destabilising mutations last | +0.529 |

The score is built so that higher means tighter, so a working score gives a **negative** Spearman and AUCs **above 0.5**.

Per structure, where at least five mutations were measured:

| structure | n | Spearman |
|---|---|---|
| 1AO7 | 49 | -0.069 |
| 1MI5 | 38 | -0.378 |
| 1OGA | 22 | -0.775 |
| 2AK4 | 13 | -0.124 |
| 2VLJ | 13 | -0.287 |
| 3QDG | 14 | +0.112 |
| 3QDJ | 16 | +0.038 |
| 5BRZ | 7 | +0.750 |
| 5BS0 | 7 | +0.750 |
| 5E9D | 8 | -0.476 |
| 6V1A | 15 | -0.011 |
| 8I5C | 6 | +0.314 |
| 8TRR | 5 | -0.800 |
| 9EJG | 11 | +0.309 |
| 9EJI | 9 | +0.083 |
| 9NMU | 5 | -0.300 |

## 2. Refitting the six weights

| metric | held out | in sample |
|---|---|---|
| n | 261 | 261 |
| Spearman | +0.125 | +0.241 |
| Pearson | +0.125 | +0.266 |
| RMSE (kcal/mol) | +0.980 | +0.944 |
| AUC, finding stabilising mutations | +0.568 | - |

The AUC row is the one that matters for design. Most of this benchmark is alanine scanning, so a model can score well overall simply by learning that deleting a large buried side chain weakens binding. Finding the rarer mutations that *improve* binding is the actual task, and it is harder.

Predicting the training mean for everything gives RMSE +0.979. This model **does not beat** that baseline.

Fitted coefficients (kcal/mol per unit of each feature):

| feature | coefficient |
|---|---|
| packing | -1.026 |
| hbond | +0.297 |
| substitution_prior | +0.243 |
| apolar_burial | +0.241 |
| electrostatic | -0.224 |
| conformational | +0.071 |
| _intercept_ | +0.549 |

## 3. Fitting the raw structural measurements

| metric | held out | in sample |
|---|---|---|
| n | 261 | 261 |
| Spearman | +0.429 | +0.550 |
| Pearson | +0.438 | +0.587 |
| RMSE (kcal/mol) | +0.884 | +0.795 |
| AUC, finding stabilising mutations | +0.718 | - |

The AUC row is the one that matters for design. Most of this benchmark is alanine scanning, so a model can score well overall simply by learning that deleting a large buried side chain weakens binding. Finding the rarer mutations that *improve* binding is the actual task, and it is harder.

Predicting the training mean for everything gives RMSE +0.979. This model **beats** that baseline.

Fitted coefficients (kcal/mol per unit of each feature):

| feature | coefficient |
|---|---|
| delta_volume | -0.458 |
| contacts_peptide | +0.423 |
| buried_fraction | +0.361 |
| relative_sasa_free | -0.274 |
| delta_sasa | -0.172 |
| wild_type_is_glycine | +0.171 |
| hydropathy_times_burial | -0.160 |
| partner_polar_atoms | -0.125 |
| growth_reaches_partner | +0.111 |
| delta_charge | +0.108 |
| volume_times_burial | -0.101 |
| mutant_is_proline | -0.086 |
| wild_type_is_proline | +0.086 |
| delta_hydropathy | +0.082 |
| partner_charge | +0.076 |
| delta_hbond_capacity | +0.075 |
| void_slots | -0.059 |
| mutant_is_glycine | +0.054 |
| blosum62 | -0.038 |
| contacts_mhc | +0.015 |
| _intercept_ | +0.323 |
