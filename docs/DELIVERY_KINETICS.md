# Mechanistic delivery-kinetics model (Phase 4)

Model: `lipidlib/kinetics.py` · Analysis: `analysis/delivery_kinetics.py` ·
Figure: `results/figures/delivery_kinetics.png`

Complements the ML potency score (Track A) with an **interpretable, dynamic**
model: given biophysical rates it produces the protein-expression *time-course* and
summary metrics (peak height, time-to-peak, total exposure AUC) that a single
potency number can't give.

## The model
A 4-compartment first-order cascade, following Müller et al. 2024 (sequential
stochastic transfer) and Mihaila et al. 2017/2019 (LNP ODE models):

```
L_ex --k_uptake--> L_endo --k_escape--> M(cytosol) --k_transl--> P(protein)
                     |  \                    |                       |
                  k_lyso  k_recycle       k_mdeg                   k_pdeg
```

**Endosomal-escape efficiency = k_escape / (k_escape + k_lyso + k_recycle)** — the
rate-limiting ~1–2% bottleneck, and the step the ionizable lipid (Track A) governs.
Because the cascade is linear, total protein exposure has a closed form
(`analytic_auc`), which the numeric integration reproduces to <1%:

    AUC_P = dose · f_escape · mRNA_per_lnp · k_transl / (k_mdeg · k_pdeg)

Literature-informed defaults (per hour): uptake 0.3, escape 0.02, lysosomal 1.0,
recycle 0.1, mRNA decay 0.058 (t½≈12 h), translation 10, protein decay 0.23
(t½≈3 h, Fluc-like). These give **escape ≈ 1.8%** and **protein peaking ≈ 13 h** —
consistent with observed mRNA-LNP expression kinetics.

## What the analysis shows (figure)
- **A — the cascade in time**: endosomal LNP spikes early and is cleared within
  hours; cytosolic mRNA peaks next; protein peaks ~13 h then decays. Textbook shape.
- **B — sensitivity**: protein AUC scales linearly (log-log slope 1) with **escape,
  mRNA-stability, and translation**, but is **flat in uptake** (all internalised LNP
  is eventually processed here). So expression is set multiplicatively by escape ×
  cargo-stability × translation.
- **C — link to Track A**: mapping a lipid's LiON predicted-delivery `z` to an
  escape rate (`k_escape = 0.02·10^z`, capped at ~8% escape) turns the ML ranking
  into expression **dynamics**. Top vs weak combo_v1 lead spans **~12× in AUC**,
  but **time-to-peak is lipid-independent (~13 h)** — the lipid sets expression
  *amplitude*, while mRNA/protein stability set the *timing/duration*.

| combo_v1 lead | z | escape | peak | AUC |
|---|---|---|---|---|
| H7+amide+C10 (top) | +1.31 | 8.3% | 13 h | 6247 |
| H3+amide+C10 (median) | +0.09 | 2.2% | 13 h | 1646 |
| H12+ester+C6 (weak) | −0.42 | 0.7% | 13 h | 516 |

## Why this matters
- **Interpretability**: decomposes "potency" into uptake / escape / stability /
  translation — you see *which* step limits a design.
- **Dynamics the ML lacks**: peak time and duration (e.g. for dosing, or comparing a
  stable vs unstable mRNA on the same lipid).
- **What-if**: e.g. a 5× more stable mRNA → 5× AUC and longer duration, independent
  of the lipid; a better-escaping lipid → higher amplitude, same timing.

## Caveats
Deliberately minimal (per-cell, first-order, no explicit protein corona / PEG
desorption sub-steps, no extracellular clearance — so uptake looks lossless).
Parameters are literature-typical, not fit to a specific dataset. It's a mechanistic
scaffold to reason about dynamics and to fuse with Track A, not a calibrated
predictor. Next refinement: fit rates to a measured expression time-course, and add
PEG-desorption / corona and extracellular-clearance sub-steps from Müller 2024.
