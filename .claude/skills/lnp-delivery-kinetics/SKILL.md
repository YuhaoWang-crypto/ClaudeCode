---
name: lnp-delivery-kinetics
description: >-
  Mechanistic (ODE) model of LNP → cell → protein delivery kinetics. Use when the
  user wants an interpretable, time-resolved estimate of mRNA/protein expression
  from biophysical rates — uptake, endosomal escape, mRNA/protein degradation,
  translation — rather than a single ML potency score. Answers "when does
  expression peak / how long does it last", "which step is rate-limiting", "what if
  the mRNA were more stable / the lipid escaped better", and links a lipid's LiON
  potency to expression dynamics. Triggers: "delivery kinetics", "endosomal escape
  model", "expression time-course", "mechanistic LNP model", "AUC of expression".
---

# LNP delivery-kinetics model (mechanistic, Phase 4)

CPU-only ODE model (no GPU). 4-compartment first-order cascade following Müller 2024
(sequential stochastic transfer) and Mihaila 2017-19 (LNP ODE):

```
L_ex --k_uptake--> L_endo --k_escape--> M(cytosol) --k_transl--> P(protein)
                     | k_lyso, k_recycle    | k_mdeg              | k_pdeg
```
Endosomal-escape efficiency = k_escape/(k_escape+k_lyso+k_recycle) is the
rate-limiting ~1–2% bottleneck and the step the ionizable lipid governs.

## Use it
```python
from lipidlib.kinetics import DeliveryParams, simulate, metrics, analytic_auc
p = DeliveryParams()                 # literature defaults (Fluc-like reporter)
df = simulate(p, t_end=96)           # time-course: t, L_ex, L_endo, M, P
m  = metrics(p)                      # escape_eff, peak_protein, t_peak_h, AUC (numeric≈analytic)
```
Full analysis + figure (cascade, sensitivity, Track-A link):
```bash
python analysis/delivery_kinetics.py   # -> results/figures/delivery_kinetics.png
```

## Key relationships
- Total protein exposure has a closed form (numeric match <1%):
  `AUC_P = dose · f_escape · mRNA_per_lnp · k_transl / (k_mdeg · k_pdeg)`.
- Expression AUC scales linearly with **escape, mRNA-stability, translation**; it is
  **flat in uptake** (all internalised LNP is processed).
- **Lipid sets amplitude, cargo sets timing**: mapping a lipid's LiON potency `z` to
  escape rate (`k_escape ≈ 0.02·10^z`, cap ~8% escape) changes expression height
  ~12× across leads but leaves time-to-peak (~13 h) unchanged.

## When to use vs the ML model
- ML (Track A / `lipid-library-screen`): ranks lipids by a potency scalar.
- This model: turns rates into **dynamics** (peak time, duration) and a
  **mechanistic decomposition** (which step limits a design) — use for dosing
  reasoning, cargo (mRNA-stability) what-ifs, and interpretability.

## Caveats
Minimal by design (per-cell, first-order, no explicit PEG-desorption/corona sub-steps
or extracellular clearance; parameters are literature-typical, not fitted). It's a
scaffold for reasoning and Track-A fusion, not a calibrated predictor. See
`docs/DELIVERY_KINETICS.md`.
