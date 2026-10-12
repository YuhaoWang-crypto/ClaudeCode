---
name: network-jacobian
description: >-
  Network-Jacobian diagnostics with the `netjac` platform: turn any SBML
  kinetic model (BioModels) plus omics/time-series data into (1) EDGE
  biomarkers (elasticity changes invisible in concentrations), (2)
  PERTURBATION-RELAXATION diagnostics (recovery time / ringing = eigenvalues,
  critical slowing), and (3) FRAGILITY profiling (which enzyme/edge/probe
  controls the critical eigenvalue). Use when asked to compute a Jacobian /
  flux-control (dFlux/dEnzyme) matrix, find network-state or resilience
  biomarkers, design a challenge test or a multiplex probe panel, infer a
  Jacobian from metabolomics covariance, analyse CGM/OGTT recovery, or extend
  netjac with a new model or dataset. Enforces ✅-rigorous vs ⚠️-hypothesis
  labeling.
---

# Network-Jacobian diagnostics (`netjac`)

`netjac/` turns a curated kinetic model into three DIFFERENT matrices. Keep
them apart, because most confusion in this area comes from mixing them up.

| Object | Definition | What it is for | Where |
|---|---|---|---|
| state Jacobian `Jz` | ∂(dx/dt)/∂x on the conservation-reduced subspace | eigenvalues → stability class, relaxation time, Hopf/saddle-node | `core.reduced_jacobian` |
| elasticities `ε` | ∂v/∂x (reactions × species) | the **edges**: enzyme r senses species k (substrate, product, allosteric) | `core.elasticities` |
| control matrix `C^J` | ∂lnJ/∂lnE (= the "dFlux_i/dEnzyme_j" matrix) | which enzyme controls which flux (MCA) | `core.control` |

They are linked exactly: `C^S = −B Jz⁻¹ BᵀN`, `C^J = I + ε C^S`. The eigenvalues of
`C^J` have **no** stability meaning. Use `Jz`.

## Quick start

```bash
pip install numpy scipy networkx matplotlib libroadrunner python-libsbml scikit-learn pandas
python3 -m netjac.run_all            # D1+D2+D3, figures → figures/netjac/
python3 -m netjac.run_all --offline  # skip the real CGMacros download
```

```python
from netjac.core import KineticModel
from netjac import spectrum as sp
m = KineticModel("BIOMD0000000064")      # any BIOMD id or SBML path
x = m.steady_state()                     # or m.fixed_point(x) for unstable states
c = m.control(x)                         # eps, Jz, C^J (scaled), theorem residuals
sp.classify(c["Jz"])                     # kind, λ_lead, τ, period
sp.total_fragility(m, "vGLT", x, solver="integrate")   # dReλ/dlnE incl. state shift
```

## The three diagnostic directions and what the demos showed

1. **Edge biomarkers** (`demos/d1_edge_biomarker.py`). A phenotype can leave
   the metabolome *exactly* unchanged (max |Δln x| = 9e-16) and still remove
   an edge. Static nodes, ratios and an 18-analyte panel cannot see it
   (AUC 0.52–0.62). A **model-predicted two-time-point readout** around a small
   standard challenge can (ΔlnF6P/ΔlnATP after ATPase ×0.8, AUC 0.91). Passive
   covariance inversion (Lyapunov equation) is exact with true C but needs
   ~10⁶ snapshots on a stiff network.
2. **Perturbation-relaxation** (`demos/d2_perturbation_relaxation.py`). A
   bolus with 1–2 observables gives λ via Hankel-DMD (init) followed by
   `relax.fit_modes` (global, unbiased). Near a Hopf, read the **damping of
   the ringing mode**, not the overall recovery, because a non-critical real
   mode can be slower. Real CGMacros CGM data: meal recovery rate falls with
   HbA1c (ρ=−0.44), but adds nothing beyond fasting glucose at n=45.
3. **Fragility** (`demos/d3_fragility.py`). Use the **total** fragility
   (local + steady-state shift); the local ∂λ/∂J formula alone mispredicted
   (r=0.16 vs 0.99). Fragile enzymes need not carry flux control (vPDC: C^J≈0,
   yet a 30 % change collapses glycolysis). The best probe is the species with
   the largest critical-mode participation |u_k|/x_k (ρ=+0.96 with AUC), not
   the most abundant one (ρ≈0).

4. **Patient-specific Jacobians** (`demos/d4_patient_jacobian.py`,
   `data/tcga.py`). TCGA-KIRC: 609 samples, each instantiated on the
   Holzhütter 2004 human glycolysis/PPP chassis (BIOMD70) by Vmax ∝
   isoform-summed TPM / median normal. The pipeline works (theorems hold),
   and it recovers Warburg shifts and an ALDOB-loss aldolase bottleneck (10 %
   of tumours have no steady state). However, every survival signal
   dissected to a single-gene proxy (λ_lead ≈ BPGM, collapse ≈ ALDOB).
   The chassis is the bottleneck.

5. **Pairwise → network attractors** (`loops.py`, `demos/d5_*`, `demos/d6_*`,
   method in `reference/loop-analysis.md`). Find every equilibrium and
   classify it. Analyse every k×k sub-system, then combine them EXACTLY
   through the Levins loop decomposition of the characteristic polynomial
   plus Routh–Hurwitz. Add bifurcation scans and a global portrait. On
   real data, use VAR → Jacobian. Pairwise stability neither implies nor
   excludes whole-network stability (Goodwin, Sturis, Topp).

## Mandatory discipline

- Label every claim **✅ rigorous** (deterministic computation, theorem check,
  analytic = finite difference, decimal-exact model, measured statistic on
  real data) or **⚠️ hypothesis** (synthetic cohort, disease analogy,
  extrapolation to patients). Virtual cohorts are always ⚠️.
- Any data-driven feature selection uses a discovery/validation split (or
  CV). Report the pre-specified, model-predicted readout separately.
- Always compare a new readout against the boring baseline (best single node,
  full static panel, and on real data the clinical standard, e.g. fasting
  glucose, including **partial** correlation given that standard).
- Report negative results (e.g. covariance inversion at realistic N, CGM
  recovery adds no information beyond fasting glucose) as findings.

## Hard-won gotchas (real bugs / traps from building this)

- Conservation laws (NAD+NADH, adenylates) give structural zero eigenvalues.
  Always use `Jz = BᵀNεB`, never the raw full Jacobian.
- Teusink has assignment-rule species (ATP/ADP/AMP from SUM_P, P). The ODE
  state is the non-rule species only, and roadrunner's built-in MCA refuses
  this model. `netjac.core` handles it generically.
- Reaction-local parameters must be promoted to globals (`vPFK_VmPFK`) to be
  addressable. `core` does this through libsbml.
- Enzyme knobs are auto-detected numerically: the parameter set that scales
  v_r linearly. Two-parameter Vmax (GAPDH f/r) is handled as a group.
- Diagnostic challenges must stay inside the basin. ATPase ×0.6 drove 30 % of
  virtual individuals into hexose-phosphate runaway; ×0.8 drove ~5 %.
- Plain DMD is biased by measurement noise. Use it only to initialise
  `fit_modes`. Closely spaced eigenvalues (−4.0 vs −7.2) remain hard: report
  the spread over noise realisations.
- Pick the RIGHT eigenvalue: `spectrum.critical_pair` / `target=` follows the
  Hopf pair even when a real mode is slower.
- Chassis modes: on a borrowed kinetic model the slowest mode may be
  chassis physiology (RBC 2,3-BPG pool), not the tissue's. Label modes by
  pathway (`d4.pathway_modes`) before calling anything "resilience".
- Every network feature that predicts an outcome must be dissected against
  its strongest single-gene correlate (`d4.dissect`). On KIRC, all of them
  reduced to one gene.
- Eigenvalue features have extreme outliers. Use rank-based inverse-normal
  covariates in Cox. A post-hoc p = 0.0003 vanished (p = 0.2) after doing so.
- GDC clinical data: follow-up time is in `follow_ups.days_to_follow_up`,
  not only in `diagnoses`. Missing it dropped 55 % of censored patients.
- BioModels GitHub mirror covers BIOMD 1–548 only.

## Where to look

- Public data sources (what is reachable from this sandbox) → `reference/data-sources.md`
- Math, definitions, theorem checks, honest limits → `reference/methodology.md`
- Add a model / dataset / readout → `reference/extending.md` + `assets/new_model_template.py`
- Full Chinese write-up with all numbers → `NETJAC_REPORT.md` (repo root)
- Sister skill for symmetry/CRNT/EFM/critical-slowing atlas → `network-biomarker`
