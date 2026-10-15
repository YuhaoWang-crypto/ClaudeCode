# The five steps in detail (`netjac/statepipe.py`)

Notation:
- x = log-scale profile (genes or proteins) of one sample.
- ref_s = balance reference of subject s.
- x0_s = the subject's first perturbed / disease sample.
- dis = disease features.
- u_s = leave-one-out cohort disease axis.

## 1 Steady state — `steady_state(D)`

1. **Balance reference.** Use the subject's own reference samples when they exist: non-lesional tissue, pre-challenge samples, the same plate. Otherwise use the cohort reference mean (control tissue, healthy volunteers).
2. **Disease features.** Test x0 − ref per feature with a one-sample t-test and BH q < 0.05, and require |log2FC| > 1. If fewer than 50 features pass, fall back to 0.585 and then 0.3. The test never sees outcome labels.
3. **Leave-one-out axis.** u_s is the mean deviation of all other subjects. Never use the subject's own deviation: its measurement noise then reappears as "restoration" (regression to the mean).
4. **Displacement in noise units.** Median on-axis displacement ÷ noise floor. Take the noise floor from repeated references within the same subject if available. Otherwise use the between-subject spread of the reference samples, which is only an upper bound.
5. **Pre-specified exclusion.** Drop a subject whose baseline on-axis displacement is < 25 % of the cohort median. Its sample carries no disease signature, so the restoration ratio is undefined. Without this rule one psoriasis patient scored R = 40.
6. **Is the disease state itself a steady state?** Check the untreated or placebo arm: its restoration over time should stay near 0. In the AD placebo arm it was −0.04 at week 4. If it drifts, model the drift (b0 in step 4).

## 2 Order parameters — `order_parameters(D, S, modules)`

- **φ (disease coordinate).** φ = ⟨x − ref, u⟩ / ⟨x0 − ref, u⟩ on the disease features. φ = 1 is the subject's own disease state and φ = 0 is balance. Projection is robust and distance is not: on low-input data the noise norm can match the signal norm (D12).
- **Module coordinates.** Use the same projection restricted to the module's disrupted features. Normalise by the cohort-median baseline, not per subject: small own baselines amplify noise.
- **Label-free check.**
  - Count collective modes of the state deviations against a column-shuffle null (`statmech.collective_modes`). Report the participation ratio.
  - Check whether a top PC aligns with φ (|r|). If PC1 aligns, the state is a collective mode, and one high-loading prior gene is a sufficient readout (D13: RSAD2 / IFI44L).
  - If no PC aligns, φ is not the dominant source of variance (e.g. AD: inter-individual heterogeneity dominates). Report it.

## 3 Minimal core network — `core_network(D, S, OP)`

- **Modules in the state.** A module enters only if ≥ 3 of its genes are disrupted (|cohort deviation| > 0.585). A module that is not disrupted has no defined restoration; housekeeping genes make the built-in negative control.
- **Couplings.** Partial correlations among module coordinates across all state samples (Ledoit-Wolf). Directed couplings come from J in step 4.
- **Core genes.** Disrupted in every replication group separately (strain, drug, centre), same sign, p < 0.01 in each. Compare with a label-permutation null (D13: 111 genes vs null median 0).
- **Perturbation screens.** With a Perturb-seq-type dataset, get directed gene-level edges from MRA / regression on knockdown responses (`statmech.mra_regression`, D9 / D10) instead.

## 4 Dynamics — `fit_dynamics(D, OP)`

Model, in module coordinates (0 = balance; y_D = the subject's untreated state):

    dy/dt = J (y − y_D) + b0 + b·u        y(t) = y* + e^{Jt}(y0 − y*),   y* = y0 − J⁻¹(b0 + b u)

- **Fit.** Least squares on whole trajectories from each subject's first sample, using the exact solution, so irregular intervals (weeks 1, 2, 4, 16) are handled correctly.
- **Stability prior.** J = −(PPᵀ + εI) + (Q − Qᵀ). The symmetric part is negative definite, so every mode decays. An unconstrained J fitted to 1–3 intervals per subject was unstable and predicted worse than "no change". Ridge on the couplings, λ by 5-fold subject-level CV.
- **b0 vs b.** b0 is natural-course drift and is identifiable only with a control arm. b is the intervention-specific direct input. Without a control arm, b absorbs regression to the mean and standard care: say so.
- **Validation.** Subject-level cross-validated RMSE against two baselines: no change from y0, and the mean trajectory of the same arm. The arm-mean baseline is strong. The model must beat no-change; matching the arm mean is acceptable, because individual information only enters through y0.
- **Read-outs.**
  - Time constants −1/Re λ.
  - The slow mode: which modules lag; this predicts the residual.
  - b: where the intervention acts directly. Positive control: dupilumab → Th2.
  - y*: predicted on-treatment steady state.

## 5 Minimal adjustment — `minimal_adjustment(dyn)`

- **Functional.** A sustained direct adjustment u_S of a module set S moves the steady state to y = y_D − J⁻¹ B_S u_S, with B_S the unit vectors of S and u_S fitted by least squares.
  - Restoration = 1 − |y| / |y_D|.
  - Greedy forward selection until restoration ≥ 0.9.
  - For each single module, also report its **propagation**: restoration of the OTHER modules. Weak couplings mean every module must be adjusted itself, so the minimal set is large.
- **Ready-made skin library.** `netjac.data.keratinocyte.lever_library()`
  (nine cytokines in primary keratinocytes, IL-4 included) with
  `rank_levers`, or the template `assets/keratinocyte_levers.py`. Columns from
  different studies are directionally comparable (each is scored against the
  disease vector on its own) but NOT magnitude-comparable, so run a multi-lever
  least-squares solve within one source; every column carries its `source`.
  Duplicated levers are taken from GSE36287 (3 donors, same platform as the
  atopic-dermatitis cohort) and the duplicate is listed in `also_in`.
- **Structural.** Driver nodes of the |J| sparsity pattern by maximum matching (`recovery.driver_nodes`). This says how many inputs FULL control would need.
- **Applied intervention as a lever.** Put b into the same formula. The model's restoration should match the observed one.
- **Predicting a lever for a drug that is not in the data.** Use J from any intervention dataset. The prediction is only as good as J's couplings; if propagation is ~0, the model cannot rank levers beyond "the module itself".

## Decision rules (from D11–D15)

| Question | Rule |
|---|---|
| Static data enough to define the state? | Yes for compositional states (genotype, cell type, HbA1c-like integrators). No for functional states (routing, susceptibility, reserve, compensated feedback): standardised perturbation needed |
| Which readout? | Whole-trajectory area or projection > single time point > static. Prior modules beat data-driven picks at small n |
| Has the object returned to balance? | Restoration ≥ 0.75 along the LOO axis, AND no increase in the off-axis deviation, AND a control arm stays ≈ 0 |
| Which component to target? | Module with the largest direct input b in a placebo-controlled fit (mechanism), and the largest propagation in J (leverage). Check that an untargeted module stays (specificity) |
| What will remain? | The slow mode of J and the modules with the lowest restoration in responders: structural (barrier, epithelium) in every tissue so far |
