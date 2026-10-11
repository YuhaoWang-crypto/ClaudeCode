# Extending netjac

## Add a kinetic model
1. Find the BIOMD id (PubMed paper → BioModels curated list / `biomodels`
   GitHub org). `KineticModel("BIOMD00000xxxxx")` fetches and caches it.
2. Sanity checks (copy `assets/new_model_template.py`):
   - `m.state_ids`, `m.n_cons` (expected conservation laws?)
   - `m.enzyme_params()`: every reaction should get its Vmax knob. Reactions
     with fixed influx (e.g. Wolf `v0`) have no elasticity, but they DO have a
     total fragility.
   - `m.control(x)`: summation and connectivity residuals < 1e-6, and the
     analytic result agrees with `control_fd`.
   - stable? `sp.classify`; oscillating? use `m.fixed_point` and
     `sp.critical_pair`.
3. Known failure modes: rate rules or non-unit compartments make
   `_check_rhs` raise (extend `core` instead of silencing it). A model with
   events needs care: the steady state may not be defined.

## Add a phenotype (disease) to a model
- Express it as parameter overrides (as in D1 `phenotypes`). To test "same
  concentrations, different network", add a compensation solve (brentq on
  one parameter to restore a marker) and report max |Δln x|.
- Always run a severity sweep so a perfect AUC at the extreme is not
  presented as the typical case.

## Add a readout
- Static: any function of the steady state, with assay noise via `cohort.assay`.
- Challenge: scale an enzyme / add a bolus. Check the challenge does not
  collapse the system (count failures per group).
- Dynamic: `relax.pulse_response` → `relax.fit_modes`.
- Score with `discriminate.auc` (bootstrap CI). Any selection goes through a
  discovery/validation split, and panels through `discriminate.cv_auc`.

## Add a real dataset
- Put a loader in `netjac/data/<name>.py` that caches under
  `figures/_data_cache/` (git-ignored, so raw data is not committed) and states
  the license and citation in its docstring.
- Large archives: copy the HTTP-range zip reader in `data/cgmacros.py`.
- Compare against the clinical standard, including the partial correlation.

## Natural next steps (not done yet)
- Patient-specific Jacobians: TCGA/CPTAC enzyme levels → Vmax scaling of a
  human kinetic model (e.g. hepatic/erythrocyte glycolysis on BioModels) →
  tumour vs normal fragility profiles.
- Group-level inverse Jacobian on a real Metabolomics Workbench study,
  restricted to a non-stiff sub-module.
- CGM: per-meal mixed model (carbs/fat/fiber covariates); a minimal-model fit
  (two eigenvalues: S_G and insulin action) instead of a single rate; test in
  larger CGM cohorts.
