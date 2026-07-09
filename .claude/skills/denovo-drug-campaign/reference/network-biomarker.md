# Stage 2b — network-biomarker layer (irreducibility + critical transition)

Independent cross-check of the target axis using dynamical-systems theory. Three
rigorous outputs; one honest hypothesis boundary.

## 1. Fibration -> irreducible core  (`input_fibers`)
Minimal balanced colouring by signed in-edge signature (Morone-Leifer-Makse 2020).
Nodes with identical signed inputs synchronize -> collapse to one fiber. Payoff:
each fiber = one pan-assay representative readout; the hub fiber names the
irreducible core target. Deterministic -> **rigorous**.

```python
fib = input_fibers(config["pathway_edges"])   # {signature: [nodes]}
# NSCLC: 5 nodes -> 3 fibers {KRAS,KEAP1} inputs / {NRF2} hub / {SLC7A11,GLS} effectors
```

## 2. CRNT deficiency -> switch capacity  (`crnt_deficiency`, `schlogl_states`)
delta = n - l - s on the EFFECTIVE feedback core (strip chemostatted species first —
this is the real bug that makes Schlogl falsely delta=0). Feinberg: delta=0 +
weakly reversible => bistability impossible; delta>=1 => permitted. ALWAYS confirm
by integrating the mass-action ODE from multiple initial conditions.

```python
delta = crnt_deficiency(n_complexes=4, linkage_classes=2, stoich_rank=1)   # =1
roots, stab = schlogl_states(k1=0.4,k2=0.03,k3=0.3,p=1.05)   # 2 stable + 1 saddle
```

## 3. Critical slowing -> early-warning biomarkers  (`critical_slowing`, `langevin_earlywarning`)
As drive p -> saddle-node fold, leading eigenvalue lam->0, recovery tau->inf, and
variance + lag-1 autocorrelation rise. `critical_slowing` gives the analytic
trajectory; `langevin_earlywarning` MEASURES var/ar1 on a simulated SDE (the
model-free statistic — the rigorous, transferable claim). Track the DISAPPEARING
(OFF/sensitive) branch, branch="lo".

```python
near = langevin_earlywarning(0.4,0.03,0.3, p=0.70, seed=3)   # near fold
far  = langevin_earlywarning(0.4,0.03,0.3, p=1.30, seed=3)   # far
fold_change = near["var"]/far["var"]   # ~5-6x variance rise toward tipping
```

## Honest boundary (label in every output)
- **Rigorous:** fibration (deterministic), delta (integer invariant), bistability
  (ODE-confirmed), critical-slowing GEOMETRY (lam->0, var/ar1 rise — analytic + measured).
- **Hypothesis:** mapping the high attractor to a specific CLINICAL state (e.g.
  therapy resistance) needs patient-data calibration; canonical Schlogl rate
  constants are illustrative, not fitted to the real pathway kinetics.

## Decimal-exact upgrade
Replace illustrative constants with a curated literature model: fetch from the
BioModels GitHub mirror (`raw.githubusercontent.com/biomodels/{id}/master/{id}/{id}.xml`)
+ libRoadRunner, read `rr.getFullJacobian()` for lam_max. See the network-biomarker
skill's data-access reference.
