# Defining a state by its perturbation-recovery behaviour (D11)

`netjac/recovery.py` is the engine and `demos/d11_recovery_state.py` holds the validated cases. Results are in `NETJAC_REPORT.md` §12.

## Protocol for a new cell, tissue or indication (transfer recipe)

1. **Fix the known positives first.** Use a defect class, outcome or gene family defined independently of the readouts. Add decoys: negatives with a large static difference.
2. **Build nested feature sets on the same subjects:**
   - M0 static;
   - M1 + one standardised stimulated endpoint;
   - M2 + a second dose;
   - M3 + recovery kinetics after the stimulus is removed.
3. **Compare them** with `incremental_classification`, or `incremental_survival` for outcomes. The comparisons are:
   - paired bootstrap of ΔAUC / ΔC-index;
   - Δ false-positive rate at a fixed sensitivity;
   - groups = donor / family for leave-group-out (transfer).

   Stop at the first step that adds nothing.
4. **Core components:** total fragility of the mode the defect changes (`spectrum.total_fragility(..., target=)`). For a measured response network, use out-/in-strength (`core_ranking`).
5. **Minimal control set:**
   - structural: `driver_nodes` (maximum matching) says how many inputs FULL control needs;
   - functional: `greedy_restoration_set` says how few knobs restore one target recovery signature.

   Include the suspected true fix as a positive control.
6. **Restoration check:** the objective must contain the DYNAMIC signature (challenge response + recovery). A static objective can be satisfied while the dynamics are still broken.

## What the validated cases say (when dynamics reduce misclassification)

| Case | Static | + stimulus / response | Verdict |
|---|---|---|---|
| Teusink, compensated feedback defect (ground truth) | AUC 0.66, FP 74 % | one challenge: AUC 0.996, FP 0.4 % | dynamics essential when compensation hides the defect |
| K562 homeostats vs co-regulated machines | AUC 0.59 (unseen family 0.06) | 0.81 (unseen family 0.63) | response transfers, static co-fluctuation does not |
| Norman 2019: which double restores the state | static co-expression 0.55 | additive single responses 0.89 | single responses predict direction; magnitude is epistatic |
| NHANES OGTT → mortality (10,748 people, 1,020 deaths) | fasting + HbA1c | + 2-h OGTT: ΔC +0.000 (all-cause), +0.010 n.s. (diabetes deaths) | no gain once the static marker saturates |
| CGMacros meals → T2D (n = 45) | fasting AUC 0.92 | recovery kinetics 0.91 | no gain; small n overfits |
| GSE75210 macrophages: which inputs route through MyD88 (real time series) | t=0 AUC 0.50 | any stimulated time point ≈ 1.0; trajectory area 1.0, replicate r 0.91, 0 errors on held-out plate | a functional state exists only under stimulation |
| GSE75210: is this macrophage Myd88⁻/⁻ | modules AUC 1.0, transcriptome 0.93 | + trajectory: 1.0 / 0.80 (worse) | genotype is compositional: static suffices |

The rule that emerges: dynamics help when the system COMPENSATES at the set point (a hidden feedback, homeostats), and do not help when the static marker already integrates the defect (HbA1c, fasting glucose).

A second distinction from D12: is the state you want to define **functional** (how the object answers an input — routing, reserve, feedback) or **compositional** (what it is made of — genotype, cell type, basal expression)? Functional states need a perturbation by definition. Compositional states are usually visible statically.

## Restoration on a real time series (D12, `demos/d12_timecourse_cell_state.py`)

1. **Target**: the WT response to the inputs whose dependence was lost (genes × time, relative to t = 0).
2. **Positive control**: each WT input scored against the mean of the others (leave-one-out). It defines "fully restored" at the data's noise level.
3. **Negative control**: the defective object's own response to those inputs.
4. **Scores**: score each candidate rescue by
   - **gain** = projection / |target|²;
   - **shape** = cosine;
   - **off-target cost** = response in a module the target does not use (here the interferon arm).
5. **Result**: bypass inputs (LPS, zymosan, Poly(I:C), MDP+LPS) reach the positive-control gain and shape but add +3 to +8 log2·h interferon. The minimal restoration set is therefore two knobs: a bypass input plus suppression of the side programme.

## Pitfalls hit

- **Euclidean restoration distance is noise-dominated on low-input time courses.** Myd88⁻/⁻ cells' "no response" vector already had a norm of 1.07 × the target. As a result, doing nothing scored better (0.30) than every bypass, and the positive controls only reached about 0.45. Use projection and cosine, and calibrate with positive and negative controls.
- **Data-driven "response genes" on 2-replicate time courses replicate at r ≈ 0.** Fix the modules from prior knowledge before looking at the data.
- **Genotype confounded with plate.** In GSE75210 each genotype has its own plates, so a basal difference (NF-κB +1.8 log2 at t = 0) cannot be separated from batch. Hold out whole plates.

- **Restoration needs multi-start.** A feedback strength has a flat objective near the defective value, so plain L-BFGS from 0 misses the true fix (`starts=`).
- **Decoys must stay inside the basin.** Teusink GLT ×1.4 sent 92 % of the decoys into runaway; use ×1.15.
- **Rare-event Cox folds return NaN coefficients** (47 diabetes deaths). Refit with BFGS (`oof_cox` does this).
- **NHANES 2005–16 has no `LBXIN` in GLU_D…I** under that name. Fasting insulin is not needed for these models.
- **Zenodo throttles to about 100 KB/s.** The Norman 2019 GEARS release on Harvard Dataverse (doi:10.7910/DVN/Q2ZV3E) downloads in seconds.
- **`pkill -f <pattern>` kills the calling shell** when the pattern appears in its own command line (exit 144). Select processes by `pgrep -x curl` plus a cmdline check.
