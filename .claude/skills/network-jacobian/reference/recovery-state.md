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

The rule that emerges: dynamics help when the system COMPENSATES at the set point (a hidden feedback, homeostats), and do not help when the static marker already integrates the defect (HbA1c, fasting glucose).

## Pitfalls hit

- **Restoration needs multi-start.** A feedback strength has a flat objective near the defective value, so plain L-BFGS from 0 misses the true fix (`starts=`).
- **Decoys must stay inside the basin.** Teusink GLT ×1.4 sent 92 % of the decoys into runaway; use ×1.15.
- **Rare-event Cox folds return NaN coefficients** (47 diabetes deaths). Refit with BFGS (`oof_cox` does this).
- **NHANES 2005–16 has no `LBXIN` in GLU_D…I** under that name. Fasting insulin is not needed for these models.
- **Zenodo throttles to about 100 KB/s.** The Norman 2019 GEARS release on Harvard Dataverse (doi:10.7910/DVN/Q2ZV3E) downloads in seconds.
- **`pkill -f <pattern>` kills the calling shell** when the pattern appears in its own command line (exit 144). Select processes by `pgrep -x curl` plus a cmdline check.
