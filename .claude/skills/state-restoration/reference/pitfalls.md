# Pitfalls that produced wrong answers before they were fixed

Each entry gives the symptom, the cause and the fix, with the case where it happened.

## Measuring restoration

| Symptom | Cause | Fix | Case |
|---|---|---|---|
| Doing nothing scores better than every rescue | Euclidean distance on low-input data: the "no response" noise norm was 1.07 × the target | Projection (gain) + cosine, calibrated with positive (WT leave-one-out) and negative (KO) controls | D12 GSE75210 |
| Housekeeping genes "restore" 0.5–0.7; placebo +0.13 | Own disease deviation used as the axis, so baseline noise regresses to the mean | Leave-one-out cohort axis → placebo −0.04 (paper: −10.5 %) | D14 GSE130588 |
| One patient restores 40× | Lesional biopsy with no disease signature (ratio of ~0) | Pre-specified exclusion: baseline projection < 25 % of cohort median | D14 GSE85034 |
| Module restoration erratic | Module genes not disrupted, or a per-subject small baseline | Use only module genes with \|cohort deviation\| > 0.585 (≥ 3); normalise modules by the cohort-median baseline | D14, D15 |
| Placebo n shrinks at late visits | Dropout and rescue therapy | Make the early visit the primary comparison (week 4: placebo n = 17; week 16: n = 7) | D14 |

## Defining the state

| Symptom | Cause | Fix | Case |
|---|---|---|---|
| 4,030 "response genes" with replicate r ≈ 0 | Data-driven gene selection on a 2-replicate time course | Fix modules from prior knowledge before looking | D12 |
| Data-driven top genes fail on a new strain (AUC 0.46) and do not beat random gene sets | n = 17; the state is a collective mode, so many genes carry it | Take prior high-loading genes (RSAD2, IFI44L: 0.85–0.90) | D13 |
| Static genotype "difference" of 1.8 log2 | Genotype = plate (confounded) | Hold out whole plates; label it as possibly batch | D12 |
| Non-lesional vs healthy skin differ in 2,661 genes | Healthy volunteers likely from a different batch | Prefer the subject's own reference; flag cross-batch references | D14 |
| Static marker "works" in one virus (0.80) and is inverted in another (0.38) | No functional information in the static state | Report the per-group values; never pool them silently | D13 |
| The perturbed state has not returned at the last sample | Sampling window ends before recovery (IFN +13 SD at 108 h) | Call it a window limit, not "non-recovery" | D13 |

## Dynamics and control

| Symptom | Cause | Fix | Case |
|---|---|---|---|
| Unstable J; CV error 5× worse than "no change" | Full J fitted by finite differences on 1–3 intervals per subject | Exact-solution trajectory fit, stable parametrisation, ridge on couplings, subject-level CV | D15 |
| Drug effect overestimated without a control arm | b absorbs natural course and regression to the mean | Fit b0 only when a control arm exists; otherwise say b = drug + natural course | D15 |
| Every module is its own best lever | Couplings ≈ 0 (J nearly diagonal) | Report propagation; a lever prediction needs measurable off-diagonal J | D15 |
| Full structural control needs about half of all nodes | Sparse directed networks (Perturb-seq: 51 % drivers) | The target is functional: restore one signature with the fewest knobs | D11 |

## Engineering

- **GEO series matrices.** Some are linear MAS5 (GSE52428, GSE16879): use log2(max(x, 16)), then quantile-normalise. Others are already log scale. Check the range before choosing.
- **BH correction.** NaN p-values from constant floor genes make every q NaN. Set them to 1 first.
- **Challenge cohorts overlap.** GSE17156's H3N2 arm is the GSE52428 DEE2 cohort. Drop it before calling a test independent.
- **`pkill -f <pattern>` kills the calling shell** if the pattern appears in its own command line, including a background task you just started. Select processes by PID.
- **Verify every PMID before writing it into a report.** Two of the citations written before checking were wrong:
  - GSE75210 is O'Connell 2016, not Tong.
  - GSE16879 is Arijs, *PLoS One* 2009, not *Gut*.
