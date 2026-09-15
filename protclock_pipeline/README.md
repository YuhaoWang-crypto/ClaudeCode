# protclock-pipeline

A runnable reimplementation of the analysis in
[*Integration of proteomic aging clocks in a phase 2a clinical trial supports
simultaneous geroprotective assessment*](https://www.nature.com/articles/s41587-026-03286-y)
(Zhavoronkov et al., Nature Biotechnology, September 2026): six proteomic
aging clocks read side by side on serum from a 12-week trial of rentosertib
in idiopathic pulmonary fibrosis.

**The paper's numbers are not reproducible from public material.** The trial
proteome is controlled-access, the reference cohort needs a UK Biobank
application, and three of the six clocks ship no usable weights. So this
package does the next best thing: it implements the whole analysis and
validates it against a synthetic cohort whose answer is known in advance,
ready to switch to the real data the moment access is granted. See
[`REPRODUCTION.md`](REPRODUCTION.md) for the full assessment.

## Run

```bash
pip install numpy pandas scipy statsmodels scikit-learn matplotlib openpyxl
python3 -m protclock_pipeline.run_all              # ~2 min, full run
python3 -m protclock_pipeline.run_all --quick      # ~30 s, skips the power sweep
python3 -m protclock_pipeline.run_all --figures    # also writes figures/protclock/
python3 -m protclock_pipeline.m4_ageaccel          # or any single module
```

With real clock coefficients on disk:

```bash
python3 -m protclock_pipeline.run_all --weight-dir data/weights
```

## Modules

| Module | Does | Key result on synthetic data |
|---|---|---|
| `m1_cohort` | builds Olink-shaped cohorts with a planted answer | 2,841 assays, 326 responsive, true enrichment OR 1.730 |
| `m2_preprocess` | LOD censoring, QC, imputation, plate centring, harmonisation | 2,841 assays in, 2,832 shared panel out |
| `m3_clocks` | the six clocks, published weights or surrogates | held-out r 0.42 to 0.79, MAE 3.8 to 4.2 yr |
| `m4_ageaccel` | age acceleration on one scale for all six | 5 or 6 of 6 track the planted offset, r 0.28 to 0.86 |
| `m5_trialstats` | within-patient change, mixed model, concordance | null false-positive rate 0.067 vs nominal 0.05 |
| `m6_enrichment` | differential abundance, aging set, Fisher test | OR 1.22 or 2.22 depending on contrast |
| `m7_pathways` | over-representation, mean shift, aging-aligned shift | senescence rejuvenated q=1e-4, immune control null |
| `real_data` | adapter for OMIX008341 and a real reference cohort | untested, the files are not obtainable |

## What the numbers mean

Every biological number the demo prints was written into the simulation by
`m1_cohort` and then measured back out. That validates the code and says
nothing about rentosertib. The pipeline's own behaviour, the calibration, the
power curve, the contrast comparison, is real and reproducible.

Four findings came out of building it, each a defect caught by the
validation rather than a design choice:

- Seeds must reach the cohort simulators or resampling replicates silently
  share one draw, which made a null cohort look like a 0.585 SD effect.
- Expressing a mortality clock's output "in years" divides by a slope near
  0.01 and inflates its noise roughly tenfold against the age clocks.
- A plain mean shift cannot test an aging pathway, because members move in
  opposite directions and cancel; project onto each protein's aging
  direction instead.
- Thresholded enrichment is biased upward, and the choice between a
  treated-vs-placebo and a paired contrast moves the altered-protein count by
  an order of magnitude from identical data.

## Getting the real data

`data/README.md` lists the three inputs, where each lives, which are
obtainable today, and the identifier-namespace trap that will silently score
a clock on a fraction of its proteins if you skip the match-count check.
