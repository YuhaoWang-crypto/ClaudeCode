# protclock-pipeline

A runnable reimplementation of the analysis in
[*Integration of proteomic aging clocks in a phase 2a clinical trial supports
simultaneous geroprotective assessment*](https://www.nature.com/articles/s41587-026-03286-y)
(Zhavoronkov et al., Nature Biotechnology, September 2026): six proteomic
aging clocks read side by side on serum from a 12-week trial of rentosertib
in idiopathic pulmonary fibrosis.

**Partly reproducible.** The trial proteome is controlled-access, so the
trial-level numbers cannot be recomputed. But the paper releases its clock
library with real weights for three of the six clocks, and its supplementary
workbook is openly downloadable and carries the per-sample clock predictions.
So this package runs the real clocks, analyses the real tables, and validates
the rest against a synthetic cohort whose answer is known in advance. See
[`REPRODUCTION.md`](REPRODUCTION.md) for the full assessment.

The most consequential real-data finding: the six clocks correlate at +0.622
on average over the 168 real trial samples, which is three effective
independent clocks, not six. The paper's central argument is that six
independently built clocks agree. Six of six agreeing gives a binomial p of
0.031 under independence and 0.25 at three effective clocks.

## Run

```bash
pip install numpy pandas scipy statsmodels scikit-learn matplotlib openpyxl
python3 -m protclock_pipeline.run_all              # full run, ~3 min
python3 -m protclock_pipeline.run_all --quick      # skips the power sweep
python3 -m protclock_pipeline.run_all --figures    # writes figures/protclock/
python3 -m protclock_pipeline.m8_supplementary     # REAL data only, ~20 s
```

For the three clocks with real published weights, install the paper's own
library first:

```bash
git clone https://github.com/Insilico-org/proteoclock
cd proteoclock && pip install -e . && pip install --upgrade scikit-posthocs
```

That last upgrade is required: `setup.py` pins `scikit-posthocs==0.11`, which
imports `multipletests` from a statsmodels location removed in 0.15, so the
package will not import without it.

## Modules

| Module | Does | Key result on synthetic data |
|---|---|---|
| `m1_cohort` | builds Olink-shaped cohorts with a planted answer | 2,923 real gene symbols, 326 responsive, true OR 1.730 |
| `m2_preprocess` | LOD censoring, QC, imputation, plate centring, harmonisation | 2,923 assays in, 2,832 shared panel out |
| `m3_clocks` | the six clocks, published weights or surrogates | held-out r 0.42 to 0.79, MAE 3.8 to 4.2 yr |
| `m4_ageaccel` | age acceleration on one scale for all six | 6 of 6 track the planted offset, r 0.48 to 0.96 |
| `m5_trialstats` | within-patient change, mixed model, concordance | null false-positive rate 0.067 vs nominal 0.05 |
| `m6_enrichment` | differential abundance, aging set, Fisher test | OR 1.22 or 2.22 depending on contrast |
| `m7_pathways` | over-representation, mean shift, aging-aligned shift | senescence rejuvenated q=1e-4, immune control null |
| `m8_supplementary` | **REAL data**: the paper's tables S2, S5, S8 | mean clock r +0.622, 3.0 effective clocks of 6 |
| `proteoclock_backend` | the paper's released library | real weights for PAC and both OrganAge clocks |
| `real_data` | adapter for OMIX008341 and a real reference cohort | untested, the files are not obtainable |

## What the numbers mean

M8 is real. Everything in M1 through M7 is simulated: those biological
numbers were written in by `m1_cohort` and measured back out, which validates
the code and says nothing about rentosertib.

Five findings came out of building it, each a defect caught by the
validation rather than a design choice:

- A namespace mismatch sent all three real clocks silently to surrogates
  while reporting success, because a two-column tab-separated symbol file was
  read a whole line at a time. Matched-protein counts are now printed on
  every clock load.
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
