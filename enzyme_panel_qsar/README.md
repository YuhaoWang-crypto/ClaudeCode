# A virtual enzyme-activity panel

In silico counterparts of BPS Bioscience biochemical screening assays, built as a
**panel** rather than one target at a time, and screened against every clinical
and marketed small molecule in ChEMBL.

This generalises the SARS-CoV-2 3CLpro virtual-assay work (BPS 79955) from one
assay to sixteen, across eleven of BPS's assay families. Going wide is not just
more of the same: BPS sells *screening and profiling*, and profiling means a
compound's activity across a panel. One library against one panel gives every
compound a predicted activity profile, which is where selectivity, and the
failure modes that masquerade as selectivity, become visible.

---

## What the panel covers

Sixteen targets, all resolved from a UniProt accession to a ChEMBL SINGLE
PROTEIN record at run time, all past the data gate:

| target | BPS family | assay readout | compounds | scaffolds | noise floor |
|---|---|---|---|---|---|
| BRD4 | Bromodomains | acetyl-lysine displacement (TR-FRET) — a *reader* | 9,107 | 3,610 | 0.66 |
| HDAC1 | Deacetylases/Sirtuins | fluorogenic deacetylation | 8,095 | 3,538 | 0.51 |
| HDAC6 | Deacetylases/Sirtuins | fluorogenic deacetylation | 6,320 | 2,658 | 0.67 |
| PARP1 | PARP, PARG, ADP-ribosylation | ADP-ribosylation / NAD⁺ | 4,479 | 1,992 | 0.71 |
| NAMPT | Metabolic Enzymes | NMN formation, coupled | 4,147 | 2,046 | 0.17 |
| KDM1A | Demethylases | FAD-dependent H3K4 demethylation | 3,841 | 1,625 | 0.34 |
| IDO1 | Immunotherapy / Metabolic | Trp → kynurenine, heme | 3,573 | 1,289 | 0.62 |
| PTP1B | Protein Tyrosine Phosphatase | pTyr dephosphorylation | 3,102 | 1,259 | 0.30 |
| PTPN11 | Protein Tyrosine Phosphatase | DiFMUP, allosteric site | 2,328 | 1,001 | 0.58 |
| PDE4B | Phosphodiesterases | cAMP hydrolysis | 2,223 | 860 | 0.70 |
| PDE5A | Phosphodiesterases | cGMP hydrolysis | 2,148 | 896 | 0.73 |
| CTSK | Proteases | fluorogenic peptide cleavage | 2,144 | 1,105 | 0.54 |
| HSP90AA1 | Heat Shock Proteins | ATP-site binding / ATPase | 1,605 | 637 | 0.71 |
| PRMT5 | Methyltransferases | arginine methylation, SAM | 1,051 | 552 | 1.27 |
| EZH2 | Methyltransferases | H3K27 methylation, SAM | 1,003 | 441 | 0.45 |
| SIRT1 | Deacetylases/Sirtuins | NAD⁺-dependent deacetylation | 616 | 301 | 0.52 |

55,782 compounds in total. Noise floor is in log units and is defined below — it
is the number model RMSE has to be judged against, and getting it right took
two attempts.

Screening library: **9,746 unique clinical or marketed small molecules**
(2,053 marketed, 890 phase 3, 6,126 phase 2, 677 phase 1; 263 withdrawn),
6,246 scaffolds, one row per distinct standardised structure.

---

## Pipeline

```
p01  curate panel      UniProt -> ChEMBL SINGLE PROTEIN, assay-confidence floor,
                       data gate, inter-assay noise floor, raw cache
p02  benchmark         ridge / RF / LightGBM bake-off, scaffold-split headline,
                       permutation null per target, model selection
p03  build library     ChEMBL max_phase >= 1, standardised identically to p01
p04  screen panel      library x every modellable target, AD tiers, PAINS,
                       reference-inhibitor recall
p05  panel report      per-target hits, cross-target selectivity, frequent hitters
```

Each stage is runnable alone and deterministic given the seeds in its header.
Raw ChEMBL activities are cached per target, so re-analysis never re-pulls.

---

## Three things this gets right that are easy to get wrong

These are documented because each one silently produced a plausible-looking
wrong answer first, and the wrong version is the one a reader would not catch.

### 1. A quality field read from the wrong endpoint empties the dataset

`confidence_score` says how securely ChEMBL ties an assay to a molecular target;
below 8 means a homologous-protein or subunit-level assignment. It lives on the
**assay** record. The `/activity` endpoint does not return it.

A first version read it off the activity rows, got `None` for all 80,000 of
them, and filtered every target in the panel to zero compounds — reported as
sixteen clean "FAILED: all_filtered" lines. Nothing errored. The fix is to fetch
the distinct assay IDs from `/assay` in batches and join; the audit now also
records how many assays came back without a score, so the same failure cannot
recur silently.

### 2. The ChEMBL "replicate" distribution is bimodal, and a median lies about it

The noise floor is what tells you whether a model is near the limit of its data.
Measuring it from raw activity rows gave **0.01 log** on PARP1 — an assay
reproducibility no biochemical assay has.

Grouping to one median per (compound, assay) first and taking the range across
distinct assays gave 0.04. Still wrong, because the distribution is bimodal:

| source documents | n | median range | q75 |
|---|---|---|---|
| 1 | 67 | 0.700 | 1.380 |
| 2 | 209 | **0.000** | 0.040 |
| 3+ | 82 | 1.010 | 2.822 |

The two-document group agreeing to *exactly* 0.000 is one primary measurement
re-reported, not two independent assays. Half of all apparent replicates sit in
that mode; a third disagree by more than 0.5 log, with a 90th percentile of 2.0.
A median over both modes understates real assay noise by roughly 20×.

The floor is therefore the median over **non-zero** ranges, reported with the
full quantile set and the exact-tie fraction alongside. The resulting values
(0.17–1.27 log, mostly 0.3–0.7) match published inter-laboratory IC50
reproducibility. PRMT5's 1.27 and NAMPT's 0.17 are both real features of those
datasets, not artifacts, and both are flagged in the per-target report.

### 3. Reliability and novelty are different axes

The 3CLpro work used a single applicability-domain tier, and compounds with
Tanimoto > 0.7 to the training set came out as `out_of_domain`. Its own
limitations section then had to warn the reader that those predictions are in
fact the *most* reliable ones and that the tier reads backwards for repurposing.

That is two questions collapsed into one label. They are separated here:

* `ad_tier` — **reliability**, by maximum Tanimoto to training: `high` ≥ 0.40,
  `borderline` ≥ 0.25, else `out_of_domain`. More similar means the model is
  interpolating. This reads in the ordinary direction.
* `novel_scaffold` — **novelty**, whether the Murcko scaffold is absent from
  training. Interesting for discovery, orthogonal to reliability.
* `in_training` — exact structure match, so memorisation is never mistaken for
  prediction.

A repurposing candidate wants high reliability; a novel-chemotype programme
wants novelty. Neither has to be explained as the inverse of a tier.

---

## Honesty machinery

* **Scaffold-grouped 5-fold CV is the headline**, repeated over 3 fold
  partitions. Random-split CV is computed alongside only to quantify leakage, and
  the ratio is reported per target as `leakage_factor`.
* **A label-permutation null per target**, through the identical scaffold-split
  pipeline. A target is called modellable only if its real Spearman clears its
  own null 95th percentile, not zero — and only then does it enter the screen.
* **Model selection is by scaffold Spearman, and the winner is whoever wins.**
  Ridge competes with the trees on every target; where the linear model wins,
  it is used and reported.
* **RMSE is quoted next to that target's noise floor**, as a ratio. A model at
  RMSE ≈ noise floor is at the ceiling the data allows, and saying so prevents
  reading a large RMSE as a bad model when the assay itself is that noisy.
* **Censored values are dropped.** A ">" IC50 means inactive up to the top dose,
  not a weak number.
* **Assay-type mixing is tracked.** IC50/Ki/Kd are pooled into one pAffinity, but
  the binding-vs-functional group is retained and the per-group mean offset is
  reported per target.
* **Reference inhibitors are recall-tested, not trained against.** Each target's
  real assay controls (olaparib for PARP1, roflumilast for PDE4B, and so on) are
  looked up in the screening output and their percentile reported, with a flag
  for whether they were in training — an in-training control validates the fit,
  not generalisation, and the two are not conflated.
* **Frequent hitters are named.** A compound in the top decile of many unrelated
  targets is more likely exploiting a fingerprint shortcut than genuinely
  polypharmacological, and it is flagged instead of being left at the top of
  every list.
* **Selectivity is computed only over in-domain targets**, and the number of
  in-domain targets is reported next to it, so a two-target profile is not read
  as a panel-wide selectivity claim.

---

## What this deliberately does not do

* **No docking.** The 3CLpro work established, on its own data, that Vina scores
  could not rank potency (Spearman −0.057, potent-vs-weak AUC 0.508) while still
  enriching binders over non-binders (AUC 0.735). Its honest conclusion was that
  QSAR is the ranking engine and docking is orthogonal binding evidence. A panel
  of sixteen targets does not change that, and 16 × 9,746 docking runs would not
  buy a ranking. Structure-based evidence belongs in a follow-up on a shortlist.
* **No FRET-interference modelling.** Fluorescence quenching, aggregation and
  thiol reactivity are real sources of false positives in these assays and are
  not simulated. PAINS alerts are attached as an advisory flag, which is the in
  silico compensation, not an equivalent.
* **No cell-level extrapolation.** These are purified-enzyme assays; permeability
  and off-target effects are out of scope by construction, which is exactly why
  this assay class is the tractable one.
* **Predictions are rankings, not calibrated affinities.** Pooled IC50/Ki/Kd, and
  a random forest does not extrapolate past its training range.

Discovery-grade throughout: every hit needs the real assay to confirm it.

---

## Licence and attribution

Training and screening data are from ChEMBL, **CC BY-SA 3.0**. Derived datasets
here keep the ChEMBL molecule and assay IDs, record the release in
`results/panel_provenance.json`, and must be shared alike.

BPS Bioscience assay families are referenced to describe what each virtual assay
is modelled on. No BPS data, protocol detail or reagent is used or reproduced —
the assay readout column is a description of the experiment being stood in for,
and all training labels come from the public literature via ChEMBL.

---

## Layout

```
configs/panel.json         one entry per target: accession, family, readout, controls
pipeline/p01..p05          one stage per file
data/                      per-target training sets, raw ChEMBL cache, library
models/                    selected model + uncertainty forest + AD reference per target
results/                   provenance, benchmark, panel matrices, hits, profiles
```
