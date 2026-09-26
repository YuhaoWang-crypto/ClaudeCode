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

## Results

### All sixteen virtual assays are modellable

Scaffold-split Spearman **+0.649 to +0.881**, every target clearing its own
permutation null by a wide margin (nulls land at +0.004 to +0.102). Model
selection went to LightGBM on 13 targets, random forest on 2, and **ridge on
PTPN11** — the linear model won there and is what PTPN11 uses.

**Leakage factors are 1.02–1.12**, far milder than the 2.7× seen on a
CDK2 dataset of 784 compounds. That is not a property of the method: these
training sets carry 616–9,107 compounds over 301–3,610 scaffolds, so holding out
a scaffold removes proportionally much less information. Leakage severity is a
function of scaffold density, which is why it has to be measured per dataset
rather than assumed from another one.

RMSE sits below the target's noise floor on 4 of 16 targets. **That is not
"better than the experiment."** Cross-validation predicts a median over pooled
sources, and folds share source-specific consistency, so the CV task is easier
than reproducing an independent assay. NAMPT's 0.17-log floor makes its ratio
the least meaningful of the set.

### Reference inhibitors come back — but that proves less than it looks

44 of the panel's reference inhibitors are present in the screening library, and
**43 reproduce their own measured value to within 1.5× the target's noise floor**
(median |residual| 0.01–0.80 log). The one miss is navoximod on IDO1
(predicted 6.54, measured 7.55).

The number that matters is different: **only 2 of those 44 are genuinely
out-of-sample.** The rest are in their target's training set, so reproducing
them is a check that labels are plumbed through correctly, not evidence of
generalisation. The two that are out-of-sample both land correctly —
resveratrol at the 5th percentile on SIRT1 (it is an activator, not a potent
inhibitor) and migoprotafib at the 100th on PTPN11 — but n = 2 is anecdote. The
generalisation evidence is the scaffold-split CV, and nothing else here.

Getting to that number required fixing three errors in the control list itself,
all of which would have read as model failures:

* **Resveratrol was listed as a SIRT1 reference inhibitor.** It is an activator.
  Corrected to `expect: weak` with the reason recorded rather than deleted.
* **Daretuzumab was listed as a NAMPT control.** It is an antibody and cannot
  appear in a small-molecule library — a modality error, not a coverage gap.
* **FK866 returned "not in library".** ChEMBL lists it under its INN,
  **daporinad**. The panel's own output surfaced this: daporinad came back with
  NAMPT as its top predicted target.

A first version scored controls by percentile band and produced boundary
artifacts — cilomilast failed at exactly the 75th percentile, and rolipram's
correct ~1 µM prediction failed a "weak" band only because a clinical library
skews weaker than that. Scoring against each control's own measured value
removed both.

### Selectivity: the panel can call it for one paralog pair and not the other

This is the test that a single-target benchmark cannot run, and the one that
decides whether the selectivity column is worth shipping.

| pair | shared compounds | ρ(Δ measured, Δ predicted) | sign agreement | verdict |
|---|---|---|---|---|
| HDAC1 / HDAC6 | 3,847 | **+0.772** (+0.710 above noise) | 86% (95% above noise) | predictable |
| PTP1B / PTPN11 | 150 | **+0.174** (+0.405 above noise) | 63% (68%) | **not reliable** |

Both PTP models are individually strong — scaffold Spearman +0.716 and +0.855 —
and their *difference* is close to a coin flip. Two good models do not make a
good selectivity call. HDAC1/HDAC6 has 25× more shared measurements and
selectivity is an explicit design objective across that series; PTPN11
inhibitors are allosteric while PTP1B inhibitors target the active site, so
there is little shared chemistry from which the gap could be learned.

**A blanket "the panel does selectivity" claim would have been false for half
the pairs tested.** The selectivity column in `compound_profiles.csv` should be
read per pair, against this table.

### The screen

9,746 clinical and marketed compounds × 16 targets = **155,936 predictions**.
2,624 compounds are inside the applicability domain of at least one target;
**1,143 in-domain hits** at predicted pAffinity ≥ 6.0 and not already in
training (101 carry a PAINS alert). 1,665 compounds (17%) are flagged frequent
hitters and excluded from the selectivity view.

The panel recovers known target assignments it was never told:

| compound | panel's top target | reality |
|---|---|---|
| tadalafil | PDE5A (gap 1.97 log) | marketed PDE5 inhibitor |
| daporinad (FK866) | NAMPT | the canonical NAMPT tool inhibitor |
| AZD2461 | PARP1 | PARP inhibitor |
| INCB-057643 | BRD4 | BET/BRD4 inhibitor |
| tulmimetostat | EZH2 | EZH2 inhibitor |
| KA-2507 | HDAC6 | selective HDAC6 inhibitor |
| OBP-801 | HDAC1 | HDAC inhibitor |

And it produces implausible ones next to them: **padimate A**, a sunscreen UV
filter, ranks as an HDAC6 hit, and **miramistin**, an antiseptic surfactant,
as a NAMPT hit. Both pass the applicability-domain filter. This is the
characteristic ligand-based failure — small or greasy molecules whose
fingerprint similarity to a training series is superficial — and it is why the
hit lists are labelled discovery-grade and why the frequent-hitter flag exists.
A shortlist from this panel needs a chemist's eye before it needs a plate.

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
