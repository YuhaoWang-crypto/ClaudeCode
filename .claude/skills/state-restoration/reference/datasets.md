# Public datasets for state-restoration studies

Every accession below was fetched from NCBI and its metadata inspected (Oct 2025).
Grouped by the question each class of data can answer. Stars = value for this
framework, not scientific merit of the original study.

## What each class can and cannot answer

| Class | Answers | Does NOT answer |
|---|---|---|
| A withdrawal / stimulus removal | did the state return to a SELF-SUSTAINING balance; who relapses | which target to hit |
| B treatment time course + outcome | is molecular recovery consistent with functional recovery; what restores last | causal upstream target |
| C perturbation library | **which target to hit**, how few knobs | long-term maintenance |
| D reference cohort | where balance is, how large the noise floor is | dynamics |

A complete answer needs C (target) + B (consistency) + A (maintenance).

## A. Withdrawal / removal — scarcest and most valuable

| Accession | Design | Tests |
|---|---|---|
| **GSE145408** ★★★ | kidney transplant, tacrolimus withdrawal vs none; pre and 3 mo post, PBMC (n=36) | does the state hold without the drug; outcome is graft stability off-drug |
| **GSE30768** ★★★ | psoriasis efalizumab: NL / LS / week 12 / **FLARE** per patient (n=14) | does relapse return to the original disease state; does an unrestored module predict flare |
| **GSE7895** ★★★ | never / current / **former** smokers, bronchial epithelium, with pack-years | which modules recover after removing a chronic stimulus and which are permanent |
| **GSE290636** ★★ | HIV, blood before analytical treatment interruption (n=74) | predicting who rebounds off therapy |

Searching repeatedly found very few human "stop therapy → serial biopsy → record
relapse" transcriptomes. Treat this as a finding: the maintenance question
mostly cannot be answered from existing public data.

## B. Treatment time course with clinical outcome

| Accession | Design | Notes |
|---|---|---|
| **GSE85034** ★★★ | psoriasis adalimumab/methotrexate, weeks 0/1/2/4/16 + own NL + PASI | used in D14/D15; week-1 molecular recovery predicts PASI75 at AUC 0.78 (n.s., n=25) |
| **GSE102725** ★★★ | psoriasis **secukinumab (anti-IL-17A)**, weeks 0/16/32, **randomised withdrawal at PASI 98** (5 continue / 5 placebo), lesional + non-affected, PMID 39565228 | used in D20: closed-loop lever test, plus the only randomised withdrawal arm in this catalogue |
| **GSE69967** ★★ | psoriasis tofacitinib (JAK), randomised phase 2, per-patient time course, placebo arm n=3 | used in D23: third mechanism; the placebo arm is too small to be a strong control |
| **GSE11903** ★★ | psoriasis etanercept (anti-TNF), weeks 0/1/2/4/12 | with the two above: three mechanisms, one disease |
| **GSE130588** / **GSE59294** ★★★ | atopic dermatitis dupilumab, **placebo-controlled** | separates drug effect from natural course |
| **GSE16879** / **GSE73661** ★★★ | IBD infliximab / vedolizumab, before-after and W0/W4-6/W12, endoscopic response | two mechanisms + hard outcome |
| **GSE91061** ★★ | melanoma nivolumab, pre/on treatment, RECIST response | "balance" must be redefined for tumours ⚠️ |
| **GSE254592** ★★ | hepatitis C before/after direct-acting antivirals (cure), PBMC | the cause is fully removed — does the state fully return |
| **GSE322560** ★★ | COVID convalescence, longitudinal, persistent interferon signature | how long slow modules take to fall back |

## C. Perturbation libraries — the only data that sets targets

| Source | Content | Notes |
|---|---|---|
| **GSE36287 + GSE305982** ★★★ | keratinocyte cytokine library, 9 levers incl. **IL-4** | the skill default (`netjac.data.keratinocyte`) |
| **GSE24767** ★★★ | keratinocyte IL-17 / TNF / **IL-17+TNF**, same donors | singles and the combination both measured |
| **GSE92742** ★★★ | LINCS L1000, ~1.3M profiles, compounds + knockdowns, many cell types incl. NHEK | extends levers from cytokines to real drugs; 978 landmark genes only, mostly cancer lines, strong batch effects |
| Replogle genome-wide / Norman doubles ★★★ | 11,258 knockdowns; 131 double perturbations | ground truth for the lever test (66 % vs 7 % chance); K562 only |
| scPerturb (Zenodo) ★★ | harmonised perturbation datasets across cell types | first place to look for a library in YOUR tissue |
| **GSE282371** ★★ | 3D epidermal equivalents, Th2 + IL-17 / IL-22 combinations | organoids are the realistic route to the right cell type |

## D. Reference cohorts

| Accession | Use |
|---|---|
| **GSE52428** / **GSE17156** ★★★ | two pre-challenge samples per subject: mean = personal balance, difference = personal noise floor |
| **GSE75214** | IBD mucosa reference ensemble (mean and covariance of controls) |
| **GSE75210** ★★ | Myd88⁻/⁻ macrophages, 12 stimuli × 8 times — a known-answer test set for the method |
| **GSE135893 / GSE141259** ★★ | alveolar repair: human IPF (10 control / 20 disease donors) and a mouse bleomycin time course (days 2–54, epithelium-enriched), both with cell-level annotations only — enough for compartment ratios, but the AT1 compartment is lost to dissociation and its recovery moves with condition (D26) |
| **GSE164949 / GSE164947** (Weber 2021) ★★★ | CAR-T exhaustion with transient rest: RNA-seq over days 7/11/15 with always-ON, always-OFF, rest from day 7 or 11, and nivolumab arms, 3 donors — a genuine dose-of-rest design. Used in D25. Caveat: the ATAC deposit has no rested arm and the ChIP matrix labels cannot be mapped, so the chromatin side of any two-timescale question is unanswerable here |
| **FrangiehIzar2021** (scPerturb) ★★★ | melanoma Perturb-CITE-seq: 218k cells x 24 surface proteins (incl. four isotype controls) x CRISPR KO x {control, IFN-γ, TIL co-culture}. Used in D24. The isotype controls make it the reference dataset for testing any within-cell covariance claim; note ~100–500 cells per KO is too few for second moments |
| **GSE7895** ★★★ | bronchial epithelium, never / current / former smokers with 1–504 months since quitting — a dose-of-time cause-removal design (used in D22; age-confounded, see pitfalls) |
| **GSE145358** ★★ | coeliac duodenum, 15 patients on a gluten-free diet then gluten-challenged, paired, plus healthy controls — cause removed AND restored in the same people (used in D22) |
| **GSE70138** (LINCS L1000 Phase 2) | 118,050 signatures, 2,171 compounds, 7 main cell lines; 978 MEASURED genes, 11,350 inferred. Used in D21: fails the psoriasis positive control; useful as a cautionary benchmark, not as a lever library for a tissue disease |

## Five studies worth doing first

1. ~~**Do slow modules predict relapse?**~~ **DONE — D19 (and D20).** Premise
   replicated (barrier unrestored at clinical clearance in a 4th and 5th drug) and
   relapse returns to the lesion position, but the predictive claim is not
   answerable in either cohort: one post-clearance biopsy makes residual and rise
   share a sample. Still open, and still needs ≥ 2 post-clearance biopsies with
   relapsers and non-relapsers in one arm.
2. ~~**Close the loop on IL-17A** (GSE102725).~~ **DONE — D20.** IL-17A ranks #1
   for secukinumab, but the ranking is 95 % explained by the disease axis and the
   anti-TNF drug's effect ranks IL-17A first too: lever ranking is indication-level,
   not drug-level. The same trial's randomised withdrawal arm reproduced the
   attractor result (+0.86 vs −0.01) and again gave no support for "residual
   predicts relapse".
3. ~~**Do three mechanisms share one path?**~~ **DONE — D23.** Same disease axis
   (cross-cohort cosines equal the within-cohort patient spread) and the same
   ordering (barrier last under all three). The 0.33 of "difference" between drug
   effects turned out to be mostly platform: cross-platform pairs 0.483 vs 0.812
   within-drug, and only 0.09 of margin survives in the one same-platform pair.
   Drug and study are perfectly confounded with one cohort per drug, so this
   question needs head-to-head sampling in one cohort.
4. ~~**Extend levers to real drugs** (LINCS).~~ **DONE — D21, and it FAILED the
   positive control.** On GSE70138 (Phase 2, 2,171 compounds) the 23 psoriasis
   therapies do not rank above 28 unrelated drugs in any of four score/gene-space
   variants, the top of the list is HSP90/mTOR inhibitors, the ranking is
   cell-line-dependent and not disease-specific, and the barrier module has zero
   measured genes. Do not spend time on GSE92742 (Phase 1, 21 GB) expecting a
   different answer: the limit is the assay (drugs on unstimulated cancer lines),
   not the compound count. The experiment that would work is drug x cytokine-stimulus
   screening in keratinocytes.
5. ~~**Does removing the cause leave a scar?**~~ **DONE — D22**, with GSE7895
   (smoking cessation) and GSE145358 (coeliac on a gluten-free diet) instead of
   GSE254592 (that one is PBMC, not tissue). Exposure-driven modules return,
   secretory/structural ones do not, but only at q ≈ 0.09 in the age-matched
   window — a cause-removal cohort is confounded by age (corr 0.79 here), and in
   coeliac bulk tissue the housekeeping control drifts as much as the modules.
   Still open, and the design that settles it is a cohort matched on age and
   exposure, or a longitudinal one sampled before and years after removal.

## Limits to state every time

- Without a perturbation library, report **bottlenecks**, never upstream targets
  (covariance-based ranking: 15 % vs 66 %).
- A library from the wrong cell type answers confidently and wrongly (fibroblasts 0/4).
- A library missing a validated target cannot answer that indication (IL-4).
- Linear models cannot see epistasis; whole-transcriptome non-additivity is mostly
  noise (synergy donor r 0.10 vs 0.54).
- Molecular recovery tracks clinical improvement at ρ 0.58–0.66 — consistent, not
  interchangeable.
- Rankings can flip between cosine and rank-based metrics; report both.
- If the sampling window is shorter than the recovery time, the long-term steady
  state is not identifiable.
