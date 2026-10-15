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
| **GSE102725** ★★★ | psoriasis **secukinumab (anti-IL-17A)**, lesional + non-affected | closes the loop: the lever D17/D18 predicted is this drug's target |
| **GSE69967** ★★ | psoriasis tofacitinib (JAK), randomised phase 2, per-patient time course | different mechanism, same disease |
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

## Five studies worth doing first

1. **Do slow modules predict relapse?** (GSE30768 + GSE85034) Prediction from D14:
   patients whose barrier/lipid module is least restored at week 12 relapse sooner.
   Small n, directional only — but nobody has tested it, and it would give a
   molecular criterion for when therapy can be stopped.
2. **Close the loop on IL-17A** (GSE102725). The method predicted IL-17A from a
   cytokine library; secukinumab blocks exactly that. Check the state moves along
   the predicted direction with the predicted gain.
3. **Do three mechanisms share one path?** (GSE11903 anti-TNF, GSE102725 anti-IL-17,
   GSE69967 JAK). Prediction: same disease axis, same fast modules, different
   extent of slow (barrier) restoration.
4. **Extend levers to real drugs** (GSE92742 LINCS). Project the disease
   displacement onto thousands of drug responses; known-effective drugs ranking
   high is the positive control.
5. **Does removing the cause leave a scar?** (GSE254592 HCV cure + GSE7895 quitting).
   Prediction: inflammatory modules fully recover, structural/epigenetic ones do not.

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
