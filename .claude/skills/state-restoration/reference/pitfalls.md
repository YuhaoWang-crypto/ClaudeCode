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

## Lever prediction

| Symptom | Cause | Fix | Case |
|---|---|---|---|
| Top lever is the slowest module, not the drug target | J from one observational treatment has no measurement of other levers | Use a measured response library (perturbation screen, cytokine/drug panel) | D15 → D16/D17 |
| Lever library gives a confident but wrong answer | Library from the wrong cell type (fibroblasts for an epidermal disease): 0/4 correct | Match the cell type to the tissue compartment the disease signature comes from | D17 |
| A big public drug library (L1000/CMap) ranks none of the known therapies | The library measures drug effects on unstimulated cancer lines, not responses to the disease's own perturbation, and has no keratinocyte | Require both the right cell type and a disease-relevant stimulus; run a machinery check on a module whose genes the platform MEASURES before believing any ranking | D21 (AUC 0.51–0.63, best therapy #166/1768) |
| The top of a connectivity ranking is HSP90 / mTOR / proteasome inhibitors | Connectivity rewards signature magnitude, so the most violently acting compounds float up | Correlate score with signature norm (here Spearman +0.29, p 1e-35) and report it; check the ranking against a second disease | D21 |
| The module you care about is "in" the platform | Inferred genes are model output, not measurement: the barrier module has 7/12 genes in L1000 and 0 measured | Count MEASURED genes per module before asking the question, and say the question is unanswerable if the count is 0 | D21 |
| A within-cell covariance / co-expression statistic looks strong | Two features measured in the same cell share that cell's total count, which manufactures positive covariance (two isotype control antibodies covaried at +0.074 on raw ADT) | Normalise within the cell (CLR for ADT) and prove it worked on a pair that cannot co-vary biologically — an isotype pair, or two unrelated spike-ins | D24 |
| A second moment is reported per group like a mean | Covariance needs far more cells than a mean: split-half agreement was 0.78–0.94 for the mean and ~0 for the covariance in the same cells | Run split-half on the statistic itself, and compute its bootstrap SE against the between-group spread; report the cells-per-group the estimate would need (here ~3,000 vs the 100–500 available) | D24 |
| A deposited matrix's sample labels do not match the series' sample titles | Figure-level matrices are often labelled by construct or internal code (D11 / HA / SI) while GEO titles use the experimental arm (rested / always-off / always-on); the two cannot be aligned by position | Do not guess: a wrong mapping inverts the result. Resolve from the paper's figure legend, or report the arm as a design gap and use only the modality whose labels are self-describing | D25 |
| RNA is used as a proxy for chromatin state | The expression of the epigenetic machinery (EZH2, DNMT1/3A/3B, TET2, KDM6B) did not differ at all between exhausted and non-exhausted CAR-T: the difference is in chromatin, not in the machinery's transcripts | Test the proxy before relying on it — if the machinery module is undefined, RNA cannot stand in for the chromatin arm | D25 |
| Ranking flips between metrics | Cosine weights large-fold-change genes; Spearman weights all genes equally | Report both, and prefer a crossed specificity contrast over top-1 | D17 (psoriasis: IL-17A by cosine, IFN-γ by rank) |
| Top-1 lever wrong although the direction is right | Lesions carry several cytokine programmes (chronic AD has Th22/Th17 besides Th2) | Test the crossed contrast: is the lever more aligned with its own indication than the other? | D17 (IL-13 → AD, p < 0.001) |
| The drug's own target ranks #1 against the clinical response — and it means less than it looks | On-treatment change ≈ reverse of the disease (cos +0.975), so ranking levers against the drug effect is ranking them against the disease (rank correlation +0.95) | Project the disease axis out of both sides before claiming the ranking is about the drug; report the raw and the residual ranking | D20 (IL-17A +0.375 raw → +0.058 residual) |
| Lever ranking read as "this drug worked through this pathway" | Mechanistically different drugs produce near-identical transcriptome reversal in the same disease (secukinumab vs etanercept cos +0.671) | Only a crossed design (≥2 drugs × the same library) can support a drug-level claim; without it, the claim is indication-level | D20 (anti-TNF effect also ranks IL-17A first; DoD −0.120 [−0.296, +0.083]) |
| "Drug A's response differs from drug B's" with one cohort per drug | Drug, study, platform and batch are perfectly confounded: "different drug" also means "different study". Within-drug pairs are within-study by construction | Split the contrast by platform first. Here cross-platform pairs sat at 0.483 vs 0.812 within-drug, while the one same-platform pair was 0.727 — most of the "specificity" was technical | D23 (p 2e-18 naively, 0.09 of margin left after the platform split) |
| A mechanism-specific module prediction is assumed rather than tested | Pathway logic ("JAK blocks interferon so the interferon module should clear first") does not survive contact with clinical transcriptomes | State the prediction before looking and report the excess over the other modules; here it was +0.05 under JAK vs +0.01 under anti-TNF | D23 |

## Combination and library design

| Symptom | Cause | Fix | Case |
|---|---|---|---|
| The indication's own lever ranks low or flat | The library does not contain the lever the clinic uses (IL-4 absent, only IL-13) | List the validated targets first, then check the panel covers them | D17 → D18 |
| Non-additivity looks enormous (65 %) | The synergy vector is a difference of differences, so it carries the noise of three measurements | Check the between-donor correlation of the synergy vector against that of the response itself (0.10 vs 0.54 here); test synergy per gene by leave-one-donor-out | D18 |
| A validated drug target is not indication-specific | Shared inflammatory inputs (TNF) drive several diseases | Do not require every validated target to win the crossed contrast; report it | D18 (TNF −0.032) |

## Relapse / residual questions

| Symptom | Cause | Fix | Case |
|---|---|---|---|
| "Residual predicts relapse" cannot be tested cleanly | With one post-treatment sample, residual = x−ref and rise = flare−x share x with opposite signs, so their cosine is biased downward | Need ≥ 2 post-clearance samples per patient, or a cohort with non-relapsers | D19 |
| A random-gene null looks "significant" | Disease genes carry signal in both vectors, which alone beats random genes — it does not show prediction | Ask for a contrast that isolates the claim (own vs other patients), and check whether the bias explains the gap | D19 |
| Everyone in the cohort has the outcome | No contrast exists; a relapser-vs-non-relapser test is impossible | Say so and switch to within-patient questions (what restored, where relapse lands) | D19 |
| A control arm seems to prove the residual points at relapse | The shared-sample bias shrinks as the rise grows, and the relapsing arm has the larger rise by construction, so the arm gap comes for free | Compare |rise| between arms first; use pairings that do not share a sample (own vs other patients) | D20 (placebo −0.18 vs continued drug −0.67, but |rise| 25.7 vs 13.5) |
| A "cause removed" cohort shows a large residue | The exposed group is older / more exposed than the reference: former smokers are 31 years older than never smokers, corr(group, age) = 0.79 | Compute corr(group, confounder) FIRST. Above ~0.7 a regression adjustment is not identifiable (it silently returns "nothing survives"); find a matched window and report its n | D22 (25 % of genes persist raw, 0 % with age in the model, 5 vs 7 in the matched window) |
| The housekeeping module moves too | In bulk tissue whose cell composition changes (epithelium replaced by immune/crypt cells), relative abundance of every gene shifts, housekeeping included | Report the control module's own drift and count only the excess above it; if the margin is small, say the cohort cannot separate a molecular residue from a compositional shift | D22 (coeliac: control +0.36 vs crypt proliferation +0.51) |

## Engineering

- **GEO series matrices.** Some are linear MAS5 (GSE52428, GSE16879): use log2(max(x, 16)), then quantile-normalise. Others are already log scale. Check the range before choosing.
- **BH correction.** NaN p-values from constant floor genes make every q NaN. Set them to 1 first.
- **Challenge cohorts overlap.** GSE17156's H3N2 arm is the GSE52428 DEE2 cohort. Drop it before calling a test independent.
- **`pkill -f <pattern>` kills the calling shell** if the pattern appears in its own command line, including a background task you just started. Select processes by PID.
- **Verify every PMID before writing it into a report.** Two of the citations written before checking were wrong:
  - GSE75210 is O'Connell 2016, not Tong.
  - GSE16879 is Arijs, *PLoS One* 2009, not *Gut*.
