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
| Cell-type proportions from scRNA are read as compartment sizes | Dissociation recovery differs per type AND per condition: AT1 was 0.84 % of recovered cells in control mouse lung and 3.3 % after bleomycin, which destroys AT1; in human IPF it moved the other way (8.9 % to 5.0 %) | Check the fragile compartment's share across conditions before using any ratio; convert the question into "what relative recovery would erase this result" (here r = 0.39) and say so | D26 |
| A flux model is fitted to a snapshot | At steady state the compartment fractions equal the residence-time ratios, so a snapshot identifies 1/k1 : 1/k2 : 1/d and nothing more — no absolute rate, and during an acute injury the steady-state assumption does not even hold | Report the identifiable combination only; use the plateau (or a chronic cohort) for any steady-state statement, and get absolute days from pulse-chase or lineage tracing | D26 |
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

## Feasibility and identifiability audits

| Symptom | Cause | Fix | Case |
|---|---|---|---|
| A proportion time course is sold as pinning down the kinetics | Proportions are invariant to anything that scales both compartments: shifting both growth rates by the same amount leaves p(t) identical (2.2e-16) while total cell number changes 16-fold | Write the observable's own ODE before promising anything. If the model is linear, p(t) depends on a strict subset of the parameters, and the FIM rank names which | D27 |
| "Proportions cannot distinguish X from Y" is taken as obvious | It is often false. For AML, r = M/I obeys dr/dt = delta − Cr, so the differentiation rate IS identifiable; it is the kill rate that lumps (gI − dI) | Derive the ratio's ODE. Parameters that survive the derivation are identifiable; the rest enter only in lumped combinations | D27 |
| Absolute counts are assumed to fix everything | They raise the AML rank from 2/4 to 3/4 only — proliferation and death still enter as gI − dI, so a third channel (EdU, dye dilution) is required | Recompute the rank with each channel added before quoting a price for the experiment | D27 |
| A pseudo-inverse reports a finite error bar for a non-estimable parameter | `pinv(F)` silently splits a null direction between the parameters that span it (gS and gR each "sd 1.516" when neither is estimable) | Test each parameter's own unit vector against the null space and print inf when it leaks; or reparameterise into the identifiable basis | D27 |
| A good chi2 is read as a good fit | The model can fit with a completely flat direction: chi2 0.49 on 2 dof while dg's 95% profile covers the entire grid and b and p0 sit on bounds | Profile every parameter that enters a claim; report which ones sit on a bound as flat directions, not estimates | D27 |
| Structural identifiability is quoted as the answer | Structural identifiability says what no amount of data fixes; it does not say what the actual design estimates. a and dg are both structurally identifiable here, yet the practical block has cond 7.1e+05 and only one combination is estimable | Report both layers with numbers. First layer alone overstates capability; second alone makes "collect more cells" look like a solution | D27 |
| "More cells" or "more timepoints" is proposed as the fix | Neither touches a structural null direction, and for a sloppy direction the return is weak: 6 → 48 timepoints only moves CV(dg) from 1.81 to 0.63 | Compute the Cramer-Rao grid over cells x measurement precision x sampling density and read off which axis actually binds | D27 |
| Replicates are counted from a GEO SuperSeries' subseries list | Submitters re-deposit the same library under several subseries: 3 of 14 advertised PC9 libraries are byte-identical after decompression | Hash the decompressed files before claiming independent replicates | D27 |
| Sequencing depth is confounded with the experimental axis | Later timepoints were sequenced deeper (corr(day, log median UMI) = +0.865) and PC1 then *is* depth (corr +0.964) | Thin every cell to a common UMI budget and report corr(PC1, log depth) before and after (+0.964 → −0.013); only then compare libraries | D27 |
| Independent libraries at the same timepoint disagree wildly | Usually depth, not biology: day-11 p was 0.870 / 0.341 / 0.643 raw and 1.000 / 0.954 / 0.971 (sd 0.023) after thinning | Use same-state libraries as the batch floor, and fix depth before blaming batch | D27 |
| "The X-positive fraction" is treated as a well-defined quantity | Different definitions of the second state give qualitatively different time courses — a label-free mixture here tracks proliferation arrest (cycle shift −0.369) and is non-monotone, while a supervised axis rises monotonically | Define the second state explicitly and report the agreement between at least two independent definitions (+0.44 to +0.86 here) | D27 |
| A supervised axis is built from the endpoints and then fitted at the endpoints | The two defining libraries are not evidence; p is forced to ~0 and ~1 there | Refit using independent libraries at the same days as a sensitivity check | D27 |
| Recovered cell counts are used as a survival read-out | Barcode counts span 24-fold within one experiment with 52–95% empty droplets, and the four xenograft arms of the same study carry 29,340–29,983 cells (CV 1.0%) — set by loading | Never infer survival from the number of cells in a deposit; it needs beads, volumetric counting or imaging | D27 |
| A deposit is assumed to carry the counts because the paper reports them | Across 120 audited GEO series, 0 carried any absolute-count field in sample characteristics, and only 1/60 carried ≥3 structured time values | Audit the structured metadata programmatically before designing around it; also expect time to live in sample titles rather than characteristics | D27 |

## Direction of perturbation (knockdown vs overexpression)

| Symptom | Cause | Fix | Case |
|---|---|---|---|
| A lever ranking from a knockdown library is used to prescribe an overexpression | The step is the linear-response identity d(-u) = -d(u), which is assumed, not measured. In K562 the matched cosine is **-0.083** against a prediction of -1, and 70/84 genes are directionally independent | Measure the antisymmetry on genes perturbed in both directions before transferring a ranking; report the corrected cosine, not the raw one | D28 |
| A loss of antisymmetry is explained by "saturation" or "the dose was different" | An odd input changes magnitude only: a 3x stronger perturbation still gives cos = -1.000000 exactly. Only a direction-specific channel can rotate the response | Rule out dose first with the planted check; then look for a program that only one direction engages | D28 |
| matched-minus-mismatched is used as the control for a global mode | The contrast collapses along with the raw value: at a global mode 4x the signal it goes from -1.00 to -0.06 | Subtract the across-perturbation mean field explicitly AND report the mismatched null; neither alone is sufficient | D28 (this refuted the planned design) |
| A cosine and a regression slope are corrected for noise the same way | A cosine is attenuated by the noise on BOTH sides, a slope only by the noise on the regressor. Here the same data gives a corrected correlation of -0.17 to -0.39 but a corrected slope of -0.87 | Correct them separately; mixing the two reads "right sign, wrong content" as either "fully right" or "fully wrong" | D28 |
| A near-zero median correlation is reported as "no relationship" | Noise fills both tails equally, so the tail COUNTS carry information the median hides: 14 significantly antisymmetric vs 0 co-directional, binomial p = 6e-5 | Classify per unit against the mismatched null and test the tail asymmetry | D28 |
| A perturbation library is assumed to have signal because the screen is large | The responses can sit on the library's own noise floor for the genes you care about: 7/84 above the 95th percentile of 585 non-targeting rows, two-TSS replicate cos +0.050 | Use non-targeting rows as the noise floor and report how many perturbations clear it before interpreting any correlation | D28 |
| Two directions are compared across two studies | Cross-study, cross-lab and cross-technology differences are perfectly confounded with the biological contrast, exactly as drug and study were in D23 | Say so; the design that settles it is both directions on the same genes in one experiment | D28 |

## Engineering

- **GEO series matrices.** Some are linear MAS5 (GSE52428, GSE16879): use log2(max(x, 16)), then quantile-normalise. Others are already log scale. Check the range before choosing.
- **BH correction.** NaN p-values from constant floor genes make every q NaN. Set them to 1 first.
- **Challenge cohorts overlap.** GSE17156's H3N2 arm is the GSE52428 DEE2 cohort. Drop it before calling a test independent.
- **`pkill -f <pattern>` kills the calling shell** if the pattern appears in its own command line, including a background task you just started. Select processes by PID.
- **Verify every PMID before writing it into a report.** Two of the citations written before checking were wrong:
  - GSE75210 is O'Connell 2016, not Tong.
  - GSE16879 is Arijs, *PLoS One* 2009, not *Gut*.
