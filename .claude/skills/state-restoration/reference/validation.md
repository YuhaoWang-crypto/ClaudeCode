# Pre-specified checks and their outcomes

Every number is reproducible with the named demo. Status:
- ✅ pass;
- ⚠️ partial, or the method's reading changed;
- ❌ fail (kept, because the failure is informative).

## Step 1 — steady state

| Check | Case | Result | Status |
|---|---|---|---|
| Disease state far from balance | D15 PSO / AD / IBD | Displacement 3.8 / 2.7 / 7.5 × noise floor (between-subject upper bound) | ✅ |
| Disease state persistent without intervention | D15 AD placebo arm | Restoration −0.05 at week 4; +0.30 at week 16 (n = 7, natural course plus selection) | ✅ early / ⚠️ late |
| Personal reference noise floor | D13 two pre-challenge samples | ISG module SD 0.23 log2 | ✅ |
| Exclusion rule needed | D14 psoriasis | Without it one patient scored R = 40 | ✅ (rule kept) |

## Step 2 — order parameters

| Check | Case | Result | Status |
|---|---|---|---|
| φ separates the known positive after perturbation, not before | D15 AD (EASI-75 wk 16) | φ AUC t=0: 0.50 → wk 4: 0.67 → wk 16: 0.85 | ✅ |
| | D15 PSO (PASI75) | 0.50 → wk 1: 0.68 → wk 4: 0.79 → wk 16: 0.73 | ✅ |
| | D15 IBD (endoscopic response) | 0.50 → wk 5: 0.86 | ✅ |
| Label-free PC aligns with φ | D15 PSO / IBD / AD | \|r(PC1, φ)\| 0.83 / 0.78 / 0.13 (AD: PC2 0.51) | ✅ / ✅ / ⚠️ heterogeneity dominates AD |
| Collective mode ⇒ one prior gene suffices | D13 influenza | PC1 (no labels) AUC 0.90; RSAD2 alone 0.85–0.90 on an unseen strain and on unseen viruses; data-driven top genes within the random-set null | ✅ |
| Functional state needs perturbation | D12, D13 | Static 0.50 → 1.00; 0.59–0.72 → 0.93 | ✅ |
| Compositional state static | D12 genotype | Static 1.00; dynamics adds 0 (transcriptome: −0.11) | ✅ |

## Step 3 — core network

| Check | Case | Result | Status |
|---|---|---|---|
| Core replicates across groups vs permutation null | D13 H3N2 / H1N1 | 111 up / 41 down vs null median 0 | ✅ |
| Core = textbook pathway | D13 | Reactome IFN α/β signalling p = 3e-15 | ✅ |
| | D12 | KEGG NF-κB p = 0.02; IRF arm kept under LPS | ✅ |
| Core replicated across drugs | D15 PSO (ADA vs MTX) | 190 of 206 disease genes | ✅ |
| Couplings biologically sensible | D15 PSO | IL-17 keratinocyte – proliferation partial r 0.48; interferon – T cell 0.42 | ✅ (descriptive) |

## Step 4 — dynamics

| Check | Case | Result | Status |
|---|---|---|---|
| Beats "no change" (5-fold subject CV RMSE) | D15 AD / PSO / IBD | 0.90 vs 1.09; 0.40 vs 0.58; 0.36 vs 0.49 | ✅ |
| Matches the same-arm mean trajectory | | 0.90 vs 0.88; 0.40 vs 0.40; 0.36 vs 0.36 | ✅ (individual info only via y0) |
| Slow mode = the module observed to restore last | D15 PSO vs D14 | Slowest τ 23 wk, barrier / lipid 0.44 of the mode. D14 observed: barrier restored 0.65 at wk 16 (lowest) | ✅ |
| | D15 IBD | Epithelial metabolism is the slow mode | ✅ |
| b points to the drug mechanism | D15 AD dupilumab | b: T cell / DC −0.55, Th2 −0.31, IL-17 keratinocyte −0.30; interferon +0.17 (not targeted) | ⚠️ Th2 second, specificity ✅ |
| Steady-state identifiability guard fires | D15 IBD (1 follow-up) | τ 129 wk vs window 5 wk → y* not reported, y(5 wk) instead | ✅ |
| | D15 AD | τ 463 wk complex pair vs 16 wk | ✅ |
| Applied intervention restoration | D15 PSO | Model 0.74 vs observed medians 0.96 (responders) and 0.59 (non-responders) in D14 | ✅ |
| | D15 AD auto modules | Model 0.70 vs observed 1.06 at wk 16 | ⚠️ under-predicts |

## Step 5 — minimal adjustment

| Check | Case | Result | Status |
|---|---|---|---|
| Ground truth: the fix is found | D11 Teusink | Greedy set = feedback strength + PFK (−99.4 %); Vmax-only −2.6 % | ✅ |
| Rescue pairs predicted from single perturbations | D11 Norman 2019 | AUC 0.89 vs static co-expression 0.55 | ✅ |
| Clinically validated upstream target is the top lever, from observational treatment data | D15 PSO (no IL-17 drug in the data; IL-17 blockade is the most effective therapy) | Top lever barrier / lipid (0.38); IL-17 keratinocyte 0.17 | ❌ levers = bottlenecks |
| | D15 AD | Top lever IL-23 / IL-17 cytokines; Th2 fourth; steady state not identifiable | ❌ / unreliable |
| Structural control on large networks | D11 Perturb-seq | 51 % driver nodes (null 48 %) → full control impractical; target functional restoration | ✅ (informative) |

### Step 5 re-tested with a causal J (D16)

| Check | Result | Status |
|---|---|---|
| Response matrix reproducible | split-half r 0.91 ± 0.04 (15 modules, Replogle genome-wide CRISPRi) | ✅ |
| Identity J⁻¹e_j = −R[:,j] | verified numerically in `statepipe.causal_J` | ✅ |
| **Lever names the module actually broken** (439 held-out knockdowns, leave-one-out, perturbed gene excluded) | causal 66 % [61–70 %]; naive 25 %; covariance 15 %; chance 7 %; permutation null 10 % (p < 0.0005) | ✅ **the D15 failure is fixed** |
| Advance diagnostic of which levers are nameable | none found: \|R_j\| Spearman +0.43 (p = 0.11), max\|cos\| −0.11 (p = 0.69), n = 15 | ❌ hypothesis not supported |
| Errors are structured | named lever resembles the true one more than random (median \|cos\| 0.55 vs 0.40, p = 0.01); only 32 % strongly collinear; confusions are causally nested (mitoribosome → complex I) | ⚠️ |
| Doubles: beat additive superposition | AUC 0.89 both ways; in linear response d_AB = d_A + d_B exactly, so epistasis is outside the model class | ⚠️ by construction |
| Doubles: pick the best partner | 54 % correct (37 perturbations, chance 21 %) | ✅ |

### Step 5 in the relevant cell type (D17, the psoriasis question itself)

Library: GSE305982, primary keratinocytes and fibroblasts, six T-cell cytokines
(IL-9, IL-13, IL-17A, IL-22, IL-26, IFN-γ) over ImmunoCult, 3 h and 24 h, one donor.
Ground truth is crossed and clinical: psoriasis → IL-17A, atopic dermatitis → IL-13.

| Check | Result | Status |
|---|---|---|
| **Psoriasis top lever (keratinocytes, 24 h)** | **IL-17A**, cos +0.231 vs IFN-γ +0.132; coefficient positive, so restoration = blockade; IL-13 negative | ✅ **the D15 answer is corrected** |
| Wrong cell type must fail | fibroblast library: 0/4 settings correct (always IFN-γ) | ✅ |
| Crossed specificity IL-13 → AD | +0.029 [0.009, 0.049] at 3 h; +0.181 [0.154, 0.207] at 24 h | ✅ |
| Crossed specificity IL-17A → psoriasis | +0.008 [−0.014, 0.031] at 3 h; +0.026 [−0.009, 0.061] at 24 h (p = 0.071) | ⚠️ right direction, n.s. |
| AD top lever | IL-22 at 24 h (IL-13 fourth); IL-13 only at 3 h | ❌ |
| Psoriasis robust to metric | survives dropping the 200 strongest disease genes, but rank-based (Spearman) ranking puts IFN-γ first | ⚠️ metric-dependent |

Limits: one donor; no IL-4 in the panel (dupilumab blocks IL-4Rα, covering IL-4
and IL-13); whole-biopsy disease vectors against a pure keratinocyte library;
RNA-seq library vs array disease data.

### IL-4 and measured combinations (D18)

| Check | Result | Status |
|---|---|---|
| Psoriasis → IL-17A replicates on independent data (GSE36287, 3 donors, GPL570) | cos +0.220 vs next +0.104; IL-4 and IL-13 negative | ✅ |
| Its crossed specificity with a better library | +0.174 [0.147, 0.201] (was +0.026, p = 0.071 in D17) | ✅ |
| AD-specific lever | **IL-4** +0.097 [0.070, 0.125]; IL-13 +0.004 n.s. — a panel without IL-4 cannot answer AD | ✅ (explains the D17 miss) |
| AD top-1 | still IFN-γ (IL-4 third) | ❌ |
| TNF specific to psoriasis | −0.032 [−0.055, −0.007]: TNF is a shared inflammatory input despite being a validated psoriasis target | ❌ |
| Measured IL-17+TNF beats linear superposition | +0.053 [0.031, 0.074] | ✅ |
| Whole-transcriptome non-additivity is real | looked like 65 %, but synergy-vector between-donor r = 0.10 (combination response r = 0.54): **mostly noise** | ❌ claim withdrawn |
| Synergy real at specific genes | leave-one-donor-out: +1.1 to +1.6 log2 at the top-20 synergy genes found in other donors vs ~0 random, p < 0.001 in all three; genes are IL19, S100A7A, DEFB4A, SERPINB4, IL23A, S100A7 | ✅ |
| Two levers beat one for this lesion | IL-17 alone +0.284 ≈ measured combination +0.282 | ❌ no gain |

### Does an unrestored slow module mark relapse? (D19)

| Check | Result | Status |
|---|---|---|
| Restoration ordering replicates in a 4th drug (efalizumab) at clinical clearance | proliferation +1.05, T cell +0.89, interferon +0.59, IL-17 keratinocyte +0.33, **barrier +0.10** | ✅ |
| Housekeeping module undefined (negative control) | fewer than 3 disrupted genes | ✅ |
| Relapse returns to the original lesion state | median φ +0.95 (range 0.71–1.55) from φ 0.18–0.65 at week 12 | ✅ new |
| **Residual at clearance predicts the flare** | **confounded**: one post-treatment sample enters the residual (+) and the rise (−), so the cosine is biased down; raw median −0.06, own-vs-other gap ≈ the bias; n=4 ordering runs the wrong way | ❌ not answerable |
| Long-resolved skin still inflamed? | 13 paired patients: all inflammatory modules n.s.; no over-response to OKT3 (resolved +1.08 vs non-lesional +1.35) | ✅ informative negative |

Design needed to answer it: ≥ 2 post-clearance biopsies per patient (to separate
noise from signal), or a withdrawal cohort where some patients do NOT relapse.

### Closing the loop with the drug that blocks the predicted lever (D20)

GSE102725 (PSORITUS secukinumab withdrawal trial, 10 patients, weeks 0/16/32,
randomised continue vs placebo at PASI 98) + GSE11903 (etanercept, 15 patients).

| Check | Result | Status |
|---|---|---|
| The drug's own target ranks #1 of 9 levers | IL-17A cos +0.375 [0.247, 0.483]; gene-shuffle p < 1e-4 | ✅ |
| ... and that ranking carries drug information | ranking vs "looks like the disease" Spearman **+0.95**; after removing the disease axis IL-17A +0.058 and IL-4 tops | ❌ |
| Multi-lever magnitude fit is usable | within one source R² 0.268 with IL-13 coefficient −1.25 (absurd sign); IL-17A alone R² 0.140 | ❌ collinear |
| **Crossed control: library distinguishes the two drugs** | etanercept also ranks IL-17A first (+0.54 > its own target TNF +0.40); DoD −0.120 [−0.296, +0.083], p 0.88; the two drug effects are cos +0.671 | ❌ decisive |
| Barrier lags at the strictest endpoint (PASI 98), 5th drug | +0.38 vs IL-17 keratinocyte +1.08, proliferation +0.98, interferon +0.77 — but only 4/12 barrier genes disrupted, sd 0.90 | ✅ weakly |
| Randomised withdrawal returns to the lesion state | week 32: placebo +0.86 vs continued drug −0.01 | ✅ (D19's attractor result, now randomised) |
| Residual at clearance predicts relapse | own-vs-other −0.33 vs −0.11; residual size vs relapse size ρ +0.10, p 0.87 | ❌ second negative |

**Rule it sets.** Lever ranking from a response library is an *indication-level*
statement ("this disease is an IL-17-type state"). A *drug-level* claim ("this
response was produced through this pathway") needs at least two mechanistically
different drugs sampled against the same library, because in one disease their
transcriptome effects are nearly the same vector.

### Can a generic drug library replace the response library? (D21)

GSE70138 (LINCS L1000 Phase 2): 118,050 level-5 signatures, 2,171 compounds,
7 main cell lines, 978 measured genes + 11,350 inferred.

| Check | Result | Status |
|---|---|---|
| Machinery works where the platform measures the genes | 16 antiproliferatives lower the proliferation module (4 landmark genes) vs unrelated drugs, AUC 0.72 | ✅ |
| **23 psoriasis therapies rank above 28 unrelated drugs** | AUC 0.634 (p 0.055) measured-genes/cosine; 0.514 measured/rank-based; 0.520 disease-genes/cosine; 0.556 disease-genes/rank-based — the two metrics agree at Spearman +0.85 | ❌ all four |
| Best psoriasis therapy's rank of 1768 compounds | #166 (dexamethasone-acetate) / #50 (methylprednisolone) | ❌ |
| Known psoriasis aggravators score below the therapies | aggravators +0.039 vs therapies +0.004 in the disease-gene space | ❌ reversed |
| Ranking survives the cell line | per-line AUC 0.31–0.79; between-line ranking Spearman +0.15 to +0.37 | ❌ |
| Score is not a magnitude artefact | Spearman with signature norm +0.29 (p 1e-35); top of the list is HSP90 / mTOR / proteasome | ❌ |
| Ranking is disease-specific | psoriasis vs atopic dermatitis Spearman +0.41 to +0.52, 9–10/20 top-20 overlap, therapies' psoriasis-minus-AD +0.027 = the library's +0.029 | ❌ |
| The barrier module can be asked at all | 0 measured genes (7/12 inferred); "best" barrier-raising compounds are CDK/proteasome inhibitors | ❌ not answerable |
| Drugs are simply inactive in these lines | median |z| by class 30–40 vs library 31 — they do perturb the cells, the direction is wrong | ⚠️ explanation withdrawn |

**Rule it sets.** A step-5 response library needs the right cell type AND
responses to a disease-relevant perturbation. Drug effects on unstimulated
cancer lines satisfy neither, however many compounds they cover.

### Remove the cause instead of the drug: what stays? (D22)

GSE7895 (bronchial epithelium: 21 never, 52 current, 31 former smokers, 1–504
months since quitting) and GSE145358 (coeliac duodenum: 15 patients on a
gluten-free diet then gluten-challenged, 6 healthy controls).

| Check | Result | Status |
|---|---|---|
| Exposure-driven modules return once the cause is gone | xenobiotic metabolism +1.07 → +0.10 (slope per doubling of time −0.076, p 0.008); NRF2 +1.03 → +0.32 (p 0.004); cilia −0.09, n.s. vs never | ✅ |
| Structural / differentiation modules lag | secretory/mucous still +0.56 after 10+ years (p 0.0005 vs never, raw) | ✅ raw |
| ... and that survives the age confound | age-matched 40–60 window (n 5 vs 7): secretory/mucous +0.93 vs +0.20 and squamous/basal +0.62 vs −0.39, both p 0.030, **BH q 0.091**; xenobiotic p 0.20, cilia p 1.0 | ⚠️ suggestive |
| Gene-level "scar" fraction | 80/323 (25 %) raw; **0 with age in the model** — corr(group, age) = 0.79, so this is non-identifiability, not absence | ❌ not estimable |
| Coeliac: residue on a gluten-free diet | crypt proliferation +0.51, epithelial absorption +0.41 (BH q 0.073) — but the **housekeeping control drifts +0.36** | ❌ not separable |
| Re-exposure returns along the disease axis | leave-one-patient-out median cos +0.52 (and the GFD residue lies on the same axis, +0.67) | ✅ third disease |

**Rules it sets.** (1) In a cause-removal cohort, compute corr(group, confounder)
before anything else; above ~0.7 the regression adjustment is not identifiable
and only a matched window answers the question. (2) In bulk tissue with large
composition change, housekeeping genes are not a negative control; report the
control module's own drift and count only the excess.

### Three mechanisms, one disease (D23)

GSE11903 (etanercept) + GSE102725 (secukinumab) + GSE69967 (tofacitinib, with a
placebo arm): three mechanisms, three platforms, per-patient effect vectors.

| Check | Result | Status |
|---|---|---|
| The three cohorts share one disease axis | cross-cohort cos +0.51 to +0.72, the same range as patient-to-patient cos inside one cohort (+0.48 to +0.66) | ✅ a disease axis transfers across studies |
| The restoration ordering repeats under every mechanism | barrier/lipid least restored in all three (+0.53 etanercept, +0.38 secukinumab, +0.77 tofacitinib) | ✅ |
| Barrier restoration differs by mechanism | tofacitinib +0.77 vs secukinumab +0.38 — but different visit (12 vs 16 weeks) and cohort | ⚠️ hypothesis |
| **A drug-specific response component** | naive test: same drug +0.812 vs different drug +0.538, p 2e-18. **Confounded**: one cohort per drug. Same-platform pair 0.727 vs 0.815 within (p 0.005); cross-platform 0.483 | ❌ not attributable |
| Mechanism-specific module (JAK clears interferon first) | excess over the other modules +0.05 under tofacitinib vs +0.01 under anti-TNF, −0.21 under anti-IL-17A | ❌ |
| Placebo arm as a negative control | n = 3, median restoration ≈ 0 | ⚠️ too small |

**Rule it sets.** With one cohort per drug, every drug-vs-drug contrast is also a
study-vs-study contrast. Estimate the technical gap from a same-platform subset
before interpreting anything as mechanism.

### Does within-cell co-occurrence carry information the means do not? (D24)

Frangieh/Izar 2021 Perturb-CITE-seq melanoma (scPerturb): 218,331 cells x 24
surface proteins including HLA-A, CD58 and four isotype controls, CRISPR KO
across control / IFN-gamma / TIL co-culture.

| Check | Result | Status |
|---|---|---|
| Depth artefact in a within-cell covariance | raw ADT: isotype x isotype +0.074, HLA-A x CD58 +0.132, r(HLA-A, depth) +0.55; CLR removes it | ✅ quantified |
| The ADT panel resolves its own knockouts | B2M KO → HLA-A −1.46 (z −56), CD274 −1.37, CD58 −0.94 (z −30) | ✅ |
| **Per-knockout covariance is estimable** | split-half r −0.005 to +0.170, and the isotype pair scores as high (+0.209); means split at r 0.78–0.94 | ❌ not estimable |
| Cells needed for a usable estimate | bootstrap SE 0.030 at n=100 vs a 0.020 total between-KO spread; ~3,000 cells/KO for SE = spread/4 | ✅ design number |
| The screen's resistance read-out is sound | STAT1 +1.25, IFNGR2 +1.11, JAK1 +0.84, CD58 +0.35 (median over all KOs +0.03) | ✅ reproduces the paper |
| Covariance predicts resistance beyond the means | added R² +0.001, permutation p 0.71; ρ(ΔCov) −0.09 n.s. — while ρ(Δmean CD58) −0.343 (p 8e-7) and ρ(Δmean HLA-A) −0.216 (p 0.002) | ❌ |
| Size of the co-occurrence term | sd(Cov) is 10–21% of sd(E[A]E[H]) | ⚠️ minor term |

**Rule it sets.** A second moment is not a mean: before interpreting any
within-cell covariance, (1) normalise away the shared per-cell total and prove
it on a pair that cannot co-vary, (2) split-half the statistic itself, and
(3) convert its bootstrap SE into a cells-per-group requirement. The hypothesis
here is under-powered by 10–30x, not refuted.

### Is there a slow mode behind "rest" in CAR-T exhaustion? (D25)

Weber 2021 (GSE164949 RNA-seq: 3 donors x HA-28z CAR-T at days 7/11/15, always
ON / always OFF / transient rest from day 7 or 11 / ON + nivolumab; GSE164947
ATAC + H3K27me3).

| Check | Result | Status |
|---|---|---|
| One exhaustion axis across days | cross-day cos +0.55 to +0.72 vs +0.51 between donors within a day | ✅ |
| Rest restores, as a dose-response | 4 d: +0.63 and +0.85 (two independent arms); 8 d: +1.07 | ✅ |
| The read-out is specific (not "any intervention restores") | nivolumab +0.05 (day 11) and +0.35 (day 15) vs rest, paired p 0.0033 / 0.040 | ✅ |
| **A slow mode inside the transcriptome** | every defined module has tau 1.2–2.6 days (ratio 2.3); effector function overshoots to +1.87 | ❌ one fast mode |
| Housekeeping module undefined | 1 disrupted gene | ✅ |
| RNA can proxy the chromatin arm | the epigenetic-machinery module is undefined — EZH2/DNMT/TET2 transcripts do not differ | ❌ |
| The cross-modality (false recovery) test | deposited ATAC has no rested arm; H3K27me3 matrix labels (D11/HA/SI) are not mappable to rested/always-off/always-on | ❌ design gap |

**Rule it sets.** Before claiming two timescales, check that the slow modality is
actually measured in the same rest/recovery arms. A transcriptome time course can
show that there is no slow mode *in RNA*, which is not the same as showing there
is none — and the expression of chromatin machinery is not a usable stand-in.

### Can a composition snapshot prove a blocked exit? (D26)

GSE135893 (human IPF, Habermann 2020) + GSE141259 (mouse bleomycin time course,
Strunz 2020), cell-level annotations only. Three-state repair model
AT2 -> transitional -> AT1.

| Check | Result | Status |
|---|---|---|
| Fractions equal residence-time ratios (1/k1 : 1/k2 : 1/d) | max relative error 4e-16 over 200 random parameter sets | ✅ exact |
| **D/T = d/k2, independent of k1 and s** | spread 1.1e-16 with k1 over 5 orders of magnitude — entry acceleration cannot move the ratio at all | ✅ sharper than the "ceiling" claim |
| k1 -> infinity fraction ceiling = healthy f_D/(f_D+f_T) | predicted 0.3636 = simulated | ✅ |
| Human: disease exceeds the fraction ceiling | ceiling from controls 0.484 vs disease f_D 0.121 | ❌ not exceeded |
| Human: D/T rises | 2.55x but [0.35, 11.33] on a **donor-level** bootstrap (10 vs 20 donors) | ❌ not significant |
| Mouse: usable baseline | PBS has 3 transitional and 20 AT1 cells; days 2 and 4 give D/T of 40 and 72 off 2–3 AT1 cells; whole-lung deposit has AT1 = 38 of 29,297 | ❌ |
| **AT1 recovery is condition-independent** | mouse 0.84 % (PBS) to 3.3 % (injury) although bleomycin destroys AT1; human 8.9 % to 5.0 %; r = 0.39 erases the human effect | ❌ decisive |

**Rule it sets.** The bound is real and parameter-free, but dissociation-based
composition cannot supply its inputs when one compartment is fragile. Ask what
relative recovery would erase the result before interpreting any ratio of cell
types, and get absolute rates from pulse-chase, not from a snapshot.

**What the step-5 results mean.** A J fitted from one observational intervention describes how the system relaxes under that intervention. Its "levers" are the modules that relax least by themselves: bottlenecks and candidate add-on targets.

Confirmed in D16: with a measured response library the same formula recovers the
right lever 66 % of the time against a 7 % chance rate, while covariance-based
ranking stays at 15 %. To rank causal upstream targets, J must come from one of these:
- perturbations of individual modules or genes (Perturb-seq MRA, CRISPR screens in the relevant cell type);
- several interventions that act on different modules;
- a mechanistic kinetic model (the network-jacobian skill).

D20 adds the ceiling on this: even with the right cell type and the right
lever in the library, the ranking is an indication-level statement. The
clinical response to a drug cannot confirm its own mechanism, because other
mechanisms in the same disease produce almost the same response vector.

## D27 — EGFR / AML feasibility audit (scenarios 4 and 5): what proportions pin down

Run as the NJ-FEAS-01 product, not a result product. `demos/d27_egfr_aml_feasibility.py`.

**Structural identifiability (algebra, verified to machine precision).**

| Claim | Check | Verdict |
|---|---|---|
| EGFR: p(t) obeys dp/dt = a(1−p) − bp + (gR−gS)p(1−p) | shifting gS and gR by the same +0.25/day moves p(t) by 2.2e-16 while N(11) changes 16-fold | ✅ the common growth level is an exact null direction |
| EGFR proportions see only (a, b, dg) | FIM rank 3/4, null vector (0, 0, +0.707, +0.707) | ✅ |
| EGFR + ln N is full rank (negative control for the rank test) | rank 4/4, cond 1.0e+04 | ✅ the deficiencies are structure, not finite-difference noise |
| "Proportions cannot tell switching from selection" | both a and dg are structurally identifiable | ❌ the common claim is wrong |
| AML: r = M/I obeys dr/dt = delta − Cr, C = gI − dI − delta + dM | two parameter sets with equal (delta, C) agree in mature fraction to 1.3e-15 (0.280 at day 7) with 67-fold different totals | ✅ |
| "Proportions cannot distinguish differentiation from killing" | **delta IS identifiable** from proportions (CV 11% at 2,000 cells x 6 days); the kill rate is not | ❌ half wrong, and the wrong half is the one that was being sold |
| Absolute counts fix the AML kill rate | rank goes 2/4 → 3/4 only; gI and dI stay lumped (null vector (+0.707, +0.707, 0, 0)) | ❌ a third channel is required |

**Practical identifiability is a separate question.** In the (a, b, dg) basis at
6 timepoints and sigma_p = 0.02, the best-determined combination is essentially
a alone (sd 0.0045, CV 6%) while b and dg are confounded (worst direction
−0.01a +0.60b +0.80dg, sd 3.78; block cond 7.1e+05). Report both layers: the
structural layer says what no amount of data can fix, the practical layer says
what this design actually estimates.

**The real deposit (GSE149383, PC9 + erlotinib Drop-seq).**

| Finding | Number | Consequence |
|---|---|---|
| Re-deposition across subseries | 3 of 14 advertised PC9 libraries are byte-identical after decompression (D0 = D0_b, D2 = D2_b, D11_ERL = D11_ERL_c) | "three independent replicates" were one library; always hash before counting replicates |
| Depth confounded with day | corr(day, log median UMI) = +0.865 (2,142 → 8,425 UMI from day 0 to day 11) | the time axis and the depth axis cannot be separated by design |
| PC1 is depth | corr(PC1, log UMI) = +0.964 raw, −0.013 after exact thinning to 1,000 UMI/cell | any unsupervised state call on raw counts is a depth call |
| Depth matching rescues reproducibility | day-11 p across three independent libraries: 0.870 / 0.341 / 0.643 raw (spread 0.53) → 1.000 / 0.954 / 0.971 (sd 0.023) | the state fraction IS reproducible; depth was moving it |
| The second state is not well defined | three definitions (label-free mixture, supervised D11−D0 axis, DTP marker score) agree pairwise at only +0.44 to +0.86; the label-free one tracks proliferation arrest (cycle shift −0.369) and is non-monotone in time | define the second state before fitting any two-state model |
| The fit pins nothing | a 95% profile [0.05, 5.74]/day; dg's profile spans the whole [−5, 5] grid; b and p0 sit on bounds. Half-conversion 0.22 day, so the design samples the plateau | good chi2 (0.49 on 2 dof) is not evidence of estimability |
| Circularity check | refitting with the independent libraries at days 0 and 11 (p 0.116 → 0.963) gives a = 1.43 and the same unconstrained profiles | the conclusion is not an artefact of the axis construction |
| Burden span | 6.0e+04-fold at day 11 for gS in [−0.5, +0.5]/day with dg fixed; 4.6e+47 with dg free | efficacy is unidentifiable by construction |
| Barcodes are not cells | deposited barcode counts span 24-fold (494–12,000) with 52–95% empty droplets; the xenograft arms of the same study carry 29,340–29,983 cells, CV 1.0% | recovered cell number is set by loading, never by survival |

**Audit of what public deposits carry** (top 60 series of each of two
reproducible GEO queries, scanning `!Sample_characteristics`):

| Query | GEO hits | audited | >=3 time values | any count field | both |
|---|---|---|---|---|---|
| EGFR drug-tolerant persister | 183 | 60 | 1/60 | **0/60** | **0/60** |
| AML differentiation therapy | 174 | 60 | 1/60 | **0/60** | **0/60** |

**Design that unblocks each scenario** (Cramer-Rao, not opinion). EGFR: without
a count, gS is not estimable at 500, 2,000 or 10,000 cells per timepoint; with
2,000 cells and a count at CV 5%, CV(a) 0.11, CV(dg) 0.48, CV(gS) 0.09. Densifying
sampling does not help dg — 6 → 48 timepoints only takes CV(dg) from 1.81 to 0.63.
AML needs three channels: proportions alone give delta at CV 0.11 and nothing
else; adding counts buys dM (0.25); only adding a proliferation readout (EdU or
dye dilution, CV 20%) makes the design full rank, and even then dI lands at CV 1.04.

**Rules it sets.**
- Split every feasibility answer into structural identifiability (what no amount
  of data fixes) and practical identifiability (what this design estimates).
  Quoting only the first overstates capability; quoting only the second makes
  "collect more cells" look like a solution when it is not.
- Before claiming replicates from a GEO SuperSeries, hash the decompressed files.
- Before any unsupervised state call, regress or thin out sequencing depth and
  report corr(PC1, log depth) both ways.
- A two-state model needs the second state defined first; report the agreement
  between at least two independent definitions of it.
- A good chi2 with a flat profile is not an estimate. Profile every parameter
  that enters a feasibility claim.
