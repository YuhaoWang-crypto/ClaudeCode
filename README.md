# grn-pipeline

A small, fully-runnable pipeline that applies four "irreducibility / symmetry"
mathematical tools to gene-regulatory and metabolic networks, on concrete
literature-grounded systems where every number is *computed*, not asserted.

| Module | Tool | System | Key result |
|---|---|---|---|
| `m1_symmetry` | graph automorphism → quotient | RTK/RAS/RAF/MEK/ERK | \|Aut\|=S₃; 9→7 nodes (3 RAS paralogues → 1 core) |
| `m2_crnt` | CRNT deficiency δ | A⇌B⇌C vs Schlögl | δ=0 monostable / δ=1 bistable switch |
| `m3_efm` | elementary flux modes | 4-metabolite network | 3 irreducible flux generators span the cone |
| `m4_dnb_lyapunov` | DNB / critical slowing / Lyapunov | 2-gene fold bifurcation | LLE→0, SD/autocorr/DNB rise at tipping point |
| `m5_kras_real` | symmetry breaking on a real target | KRAS G12C + covalent drugs (ChEMBL/Boltz/Inductive Bio) | covalent G12C drug breaks paralog symmetry S₃(6)→S₂(2) |
| `m6_integrate` | binding → network stability | sotorasib vs adagrasib | real ChEMBL+Boltz binding → engagement → DNB biomarker |
| `m7_screen` | Boltz-2.1 library screen | 10 G12C ligands (ChEMBL) | ranked by binding; 4 analogues out-rank sotorasib |
| `m8_clinical` | biomarker → trial endpoints | CodeBreaK 100 (NCT03600883) | layers mapped to ORR/DOR, PFS/OS, Cmax/AUC, QTc |
| `m9_occupancy` | PK occupancy → μ calibration | sotorasib (IC50=30 nM) | 98% occupancy at approved dose → network near tipping |
| `m10_validate` | Boltz ranking vs ChEMBL truth | 5 G12C ligands w/ measured IC50 | opt_score tracks potency (ρ=+0.6); binding_confidence doesn't (ρ=−0.2) |
| `m11_fibration` | input-tree fibration (Morone) | expanded MAPK paralogue graph | 27→11 fibers; generalises M1 automorphism to fiber representatives |
| `m12_dualphos` | real ERK double-phospho core | Markevich-style mass-action | CRNT deficiency δ=2; bistable ERK switch; EFM = 2 futile cycles |
| `m13_fim_sloppy` | FIM / sloppy / stiff axes | ERK dual-phospho ODE | sloppy spectrum (38 orders); flux-ratio observables load best |
| `m14_atlas` | 18-pathway systematic atlas | JAK-STAT…mevalonate | fibration compression + biomarker class per pathway; JAK-STAT top (3.0×) |
| `m15_markevich_mm` | exact Markevich 2004 MM ERK cycle | published parameters (JCB 2004) | reproduces bistable window [39.25, 57.38] nM + 3-state table to the decimal |
| `m16_erk_dnb` | DNB / critical slowing on the real switch | M15 saddle-nodes 39.25/57.38 nM | λ_max→0, τ≈4720 s at boundaries; SD/autocorr/DNB rise (early warning) |
| `m17_realdata` | validate M16 on real single-cell ERK imaging | Pertz-lab EKAR traces (FGF pulses + EGF dose) | lag-1 autocorr rises before ERK pulses (p≈0.004); variance flat; EGF all supra-threshold — partial validation |
| `m18_titration_benchmark` | positive control: MEKi titration across the real bifurcation | simulated from M15 Markevich switch | variance & lag-1 autocorr PEAK near threshold (≈5×), tracking τ — pipeline is sensitive, not blind |
| `m19_switch_library` | migrate the critical-slowing engine to many pathways | 10 canonical bistable switches (MAPK, Rb-E2F, apoptosis, Cdc2, CaMKII, Wnt, Cdc42, lac, Schlögl, master-TF) | 10/10 show variance+autocorr rising to their saddle-node — biomarker is universal to the bifurcation |
| `m20_literature_bistable` | multi-variable literature switches + hysteresis | Rb-E2F (Yao 2008), apoptosis (Eissing 2004) topology | bistability + hysteresis loops reproduced; eigenvalue→0 at folds (exact params not open-access-fetchable) |
| `m21_oscillators` | extend framework to oscillatory (Hopf) pathways | Goodwin (circadian), p53-Mdm2, Brusselator (glycolytic) | approaching Hopf: variance rises AND a spectral peak sharpens at the intrinsic frequency — distinct from saddle-node |
| `m20b_biomodels_exact` | fetch + simulate EXACT curated models (fills M20 gap) | Markevich2004 (BIOMD27), Legewie2006 apoptosis (BIOMD102) | download method = biomodels GitHub mirror + libRoadRunner; official Km5=78 confirms hand-coded M15 (states to the decimal); Legewie caspase switch bistable in XIAP synthesis |
| `m22_snic_mixed` | mixed bifurcation: saddle-node ON a limit cycle (SNIC) | θ / Ermentrout-Kopell normal form (cell-cycle / excitable) | finite-amplitude spikes whose period diverges (T~π/√I, log-log slope −0.50; frequency→0) — signature distinct from both Hopf and pure saddle-node; ISI mean+CV both grow |

## Run

```bash
pip install numpy scipy networkx matplotlib
python3 -m grn_pipeline.run_all       # full pipeline + figures
python3 -m grn_pipeline.m1_symmetry   # or any single module
```

Figures are written to `figures/`. A full write-up with numbers, rigour
labels, and the interpretation (including the Lyapunov-exponent biomarker
question) is in [`REPORT.md`](REPORT.md).

---

# recorder-pipeline

A second, IVD-facing package. Where `grn_pipeline` asks *which network node is
irreducible or near a tipping point*, `recorder_pipeline` asks the assay
question: **given that a disease process happened, where was it written, how
long does the writing survive, and what should it be divided by?**

It formalises the "Disease Recorder" framework (biomarkers as integrals of past
pathology — HbA1c, not glucose) into a five-tuple `⟨Writer, Carrier, Kernel,
Detector, Normalizer⟩`, and turns three of its qualitative claims into closed
form.

| Module | Tool | Key result |
|---|---|---|
| `r1_kernel` | renewal theory on carrier age | the persistence kernel is the carrier **age** survivor, not its decay curve. One input (RBC lifespan 120 d) reproduces HbA1c's clinical properties: 43.8% of signal from the last 30 d, 90% horizon 11.7 weeks, mean lag exactly L/3 |
| `r2_pairing` | Fisher separation in log space | a fixed ratio A/B beats the single marker **iff ρ > κ/2** (κ = denominator noise ratio) — so a denominator >2× noisier can never help; the fitted residual gains `1/√(1−ρ²)` and can never hurt. Derives the framework's own L1>L3>L0 denominator hierarchy |
| `r3_individuality` | biological-variation screen | 11/15 routine analytes have index of individuality < 0.6. The four that already drifted to personal-baseline rules in practice (creatinine/KDIGO, hs-troponin delta, CA125 ROCA, PSA velocity) are exactly the lowest-II ones |
| `r4_catalog` | taxonomy as data | 29-entry machine-readable recorder catalog (`recorder_catalog.tsv`) across 10 classes, scored by an explicit six-gate design rubric |
| `r5_datacases` | published clinical values → feasibility | effect sizes span **two orders of magnitude** (GFAP 51.8× vs PRO-C3 1.7×). Extends r2 with measurement error: a ratio helps iff **ρ > (κ²+a²)/(2κ)**, so PRO-C3's 11% inter-assay CV raises its required correlation from 0.50 to 0.62. Carbamylated albumin cannot reach AUC 0.70 for its own intended use at any sample size |

```bash
python3 -m recorder_pipeline.run_all      # all five modules + figures
```

Three write-ups:

- [`RECORDER_FRAMEWORK.md`](RECORDER_FRAMEWORK.md) — the formalisation, the
  three derived laws, and 10 recorder classes (5 added beyond the source
  framework: glycoform, failed-writing, epigenetic, autoantibody-amplified,
  clonal).
- [`RESCUE_MINING.md`](RESCUE_MINING.md) — six failure modes mapped to six
  transforms, then 12 worked case studies re-mining clinically failed
  biomarkers (NGAL, NT-proBNP, CA-125, creatinine, CRP, serum HER2, MMP-9,
  anti-dsDNA, procalcitonin, plasma Aβ42, AFP, total oxidative markers).
- [`CASE_STUDIES.md`](CASE_STUDIES.md) — the same candidates tested against
  **real published normal-vs-disease values**: data-space table, achievable
  AUC, analytical error budget, Gate-2 sample sizes, and a speed-to-implement
  ranking. Headline: the fastest new readout available (carrier-age-corrected
  cell-bound C4d in SLE) needs **no new assay at all** — EC4d/BC4d are
  commercial and reticulocyte/immature-platelet fractions already print on
  every routine CBC.
