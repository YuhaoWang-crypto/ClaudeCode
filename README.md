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

# mpro-pipeline

Enzymatic QSAR for SARS-CoV-2 Mpro (3CLpro), i.e. the in-silico counterpart of
a purified-enzyme FRET protease assay (BPS Bioscience #79955 type: recombinant
3CLpro + DABCYL-KTSAVLQ|SGFRKME-EDANS, Ex360/Em460, GC376 control). Buffer and
protein only — **no cell**. Built to four standing rules: measured pIC50 as the
only label, pool-then-quarantine instead of stratify, preincubation as a tier,
and a null model before any reported statistic.

| Module | Job | Key measured result |
|---|---|---|
| `m1_labels` | label QC on 3,629 enzymatic IC50 / 2,858 compounds / 507 assays | per-assay stratification impossible (median 2 records, largest 77); pooling sound for the 59.2% agreeing within 0.5 log; 2,787 compounds retained after quarantining the >1 log tail |
| `m2_null_model` | descriptor bar every score must clear | strongest null \|ρ\|=0.486 (cLogP); Boltz `ptm`/`plddt` score *below* it |
| `m3_boltz_calib` | pre-registered Boltz-2 calibration, Mpro **dimer** + monomer control | `optimization_score` ρ=**+0.779** (p=0.005, CI [+0.33,+0.94], +0.293 over null) → trust for ranking; `iptm` +0.519 but only +0.033 over null → useless |
| `m4_qsar` | scaffold-split QSAR, null-gated, tier-reported | RF ρ=**+0.718**, R²=+0.516, RMSE 0.65 log (4.5× in IC50), +0.420 over null; tier A ρ=+0.764 vs tier B ρ=+0.568 |
| `m5_warhead_mmp` | covalent-warhead matched pairs mined from the measured data (29 clusters found, 4 co-folded, 44 compounds) | within-series pairwise: warhead swaps **70.1%** (61/87), recognition swaps **71.8%** (74/103), large effects (≥1.5 log) **79.5%**; blocked within-cluster ρ=+0.536 (permutation p=0.0001) |

## What this pipeline REFUTES

**"opt_score is warhead-blind, so covalent docking can be skipped."** This was
our own hypothesis after M3 — a noncovalent model reaching ρ=+0.78 on covalent
inhibitors invited it. `m5` tested it on matched pairs differing only in
warhead and it does not hold: warhead swaps (70.1%) and recognition swaps
(71.8%) are called at the same rate. The model sees both, weakly.

What *is* licensed: noncovalent co-folding calls ~70% of within-series pairs
correctly, ~80% when the true gap is ≥1.5 log. That is the bar a covalent or
QM/MM method must clear to be worth its cost — not a reason to skip it.

`m5` also documents a power trap worth remembering: the pre-registered
per-cluster test returned "blind" for all four clusters (every bootstrap CI
spanning 0, n=6–17). A blocked within-cluster analysis on the same data gives
ρ=+0.536 at permutation p=0.0001. **Underpowered is not the same as negative.**
The blocked test is labelled post-hoc in the module.

## What this pipeline does NOT establish

- **Dimer > monomer.** Direction is consistent on both discriminative metrics
  (Δρ=+0.123 each) but the bootstrap CI on Δ spans 0 (P(Δ≤0)≈0.17). The dimer
  is required by structural biology (the N-finger of one protomer completes the
  other's S1 pocket); n=11 does not demonstrate it.
- **Stereo-SAR.** Boltz scores nirmatrelvir and its γ-lactam epimer within
  0.008 `optimization_score` and identically on `binding_confidence` (0.999).
- **Anything cellular.** Enzymatic pIC50 explains ~14% of cellular antiviral
  variance (ρ=+0.43, n=615 paired); the rest is permeability / efflux /
  glutathione / esterase, invisible here. Different axis, out of scope.

## The label ceiling

Nirmatrelvir's own measured pIC50 spans **6.12–9.10 (2.98 log)** across 35
assays; ebselen's spans 4.99–8.00. RMSE near 0.65–0.9 log is the noise floor of
the labels, not a modelling failure.

## Two traps this code records

1. ChEMBL `molecule/search` returns an **unnamed duplicate** of nirmatrelvir
   (`CHEMBL5201264`, `pref_name: None`, 2 records) *ahead of* the curated entry
   (`CHEMBL4802135`, 35 records); they differ at one γ-lactam stereocentre.
   Taking result `[0]` submits the wrong epimer. Check `pref_name` and record
   count, never search order.
2. Boltz SMARTS filtering must be **disabled** for a calibration set — at the
   default `recommended` level ebselen (Se) and disulfiram (thiuram disulfide)
   are dropped before prediction and the weak end of the potency range
   disappears.

## Run

```bash
pip install numpy scipy rdkit scikit-learn
python3 -m mpro_pipeline.run_all       # full pipeline
python3 -m mpro_pipeline.m1_labels     # or any single module
python3 -m mpro_pipeline.m5_warhead_mmp   # the warhead matched-pair test
```

Cached inputs are in `mpro_pipeline/data/` (ChEMBL labels, the exact SMILES
submitted to Boltz, and the returned metrics for all three target
configurations), so every number above reruns offline.
