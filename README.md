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

## netjac — network-Jacobian diagnostics platform

A second package, [`netjac/`](netjac), builds on exact BioModels kinetics (SBML
via libRoadRunner). It keeps three objects separate: the state Jacobian,
which gives eigenvalues and stability; the elasticities, which are the network
edges; and the flux-control matrix dFlux/dEnzyme. It then demos three
diagnostic directions:

| Demo | System | Key result |
|---|---|---|
| `d1_edge_biomarker` | Teusink 2000 glycolysis (BIOMD64) | loss of PFK–ATP feedback leaves the metabolome identical (max \|Δln x\| 9e-16). Static panel AUC 0.62 vs a model-predicted challenge edge readout AUC 0.91. Passive covariance inversion needs ~1e6 samples (stiff) |
| `d2_perturbation_relaxation` | Teusink, Wolf 2000 oscillations (BIOMD206), **real CGMacros CGM (n=45)** | eigenvalues recovered from 1–2-analyte pulses; Hopf read as ringing damping from NADH alone. Real data: meal recovery rate vs HbA1c ρ=−0.44 (p=0.002), but no information beyond fasting glucose |
| `d4_patient_jacobian` | Holzhütter 2004 human glycolysis/PPP (BIOMD70) × **TCGA-KIRC (537 tumours + 72 normals)** | 609 patient-specific Jacobians. Recovers Warburg and an ALDOB-loss aldolase bottleneck (10 % of tumours lose the steady state, OS HR 1.69). Every survival signal reduces to a single-gene proxy (BPGM, ALDOB). The chassis is the bottleneck |
| `d5_pairwise_to_network` | user's tumour–immune & macro models, Goodwin, random nets, **curated physiology atlas** (Topp, Sturis, Tolić, Lemaire, Komarova, Miao, Borghans) | Every equilibrium classified; 2×2 sub-systems combined exactly through the Levins loop decomposition + Hurwitz. All pairs stable ≠ whole stable (Sturis: 0/15 unstable pairs, 110-min Hopf cycle via the 5-step delay loop). Topp diabetes saddle-node at d0 = 0.0735 (= analytic). A linear centre is not a limit cycle. May bound |
| `d6_baroreflex_jacobian` | **real** beat-to-beat ECG + BP, PhysioNet Autonomic Aging | 286 subjects. BRS positive control: ρ=−0.54 with age (p=9e-18). With age the delay-embedded loop mode decays more slowly (ρ=−0.23, beyond BRS p=6e-4) and the VAR spectral radius moves toward 1. The pairwise 2×2 Jacobian alone cannot produce the Mayer oscillation |
| `d3_fragility` | Teusink, Wolf | total fragility predicts Hopf crossing 100 % (local-only formula r=0.16). Zero-flux-control vPDC is a collapse node. Best probe = critical-mode participation (ρ=+0.96), not abundance (ρ≈0) |

```bash
pip install numpy scipy networkx matplotlib libroadrunner python-libsbml scikit-learn pandas
python3 -m netjac.run_all            # figures -> figures/netjac/ ; --offline skips real data
python3 -m netjac.run_all --tcga     # + D4 patient-specific Jacobians (~15 min first run)
```

Chinese write-up: [`NETJAC_REPORT.md`](NETJAC_REPORT.md). Skill:
`.claude/skills/network-jacobian/`.

## Run

```bash
pip install numpy scipy networkx matplotlib
python3 -m grn_pipeline.run_all       # full pipeline + figures
python3 -m grn_pipeline.m1_symmetry   # or any single module
```

Figures are written to `figures/`. A full write-up with numbers, rigour
labels, and the interpretation (including the Lyapunov-exponent biomarker
question) is in [`REPORT.md`](REPORT.md).
