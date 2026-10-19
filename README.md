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

# landscape_flux — 势能–通量景观诊断设想的公开数据实测

Seven runnable demos that test the "landscape–flux" diagnostic proposal
(entropy / potential / flux / macroscopic physical readouts) on public data
or PMID-referenced literature measurements. Full write-up (Chinese):
[`LANDSCAPE_FLUX_REPORT.md`](LANDSCAPE_FLUX_REPORT.md); cross-check against an independent
analysis: [`COMPARISON_WITH_SYNTHESIS.md`](COMPARISON_WITH_SYNTHESIS.md).

| Module | Claim tested | Data | Key result |
|---|---|---|---|
| `lf1_maxent` | macro entropy → P(X); MaxEnt residual exposes hidden interactions | pancreas scRNA, PBMC-3k | equal entropy (ΔH=0.022) yet AUC 0.51 vs 0.97 for one gene; 79% of MaxEnt residual = hidden cell identity |
| `lf2_landscape_flux` | curl flux & entropy production in gene circuits | exact Fokker–Planck on grid | symmetric A⇄B σ≈5e-4 (near detailed balance) vs A→B⊣A σ=7.7 |
| `lf3_flux_realdata` | measuring flux without / with time series | GEO GSE3431, pancreas RNA velocity | YMC rotation p=2e-4; 81% of cycling cells rotate G1→S→G2M (p=0.001); Hodge rotational 0.46 (cycle) vs 0.06 (differentiation) |
| `lf4_hamiltonian` | U + K = const ("potency + metabolism = const") | model, pancreas, Paul15 | refuted (corr +0.91 / +0.26 / ≈0 after depth control); rotational flux runs along iso-U contours |
| `lf5_nuclear_texture` | nuclear fractal dimension / texture entropy as chromatin readout | BBBC021 MCF-7 (7923 nuclei) | taxol texture entropy AUC 0.71 vs area 0.52, but 4 trivial scalars incl. sharpness already 0.75 (→0.77 with entropy/FD); box-counting FD fails synthetic check |
| `lf6_dielectric` | membrane capacitance tracks EMT / T-cell state | single-shell model + Lannin 2016, Pethig 2002/2007 | model reproduces crossovers ≤5%; EMT C_m unchanged (20.2→20.1); T-cell activation single-cell AUC≈0.66 |
| `lf7_calorimetry` | heat output reads out state entropy production | published pW/cell + TAM limits | 2.32 pW = 5.4e8 k_BT/s; circuit ≈0.2% of heat; single-cell fluctuations averaged to 0.09% |
| `lf8_jacobian_slowing` | critical slowing at the Pre-endocrine fork (replication of an independent report) | pancreas RNA velocity | raw τ 1.16x vs shuffled-velocity null 1.17x across 2–10 PCs: not replicated |

```bash
pip install numpy scipy pandas matplotlib scikit-learn scikit-image anndata h5py scanpy pooch tifffile
python3 -m landscape_flux.run_all
```
