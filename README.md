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
| `m23_virtual_tcell` | a "virtual T cell": kinetic proofreading → ERK switch → relative IL-2 | McKeithan 1995 chain + Altan-Bonnet/Germain 2005 feedback topology | proofreading exponent → N+1 (rigorous, verified against the closed form); a 33× ligand-dose change moves the quality threshold only 2.58×; ERK hysteresis window τ∈[0.02, 1.03] s. **Ligand ranking only, never pg/mL** |
| `m24_rna_secretion_coupling` | does cytokine mRNA predict what the same cell secreted? | TRAPS-seq, GEO GSE200690 (4,875 + 3,223 cells) | own-mRNA explains R²=0.25 (IFN-γ), 0.12 (TNF), **0.02 (IL-2)** of per-cell secretion; a cell's earlier secretion predicts its later secretion far better (ρ 0.55–0.80) and adding mRNA buys +0.002…+0.033 R² |
| `m26_secretion_model` | audit the secretion predictor, then beat it | TRAPS-seq, depth-residual target, leave-one-hashtag-out | raw-count targets score up to **0.85 from capture depth alone**; under the clean protocol the best models reach 0.841 (IFN-γ), 0.776 (TNF), 0.696 (IL-2), and a cell's own earlier secretion still beats its transcriptome |
| `m27_aim_response_audit` | antigen response on peptide-pool-stimulated T cells | AIM-seq, 43,222 cells, 5 donors, leave-one-donor-out | AIM+ vs AIM− reaches **0.940** (vs 0.882 previously reported), but 0.864 of it is RNA content and the hashtag total alone scores 0.907 — a staining artifact, flagged not used |
| `m25_peptide_to_cytokine` | the full chain: peptide → signalling → IL2 mRNA → secreted IL-2 | OT-I altered-peptide series + live NetMHCpan via IEDB | presentation cannot order MHC-matched APLs (Spearman **+0.36** overall, **−0.54** on the six stimulatory ones); above threshold the chain compresses a 700× potency range into ~1.2× output |

## Can a virtual cell replace a functional T-cell assay?

Short answer: no, and [`VIRTUAL_TCELL_REPORT.md`](VIRTUAL_TCELL_REPORT.md) says
exactly why, with the numbers computed rather than asserted — including a full
pass over the IEDB T-cell assay table (577,789 records, of which 1.0% carry any
quantitative value) and over the CELLxGENE index (168M T cells, single-digit
stimulation datasets). `m23_virtual_tcell` is the model that report is about, and `m25_peptide_to_cytokine`
is the coupled peptide→pathway→RNA→cytokine chain, built so that each joint is
tested rather than assumed.

## Run

```bash
pip install numpy scipy networkx matplotlib
python3 -m grn_pipeline.run_all       # full pipeline + figures
python3 -m grn_pipeline.m1_symmetry   # or any single module
```

Figures are written to `figures/`. A full write-up with numbers, rigour
labels, and the interpretation (including the Lyapunov-exponent biomarker
question) is in [`REPORT.md`](REPORT.md).
