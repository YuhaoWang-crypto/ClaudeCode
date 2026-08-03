# Repository contents

Two independent, fully-runnable pipelines:

1. **[`grn_pipeline/`](#grn-pipeline)** — irreducibility / symmetry / dynamical-systems
   analysis of gene-regulatory and metabolic networks.
2. **[`tissue_evolution/`](#tissue-evolution)** — do different tissues evolve at
   different speeds, and can that be used to reconstruct evolutionary conditions?

---

# tissue-evolution

Everything is computed from primary data (GTEx v10, Ensembl 116, Cardoso-Moreira
2019 via Expression Atlas, AlphaFold DB, ESM-2). Nothing is quoted from a paper.

```bash
pip install numpy scipy pandas matplotlib statsmodels biopython
python3 -m tissue_evolution.run_all            # full pipeline (t12 downloads 1.2 GB of TOGA alignments)
python3 -m tissue_evolution.t5_confounders     # or any single module
```

| Module | Question | Key result |
|---|---|---|
| `t1_expression` | tissue-specificity τ over GTEx | 32 organs, 46325 genes; brain sub-regions collapsed first, or τ is an artefact of column count |
| `t2_orthologs` | 1:1 orthologues + CDS | 16.4k (macaque) → 12.7k (chicken) pairs |
| `t3_dnds` | dN/dS **from scratch** | Ensembl no longer ships dN/dS; vectorised NG86 matches Biopython to <1e-13 |
| `t4_tissue_rates` | 3 definitions of "tissue rate" | identity genes span **4.10×**, whole transcriptome only **1.13×** |
| `t5_confounders` | composition vs intrinsic | only immune/barrier + testis (fast) and CNS + muscle (slow) survive τ/expression matching |
| `t6_depth` | 29 → 319 Mya | organ dN ranking stable (ρ=0.87); ω degrades as dS saturates (74% dS>1 at chicken) |
| `t7_esm` | constraint without an orthologue | ESM reaches the 11–59% of expression dN/dS structurally cannot |
| `t8_structure` | Foldseek TM-align vs sequence | **negative result**: TM-score tracks AlphaFold pLDDT (ρ=+0.84), not divergence (ρ=−0.16) |
| `t9_lineage` | which *lineage* sped up? | **negative result**: no organ survives BH + matching; generation-time effect dominates |
| `t10_expression_divergence` | is the tissue profile itself conserved? | median human–mouse profile r=+0.80; coupling to ω is real but stage-pooling-dependent |
| `t11_drivers` | few dominant proteins, or the bulk? | **the bulk**: ranking survives deleting the fastest 25% (ρ=+0.97); 78% of testis genes sit above the genomic median |
| `t12_rer` | real RERconverge over Zoonomia/TOGA | 4195 genes x 240 branches x 121 mammals. **230 hits at q<0.05 — but a Brownian-motion null on the same tree gives a median of 194.** Empirical p=0.395 |
| `t13_compartment` | is organ rate just secreted-vs-intracellular mix? | **hypothesis falsified**: interface proteins are 1.26× faster (p=4.6e-60) and organs differ 0.22–0.80 in composition, yet compartment absorbs −0.2% of the organ term |

Full write-up with every number, rigour label and caveat:
**[`REPORT_tissue_evolution.md`](REPORT_tissue_evolution.md)**.

---

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
