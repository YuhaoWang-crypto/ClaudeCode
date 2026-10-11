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

## msm_pipeline — 构象通量 (MSM / TPT) on real protein MD

A second, independent pipeline in this repo: the same "find the irreducible
core and the measurable signal" philosophy, applied to **protein conformational
dynamics** instead of regulatory networks.

**Question:** can MSM + transition-path-theory find an allosteric site that
RMSD/RMSF analysis of the *same trajectories* would miss?

**System:** TEM β-lactamase, 48 × 500 ns (24 μs) from
[Zenodo 6565237](https://zenodo.org/records/6565237), the data behind
Avery, Baker & Jacobs, *Entropy* 2022 ([doi:10.3390/e24050729](https://doi.org/10.3390/e24050729)).

| Module | What it does | Rigor |
|---|---|---|
| `m1_dataset` | provenance + sampling budget | ✅ |
| `m2_cryptic_pocket` | opening axis fixed by **crystals** (1BTL→1PZO), not by the trajectories | ✅ |
| `m3_msm_diagnostics` | implied-timescale test — **it fails**, and why (trajectory-identity variance 0.6–0.86) | ✅ negative result |
| `m4_tpt_flux` | core-set MSM + TPT: committor, flux, bottleneck | ⚠️ topology only |
| `m5_variant_allostery` | **main result:** TEM-10 (R164S/E240K) opens a pocket 30 Å away **282×** more than WT, paired Wilcoxon p=0.0156 | ✅ |
| `m6_traditional_contrast` | RMSD/RMSF/PCA shift 1.1× on the same data | ✅ |
| `m7_druggability` | open-state export, docking template, chemical starting point | ✅ / ⚠️ |
| `m8_positive_control` | alanine dipeptide — proves the ITS test **does** pass when it should (1.12×) | ✅ |
| `m9_adaptive_benchmark` | controlled test with known ground truth: adaptive reseeding cuts median timescale error **29.1 % → 8.5 %** at matched cost; but the seeding heuristics are **statistically indistinguishable** (0/6 comparisons p<0.05) | ✅ |
| `m10_tem_projection` | measured throughput on this box (5.3 ns/day) → honest TEM cost: 8 GPUs ≈ 3–8 days, this container ≈ 8 years | ✅ |
| `adaptive_md.py` | production round-based OpenMM driver, checkpointed and resumable; validated end-to-end on real MD | ✅ |

**Headline:** a distal mutation pair shifts the population of a known cryptic
allosteric pocket by two orders of magnitude, giving a falsifiable prediction
(non-competitive 1PZO-chemotype inhibitors should be more potent against
TEM-10 than TEM-1). The conventional metrics see a 1.1× change on the
identical trajectories — significant, but not actionable. Full write-up with
the caveats in [`MSM_REPORT.md`](MSM_REPORT.md).

```bash
pip install deeptime mdtraj scipy matplotlib numpy
python -m msm_pipeline.run_all            # reproduce (downloads ~660 MB)

# ...or point it at your own trajectories, on any disk:
python -m msm_pipeline.run_gromacs \
    --top /Volumes/T7/run/md.gro --traj "/Volumes/T7/run/*/md.xtc" \
    --sel "name CA" --a-pair 220 276 --a-max 7.0 --b-min 9.5

# ...or run an adaptive-sampling campaign that reseeds from the barrier:
python -m msm_pipeline.adaptive_md \
    --pdb 1BTL.pdb --out /Volumes/T7/tem_adaptive \
    --rounds 60 --walkers 16 --ns-per-walker 20 \
    --sel "name CA and (resid 214 to 231 or resid 269 to 287)" \
    --strategy hybrid --platform CUDA
```

**On adaptive sampling** (M9): reseeding in rounds is what fixes the
connectivity failure M3 diagnosed — median error on the slowest timescale
drops from 29.1 % to 8.5 %, and "few long trajectories" failed to connect the
two basins in 3/24 repeats while adaptive never did. Targeting the committor
bottleneck specifically was **not** measurably better than simply reseeding
the least-visited states, so run rounds and don't over-engineer the rule.
