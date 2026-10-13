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
| `m23_msm_tpt` | conformational flux: MSM + Transition Path Theory (deeptime) | alanine dipeptide, 4×50 ns OpenMM (amber99sb-ildn/OBC, 300 K) | β→αR: MSM-TPT MFPT 52 ps vs 48 ps counted directly (2346 events, no Markov assumption); 94% of reactive flux crosses ψ≈+60…120, 6% wraps through ψ=±180; CK error 0.003; flux through q⁺=½ = 1.0000·F. β→αL only 4 events → flagged order-of-magnitude. αL vs αR overlap 31% in RMSD |
| `m24_msm_validation` | reliability of M23: exact ground truth + positive/negative controls + perturbation screening | triple-well SqRA jump process (exact rates/committor); alanine dipeptide with ψ-bias (positive) / methyl-bias (negative) | rate within 0.4% of exact (bootstrap CI covers it); committor r=0.999; shuffled / blind-coordinate negatives collapse the slow process; 56 virtual ligands ranked from ONE MSM with Spearman ρ=0.989 (static baseline 0.03), 100% correct sign, decoys ≈0; real MD: population shift predicted 0.73 vs 0.72, kinetics right direction but ~1.5× off; hyper-parameter spread −14…+11% |
| `m25_allosteric_retro` | retrospective test on proteins with known allosteric drugs | glucokinase (activator, 1V4T/1V4S) and SHP2 (inhibitor SHP099, 4DGP/6CRF; 5EHR only for scoring); dual-basin Cα Gō model + adaptive sampling | pre-registered site scoring: GK MSM ΔlnK AUC 0.75 vs static-structure 0.72 (top-10 hits 10 Å from the activator vs 18.5 Å random); SHP2 AUC 0.50 vs 0.68 with a weaker model (CK 0.17, failed control) — MSM scan does not beat a two-structure comparison in this CG setting |
| `m26_texh_jacobian` | order parameter -> Jacobian on real omics: does the "reduce, perturb, find fragile points" paradigm work? | GSE327731 in vivo Perturb-seq, 73-gene CRISPR KO in CD8 T cells, mouse melanoma (25,423 cells) | exhaustion order parameter reproducible across 3 definitions (rho 0.50/0.56, participation ratio 17); 30/42 KOs shift it at FDR<0.1 with a passing cutting-vs-non-cutting control; but MRA Jacobian recovery is at chance (sign consistency 51%) because CRISPR KO leaves target mRNA intact (diag/off-diag 2.2) — calibrated against synthetic networks; and eigenvalues are unidentifiable from steady state anyway (tau spans 154x under unknown turnover); CellOracle-style in-silico KO has no predictive power (rho -0.10, p=0.53) |
| `m27_crispri_mra` | does a perturbation that moves the target's own transcript rescue network inference? (follow-up to M26) | GSE190604 CRISPRa Perturb-seq, 73 targets in primary human T cells (61,041 cells, 71-gene closed module) | diagonal fully restored: 100% of targets self-activate, diag/off-diag 2.2 -> 113, condition number 6356 -> 53; MRA leaves chance (bootstrap 51% -> 79%) but split-half (independent cells) gives only 50% -> 59%, with 79% for the strongest edges — bootstrap overstates reliability by ~20 points. GSE249595 (the requested CRISPRi set) was unusable: its GEX and guide libraries use 10x's translated barcode whitelists (zero overlap verified) and pairing needs CellRanger's license-gated table |
| `m28_identifiability_tiers` | when is the "order parameter -> Jacobian -> fragile point" paradigm actually identifiable? | a candidate 3x3 GRN Jacobian (structural check) + Chassagnole E. coli central carbon metabolism (BIOMD0000000051, exact Jacobian) | structural: Thomas' rules read off a proposed matrix decide which claims are possible before any fitting (the example matrix has only positive loops → can be bistable, can NEVER oscillate; its saddle comes from the BLIMP1–XBP1 positive loop, not the mutual inhibition). Known kinetics: real timescales (τ 0.31 ms … 112 s, 4e5× span), slow mode localised on 4 of 18 metabolites (probe selection), and enzyme fragility uncorrelated with flux (ρ=−0.03) — the fragile point is genuinely hidden from flux data |
| `m29_toggle_case_study` | the full paradigm worked end-to-end on a real published switch | Gardner, Cantor & Collins 2000 genetic toggle (BIOMD0000000507, deposited parameters; my implementation reproduces the curated SBML's two attractors to 3 decimals) | bistable with an explicit saddle at u=1.32: the quiescent state sits 0.98 from the threshold but the committed state is 154 away (157× asymmetry = irreversible commitment); committed state relaxes 4.3× faster (τ 1.10 vs 4.74 h); a *transient* pulse gives a permanent switch with amplitude ~ duration^−0.73, so dose is NOT conserved (3.5× range) and short strong pulses are more efficient; 0 complex eigenvalue pairs in 666 fixed points over 400 damping combinations — mutual repression is a positive loop and cannot oscillate; threshold sensitivity is dominated by the Hill coefficient (d ln T/d ln β = +8.1), not by the abundances |

## Run

```bash
pip install numpy scipy networkx matplotlib
python3 -m grn_pipeline.run_all       # full pipeline + figures
python3 -m grn_pipeline.m1_symmetry   # or any single module
```

M23 (MSM + TPT) additionally needs `pip install openmm mdtraj deeptime`. With no
arguments it simulates/analyses the alanine-dipeptide demo (trajectories are
cached in `figures/_md_cache/`, ~15 min to regenerate on 4 cores). On your own
GROMACS/AMBER trajectories:

```bash
gmx trjconv -s md.tpr -f md.xtc -o md_whole.xtc -pbc mol -center   # whole molecules
gmx editconf -f md.tpr -o conf.gro                                  # mdtraj cannot read .tpr
python3 -m grn_pipeline.m23_msm_tpt --top conf.gro --traj md_whole.xtc \
    --dt-ps 10 --features both --n-macro 4 --source 0 --target 2 --out flux.png
```

`python3 -m grn_pipeline.m24_msm_validation` runs the reliability suite and
`python3 -m grn_pipeline.msm_html_report` builds `figures/msm_report.html`.

Run once without `--source/--target` to see the PCCA+ macro-states, then pick
A and B. TPT uses PCCA+ core sets (membership >= 0.9); a warning is printed
when fewer than 10 A->B events were sampled.

Figures are written to `figures/`. A full write-up with numbers, rigour
labels, and the interpretation (including the Lyapunov-exponent biomarker
question) is in [`REPORT.md`](REPORT.md).
