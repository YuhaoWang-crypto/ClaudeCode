# ClaudeCode — two runnable, honesty-labelled pipelines

| directory | topic | report |
|---|---|---|
| `electrolyte_pipeline/` | battery-electrolyte computation (MD / QC / interface / ML), reproducing the methods digest of Yao et al., *Chem. Rev.* 2022, 122, 10970 | [`REPORT_ELECTROLYTE.md`](REPORT_ELECTROLYTE.md) |
| `grn_pipeline/` | gene-regulatory-network irreducibility / tipping-point biomarkers | [`REPORT.md`](REPORT.md) |

## electrolyte_pipeline

One module per page of the digest; every number carries ✅ converged / ⚠️ demo / ❌ not-run.

| Module | Digest page | What it computes | Where it runs |
|---|---|---|---|
| `e0_systems` | p4 | compositions (counts, not "1 M"), force-field & dynamics spec, experimental reference table | – |
| `e1_build` / `e1_run` | p4 | packmol → OpenFF Sage 2.2.1 + NAGL charges (ions ×0.8) → OpenMM; minimise / NPT / production with block-convergence check | Modal A10G (`modal_run.py`) |
| `e2_rdf_cn` | p5 | Li–O(DME)/O(TFSI)/N/F RDF, r_min, density-weighted N(r), direct-count CN, P(n) at atom & molecule level, representative shell PDB | local |
| `e3_clusters` | p6 | contact graph → no-contact / CIP / AGG fractions, cluster sizes, anion bridging; criteria written next to the numbers | local |
| `e4_properties` / `e4_nemd_viscosity` | p8 | density, Einstein D with log-log slope check, Nernst–Einstein vs Einstein–Helfand conductivity, solvent dielectric, periodic-perturbation NEMD viscosity, all vs experiment | local / Modal |
| `e5_qc_clusters` | p2 | Li⁺–EC / Li⁺–DME from several starting placements, B3LYP/def2-TZVP//def2-SVP, vertical / CP-corrected / relaxed ΔE, MD-shell cluster | local CPU (PySCF) |
| `e6_desolvation` | p11 | same Li…O coordinate as (2) gas-phase electronic scan and (3) liquid PMF from g(r) — shown side by side, not interchanged | local |
| `e7_interface` | p9 | rigid uncharged graphite(0001) + LiTFSI/DME, z-resolved number densities from the top C plane, DME orientation, Li PMF along z | Modal |
| `e8_ml` | p12 | (A) GFN2-xTB vs DFT forces/energies on Li⁺ shells cut from the liquid; (B) property regressor under random vs leave-one-concentration-out split | local |
| `e9_reactive` | p10 | why reactive MD is not run; bond-topology event counter returns 0 on a fixed-topology FF | local |

```bash
# environment (conda-forge): openmm openff-toolkit openff-interchange openff-nagl openff-nagl-models packmol
#                            rdkit mdanalysis pyscf geometric xtb-python scikit-learn networkx  + pip 'modal[api-proxy-support]'
modal run electrolyte_pipeline/modal_run.py                    # E1 bulk series (6 systems, 2+10 ns)
modal run electrolyte_pipeline/modal_run.py --stage nemd       # E4 viscosity
modal run electrolyte_pipeline/modal_run.py --stage interface  # E7
modal volume get elyte-work / elyte_work/
python -m electrolyte_pipeline.run_all --qc                    # all analyses + CPU quantum chemistry
```

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
