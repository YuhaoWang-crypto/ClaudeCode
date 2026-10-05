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

# fftplsr

A second, independent pipeline in this repo: **combinatorial enzyme design with
FFT-PLSR**. Given measured activities for single mutants of one enzyme, rank the
2^k recombinations and emit an order list — plus the baselines and numerical
guardrails that decide whether the ranking deserves trust.

Method source: Hu *et al.*, "Machine learning-guided evolution of pyrrolysyl-tRNA
synthetase for improved incorporation efficiency of diverse noncanonical amino
acids", *Nat Commun* **16** (2025),
doi:[10.1038/s41467-025-61952-2](https://doi.org/10.1038/s41467-025-61952-2);
reference code [zjuhaoran/FPFORCOM](https://github.com/zjuhaoran/FPFORCOM) (MIT).
The paper's data is vendored under `fftplsr/data/` — see
[`PROVENANCE.md`](fftplsr/data/PROVENANCE.md).

| Module | What it establishes | Key result |
|---|---|---|
| `m1_reproduce` | re-runs the paper's prospective rounds | round-2 held-out R² = **0.833** (paper: 0.835); round-3 descriptor triple reproduces **exactly**; round-1 does **not** reproduce |
| `m2_dc_artifact` | why FFT bin 0 must be dropped | bin 0 is round-off, but `scale=True` amplifies it to O(1); cvMSE at k=10 swings **0.29–8.05** across equally valid round-off draws, vs **7.7e-14** with it dropped |
| `m3_baselines` | is the FFT worth it? | a plain one-hot PLS **beats** FFT-PLSR on the 38→64 split (R² 0.827 vs 0.750) and picks better variants on both splits |
| `m4_design` | a worked design round | ranks a ~1.5k-variant space from measured singles; reports where the paper's final winner lands |

```bash
pip install numpy scipy pandas scikit-learn aaindex joblib
python3 -m fftplsr.m1_reproduce     # or m2_dc_artifact / m3_baselines / m4_design
python3 -m pytest tests/test_fftplsr.py -q
```

Design on your own enzyme:

```python
from fftplsr import design
report = design.design_round(parent=seq, measured={"WT": 1.0, "D2N": 3.62, ...}, pick=8)
print(report.summary())    # headline + baselines + pick-list + caveats
```

Three findings worth carrying away before using this method anywhere:

1. **Drop FFT bin 0.** Mean-centering zeroes the DC term, so bin 0 is pure
   floating-point round-off — and standardizing it turns that round-off into a
   fitted predictor. This is what makes the published round-1 descriptor choice
   irreproducible across numerical stacks.
2. **Report out-of-fold scores.** The reference implementation's `R2` column is
   in-sample; at 13 samples and 227 features it reads 0.999 and means nothing.
3. **Beat the additive model or don't bother.** Only a dozen positions vary, so
   the spectrum is a deterministic function of a dozen bits that a one-hot
   regression handles directly — and usually better once you have ~40 measurements.

The methodology, traps and composition rules are packaged as the
[`enzyme-combinatorial-design`](.claude/skills/enzyme-combinatorial-design/SKILL.md)
skill.
