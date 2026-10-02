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

## `bioif/` — composing heterogeneous models into longer causal chains

A second, independent piece of the repo. `grn_pipeline/` goes *deep* on one
mechanistic layer; `bioif/` addresses the orthogonal problem: models at
different layers (small molecule, protein, transcript, variant, epigenome)
are each fitted on their own data for their own question, and do not compose
into a longer causal chain.

The thesis is that they fail at the **semantics of the seam**, not the file
format — so the fix is a thin typed contract (entity + quantity + context +
estimate + domain + provenance), not a universal model.

```bash
python3 -m bioif.demo             # six scenarios, stub adapters, no dependencies
python3 -m bioif.demo_real        # same chain, REAL ChEMBL affinity source
python3 -m bioif.demo_conformal   # a calibrated interval, with its coverage checked
python3 -m bioif.demo_routing     # 12 models for one hop; the registry chooses
python3 -m bioif.demo_pairings    # 3 more interface-map pairings, on real labels
python3 -m bioif.map21            # all 21 pairings audited against what is buildable
python3 -m bioif.datasheet        # per-pairing inputs/outputs/throughput + chain verdict
python3 -m bioif.real.heterogeneity  # the measurement that justifies the contract
python3 -m bioif.selftest         # 53 guarantees, ~70 s
```

Replacing the stub affinity source with measured ChEMBL bioactivity kept the
suite green — but only after fixing three defects in the contract that stub
data structurally could not reveal (readout-type pooling, single-sample
estimates collapsing the Monte Carlo, and an identity resolver that always
answered). Details and two real measured numbers in [`INTEROP.md`](INTEROP.md) §5.

One adapter's error bar is then **checked rather than asserted**: split
conformal prediction on the assay-transfer hop (may a potency measured in
assay A be reused where assay B is needed?). Within a pair it reaches nominal
coverage where the usual Gaussian ±1.28·sd under-covers (0.830 vs 0.671 at a
nominal 0.80); across assay pairs the guarantee collapses to 0.612 with a
worst case of 0.002, which is why the adapter refuses uncalibrated pairs.
[`INTEROP.md`](INTEROP.md) §6.

Finally the same hop is given **twelve competing adapters** and the registry
picks per claim, on declared properties (in-domain, then narrowest calibrated
interval). Running the alternatives as well as the winner prices the
uncertainty a conformal interval cannot see — the choice of model itself. It
also caught a bug in this repo that had inflated that figure ~6×.
[`INTEROP.md`](INTEROP.md) §7.

Then the whole 21-pairing interface map is audited against what can actually
be built (4 built, 3 partial, 9 blocked, 5 correctly refused), three more
pairings are built on measured public labels (Tox21 SR-p53, Hansen Ames, and
the p53→mutagenicity edge a blueprint had marked *blocked* — RR 1.81,
p=5e-9, sensitivity 0.13), and the long-chain question is settled on the
1,908 compounds that have labels at **both** ends: a direct one-hop model
beats the two-hop chain by +0.223 AP, and a *measured* intermediate adds
−0.000 on top of structure. [`INTEROP.md`](INTEROP.md) §8.

| | |
|---|---|
| Bottleneck taxonomy (12 cross-cutting axes + 19 layer crossings), ranked plan, and what the demo does *not* show | [`INTEROP.md`](INTEROP.md) |
| The contract | `bioif/core.py` |
| Adapters, domain verdicts, refusal | `bioif/adapter.py` |
| Automatic routing between quantities | `bioif/registry.py` |
| Distribution propagation + variance attribution → experiment ranking | `bioif/chain.py` |
| Real ChEMBL source, identity resolution that refuses, committed snapshot | `bioif/real/` |
| Conformal calibration, competing models, the assay-transfer adapter | `bioif/real/conformal.py`, `transfer_models.py`, `transfer_adapter.py` |
| Running every route and pricing the model choice | `bioif/ensemble.py` |
| The 21-pairing audit | `bioif/map21.py` |
| Deployment datasheet: inputs, outputs, throughput, chain verdict | `bioif/datasheet.py` |
| Tox21 / Ames pairings, QSAR + classification conformal | `bioif/real/tox.py`, `qsar.py`, `tox_adapters.py` |
| Chain vs direct, with labels at both ends | `bioif/real/chain_vs_direct.py` |
| The cheap joint model (copula) behind the conditional C→F edge | `bioif/real/copula.py` |

Dependencies: the contract layer (`core`, `adapter`, `registry`, `chain`,
`ensemble`) is **stdlib-only**. The real-data pairings need `rdkit`, `numpy`
and `scikit-learn`; without them those tests report SKIP rather than fail.

⚠️ Every numeric constant in `bioif/adapters_demo.py` is an illustrative
stub. The demo exercises interface behaviour — typed seams, automatic
routing, honest widening, refusal, variance attribution — and makes no
prediction about any real gene, compound or cell line.
