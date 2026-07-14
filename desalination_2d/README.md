# 2D-Material Desalination: DFT → MLP → NEMD pipeline

An open-source, reproducible scaffold for studying reverse-osmosis (RO)
desalination through **nanoporous 2D membranes** — **graphene first, then
MoS₂** — using the standard multiscale route:

```
  Stage 1            Stage 2                 Stage 3
  DFT / AIMD   ──▶   ML potential (MLP)  ──▶  NEMD (ns, pressure-driven)
  (QE, CP2K)         (DeePMD + DP-GEN)        (LAMMPS + DeePMD)
  gold-standard      DFT accuracy,            water flux &
  energies/barriers  MD speed                 salt rejection
```

Everything here runs on **fully open-source engines** (Quantum ESPRESSO /
CP2K / DeePMD-kit / LAMMPS / ASE / DP-GEN). No commercial code is required.

---

## What actually runs today (validated in this repo)

The geometry building and the **entire analysis + validation pipeline** run
with nothing more than `numpy + ase + matplotlib`:

```bash
pip install -r requirements.txt

# Stage 1 geometry
python 01_build/build_graphene_pore.py --pore-radius 3.5 --out figures/graphene_pore
python 01_build/build_mos2_pore.py     --pore-radius 3.5 --out figures/mos2_pore
python 01_build/build_saltwater_box.py --membrane figures/graphene_pore.xyz \
       --n-water-feed 300 --n-water-perm 300 --n-nacl 10 --out figures/nemd_cell

# Stage 3 analysis self-check
python 04_nemd/analyze_flux_rejection.py --demo

# End-to-end TOY validation (physics trends + analyzer correctness) + figure
python 05_toy_validation/toy_nemd.py --plot figures/toy_nemd.png
```

The toy validation asserts five things and **all pass**:
1. the analyzer's crossing count matches an independent ground truth (0 % error);
2. water throughput increases monotonically with pressure;
3. salt rejection increases monotonically with the ion dehydration barrier;
4. the membrane is selective (water permeates, ions are ~97 % rejected);
5. molecular flux → experimental permeance conversion is finite/sane.

![toy validation](figures/toy_nemd.png)

> The toy is a kinetic pore-gate model — it validates the *workflow and the
> analysis code*, **not** a quantitative prediction. Real numbers require
> Stages 1–3 on real compute (below).

---

## What needs real compute (and where to run it)

| Stage | Workload | Bound by | Cloud (Modal) | HPC (Slurm) |
|---|---|---|---|---|
| DFT / AIMD | Quantum ESPRESSO | **CPU** (plane-wave/MPI) | `02_dft/modal_run_dft.py` | `slurm/qe_scf.slurm` |
| MLP training | DeePMD-kit | **GPU** | `03_mlp/modal_train_deepmd.py` | `slurm/deepmd_train.slurm` |
| NEMD production | LAMMPS+DeePMD | **GPU** | (mirror the training launcher) | `slurm/lammps_nemd.slurm` |

### Can I run this on modal.com (GPU cloud)?
**Yes — and it's a good fit.** Two ready templates are included:
- `02_dft/modal_run_dft.py` runs QE on a **CPU** box. Plane-wave DFT is
  CPU/MPI-bound, so a GPU would be wasted here — Modal lets you rent a fat
  16-core box per-second instead.
- `03_mlp/modal_train_deepmd.py` trains the MLP on a **GPU** (A10G/A100), which
  *is* GPU-bound — the right place to spend GPU dollars, along with the LAMMPS
  NEMD inference.

Both are runnable once you `pip install modal && modal token new`. They define
the whole software environment in code (QE / DeePMD images), so runs are
reproducible with no cluster module juggling.

### Is there a ready-made "materials-science / DFT skill" for this?
Not in this environment. The available skills are Claude-Code / bio-oriented,
and the connected data servers (PubMed, ChEMBL, bioRxiv, …) are life-sciences,
**not** materials tools. So there's no push-button DFT skill — but the pipeline
is fully covered by the open-source engines above, and this repo is the
scaffold that wires them together. Literature search for the *materials* side is
best done on arXiv `cond-mat` / Materials Cloud rather than the bio servers here.

---

## MLP vs ReaxFF — recommendation

For **physical** (non-reactive) RO desalination — water and ions passing a pore
without breaking bonds — **use the MLP.** It captures polarization and the
complex charge distribution in the confined pore at near-classical cost.

**ReaxFF is a "sledgehammer" here** (10–100× cost, painful fitting) and is only
worth it if you specifically study *chemistry*: pore-edge protonation/
deprotonation, membrane chemical degradation under strong flow, or fouling.
Recommendation: ship the MLP path first; treat ReaxFF as an optional Stage-2b
only if a chemical question demands it. Don't develop both in parallel.

---

## Repository layout

```
desalination_2d/
├── 01_build/          # ASE geometry builders (graphene, MoS2, saltwater cell)
├── 02_dft/            # QE SCF + NEB templates, Modal CPU runner   (Stage 1)
├── 03_mlp/            # DeePMD + DP-GEN configs, Modal GPU trainer  (Stage 2)
├── 04_nemd/           # LAMMPS NEMD input + flux/rejection analyzer (Stage 3)
├── 05_toy_validation/ # numpy-only end-to-end pipeline validation
├── slurm/             # HPC batch scripts (CPU DFT, GPU train, GPU NEMD)
├── figures/           # generated geometries + validation figure
└── requirements.txt
```

## Suggested first milestone
1. `build_graphene_pore.py` → a converged QE SCF on Modal (Stage 1 sanity).
2. One CI-NEB Na⁺ barrier → your first physical selectivity number.
3. A small DeePMD model on ~a few thousand DFT frames + DP-GEN loop.
4. A 1 ns LAMMPS NEMD at ~100 MPa → first flux/rejection point.
5. Scale up pressure/pore-size scans; then repeat the whole path for MoS₂.
