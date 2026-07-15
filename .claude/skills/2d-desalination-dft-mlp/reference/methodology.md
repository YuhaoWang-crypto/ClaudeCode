# Methodology reference — 2D-material desalination DFT → MLP → NEMD

## 1. Why multiscale
DFT is the gold standard for interaction energies, charge transfer, and
pore-crossing barriers, but it reaches only ps of AIMD. Water flux and salt
rejection are ns-scale non-equilibrium phenomena. The bridge is a
**machine-learning potential (MLP)** trained on DFT: DFT accuracy at MD speed.

```
DFT/AIMD (ps, QE/CP2K) → MLP (DeePMD/MACE/NequIP + DP-GEN) → NEMD (ns, LAMMPS)
```

## 2. Stage 1 — DFT
- **Functional:** water/ion physisorption needs dispersion — use **vdW-DF2** or
  PBE+D3. Plain PBE underbinds water and gives wrong barriers.
- **Charged ions:** set `tot_charge = +1` (Na⁺) / `-1` (Cl⁻) in a compensating
  background; converge box size to limit spurious electrostatics.
- **Convergence:** always converge `ecutwfc` and the k-mesh (a porous supercell
  is large in-plane; 2×2×1–3×3×1 often suffices). Tool:
  `make_convergence_tests.py`.
- **Barriers:** CI-NEB (`neb.x`) gives the saddle. The **forward barrier**
  `E_a = E_saddle − E_feed` is the selectivity number. Compare E_a(water) vs
  E_a(hydrated ion).
- **CPU-bound:** plane-wave DFT is MPI/CPU work — never a GPU job.

## 3. The dehydration rule (most important)
Real RO rejection is dominated by the **free-energy cost of stripping the ion's
hydration shell** to fit a sub-nm pore, not by bare-ion sterics. Consequences:
- Model **hydrated ions** (ion + first shell; Na⁺·(H₂O)₆, Cl⁻·(H₂O)₆).
- A bare-ion NEB/scan will *underestimate* rejection dramatically (bare ions can
  show ~0 barrier — see `06_dft_surrogate_study/REPORT.md` §3).
- For a quantitative free energy, follow the DFT barrier with **PMF / umbrella
  sampling / metadynamics** in MD using the trained MLP.

## 4. Stage 2 — MLP
- **DeePMD-kit** (DeepPot-SE) is the aqueous-interface workhorse; **MACE** and
  **NequIP/Allegro** are equivariant and usually more data-efficient. Train all
  three on identical data and pick by force RMSE + NEB-barrier reproduction.
- **Data formats:** DeePMD reads `npy` systems; MACE/NequIP read extended-XYZ
  with `REF_energy`/`REF_forces`. Convert with
  `prepare_deepmd_data.py --to-extxyz`.
- **Cover the right states:** include confined, high-pressure, **hydrated-ion**
  frames — not just bulk water — or the MLP extrapolates and explodes under the
  NEMD piston. **DP-GEN** active learning (`dpgen_param.json`) finds
  high-uncertainty frames to label, minimizing DFT cost.
- **GPU-bound:** training and LAMMPS+MLP inference belong on a GPU.

## 5. Stage 3 — NEMD
- Rigid (or stiffly tethered) membrane; a feed **piston** applies the RO
  pressure. Compute the per-atom piston force from `P = F·N/A` with
  `piston_force.py` (metal units: eV/Å).
- **Observables:** water flux/permeance and rejection `R = 1 − c_perm/c_feed`.
  `analyze_flux_rejection.py` counts net feed→permeate crossings with a
  boolean-side test (no double-count on the membrane plane) and converts flux to
  L·cm⁻²·day⁻¹·MPa⁻¹.
- Validate the analysis code with `05_toy_validation/toy_nemd.py` before trusting
  a production run.

## 6. Surrogate vs converged DFT (honesty boundary)
A universal MLIP (e.g. CHGNet, MACE-MP) is a legitimate way to **run the whole
workflow locally and get order-of-magnitude numbers and correct trends**
(pore-size scaling, water-vs-hydrated-ion ordering, material contrast). It is
**not** a substitute for the user's own DFT for:
- absolute barrier values,
- charged-ion / dehydration energetics,
- vdW-corrected water structure.
Always tag surrogate results ⚠️ and converged DFT ✅.

## 7. MLP vs ReaxFF
| | MLP (recommended) | ReaxFF (special cases) |
|---|---|---|
| Cost | ~classical MD | 10–100× |
| Fitting | DFT + DP-GEN, automatable | manual, painful |
| Captures | polarization, confined charge, dehydration | bond breaking/forming |
| Use for | physical RO transport (default) | pore-edge (de)protonation, degradation, fouling |

## 8. Common pitfalls
- Passivate only true **pore-edge** atoms (use minimum-image coordination; the
  builders already do this) — not periodic-boundary atoms.
- Do not measure flux as a polyfit slope on a **saturating** (finite-feed) curve;
  in steady state the cumulative is linear and the slope is valid.
- Keep `r_max`/`rcut` identical across DeePMD/MACE/NequIP for a fair bake-off.
- Recompute the piston force whenever the box or piston-atom count changes.
