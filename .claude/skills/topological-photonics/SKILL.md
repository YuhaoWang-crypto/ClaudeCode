---
name: topological-photonics
description: >-
  Compute geometry-induced topological photonic/phononic band structures and
  their invariants on the Wu-Hu C6 lattice (and similar), end-to-end in
  open-source NumPy/SciPy: plane-wave-expansion bands, Dirac-cone and
  bandgap-opening detection, C6 band-inversion diagnosis, Berry curvature and
  spin-Chern number, real-rod FDFD edge states at a domain wall, 2D FDTD wave
  propagation (guiding / sharp bend / defect bypass), a scalar-acoustic phononic
  twin, the real-space spin Bott index under disorder (incl. topological Anderson
  insulator), and an absorption loss budget. Use when asked to model a topological
  photonic crystal / phononic lattice, judge whether a Dirac gap is topologically
  non-trivial, compute Chern/spin-Chern/Z2/valley or a disorder-robust invariant,
  optimise lattice geometry, simulate edge-state transport and backscattering
  immunity, or estimate ultralow-loss propagation. Enforces ✅-rigorous vs
  ⚠️-effective-model / ⚠️-illustrative labeling on every claim.
---

# Topological photonics & phononics on a Wu-Hu lattice

A reusable methodology (and a working `topological_photonics/` package) that
computes geometry-induced topology for classical waves — light and sound — and
always labels what is rigorous (ab-initio Maxwell/acoustic) vs. effective-model
or illustrative. No proprietary software (no Lumerical/COMSOL/MEEP); pure
NumPy/SciPy/Matplotlib.

## The core question this answers

Given a lattice geometry, decide whether it is topologically non-trivial for
photons or phonons, compute the right invariant **for its symmetry class**, and
show the physical consequence (edge transport, backscattering immunity) — without
overclaiming. A Dirac cone opening a gap is **not** proof of topology; the
deliverable is: invariant + interface states + edge dispersion + transport +
(if claimed) a loss budget.

## The model system (Wu-Hu photonic/phononic QSH)

Triangular lattice, hexagonal cluster of 6 rods per cell; one knob `R/a`:
`R=a/3` → double Dirac cone at Γ; `R<a/3` trivial gap; `R>a/3` band-inverted
(topological) gap where the p(l=±1) and d(l=±2) C6 doublets swap. Time-reversal
symmetric ⇒ global Chern = 0; the invariant is the **spin-Chern / spin-Bott**.

## The engines + when each applies

| Engine | Question it answers | Module | Rigor |
|---|---|---|---|
| Plane-wave expansion (TM) | photonic bands; Dirac cone; gap opening | `pwe.py`, `bands.py` | ✅ ab-initio Maxwell |
| C6 angular-momentum indicator | is the gap band-inverted (topological)? | `symmetry.py` | ✅ from real modes |
| Berry curvature + FHS spin-Chern | the invariant number (Chern total 0, spin ±1) | `bhz.py`, `topology.py` | ⚠️ effective BHZ (sign from PWE) |
| Geometry scan | phase boundary; max protected gap | `param_scan.py` | ✅ PWE |
| FDFD supercell | real-rod edge states at a topological\|trivial wall | `fdfd.py`, `fdfd_edge.py` | ✅ ab-initio Maxwell |
| 2D FDTD | guiding / sharp bend / defect bypass; energy delivery | `fdtd.py`, `fdtd_transport.py` | ✅ transport (qualitative) |
| Acoustic PWE | phononic twin; two-system comparison | `acoustic.py`, `phononic_compare.py` | ✅ ab-initio (scalar) |
| Real-space spin Bott index | disorder-robust invariant (no BZ); TAI | `bott.py`, `bott_disorder.py` | ✅ rigorous |
| Confinement-factor loss | absorption floor; L_p, dB/cm | `loss_budget.py` | ✅ method; ⚠️ illustrative ε″ |

## Choose the RIGHT invariant for the symmetry class

The single most important operational rule — do not hard-code "Chern number".

- **Broken time-reversal (gyromagnetic / Haldane):** genuine **Chern number**,
  chiral edge states. `bhz.chern_fhs` on one block.
- **Time-reversal symmetric, sz-conserving pseudospin (Wu-Hu QSH):** global
  Chern = 0; use **spin-Chern** (`bhz.spin_chern`) in the clean limit and the
  **spin Bott index** (`bott.spin_bott`) under disorder.
- **Valley-Hall (broken inversion, T-symmetric):** global Chern = 0; use the
  **valley-Chern** (integrate Berry curvature over one valley), not total Chern.
- **Z2 class:** compute via Wilson-loop / Wannier-centre flow (extension).

## Run the pipeline

```bash
pip install -r topological_photonics/requirements.txt   # numpy scipy matplotlib
cd topological_photonics/src
python run_all.py          # basic package  (~90 s)  -> fig1..fig4
python run_all.py --full   # + complete pkg (~4 min) -> fig5..fig10 + GIF
# or any single stage:
python bands.py | topology.py | fdfd_edge.py | fdtd_transport.py \
     | phononic_compare.py | bott_disorder.py | loss_budget.py
```

Every stage writes a figure to `figures/` and prints its computed numbers.
Build the standalone HTML report: `python topological_photonics/build_report.py`.
See `README.md` (deliverable tables) and `REPORT.md` (full write-up incl. §0b
scientific scope and the honesty ledger).

## The non-negotiable discipline: honesty labeling

Every result carries one of:

- **✅ rigorous / ab-initio** — solved directly from Maxwell (PWE/FDFD/FDTD) or
  the acoustic wave equation, or a symmetry indicator read from the real modes,
  or a real-space invariant (Bott index).
- **⚠️ effective model** — computed on the BHZ/tight-binding model the real band
  inversion maps onto (the spin-Chern *number*, the basic-package ribbon edges).
- **⚠️ illustrative** — depends on representative material constants pending
  measured data (the loss-budget ε″ values).

Never blur them. State the topological *classification* (band inversion, rigorous)
separately from the invariant *number* (effective, until Wilson loops are added).
When "ultralow-loss" is claimed, always attach the loss budget: topology removes
backscattering, **not** absorption/radiation.

## Common tasks → where to look

- **New lattice geometry** → give `pwe.PhotonicCrystal2D` a new rod list + lattice;
  reuse `bands.gamma_gap`, `symmetry.diagnose`. The double Dirac cone is at
  different band indices per system (see gotchas).
- **Phononic version** → `acoustic.AcousticCrystal` (set ρ,B contrast). Same
  `symmetry.c6_angular_momentum` works (scalar field).
- **Disorder robustness** → `bott.sweep_disorder` (real space). Do NOT reuse
  k-space Berry curvature once periodicity is broken.
- **Transport / bends / defects** → `fdtd.build_domain(interface=...)` with a
  path-based topological/trivial region; `drop_rod` for a defect.
- **Loss / Q / propagation length** → `loss_budget.edge_mode_confinement` gives Γ
  from the real FDFD mode; feed `loss_budget.modal_loss` with measured ε″.

## Hard-won gotchas (these were real bugs)

- **Double Dirac cone lives on different bands per system.** Photonic (TM Si
  rods): the p/d doublets are bands 1–2 / 3–4. Acoustic: there is an extra low
  acoustic branch, so they are bands 3–4 / 5–6 — reading bands 2–3 gives a
  meaningless gap. Always locate the four-fold degeneracy at R=a/3 first.
- **Spin partner must flip ONLY σx (or only σy), not both.** Flipping both dx,dy
  is a rotation and does NOT change the Chern sign → spin-Bott came out 0. Flip a
  single component (a mirror) so C↓ = −C↑.
- **FDFD needs sparse shift-invert** (`eigsh(L, M=Eps, sigma=…)`) targeted at the
  gap frequency; dense diagonalisation of the supercell is infeasible.
- **Crude FDTD transmission is only qualitative.** A soft point source + graded
  sponge (not UPML) gives Fabry-Perot noise; single-point monitors sit on
  standing-wave nodes. Use line-integrated flux or steady-state energy delivery,
  and report the bend qualitatively (snapshots/GIF) + defect/straight/trivial as
  the quantitative bars.
- **Don't fit L_p from the lossy FDTD envelope** — under strong absorption the
  wave barely enters the lossy region and the fit is non-monotonic. Use the
  perturbative confinement-factor formula with Γ from the real FDFD mode.
- **The topological side is not always R>a/3.** For the photonic crystal it is
  R>a/3; for the scalar-acoustic twin (these materials/bands) it is R<a/3.
  Read it from the C6 indicator, don't assume.
- **`matplotlib.cm.get_cmap` was removed** — use `matplotlib.colormaps["name"]`.
- **`np.cross` on 2-vectors is deprecated** — fine for a scalar z-component but
  silence/upgrade if it matters.
