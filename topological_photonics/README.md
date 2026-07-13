# Topological Photonics — Basic Computation Package

A small, fully-runnable package for the **computational** half of the proposal
*"Fabrication and Study of Topological Insulators with Silica-Coated Bi₂Se₃ for
Ultralow-loss Phononic Waveguides and Photonic Chips."*

It delivers every item of the **基础计算包 (basic package)** on one coherent,
literature-grounded system — the **Wu-Hu photonic topological insulator**
(Wu & Hu, *PRL* **114**, 223901, 2015), the photonic/phononic analogue of the
quantum spin Hall effect that underlies references [2]–[4] of the proposal.
Everything is computed with open-source NumPy/SciPy — **no Lumerical/COMSOL**.

| # | Deliverable (from the proposal) | Module | Figure | Status |
|---|---|---|---|---|
| 1 | A definite lattice; photonic band structure | `bands.py` | `fig1` | ✅ rigorous (PWE) |
| 2 | Identify Dirac cone; confirm bandgap opening after symmetry breaking | `bands.py` | `fig1` | ✅ rigorous (PWE) |
| 3 | Berry curvature / Chern / **spin-Chern** post-processing | `topology.py` | `fig2` | ✅ / ⚠️ (see below) |
| 4 | Basic geometry parameter scan / optimisation | `param_scan.py` | `fig4` | ✅ rigorous (PWE) |
| 5 | One set of edge-state / wave-propagation simulations | `edge_states.py` | `fig3` | ⚠️ effective model |
| 6 | Basic result figures + explanation | this README + `REPORT.md` | all | ✅ |

## The physics in one paragraph

A triangular lattice with a hexagonal cluster of 6 silicon rods (ε=11.7) has a
single geometric knob: the cluster radius `R/a`. At `R = a/3` the rods form a
perfect honeycomb and the two graphene Dirac cones fold onto Γ into a **double
Dirac cone**. Shrinking (`R<a/3`) opens a **trivial** gap; expanding (`R>a/3`)
opens a **topological** gap in which the p-like (l=±1) and d-like (l=±2) modes
**invert**. The inverted crystal carries a **spin-Chern number ±1** and, at an
interface with a trivial crystal, a **Kramers pair of helical edge states** that
transport light without backscattering — the "ultralow-loss, backscatter-immune"
channel the proposal targets.

## Rigorous vs. effective-model — honest labelling

* **✅ rigorous (ab-initio photonics):** the band structures, Dirac cone, gap
  opening, geometry scan, and the **band-inversion diagnosis** (C6 angular-
  momentum of the real PWE modes at Γ) are computed directly from Maxwell's
  equations for the actual rod crystal.
* **⚠️ effective model:** the **spin-Chern number** and the **helical edge-state
  spectrum + disorder robustness** are computed on the 2-copy BHZ Hamiltonian
  that the C6 band inversion maps onto, with its mass sign fixed by the PWE
  result. This is the standard description of the photonic QSH crystal; the
  *number* (Cs=±1) and the *helical crossing* are model-level, the *band
  inversion that forces them* is rigorous. Promoting the edge calculation to a
  full PWE supercell of the real rods is a **complete-package** item.

## Run it

```bash
pip install -r requirements.txt
cd src && python run_all.py      # ~90 s, writes all four figures to ../figures/
```

Individual stages: `python bands.py`, `python topology.py`,
`python edge_states.py`, `python param_scan.py`.

## Key computed results

* Double Dirac cone at Γ at `R=a/3` (four-fold degeneracy, verified).
* Trivial gap `Δ≈0.048` and topological gap `Δ≈0.046` (ωa/2πc) at `R=0.30a` / `0.36a`.
* Band inversion: lower doublet |l| goes 1→2 across the transition (C6 indicator).
* Spin-Chern: `Cs=±1` (topological) vs `0` (trivial); total Chern `=0` (T-symmetric).
* Helical edge states cross the gap only in the topological ribbon; the mid-gap
  crossing survives in **20/20** random disorder realisations (W=0.6).
* Optimal operating point in the scanned window: `R/a≈0.385`, gap `≈0.099`.

See `REPORT.md` for the full write-up, method details, and caveats.
