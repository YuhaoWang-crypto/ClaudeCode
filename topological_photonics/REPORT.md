# Basic Computation Package — Technical Report

**Project:** Topological insulators with silica-coated Bi₂Se₃ for ultralow-loss
phononic waveguides and photonic chips (computational half).
**System studied:** Wu-Hu photonic topological insulator (a C6-symmetric
honeycomb-cluster photonic crystal), the geometry-induced topological platform
that underlies the proposal's references [2] (programmable topological photonic
chip), [3] (photonic topological Anderson insulator), and [4] (GHz topological
phononic circuits — the same lattice physics for elastic waves).

---

## 0. Why this system (and a note on the proposal's framing)

The proposal's title centres on Bi₂Se₃ (an *electronic* topological insulator)
and silica coating. However, its **computational methodology** — photonic/
phononic band structure, Dirac cones, Berry curvature, (spin-)Chern numbers,
edge states — is the physics of **geometry-induced (classical) topology**, not
electronic band topology. That is also what all four cited papers do. Silica's
role is as a **low-loss host/substrate material**, not as an electronic TI
(pure silica is a trivial wide-gap insulator, as the proposal itself notes).

We therefore build the package on the canonical geometry-induced photonic TI.
The identical lattice, solved for elastic waves, gives the phononic waveguide of
ref [4]; swapping the solver kernel (Maxwell → elastodynamics) is a
complete-package extension.

---

## 0b. Scientific scope & honest caveats (materials / claims)

These bound what the package does and does **not** assert; they correct several
framings in the original proposal.

1. **Silica coating does not create topology.** SiO₂ is a trivial wide-gap
   dielectric. Ref [1] (Bi₂Se₃@SiO₂) demonstrates colloidal stability,
   biocompatibility and *retained optical properties* — **not** a topological
   photonic/phononic band structure. In this package topology is **geometry-
   induced** (the Wu-Hu C6 lattice); SiO₂/Si enter only as low-loss host media.
2. **The cited chip [2] is a different platform.** It is a programmable Si
   waveguide/microring/MZI lattice near 1525 nm — a valid *methodology* reference
   for topological transport, **not** a Bi₂Se₃@SiO₂ material model.
3. **"silicon" ≠ "silica".** The six-petal holey lattice [4] is *silicon* (a
   stiff elastic solid), not SiO₂. Any elastic model must use the correct
   density and elastic constants and the correct fabrication (unsuspended Si).
   The phononic module here is a **scalar-acoustic** analogue, not the ref-[4]
   full-vector elastic plate.
4. **A Dirac gap is not proof of non-trivial topology.** We therefore also
   provide: the C6 band-inversion indicator (§1.3), a topological invariant
   (§1.4, §4.4), interface states between two *distinct* phases (§4.1, an
   expanded|shrunken wall — both gapped but opposite character), edge-state
   dispersion and field localisation (fig5), and propagation-direction transport
   (fig6).
5. **The invariant is symmetry-dependent — not hard-coded "Chern".** These
   crystals are time-reversal symmetric, so the global Chern number is **0**
   (verified: `C_total=0`). The reported invariant is the **spin-Chern / spin
   Bott** index; a Z₂ or valley-Chern reading applies for the corresponding
   symmetry class.
6. **Under disorder, k-space Berry curvature is invalid.** §4.4 uses the
   **real-space spin Bott index** instead (fig9), which needs no Brillouin zone.
7. **"Ultralow-loss" is not implied by topology.** Topology removes
   *backscattering* (fig7), but absorption/radiation/interface/roughness losses
   set the floor. §4.5 quantifies the **absorption budget** from the real mode
   profile (fig10): with Bi₂Se₃ in the mode, loss is catastrophic regardless of
   topology.

## 1. Method

### 1.1 Photonic band structure — Plane-Wave Expansion (PWE), ✅ rigorous
TM modes (E_z) of a 2D photonic crystal solve the Maxwell master equation

    |k+G|² c_G = (ω/c)² Σ_G' ε(G−G') c_{G'} ,

a generalized Hermitian eigenproblem in a plane-wave basis (`src/pwe.py`).
ε(G) is the analytic Fourier transform of the circular silicon rods
(ε=11.7, radius 0.12a) on an air background. We use |m|,|n| ≤ 9
(361 plane waves); frequencies are reported in ωa/2πc. The solver is validated
against the empty-lattice limit and reproduces the honeycomb Dirac cone.

### 1.2 Wu-Hu geometry (`src/geometry.py`)
Triangular Bravais lattice; 6 rods per cell on a hexagon of radius `R`.
`R=a/3` → honeycomb → double Dirac cone at Γ. `R≠a/3` → gap.

### 1.3 Band-inversion diagnosis — C6 symmetry indicator, ✅ rigorous
At Γ the gap-edge modes form p (|l|=1) and d (|l|=2) doublets. We build the C6
rotation operator in the plane-wave basis (permutation G→R₆₀G) and read the
angular momentum of each doublet directly from the **real PWE eigenvectors**
(`src/symmetry.py`). Their ordering is the rigorous topological fingerprint.

### 1.4 Spin-Chern number — BHZ effective model, ⚠️ effective
The C6 band inversion maps near Γ onto a two-copy (Kramers) BHZ Hamiltonian
(`src/bhz.py`). Its mass sign is fixed by the PWE band-inversion result. Berry
curvature and Chern number per spin block are computed by the
Fukui-Hatsugai-Suzuki lattice method (gauge-invariant, exact on a discrete BZ).
`C_total = 0` (time-reversal-symmetric photonics) and `C_spin = ±1` topological.

### 1.5 Edge states & robustness — BHZ ribbon, ⚠️ effective
A ribbon (periodic in x, `Ny=40` sites in y) of the same BHZ model gives the
projected spectrum. Topological ribbons show a helical Kramers pair crossing the
gap; trivial ribbons are fully gapped. On-site disorder tests backscattering
immunity.

---

## 2. Results

### Figure 1 — `fig1_bands_dirac_gap.png`
Bands along Γ–M–K–Γ for shrunken / critical / expanded R.
- **Critical R=a/3:** four-fold degeneracy at Γ = **double Dirac cone**.
- **Shrunken R=0.30a:** trivial gap Δ ≈ 0.048; lower doublet |l|=1 (p), upper |l|=2 (d).
- **Expanded R=0.36a:** topological gap Δ ≈ 0.046; lower |l|=2 (d), upper |l|=1 (p) — **inverted**.

### Figure 2 — `fig2_berry_spin_chern.png`
Berry curvature Ω(k) concentrated near Γ. Integrated:
`C↑=−1, C↓=+1 → C_total=0, C_spin=−1` (topological) vs all-zero (trivial).
Mass sweep shows the quantized jump of C_spin at the gap-closing point.

### Figure 3 — `fig3_edge_states_robustness.png`
- Trivial ribbon: full gap, no edge channel.
- Topological ribbon: **helical edge states cross the gap** (red X).
- Topological + disorder (W=0.6): the mid-gap crossing **survives** — present in
  **20/20** random realisations. This is the backscattering-immune transport
  channel (the "ultralow-loss waveguide" of the proposal).

### Figure 4 — `fig4_param_scan_phase.png`
Geometry scan of the Γ gap vs R/a: the gap **closes at R=a/3** (Dirac) and
**reopens topological** for R>a/3. Locates the phase boundary and the R that
maximises the protected gap (largest robust bandwidth). Optimal in-window point:
R/a≈0.385, gap≈0.099.

---

## 3. What is rigorous vs. model-level (honesty ledger)

| Claim | Basis | Label |
|---|---|---|
| Double Dirac cone at R=a/3 | PWE eigenvalues of real rod crystal | ✅ rigorous |
| Gap opens on symmetry breaking | PWE | ✅ rigorous |
| Band inversion (p↔d) | C6 indicator on PWE modes | ✅ rigorous |
| Gap-vs-R phase boundary | PWE scan | ✅ rigorous |
| C_spin = ±1 quantized | BHZ + FHS (mass sign from PWE) | ⚠️ effective model |
| Helical edge states + disorder immunity (basic) | BHZ ribbon | ⚠️ effective model |
| Real-rod edge states in the gap (complete) | FDFD supercell | ✅ ab-initio |
| Guiding / bend / defect bypass (complete) | FDTD | ✅ (qualitative transport) |
| Phononic double Dirac cone + gap + inversion | acoustic PWE | ✅ ab-initio |
| Disorder invariant + TAI (real space, no BZ) | spin Bott index | ✅ rigorous |
| Absorption-limited propagation length / dB-cm | Γ (FDFD) + perturbation | ✅ method; ⚠️ illustrative ε'' |

The topological *classification* (band inversion) is rigorous; in the basic
package the invariant *number* and edge spectrum use the standard effective
model, while the complete package computes the edge states ab-initio (FDFD) and
demonstrates transport directly (FDTD). This mirrors how the Wu-Hu result is
established in the literature.

---

## 4. Complete-package stages (implemented)

### 4.1 Real-rod FDFD edge states — `fig5` (✅ ab-initio, ⚠️ removed)
A finite-difference frequency-domain supercell of the actual silicon-rod crystal
(orthorhombic a×√3a cells stacked with a topological|trivial domain wall,
periodic in y, Bloch phase in x) is diagonalised with sparse shift-invert
(`scipy.sparse.linalg.eigsh`). The bulk gap [0.436, 0.487] reproduces the PWE
value, and edge states appear **inside** it, localised at the domain wall
(localisation fraction 0.7–0.94 in the central quarter). This replaces the BHZ
ribbon of §1.5 — the edge states are now computed directly from Maxwell's
equations for the real rods.

### 4.2 FDTD wave propagation — `fig6`, `fig7`, `fdtd_bend.gif` (✅ FDTD)
A pure-NumPy 2D TM Yee-grid FDTD (leapfrog, graded absorbing sponge) launches a
CW source at f=0.46 c/a on the domain wall of the real rod crystal:
- **straight** guide: the wave follows the wall to the far end;
- **sharp double-bend** (two 90° corners): the wave turns both corners and
  continues (see the GIF);
- **point defect** (a rod removed on the wall): the wave passes essentially
  undisturbed.
Quantitatively (steady-state energy delivered ÷ input, common output plane):
straight `0.70`, defect `0.61` (barely reduced — topological immunity), trivial
reference with no wall `0.08`. The sharp bend is shown qualitatively (fig6 +
GIF); its longer bent path is not directly comparable on the same plane in this
coarse solver. **Honesty note:** the absorbing sponge is not a true UPML and the
source is a soft point source, so the transmission is qualitative; a production
run would use UPML + an eigenmode source for calibrated dB numbers.

### 4.3 Phononic twin — `acoustic.py`, `fig8` (✅ acoustic PWE)
The identical Wu-Hu lattice is solved with a genuine scalar-acoustic PWE kernel
(generalised eigenproblem with both η=1/ρ and β=1/B varying; rods ρ=6, B=2 vs
air-like background). It hosts a **double Dirac cone at R=a/3** (four-fold at
f≈1.09, on bands 3–6) and a C6-inverted gap (≈0.067). The band-inversion
diagnosis (same C6 indicator) shows the **topological side is R<a/3**, opposite
to the photonic case — a genuine, material-dependent difference, not a bug. This
provides the second system for the photonic-vs-phononic comparison (the elastic
GHz waveguide of ref [4] is the full-vector version of this scalar model).

### 4.4 Real-space spin Bott index under disorder — `fig9` (✅ rigorous)
Because disorder breaks periodicity, k-space Berry curvature no longer applies.
We put the two BHZ Kramers blocks on an L=14 torus, add Anderson on-site
disorder, and compute the Loring-Hastings Bott index of each block's occupied
projector; the spin Bott index `B_s=(B↑−B↓)/2` is a real-space invariant with no
Brillouin zone. Ensemble-averaged sweeps show two regimes:
- clean topological (M=+1): `B_s=1` stays quantised, then **collapses** at strong
  disorder (topology destroyed);
- clean trivial (M=−0.6): `B_s` jumps **0→1** at intermediate disorder — a
  **Topological Anderson Insulator** (disorder-*induced* topology), the mechanism
  of the cited photonic topological-Anderson insulator [3].
This replaces the naive "Berry curvature survives disorder" reasoning. (A full
study would add the local Chern marker map and bulk/edge scattering transmission.)

### 4.5 Loss budget — `fig10` (✅ perturbative, ⚠️ illustrative constants)
Topology removes backscattering (fig7) but not absorption. From the **actual FDFD
edge mode** we compute the confinement factor `Γ_rod=0.93` (the TM mode sits
almost entirely in the rods) and feed the confinement-factor modal-loss formula
`α = (2π/λ0)·Γ·κ`, `κ=ε''/(2√ε')`, at a≈0.71 µm (λ0≈1.55 µm):
- low-loss Si/SiO₂ (ε''≈1e-3): α≈48 dB/cm, L_p≈1.8 mm — a real ultralow-loss guide;
- SiO₂-shelled Bi₂Se₃ (ε''≈0.5): α≈2.4×10⁴ dB/cm, L_p≈4 µm;
- bare Bi₂Se₃ (ε''≈8): α≈3.8×10⁵ dB/cm, L_p≈0.2 µm — catastrophic.
**Conclusion:** "ultralow-loss" and "Bi₂Se₃ carries the mode" are in direct
tension; the absorption floor is set by the Bi₂Se₃ fraction, not by topology.
Optical constants are representative — dropping in measured Bi₂Se₃/SiO₂ values is
a data step, not a code change.

## 5. Remaining extensions (not yet done)

1. **Full-vector elastodynamics** (Lamb/plate waves) for the true ref-[4] system,
   replacing the scalar-acoustic approximation.
2. **UPML + eigenmode-source FDTD** for calibrated transmission/​reflection in dB,
   and a full disorder-strength sweep (ensemble-averaged edge vs bulk).
3. **Spin-Chern from real modes** via Wilson loops / Wannier-centre flow on the
   PWE/FDFD bands, removing the effective-model caveat on the invariant *number*.
4. **Local Chern marker / spin Bott map** in real space + scattering-matrix
   transmission, and **measured** Bi₂Se₃/SiO₂ optical constants for the loss
   budget (radiation and surface-roughness terms included).

---

## 5. Reproduce

```bash
pip install -r requirements.txt
cd src && python run_all.py
```
Runtime ≈ 90 s on a laptop core; deterministic (fixed disorder seeds).
