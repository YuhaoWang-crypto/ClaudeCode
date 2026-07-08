# Methodology — definitions, theorems, gotchas

The math behind each tool, at the depth needed to apply it correctly and to
know where the honest boundaries are.

## 1. Symmetry: automorphism → quotient, and fibration

- **Graph automorphism group** Aut(G): permutations of nodes preserving the
  (signed, directed) edge set. A non-trivial orbit = a set of interchangeable
  nodes (e.g. HRAS/KRAS/NRAS wired identically → Aut = S₃). Quotient by the
  orbit collapses redundant nodes to one effective node (`RAS*`). Deterministic
  graph theory — fully rigorous.
- **Input-tree fibration** (Morone–Leifer–Makse 2020) is the *generalization*:
  nodes with input-isomorphic trees synchronize even without a global
  automorphism. Implement as **minimal balanced colouring** (seed by role,
  refine by in-edges to a fixed point). Compresses more than automorphism (MAPK
  27→11 fibers vs. the strict S₃ orbit only).
- Biomarker payoff: each fiber → one pan-assay representative readout
  (MAPK1/MAPK3 → pERK1/2). The irreducible-core node identity is a structural,
  a-priori biomarker.

## 2. CRNT deficiency δ = n − ℓ − s

- n = number of complexes, ℓ = number of linkage classes, s = rank of the
  stoichiometric matrix.
- **Feinberg deficiency-zero theorem:** δ=0 + weakly reversible ⇒ for *any* rate
  constants there is a unique, stable steady state within each compatibility
  class ⇒ **bistability is impossible**. A switch therefore *requires* δ≥1.
- δ is an integer topological invariant that overrides all concentrations/rates.
- **GOTCHA (real bug):** compute δ on the **effective** network. Species held at
  constant concentration (chemostats, e.g. A, B in Schlögl) must be stripped
  from the complexes before counting — otherwise δ comes out wrong (Schlögl
  falsely δ=0 instead of δ=1). Always confirm predicted attractor count by
  integrating the mass-action ODE from multiple initial conditions.

## 3. Elementary flux modes (EFM)

- Steady-state fluxes satisfy S·v = 0, v ≥ 0 — a polyhedral cone. Its **extreme
  rays are the EFMs**: support-minimal, non-decomposable pathways. The nearest
  rigorous analogue to "prime factorization / irreducible generators."
- Verify: any feasible flux reconstructs as a non-negative combination of EFMs
  (nnls residual ~1e-16), and max|S·v| ~1e-15 confirms the null-space.

## 4. Critical slowing, Lyapunov, DNB — the biomarker layer

- **Largest Lyapunov exponent (LLE):** Benettin algorithm. VALIDATE the code on
  a system with a known positive exponent first — Rössler gives LLE ≈ +0.0737
  (literature ≈ +0.071). This proves "a positive λ can be computed and
  corresponds to chaos" before using λ→0 as a tipping signal.
- At a fixed point, the **leading Jacobian eigenvalue = local LLE**. Approaching
  a bifurcation it → 0 (critical slowing); **recovery time τ = −1/λ_max → ∞**.
- **Model-free early-warning signals** (need only a time series): rising
  variance / SD, rising lag-1 autocorrelation, rising DNB index (Chen–Aihara:
  within-module SD × correlation). These are the deployable biomarkers.
- **Honest boundary:** in a 2-D fold system λ rises from negative *toward* 0
  (loss of stability = state collapse = fate flip). A strictly *positive* λ
  (true chaos) needs ≥3-D dynamics. Molecular-dynamics (ns-scale conformational)
  Lyapunov and network/cell-scale Lyapunov are **different scales**; the same
  Benettin algorithm applies to both, but **"drug binding pushes λ positive ⇒
  toxicity" is a testable HYPOTHESIS**, not a proven theorem. What is rigorously
  deliverable: λ changes sign at the bifurcation as an in-model, computable
  stability biomarker.

## 5. FIM / sloppy models / stiff axes (biomarker selection)

- Build the Fisher Information Matrix FIM = SᵀS from sensitivities of outputs to
  log-parameters; eigen-decompose. Sloppy spectrum spans ~tens of orders; the
  few **stiff axes** carry ~all identifiable information.
- **GOTCHA (real bug):** rank observables by **dimensionless relative**
  sensitivity, not absolute — otherwise a single concentration wins on its unit
  (nM) magnitude alone. With relative sensitivity, **ratio biomarkers**
  (e.g. Mpp/M, ppERK/(pERK+ppERK)) load best on the stiff axis — more
  identifiable than any single species, matching the Gutenkunst–Sethna result.

## 6. Three bifurcation classes (see SKILL.md taxonomy for the read-out rule)

- **Saddle-node** — S-shaped hysteresis, two stable branches + a saddle; locate
  the fold edges by bisection on the number of stable roots.
- **Hopf** — a complex eigenvalue pair crosses the imaginary axis: Re(λ)→0 with
  Im(λ)≠0. Brusselator has the exact analytic Hopf at B = 1 + A². Signature is
  spectral, not autocorrelative.
- **SNIC** — saddle-node ON an invariant circle (θ / Ermentrout–Kopell normal
  form dθ/dt = (1−cosθ) + (1+cosθ)I, SNIC at I=0). Period T ~ π/√I diverges;
  log-log slope of T vs I = −0.50 is the exact fingerprint.
