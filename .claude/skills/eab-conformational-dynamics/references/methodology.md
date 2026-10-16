# Methodology

## Contents
1. Electron-transfer model and current
2. Generator convention and the exact propagator
3. Multi-species chemistry (C0, C1, C2)
4. Routes R0 / R1 / R2 / static control
5. R2 estimation and gates
6. Fits, folds, objective, metrics
7. Diagnostics

## 1. Electron-transfer model and current
- Surface-confined reporter, Butler–Volmer per electron-transfer edge:
  k_red = k·w·exp(−α n f η), k_ox = k·w·exp((1−α) n f η), η = E − E°, f = F/RT.
- Distance factor per conformational state i: w_i = ⟨exp(−β(z − z_min))⟩ over the samples in state i
  (z = reporter-centre height above the monolayer plane, z_min = reporter radius = contact). β is an
  assumption (main 5 nm⁻¹; 2 and 10 nm⁻¹ as sensitivity) and is not separately identifiable from k in R0.
- Current comes only from redox jumps; conformational moves carry no faradaic current.
  I = (Q/2) Σ_edges n_e Σ_i (k_red p_ox − k_ox p_red), Q = 2FAΓ, cathodic positive (state the sign
  convention and reference electrode explicitly).
- SWV (Osteryoung): staircase E_k from E_start in steps ΔE towards the reduction direction; forward
  half-pulse at E_k − E_sw, backward at E_k + E_sw, each 1/(2f); sample at the end of each half;
  net = forward − backward. State is carried through every half-pulse and step; start at Nernst
  equilibrium of E_start.

## 2. Generator convention and the exact propagator
- Column generator: L[i, j] = rate j → i, columns sum to zero, reversible: L_ij π_j = L_ji π_i.
- Shared landscape assumption (O, R and intermediates share π and L) makes the joint chain
  reversible with respect to P(E) = π ⊗ q(E) (q = Nernst species fractions).
- Propagate the deviation u = p − P(E) (homogeneous: du/dt = G u). On a potential change p is
  continuous, so u jumps by P_old − P_new. Current per edge = k_ab u_a − k_ba u_b: no cancellation.
- Symmetrize v = P^{-1/2} u; −A = GᵀG with one row per edge:
  conformational edge (i,j), species s: √F_ij (e_is/√π_i − e_js/√π_j);
  reaction edge a→b in state i: √k_ab(i) e_ia − √k_ba(i) e_ib (no q appears).
  Take the spectrum from the SVD of G: the slow reactive eigenvalue keeps relative accuracy
  ~ε√(‖A‖/|λ|) even when exchange is 10⁸–10¹³× faster than reaction. (eigh of A loses it: ε‖A‖/|λ|.)
- Half-pulse: v(h) = W diag(exp(−σ²h)) Wᵀ v. Handle q in log space.
- Verification set (all in tests/): dense expm per half-pulse; stiff Radau over a whole scan; s → ∞
  reproduces R0; s → 0 reproduces static disorder; reversible (k → ∞) and irreversible limits vanish;
  long hold reaches Nernst; per-conformation population conserved; known CTMC rates recovered.

## 3. Multi-species chemistry
| name | scheme | parameters (+Q) |
|---|---|---|
| C0 | O + 2e⁻ ⇌ R (apparent concerted) | k, α, E° |
| C1 EE | O + e⁻ ⇌ I, I + e⁻ ⇌ R (protonation fast, folded into apparent potentials) | k₁, k₂, α, E_m, ΔE = E₁ − E₂ |
| C2 E–C–E | O + e⁻ ⇌ I; I ⇌ IH (k_p, K_p fixed); IH + e⁻ ⇌ R | k₁, k₂, α, E_m (apparent overall), ΔE, k_p |
C2: E₁ + E₂′ = 2E_m − (RT/F) ln K_p. Both electron steps take w; protonation does not. Bounds used:
log k −3…8, α 0.1–0.9, E_m within the scan, ΔE ±0.3 V, log k_p 0…8.

## 4. Routes
- R0: one conformational state with w = ⟨w⟩ (fast equilibrium). Exactly the classic surface model.
- R1: finite-volume Smoluchowski on z, square-root approximation rates
  L_{i±1,i} = (D/Δz²)√(π_{i±1}/π_i), no-flux walls, π from the sampled histogram (+0.5 pseudo-count),
  w per cell from samples. D = k_BT/(6πηa) (free-draining proxy) × s. Compare 40 vs 80 cells
  (20 cells acceptable for fitting if re-checked at 40/80 at each optimum).
- R2: states from trajectories, reversible generator (section 5), × s.
- Static control: two non-exchanging populations, weights ½, rates k·e^{±δ}, δ ∈ [0, ln 10⁴].
- Historical/literature models are listed, never ranked, when their model class differs.

## 5. R2 estimation and gates
- Descriptors (pre-declared): reporter height, DNA radius of gyration, anchor–reporter distance,
  reporter–DNA non-adjacent contacts; standardized; k-means with candidate counts 4/8/12.
- π fixed to trajectory weights; parameters θ = log F_ij (symmetric equilibrium fluxes ≥ 0);
  L_ij = F_ij/π_j; maximize Σ C(j→i) log[exp(Lτ)]_ij with the analytic Daleckii–Krein gradient.
- Gates: (a) slowest two resolvable timescales (t > τ) change ≤20 % from τ to 2τ and 4τ;
  (b) CK / direct propagation at 2τ, 4τ, 8τ: max |P_pred − P_emp| ≤ 0.10 over states with π ≥ 0.05 on
  the held-out seed, and compatible with block-bootstrap SE (≤3 SE or ≤0.03); (c) states carrying 95 %
  of Σπ_i w_i: ESS ≥ 200, ≥ 50 committed entries (residence ≥ τ), entries from ≥2 seeds;
  (d) seed × half blocks: equilibrium weights, contact occupancy, correlation times compatible;
  PMF difference ≤ 1 k_BT over the important region.
- Re-estimate on all seeds after validation, keep both records.
- Never copy short trajectories to fake long sampling; never slow the dynamics to manufacture a difference.

## 6. Fits, folds, objective, metrics
- Objective: mean over training curves of (1/3) Σ_components SSE/(observed peak-to-valley)²/N_points;
  Q profiled analytically, bounded 10⁻¹²–10⁻⁵ C.
- Fit A (shared θ, s = 1, routes equal weight) and Fit B (per route + global s).
- Folds: leave-one-frequency-out (primary), leave-one-amplitude-out (secondary); nested selection only
  inside training folds.
- Report parameter counts, train vs held-out gap, all start optima, parameters at bounds, profiles
  (k and s), and R1 grid differences.
- Acceptance per held-out curve: NRMSE ≤ 30 % and |peak error| ≤ 30 %; kinetic gain ≥ 10 % mean
  reduction vs R0 within the same chemistry.

## 7. Diagnostics
- Decomposition per route solution: R0 at the same chemistry point; route at s = 1.
- R0*: best training basin across all starts of all routes, evaluated as R0.
- Identifiability: held-out range among points within 5 % of the best training cost.
- Static control at δ = 0 with the same route (not via R0: its rate is not ⟨w⟩-scaled).
- Profiles: flat s profile ⇒ the data carry no information on conformational speed.
