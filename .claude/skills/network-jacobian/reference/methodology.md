# Methodology: definitions, checks, honest limits

## Coordinates
- `x` = floating species not defined by assignment rules; `dx/dt = N v(x)`
  (N divided by compartment volume; `core` verifies ‖f − Nv‖ at load).
- `B` = orthonormal basis of range(N). Conservation laws keep x on x* + B z.
- `Jz = Bᵀ N ε B` (m × m, m = n − #conservation laws). Its eigenvalues are the
  dynamical ones. The full Jacobian adds `#conservation` structural zeros.

## Metabolic control analysis (✅ checks run every time)
- `C^S = −B Jz⁻¹ BᵀN`, `C^J = I + ε C^S` (unscaled); scaled `C^J_ij = C^J_ij v_j / v_i`.
- Summation: `C^J K = K` for K = null(N) ⇒ Σ_j scaled C^J_ij = 1.
- Connectivity: `C^J ε B = 0`.
- Brute-force check `control_fd`: perturb each enzyme ±0.1 % and re-solve the
  steady state. Teusink: max |analytic − FD| = 1.6e-6.

## Stability classification (what eigenvalues can and cannot say)
- At a hyperbolic fixed point, Re λ_lead < 0 ⇒ locally stable (node if real,
  focus if complex). Re λ_lead > 0 with Im ≠ 0 ⇒ unstable focus; a limit
  cycle **must be confirmed by simulation** (or Floquet multipliers).
- Chaos cannot be read from fixed-point eigenvalues. It needs the largest
  Lyapunov exponent along the attractor (validated Benettin code:
  `grn_pipeline.m4_dnb_lyapunov.benettin_lle`).
- "Eigenvalue near 0 ⇒ low resilience" is ✅ as a statement about the
  linearised model (recovery time τ = 1/|Re λ|). Mapping it to ageing or
  disease susceptibility in people is ⚠️.

## Fragility
- First order: `dλ/dε_rk = (wᴴ BᵀN)_r (B u)_k / (wᴴ u)`; scaled edge fragility
  `F_rk = ε_rk dλ/dε_rk`; local enzyme fragility `F_r = Σ_k F_rk`.
- **Total** fragility = d Re λ / d ln E including the shift of the steady
  state (`spectrum.total_fragility`, finite difference with re-solve). This is
  the quantity a real perturbation sees. In Wolf at 98 % of the Hopf, total
  predicted ±10 % crossings with r = 0.989 and local-only with r = 0.16.
- Probe placement: right eigenvector of the critical eigenvalue, relative
  participation |u_k|/x_k (assays have relative CV).

## Fluctuation ↔ Jacobian
- LNA: `Jz Cz + Cz Jzᵀ + Dz = 0`, D = N diag(v) Nᵀ/Ω. (Sun 2015 writes
  JC + CJᵀ = −2D: same equation with D → D/2.)
- Inverse: linear in ε on a sparsity mask; `fluctuation.inverse_elasticities`
  (column-scaled LS, optional Tikhonov prior toward a reference model).
- Honest limit (D1-D): exact C → 0.6 % error. Stiffness |λ|max/|λ|min ≈ 3e5
  ⇒ a 1e-4 relative error in C gives ~100 % error in ε; sampled C needs
  N ≈ 1e6 for <10 % error. A prior does not rescue N ≤ 3000. Use challenges
  (D1-C) or time series (D2) instead, or shrink to a non-stiff sub-module.

## Relaxation estimation
- `relax.hankel_dmd`: delay-embedded DMD; exact on noise-free data, biased by noise.
- `relax.fit_modes`: global variable-projection fit with rates shared
  across observables and amplitudes per observable, plus optional damped
  oscillations. Initialise from DMD.
- `relax.fit_decay`: single exponential + offset on a post-peak segment (CGM).

## Real-data CGM protocol (D2-C)
- Dexcom trace resampled at 5 min. Isolated meals only (≥ 20 g carbs, no meal
  within the next 3 h). Peak ≤ 120 min, excursion ≥ 15 mg/dL. Decay fit from
  peak to +3 h with R² ≥ 0.7. Per-person median.
- Groups by HbA1c (ADA: < 5.7, 5.7–6.4, ≥ 6.5). Because the groups are
  defined by HbA1c, fasting glucose is a near-circular benchmark. Always also
  report the partial Spearman given fasting glucose.
- Results: k medians 0.023 / 0.014 / 0.010 min⁻¹ (τ ≈ 43 / 73 / 100 min);
  Spearman vs HbA1c −0.44 (p = 0.002); partial given fasting glucose −0.09
  (p = 0.6). Night AR1 +0.36 (p = 0.02), partial +0.07.
