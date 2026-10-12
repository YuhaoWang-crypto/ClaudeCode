# Pairwise → network: attractor analysis of interacting indicators

`netjac/loops.py` is the engine. `demos/d5_pairwise_to_network.py` covers
models and theory, and `demos/d6_baroreflex_jacobian.py` covers real data.

## The workflow for any set of interacting indicators

1. **Write or obtain the equations.** Use a curated BioModels model when one
   exists (`loops.sbml_system(id)` handles reaction-based and rate-rule ODE
   models). Otherwise write a `loops.System(name, f, names, params)`. With
   only time series, estimate the Jacobian from data (step 7).
2. **Find every equilibrium.** Use `System.equilibria(box)` (multi-start).
   One model can hold several states, e.g. healthy, threshold and disease in
   Topp 2000.
3. **Classify each one.** `analyse_point(J)` returns stable node / stable
   focus / saddle / unstable focus pair / centre / non-hyperbolic.
   Linearisation is **local**. `global_portrait` simulates many starting
   states and reports which attractor each one reaches, or oscillation or
   divergence.
4. **Look at every sub-system.** `principal_subsystems(J, k)` gives every
   k×k block with the other variables clamped. This is the "pairwise
   matrix" view.
5. **Combine them exactly.** The whole is governed by
   det(λI − J) = λⁿ + a₁λⁿ⁻¹ + … + aₙ, and each a_k is a signed sum over sets
   of disjoint feedback loops (Levins). `loop_decomposition` lists each
   loop's contribution, and `hurwitz` decides stability.
   - Necessary: every Levins feedback F_k = −a_k < 0.
   - A long loop enters only the high-order a_k. That is how a 3-step
     negative loop destabilises a network whose pairs are all stable
     (Goodwin: unstable once the loop gain exceeds sec(π/3)³ = 8).
6. **Map the parameter dependence.** `bifurcation_scan` reports:
   - Hopf → oscillation;
   - saddle-node → the state disappears (tipping point);
   - transcritical → two states exchange stability.
   The leading eigenvalue → 0 before a saddle-node (critical slowing).
7. **Data-driven Jacobian (D6).** Fit a VAR on multichannel time series.
   - VAR(1): J = logm(A)/Δt is the direct pairwise Jacobian.
   - VAR(p): the companion-matrix eigenvalues ln(μ)/Δt include the
     delay-embedded modes. A real physiological loop with delays usually
     needs p > 1.

## Rules of thumb (all demonstrated)

| Statement | Evidence |
|---|---|
| All pairs stable does **not** imply the whole is stable | Goodwin; random 5×5: 21 % of "all pairs stable" networks are unstable; Sturis/Tolić: 0/15 pairs unstable, yet the whole oscillates (110-min ultradian cycle) |
| A pair unstable on its own does **not** imply the whole is unstable | Topp: the {G,B} and {I,B} pairs run away when insulin/glucose is clamped; the fast G–I loop holds them |
| A linear centre is **not** a limit cycle | Orbit size = shock size; ±0.01 of damping flips it. A real cycle needs nonlinearity (fixed amplitude whatever the shock) |
| Sign pattern alone can decide | Prey self-damping + predator–prey pair: `sign_pattern_stability` = 1.0 |
| Diagonal self-damping = "buffering" (γ, δ) | Degradation, dilution, product inhibition, saturation. Gershgorin: enough self-damping beats any coupling |
| More variables is **not** automatically better | May 1972: random networks lose stability past n* = (d/σ)²/C. Feedback-paired structure keeps them stable |
| A pairwise 2×2 fit to delayed real data misses the oscillation | D6 (286 real subjects): the 2×2 VAR(1) Jacobian's eigenvalues lie mostly on or near the real axis; the Mayer mode appears only with delay embedding (p ≈ 12 samples = 6 s) |
| A real feedback loop loses resilience with age, beyond its gain | D6: baroreflex mode decay rate ρ=−0.23 with age, still −0.23 (p=6e-4) after adjusting for BRS; BRS itself ρ=−0.54 (positive control) |

## Pitfalls hit while building this

- A solver "failure" can sit exactly on the solution (scipy `hybr`: "no
  progress"). Judge roots by their residual, with a tolerance relative to
  the model's own rate scale.
- A roadrunner instance is stateful: a parameter set during a scan leaks
  into later calls unless it is reset (`sbml_system` now does this).
- Branch continuation can hop branches at a transcritical point. Detect a
  zero eigenvalue as well as sign changes.
- Do not run `grep -v "^\["` on numpy output: it deletes every array line.

## Beyond eigenvalues: `netjac/statespace.py` (D7)

Each tool answers a different, testable question about the same J:

| Tool | Question | Biology use |
|---|---|---|
| `phase_map` / `classify_2x2` | where do parameter pairs sit (node / focus / centre / saddle) | bifurcation map of a 2-variable loop |
| `quadratic_energy` | a Lyapunov energy; with a diagonal P, the diagonal of J is the only dissipation and the off-diagonal loop only exchanges energy | degradation / dilution = dissipation; interaction loops = oscillation |
| `reactivity`, `max_amplification` | does a stable system first amplify a shock, how much, from which direction | rebound / overshoot after a dose |
| `noise_response` | quasi-cycles under random shocks (spectral peak), lead-lag, contemporaneous correlation | Mayer waves; why static correlations miss strong causal coupling |
| `critical_delay` + `simulate_dde` | the delay on one coupling that creates sustained oscillation | hepatic delay (Sturis), baroreflex delay |
| `gramians` | controllability (which input steers) / observability (which single readout suffices) | drug-target choice, minimal biomarker panel |
| `pole_placement`, `lqr` | the feedback gains for a chosen damping | closed-loop dosing |
| `coupling_sensitivity` | dλ/dJ_ij: which entry sets damping vs frequency | where to intervene |
| `var_jacobian` | J = logm(A)/Δt from time series, with an equilibrium | data-driven state space, era by era |

D7 results on the user's macro model: V = x² + 0.625y², dV/dt = −0.4x² − 0.125y². Reactive iff γδ < 0.0225. The quasi-cycle exists for a stable focus. Contemporaneous corr is +0.02, but lagged corr is +0.72. Critical policy lag τ* = 0.78.

D7 on real FRED data (x = UNRATE − NROU, y = CPI YoY − 2):
- 1960–83 has the user's sign pattern with δ ≈ 0 (near-centre, scenario 2).
- 1984–2007 has δ = 0.56 and 5× faster damping (soft landing).
- 2008–19 has a Phillips coupling indistinguishable from 0.
