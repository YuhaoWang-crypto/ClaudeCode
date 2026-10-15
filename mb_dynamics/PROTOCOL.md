# Pre-registered protocol — three-route comparison of conformational dynamics in pulsed electron transfer

Frozen and committed **before** any route was fitted to the 2017 curves. Changes after this commit are
listed in the report under "deviations".

## Question
With the same data and the same electron-transfer assumption, does keeping conformational exchange
predict SWV current better than a static / fast-equilibrium average?

## Shared assumptions (all routes)
- Surface-confined MB, apparent simultaneous 2e⁻ Butler–Volmer: k_red = k_loc·exp(−α n f η),
  k_ox = k_loc·exp((1−α) n f η), n = 2, f = F/RT, T = 298.15 K, η = E − E⁰.
- k_loc(z) = k_c·exp[−β(z − z_min)], z = MB-centre height above the SAM plane, z_min = MB radius
  (0.35 nm, i.e. contact). Main β = 5 nm⁻¹; 2 and 10 nm⁻¹ are pre-declared sensitivity cases.
  Modelling assumptions, not measured decay constants. Exponential averaged **within** each
  conformational state, never applied to a mean height. MB treated as an isotropic proxy.
- Column-generator convention: dp_R/dt = L p_R + K_red p_O − K_ox p_R; dp_O/dt = L p_O − K_red p_O + K_ox p_R.
  O and R share the conformational landscape and L (main analysis). Conformational moves carry no
  faradaic current. I = Q Σ_i [k_red,i p_O,i − k_ox,i p_R,i] (cathodic positive, as in the figure);
  Q = nFAΓ fitted as one number per electrode-chemistry group.
- Exact SWV waveform: staircase from 0 V towards −0.45 V in 1 mV steps; each step = forward half-pulse at
  E_k − E_sw then backward half-pulse at E_k + E_sw, each 1/(2f) long; currents sampled at the end of each
  half-pulse; net = forward − backward. State carried across every half-pulse and step (never reset);
  initial state = Nernst equilibrium at 0 V.
- Constant-potential propagation exact (matrix exponential via symmetric eigendecomposition; checked
  against scipy `expm` on the full 2M system and against a stiff ODE integrator).

## Routes
- **R0 static / fast equilibrium**: conformations always at equilibrium; single 2-state O/R system with
  k_eff = Σ_i π_i k_loc,i over the full height distribution of the flexible model.
- **R1 1-D Smoluchowski**: F(z) = −k_BT ln P(z) from the same flexible-model trajectories; finite-volume
  generator (square-root approximation, non-negative, discrete detailed balance, no-flux walls); 40 vs 80
  cells. D from Stokes–Einstein for a 0.35 nm sphere in water (free-draining proxy, not a measured
  tether D).
- **R2 MD-constrained multi-state**: states from the unbiased flexible BD trajectories; reversible
  continuous-time generator with π fixed to trajectory weights, parameterized by non-negative symmetric
  fluxes, fitted by the lag-τ transition likelihood of exp(Lτ). Never a matrix log, never edge rates fitted
  to current.
- **Historical reference**: the authors' published k⁰ = 0.4 cm s⁻¹, α = 0.37 (different model; listed, not ranked).
- **Non-MD control (5 parameters)**: two non-exchanging static populations, weights ½/½, rates k·e^{±δ},
  δ ∈ [0, ln 10⁴], plus α, E⁰, Q.

## Flexible molecular model (default)
OpenMM CPU Brownian integrator, one bead per nucleotide, 20 nt; bond 0.6 nm, 100 k_BT nm⁻²; l_p 1 nm;
DNA bead diameter 0.5 nm; MB radius 0.35 nm; C6 arm = 2 × 0.3 nm, MB linker = 3 × 0.3 nm; harmonic bonds,
cosine bending, WCA excluded volume, repulsive plane at z = 0 (SAM top), anchor fixed. No attraction to the
surface, no base pairing/stacking. Free-draining Stokes friction (water, 298.15 K). ≥3 seeds,
extended and collapsed starts. Time step chosen from {0.1, 0.25, 0.5, 1} ps: the fastest one whose halving
changes equilibrium observables and correlation times by ≤5 %.
Pre-declared one-factor sensitivity: MB linker total 0.6 / 0.9 / 1.8 nm; linker l_p 0.7 / 1 / 2 nm.

## R2 gates (molecular-trajectory validation only; current data never used)
Descriptors: z_MB, DNA R_g, anchor–MB distance, MB–DNA non-adjacent contacts. Candidate state counts
4, 8, 12; choose the smallest passing. Seeds 1–2 build, seed 3 validates. Shortest passing lag τ with:
slowest two resolvable timescales change ≤20 % across τ, 2τ, 4τ; CK / direct propagation at 2τ, 4τ, 8τ
with main-state probability error ≤0.10 and consistent with seed-3 block uncertainty. States covering 95 %
of the current contribution: decorrelated ESS ≥200 and ≥50 entries/exits, from ≥2 seeds. Block (seed
and half) agreement of equilibrium weights / contact occupancy / correlation times; block PMF difference
≤1 k_BT over the equilibrium-important region. Failing R2 is reported as a conditional prediction, not
ranked.

## Data and folds
Primary benchmark: 9 digitized DNA20–MB curves (see `data/dauphin2017/MANIFEST.md`); retrospective,
already-seen data.
- **Primary CV**: leave-one-frequency-out (hold out all 3 curves at 20, 100 or 500 Hz).
- **Secondary CV**: leave-one-amplitude-out (25, 50, 75 mV).
Same folds, preprocessing and readings for all routes.

## Fits
- **A (shared θ)**: one θ = (k_c, α, E⁰, Q) per training fold, jointly over R0, R1, R2 (routes equal
  weight, training curves equal weight, normalized squared error); s = 1 for R1/R2.
- **B (per-route calibration)**: each route fits (k_c, α, E⁰, Q); R1, R2 add one global motion multiplier
  s ∈ [10⁻³, 10³] scaling D or L (not per frequency/amplitude/test case).
- Bounds: k_c 10⁻³–10⁸ s⁻¹, α 0.1–0.9, E⁰ −0.45–0 V (vs Ag/AgCl), Q 10⁻¹²–10⁻⁵ C. Boundary hits and
  flat profiles are reported.
- Objective per curve: Σ over the three components (forward, backward, net) of SSE / (peak-to-valley of
  that observed component)², averaged over components; curves weighted equally.
- Multi-start (≥5 starts) per fit; profiles of k_c and s on the full data set.
- β and geometry sensitivities are reported, never chosen by test error.

## Metrics and acceptance
Per held-out curve: **NRMSE of the net voltammogram** = RMSE / observed net peak-to-valley (primary);
forward/backward NRMSE (secondary); relative error of net peak height. Target NRMSE ≤30 %, peak error ≤30 %.
Kinetic gain claim requires a ≥10 % relative reduction of the held-out NRMSE vs R0, also compared with R1
vs R2 and with the 5-parameter non-MD control. Digitization floor ≈3 % NRMSE. Nine curves on presumably one
electrode: no statistical claim of superiority. Labels kept separate: numerically correct / molecular
dynamics internally gated / time scale experimentally calibrated / held-out target met / target-system transfer.
