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

---
## Deviations recorded after the freeze (none used current data to choose anything)
1. Excluded volume and wall: WCA replaced by truncated harmonic repulsion (1000 k_BT nm⁻², same contact
   distances) because WCA on 0.15 nm linker beads is unstable for every allowed Brownian time step.
2. 64 non-interacting replicas per OpenMM context; a "seed" is one context.
3. R2 lag candidates extended to 5 and 10 ns; CK compatibility = deviation ≤ 3 SE or ≤ 0.03.
4. Buffer not stated in the source; water viscosity assumed.
5. Fit objective uses forward, backward and net components; primary metric is net NRMSE (as stated above).

---
## Round 2 addendum — stepwise MB electron/proton chemistry (frozen before any round-2 fit)

**Question.** Repeat the three-route comparison with MB reduction modelled stepwise instead of as a concerted
2e⁻ Butler–Volmer step, and report separately whether the chemistry change reduces the shared misfit.
Routes are compared only within a chemistry; chemistries are compared with parameter counts shown and are
not mixed into the route ranking.

**Chemistries** (Q = 2FAΓ in all; current = (Q/2)·Σ electrons transferred; cathodic positive):
- **C0** concerted apparent 2e⁻ BV (round 1): k, α, E⁰ (+Q) — 4 parameters.
- **C1 EE**: O + e⁻ ⇌ I (k₁, E₁), I + e⁻ ⇌ R (k₂, E₂), n = 1 each, shared α; protonation fast and folded
  into apparent potentials at pH 7.2. Parameterized as log k₁, log k₂, α, E_m = (E₁+E₂)/2, ΔE = E₁ − E₂
  (+Q) — 6 parameters.
- **C2 E–C–E**: O + e⁻ ⇌ I (k₁, E₁); I ⇌ IH semiquinone protonation, pseudo-first-order k_p at pH 7.2,
  equilibrium constant K_p fixed at 10 (sensitivity 1 and 100, R0 only); IH + e⁻ ⇌ R (k₂, E₂′). E_m is the
  apparent overall potential, E₁ + E₂′ = 2E_m − (RT/F) ln K_p, ΔE = E₁ − E₂′. Parameters log k₁, log k₂, α,
  E_m, ΔE, log k_p (+Q) — 7.
- Both electron steps carry the same distance factor w(z) (rates k·w); protonation is solution-side and not
  distance-dependent; O, I, IH, R share the conformational landscape and L (assumption, as in round 1).
- Bounds: log k 10⁻³–10⁸ s⁻¹, α 0.1–0.9, E_m −0.45–0 V, ΔE −0.3–+0.3 V, k_p 1–10⁸ s⁻¹, s 10⁻³–10³.

**Routes per chemistry**: R0 (fast equilibrium, ⟨w⟩), R1 (40-cell Smoluchowski, + global s), R2 (12-state
generator, conditional only — it failed its dynamic gates in round 1), and a static two-population control
(weights ½, both electron-transfer rates scaled by e^{±δ}, δ ∈ [0, ln 10⁴]), giving it the same parameter
count as R1/R2. β = 5 nm⁻¹ only (round 1 showed β changes no held-out metric by more than 0.1 %).

**Numerics**: the joint chain (conformation × redox species) is reversible with respect to π⊗q(E) because every
step obeys detailed balance and the landscape is shared, so the round-1 exact symmetric square-root/SVD
propagator generalizes. Tests before fitting: identity with the round-1 backend for C0; full-system expm per
half-pulse and stiff Radau over a scan for C1/C2 with random reversible conformational models; fast/slow
exchange limits; conservation of per-conformation population.

**Fits and evaluation**: identical data, preprocessing, folds (primary leave-one-frequency-out, secondary
leave-one-amplitude-out), objective, Q profiling and metrics as round 1. Fit B for every route and chemistry on
all folds; fit A (shared θ, s = 1) on the full data and the frequency folds. Starts: 4 per fit (3 Latin-hypercube +
1 from the R0 solution of the same chemistry and the same training fold — never from a fit that saw the
held-out curves). Kinetic-gain criterion unchanged (≥10 % held-out NRMSE reduction vs R0 within a chemistry).
Chemistry comparison: held-out NRMSE and the count of curves meeting NRMSE ≤ 30 % and |peak error| ≤ 30 %,
with parameter counts; no chemistry is declared "validated" from 9 retrospective curves on one electrode.

**Round-2 compute amendment (before any round-2 fit).** R1 under C2 costs 4.4 s per evaluation at 40 cells,
which would exceed the 24 h cap. R1 is therefore fitted on 20 cells and every R1 optimum is re-evaluated on
40 and 80 cells (pre-fit check: max curve difference 20 vs 80 cells ≤ 0.6 % of range even at s = 10⁻³ with
fast contact rates; 40 vs 80 ≤ 0.3 %). Fit A uses 2 starts (the R0 warm start + 1 Latin-hypercube point).
C0 results are taken from round 1 unchanged (6 starts, 40 cells).
