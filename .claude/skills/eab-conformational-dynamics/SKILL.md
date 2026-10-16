---
name: eab-conformational-dynamics
description: >-
  Test whether conformational dynamics of a surface-tethered redox-reporter probe (E-DNA or E-AB
  aptamer sensors; methylene blue or ferrocene on DNA, peptide or PNA) predict pulsed
  electrochemical current (square-wave, AC, cyclic voltammetry) better than a static rate
  average. Bundles a tested exact current backend for conformation x redox master equations
  (concerted 2e-, stepwise EE and E-C-E proton-coupled chemistry), OpenMM coarse-grained
  Brownian dynamics, Smoluchowski and Markov-state generators with validation gates,
  pre-registered grouped cross-validation, and checks that separate real dynamics gains from
  optimizer artifacts. Use it whenever someone models electron-transfer kinetics of
  electrode-bound DNA/aptamer probes, asks whether probe flexibility or chain dynamics explains
  SWV frequency/amplitude behaviour, fits k0/alpha/E0 to voltammograms with a structural model,
  or wants MD/MSM-informed sensor signal prediction - even without naming the method.
---

# Conformational dynamics vs static averaging in pulsed electron transfer

A reusable, pre-registration-first workflow for one question:

> With the same data and the same electron-transfer assumption, does keeping conformational
> exchange predict the current better than a static / fast-equilibrium average?

It was built and validated on a real benchmark (20-nt DNA–methylene-blue, 9 SWV curves; see
`references/case_mb_dna20.md`). The answer there was **no** — half-pulses were 10⁵–10⁶× longer than
conformational relaxation — while a stepwise MB chemistry improved fits. Expect that outcome often,
and design the analysis so a negative result is clean and credible rather than something to tune away.

Report to the user in their own language. Keep every claim labelled by what supports it (see
"Status labels").

## Scripts (in `scripts/`, run from that directory or add it to `sys.path`)

| file | role |
|---|---|
| `swv_backend.py` | Waveform; `Kinetics` (π, w, reversible generator L); exact SWV propagation for concerted 2e⁻ BV; `fast_eq_kinetics` (R0), `smoluchowski_kinetics` (R1), `two_population_kinetics` (static control) |
| `chem_backend.py` | Same propagator for conformation × multi-species chemistry: `Chemistry("E1"|"C0"|"C1"|"C2", K_p=...)` = one-electron (ferrocene) / concerted 2e⁻ (MB apparent) / EE / E–C–E |
| `cg_md.py` | OpenMM CPU Brownian dynamics of a tethered probe (1 bead/nt, linkers, wall), many independent replicas per context; `--n_dna`, `--r_reporter`, `--n_link`, `--lp_link`; writes height, Rg, end distance, contacts |
| `msm.py` | k-means states, reversible continuous-time generator by lag-τ likelihood (analytic gradient), CK test, ESS / committed entries |
| `fit_tools.py` | data loader, frozen baseline step, Q profiling, objective, metrics, grouped folds, LHS starts, multi-start fit, best-basin reference R0* |
| `timescale_check.py` | **triage**: relaxation time vs half-pulse, and \|R1 − R0\| at s = 1 and 10⁻³ |
| `tests/` | 11 numerical tests (expm, stiff Radau, exchange limits, detailed balance, reversible limits, CTMC recovery). Run `python -m pytest tests -q` before trusting any change |

Dependencies: numpy, scipy, matplotlib; openmm only for `cg_md.py`. Pin BLAS to one thread
(`OPENBLAS_NUM_THREADS=1`) when running parallel fits.

## Workflow

### 0. Scope and pre-register (before touching data)
Write `PROTOCOL.md` from `references/protocol_template.md` and commit it before any fit: routes,
chemistry, assumptions (distance decay β, contact height), folds, objective, metrics, acceptance
thresholds, what counts as a kinetic gain, and which results are retrospective. Every later change
goes into a dated "deviations" section with the reason and the confirmation that no current data
drove it. This is what makes a negative or surprising result believable.

### 1. Data
- Prefer raw instrument tables. If only figures exist, digitize the **experimental** traces only
  (never the authors' fitted curves), keep the source image, the axis-calibration map, a per-point
  reading error, and flags for occluded/imputed points. Check internal consistency
  (net ≈ forward − backward) — it gives the digitization floor; differences below it are not
  interpretable.
- Record conditions exactly (frequency, amplitude, step, scan direction, reference electrode, buffer,
  temperature). Note when curves are one electrode (no replicate statistics).
- Freeze one preprocessing operator (e.g. linear edge baseline) and apply it identically to data
  **and** model curves.
- Paywalled/closed data: say "insufficient", never substitute author simulations or another reporter.

### 2. Triage — run `timescale_check.py` first
It needs an equilibrium height distribution of the reporter (from step 3, or a quick short CG run).
If half-pulses exceed conformational relaxation by ≫10³ and |R1 − R0| at s = 1 is below the noise
floor, the routes will converge. Tell the user that prediction up front, then still run the
comparison (it is cheap) but put the effort into chemistry and identifiability.

### 3. Molecular model (`cg_md.py`)
- Defaults (all assumptions to declare): 1 bead/nt, bond 0.6 nm, 100 k_BT nm⁻², l_p 1 nm,
  DNA bead radius 0.25 nm, reporter radius 0.35 nm, C6 arm 2×0.3 nm, reporter linker 3×0.3 nm,
  truncated-harmonic excluded volume and wall (WCA on small linker beads is unstable at Brownian
  steps ≥0.1 ps), free-draining Stokes friction, no surface attraction, no base pairing.
- Use many replicas per context (`n_rep=64`): single small systems run at only ~1.5×10⁴ steps/s on
  CPU, replicas give ~6× throughput. ≥3 seeds (contexts) with mixed extended/collapsed starts;
  seeds 1–2 build, seed 3 validates.
- Time-step rule: fastest dt whose halving changes every equilibrium observable and correlation time
  by ≤5 %. Size runs so the standard error is well below 5 % (contact-weighted rates need ~10²–10³ µs
  of pooled sampling); otherwise the rule cannot be decided.
- Pre-declared one-factor geometry sensitivity (linker length, linker l_p); never pick by test error.

### 4. Routes (all share the same equilibrium distribution and the same current backend)
- **R0** fast equilibrium: `fast_eq_kinetics(pi, w)` with w = ⟨exp(−β(z − z_min))⟩ over the full
  distribution — never the exponential of a mean height.
- **R1** 1-D Smoluchowski on z with D from Stokes–Einstein × global multiplier s.
- **R2** Markov-state generator from unbiased trajectories (`msm.py`): pre-declare descriptors and
  state counts, choose the smallest passing count and shortest passing lag by trajectory gates only
  (timescales stable within 20 % over τ, 2τ, 4τ; CK error ≤ 0.10 on a held-out seed; ESS ≥ 200 and
  ≥ 50 committed entries in the states carrying 95 % of the current). If none pass, R2 is reported as
  a conditional prediction and not ranked. Never use a matrix logarithm or fit edge rates to current.
- **Static control** with the same parameter count (two non-exchanging populations, rates k·e^{±δ}).
- Chemistry is chosen independently of routes; compare routes only within a chemistry.

### 5. Fits and cross-validation (`fit_tools.py`)
- Fit A: one shared θ across routes (s = 1) — tests the information in dynamics under common calibration.
- Fit B: each route calibrated separately, plus one global s for R1/R2 (never per condition).
- Leave-one-group-out folds (frequency primary, amplitude secondary); Q profiled on training curves;
  ≥4 starts, including a warm start from the same-fold R0 solution (never from a fit that saw the
  held-out curves).
- Metrics per held-out curve: NRMSE of net (RMSE / observed peak-to-valley), forward/backward NRMSE,
  net peak error; target ≤30 %; kinetic gain = ≥10 % relative reduction vs R0.

### 6. Diagnostics that must accompany any route "gain"
Apparent gains are often optimizer artifacts, because extra routes bring extra starts into basins
that R0's own starts missed. For every R1/R2 solution:
1. Re-evaluate R0 at the **same chemistry parameters** (drop s) and the route at s = 1. If either
   reproduces the route's held-out error, the gain is a basin effect, not dynamics.
2. Compute **R0\*** with `best_basin_reference`: the lowest-training-cost chemistry point found by
   any start of any route, evaluated as R0. Compare routes against R0\*.
3. Report the held-out range among near-optimal points (within 5 % of best training cost). A wide
   range means the held-out prediction is not identified by the training data.
4. Flag s at its bound; a dynamics effect that appears only at s = 10⁻³ means motion 1000× slower
   than Stokes–Einstein.
Label these diagnostics as post hoc if they were not pre-registered.

### 7. Chemistry alternatives (`chem_backend.py`)
Reporter chemistry often dominates the misfit. C1 (EE: two 1e⁻ steps, potentials E_m ± ΔE/2) and C2
(E–C–E: semiquinone protonation with rate k_p and fixed K_p, with sensitivity) are drop-in. Both
electron steps take the distance factor; protonation does not. Show parameter counts beside every
comparison; an unidentified constant (e.g. K_p) makes a chemistry a candidate, not a validated mechanism.

### 8. Report
- Lead with the answer, then the time-scale argument, then tables per fold and route, then the
  decomposition, then gates and limitations.
- Status labels, kept separate: numerically correct · molecular equilibrium checks · molecular
  dynamics gates (R2) · time scale experimentally calibrated · held-out target met · target-system
  transfer. Use "not executed" / "insufficient data" rather than filling gaps.
- For a shareable page, publish an HTML artifact; `assets/report_page_example.html` is a worked
  example (interactive held-out chart, decomposition table, status pills).

## References
- `references/methodology.md` — equations, generator conventions, the exact propagator, gates, objective.
- `references/pitfalls.md` — failure modes met in practice and their fixes (numerics, MD throughput,
  statistics, optimization, data access, long-running jobs). Read before running anything long.
- `references/protocol_template.md` — pre-registration skeleton.
- `references/case_mb_dna20.md` — the worked benchmark with numbers, for calibration of expectations.
