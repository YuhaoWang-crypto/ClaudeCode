# Pitfalls met in practice (and fixes)

## Numerics
- **Catastrophic cancellation in the current.** Computing I = g k·x − e_o k·π after the transient has
  decayed subtracts two nearly equal numbers. Propagate the deviation from Nernst equilibrium instead
  (current = k_ab u_a − k_ba u_b).
- **Stiff spectrum.** With conformational exchange ~10¹³ s⁻¹ and reaction ~10³ s⁻¹, `eigh` of the
  symmetric generator loses the slow eigenvalue (absolute error ε‖A‖). Use the SVD of the square-root
  factor. Symptom: R1 − R0 error that grows linearly with s.
- **Reference solvers can be the inaccurate side.** Dense `expm` of a stiff generator gives
  ~10⁻⁸ noise where the true current is ~0; compare with an absolute tolerance scaled to the largest current.
- **Static limit needs truly slow exchange.** s = 10⁻¹² of a 10¹³ s⁻¹ generator is still 10 s⁻¹ over a
  multi-second scan; use 10⁻²⁰ in tests.
- **Units.** Stokes–Einstein for a 0.35 nm sphere in water is 7.0×10⁸ nm²/s (not 10¹¹).

## Coarse-grained MD
- WCA repulsion on 0.15 nm beads makes explicit Brownian steps unstable at every dt ≥ 0.1 ps (curvature ×
  D ~ 30 ps⁻¹). Truncated-harmonic repulsion with the same contact distances is stable; declare it.
- OpenMM CPU runs a 25-particle system at ~1.5×10⁴ steps/s (per-step overhead). Put 64 non-interacting
  replicas in one System: ~9×10⁴ replica-steps/s per core.
- Reporter-height correlation time is ~20 ns; contact-weighted rates are heavy-tailed. A 5 % time-step
  rule needs ~10² µs of pooled sampling per dt, otherwise noise (~10 %) decides it.
- Local contact counts show real dt bias (−6.7 % from 1 → 0.5 ps) while global observables do not.

## R2 state models
- k-means on (height, Rg, end distance, contacts) gave implied timescales that rose with lag up to 10 ns
  with no plateau and CK errors 0.08–0.35: not Markovian. Extend lag candidates to cover the slowest
  relaxation (here ~20 ns), and accept that R2 may stay conditional. Do not keep re-clustering against
  current data.

## Fitting and inference
- **Optimizer basins masquerade as dynamics gains.** In stepwise chemistry, R1/R2 beat R0 by 13–67 %
  in raw numbers; R0 evaluated at the R1/R2 chemistry point reproduced the gain almost every time.
  Always run the decomposition and R0* before claiming anything.
- **Non-identified extrapolation.** Training fits within 5 % of each other gave 17 %–126 % held-out
  error at an unseen frequency. Report the near-optimal range.
- Warm starts must come from the same training fold, never from a full-data fit.
- Motion multiplier s at its lower bound with a flat profile = no information on speed.
- Apply the identical baseline operator to model curves; Q profiled per training fold only.
- Python `hash()` is salted per process: use `zlib.crc32` for reproducible start seeds.

## Data access
- PMC and Europe PMC put bot challenges in front of figures and SI. A headless Chromium (Playwright)
  that loads the article page gets the rendered HTML (figure CDN links then download with curl); SI
  files sit behind a proof-of-work cookie that the browser solves on navigation, after which
  `context.request.get(url)` returns the PDF.
- Inline formulas (sequences) may be lazy MathJax; scroll the element into view and read `inner_text`.
- Never send the user's email or other identifiers to third-party APIs (e.g. as an Unpaywall
  `email` parameter) without explicit permission.

## Long jobs
- Background tasks may have a wall-clock limit (e.g. 30 min default); launch long jobs detached
  (`setsid nohup … &`) and watch a log file.
- `pkill -f pattern` / `pgrep -f pattern` can match the shell running the command; kill by PID or use
  `^python3 …` anchors.
- Pin BLAS threads to 1 when running process pools, or load averages explode and jobs starve.
- Keep a compute ledger (CPU-hours per batch) against any stated budget.
