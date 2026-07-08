# Adding a new pathway / network

The pipeline is designed to migrate. To add a pathway, pick its dynamical class,
reuse the matching generic engine, and follow the module conventions. Do **not**
re-derive Lyapunov / early-warning code — it already exists and is validated.

## Module conventions (every `mXX_*.py` follows these)

1. A top docstring stating the system, what's rigorous, and what's hypothesis.
2. A `report(fig_path=...)` function that:
   - prints a self-contained, human-readable summary,
   - computes every number (no asserted constants),
   - writes one figure to `figures/`,
   - returns a `dict` of key results.
3. A `_figure(...)` helper (matplotlib, `matplotlib.use("Agg")`).
4. Wire it into `grn_pipeline/run_all.py` (guard with try/except if it needs
   network or an optional dependency like libroadrunner — see the m17 / m20b
   blocks).
5. Add a row to `README.md`'s module table and a section to `REPORT.md`.

## Pick the class → reuse the engine

### Bistable switch (saddle-node)
Write the 1-D rate law `f(x, p)` with a control parameter `p`, then reuse
`m19_switch_library`:
- `bistable_window(f, prange, xmax)` → the fold edges,
- `titration(model, ...)` → Langevin traces + variance/AR1 across the window,
- `n_stable(f, p, xmax)` / `_roots(...)` → attractor count.
Expect: variance + lag-1 autocorr peak at the saddle-node. Track the
**disappearing** branch (max λ).

### Oscillator (Hopf)
Write the ODE system `f(state, p)` and reuse `m21_oscillators`:
- `find_hopf(model)` → the parameter where Re(λ)→0, Im(λ)≠0,
- `power_spectrum(model, p, ...)` → the spectral peak (use `scipy.signal.welch`),
- `_lead_eig(f, x0, p)` → leading eigenvalue.
Expect: variance ↑ **and** a sharpening spectral peak at the intrinsic
frequency. Read the spectrum, not autocorrelation.

### Mixed / excitable (SNIC)
Use the θ normal form in `m22_snic_mixed` as the template: `period(I)` (spike
period from θ crossing 2πk), `spikes_isi(I, sigma)` (noisy ISI). Expect period →
∞ (log-log slope −0.50) and ISI mean + CV growing toward the bifurcation.

## Shared validated helpers (import, don't reinvent)

- `m4_dnb_lyapunov.benettin_lle(rhs, jac, x0)` — LLE (Rössler-validated).
- `m4_dnb_lyapunov.stochastic_signals(mu, x_star, ...)` — SD / AR1 / DNB.
- `m17_realdata._ar1(series)` and `._rollmean(x, w)` — early-warning stats used
  on real data; reuse them so simulated and real analyses share one code path
  (this is exactly what made M18 a valid positive control for M17).
- For an exact literature model, use `m20b_biomodels_exact.fetch_sbml(id)` +
  `_rr(path)` and read `rr.getFullJacobian()` for λ_max — see `data-access.md`.

## Rigor labeling for the new module

State, in the docstring and the printed report, which parts are:
- ✅ deterministic computation / decimal-exact reproduction / measured statistic;
- ⚠️ mechanistic-but-illustrative parameterization or an untested hypothesis.

If you use a canonical mechanistic form (positive feedback / covalent cycle /
substrate depletion) rather than fitted literature rate constants, say so — the
shared *biomarker geometry* is the rigorous claim, not the exact parameters.
To make a module decimal-exact, fetch the curated model (`data-access.md`).
