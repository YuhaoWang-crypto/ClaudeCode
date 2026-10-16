# Pre-registered protocol — <probe / dataset>

Frozen in commit <hash> before any route is fitted. Later changes go to "Deviations" with date and reason.

## Question
With the same data and the same electron-transfer assumption, does keeping conformational exchange
predict <technique> current better than a static / fast-equilibrium average?

## Data
- Source(s) with DOI; raw vs digitized; conditions (f, E_sw, step, scan direction, reference, buffer, T).
- Independent units (electrodes, sessions) — and what is NOT independent.
- Frozen preprocessing (baseline windows, smoothing = none, grid).
- Data that are already seen / used before (retrospective) vs genuinely held out.

## Shared assumptions
- Chemistry: <C0 / C1 / C2> with parameterization and bounds.
- Distance factor: k·exp(−β(z − z_min)), β = <main> (sensitivity <values>); z_min = <reporter radius>.
- Shared landscape for all redox species; no double-layer potential dependence; isotropic reporter.
- Current sign and reference electrode.

## Molecular model
- Sampler, force field / CG parameters (each an assumption), friction model, temperature, viscosity.
- Seeds and starts; time-step rule (≤5 % on halving); sensitivity factors (pre-declared values).

## Routes
R0 · R1 (cells, D source, s bounds) · R2 (descriptors, state counts, lag candidates, gates) · static control.

## Fits and evaluation
- Fit A / Fit B definitions; starts per fit; warm-start policy.
- Folds: primary <leave-one-frequency-out>, secondary <leave-one-amplitude-out>.
- Objective; metrics; acceptance (NRMSE ≤ 30 %, |peak| ≤ 30 %); kinetic gain (≥10 % vs R0, same chemistry).
- Diagnostics to report (decomposition, R0*, identifiability) — state whether pre-registered.
- Compute budget and paid-resource policy.

## Labels reported separately
numerically correct · molecular equilibrium checks · molecular dynamics gates · time scale calibrated ·
held-out target met · target-system transfer.

## Deviations
(date — change — reason — confirmation that no current data drove it)
