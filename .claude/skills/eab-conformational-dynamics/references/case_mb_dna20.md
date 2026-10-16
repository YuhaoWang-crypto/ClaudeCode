# Worked case: 20-nt DNA–methylene blue, SWV (Dauphin-Ducharme et al., Langmuir 2017)

Source: doi:10.1021/acs.langmuir.7b00359, SI Fig. S-7. 5′-HS-(CH₂)₆-ATTATTTTTTATTTATTTTT-MB-3′ on a
6-mercaptohexanol SAM, 2 mm Au, 0 → −0.45 V vs Ag/AgCl, 1 mV step, f = 20/100/500 Hz × E_sw =
25/50/75 mV. Nine digitized curves (forward/backward/net), presumably one electrode. Digitization floor
≈ 3 % NRMSE (net vs forward − backward consistency 1.7–4.4 % RMS).

## Molecular model (CG Brownian dynamics, dt 0.5 ps, 4 seeds × 64 replicas × 1 µs)
- ⟨exp(−5(z − 0.35))⟩ = 0.022; P(contact within 0.2 nm) ≈ 2 %; τ_z ≈ 20 ns.
- R1 slowest relaxation 5.3 ns (40 and 80 cells agree); R2 slowest implied timescale ~15 ns.
- Half-pulses 1–25 ms → separation 10⁵–10⁶.
- R2 failed dynamic gates for 4/8/12 states, lags 10 ps–10 ns; equilibrium blocks passed (PMF ≤ 0.22 k_BT).
- Linker sensitivity (0.6/1.8 nm, l_p 0.7/2 nm): ⟨w⟩ 0.020–0.026, relaxation 5.1–5.9 ns.

## Round 1 — concerted 2e⁻ BV (4 parameters)
- R0/R1/R2 identical within 0.1 % in every fold, fit A/B, β = 2/5/10.
- Leave-one-frequency-out NRMSE 23.6 % (2/9 curves meet target); leave-one-amplitude-out 16.7 % (4/9).
- k_eff ≈ 6 s⁻¹, α = 0.35, E° = −0.257 V; s profile flat over 10⁻¹–10³.
- Single-rate BV cannot match the anodic backward current at 500 Hz and the large currents at 20 Hz.

## Round 2 — stepwise chemistry
| chemistry | params | train cost | LOFO NRMSE (R0 / R0*) | LOAO NRMSE (R0 / R0*) | curves OK (LOFO / LOAO) |
|---|---|---|---|---|---|
| C0 | 4 | 0.0245 | 23.6 % | 16.7 % | 2/9 · 4/9 |
| C1 EE | 6 | 0.0094 | 61.6 % / 57.0 % | 9.0 % / 9.1 % | 1/9 · 8/9 |
| C2 E–C–E (K_p 10) | 7 | 0.0049 | 19.6 % / 17.0 % | 10.6 % / 6.9 % | 2/9 · 8/9 |
- C2 parameters: k₁⟨w⟩ ≈ 70–280 s⁻¹, k₂⟨w⟩ ≈ 360 s⁻¹, ΔE ≈ −0.13…−0.17 V (potential inversion),
  k_p ≈ 600–1700 s⁻¹. K_p not identified (1/10/100 → 24/20/25 % LOFO).
- Raw R1/R2 "gains" (13–67 %) were basin effects: R0 at the same chemistry point reproduced them;
  against R0* no consistent gain. Real dynamics differences only with s near 10⁻³, mixed sign.
- C1 at held-out 500 Hz: near-optimal training points gave 17 %–126 % held-out error (not identified).
- Static two-population control under C2: 15.6 % LOFO (5/9), 6.1 % LOAO — heterogeneity is used by the
  data, but with one extra parameter and nine curves it is not a validated mechanism.
- Held-out 20 Hz under-predicted by 43–73 % peak in every chemistry: frequency extrapolation unsolved.

Compute: ~16 CPU-h total on 4 cores.
