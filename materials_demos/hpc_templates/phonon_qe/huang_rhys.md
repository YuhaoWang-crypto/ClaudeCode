# Huang-Rhys factor & configuration-coordinate diagram (STE / scintillator)

The single most cost-effective answer to the scintillator reviewers: it turns
"strong electron-phonon coupling" (R3/R7) and "STE formation" (R2/R6) from
qualitative claims into numbers, with only ground- and excited-state geometry
optimisations — no full EPW run required.

## Physical picture
A localised excited state (self-trapped exciton) relaxes to a distorted
geometry. The lattice distortion ΔQ between ground- and excited-state minima
sets the Huang-Rhys factor:

    S = (1/2) * (ΔQ)^2 * ω / ħ

- S ≫ 1  ⇒ strong electron-phonon coupling / strong self-trapping (large Stokes shift)
- total reorganisation energy  λ = S · ħω
- emission energy  E_em = E_ZPL − S·ħω  (Stokes shift ≈ 2·S·ħω)

## Recipe
1. Optimise the **ground state** geometry  → Q_g, energy E_g(Q_g).
2. Optimise the **excited state** (ΔSCF / constrained-occupation DFT in VASP or
   CP2K; or ΔSCF in a cluster with ORCA)  → Q_e, energy.
3. ΔQ from mass-weighted displacement between Q_g and Q_e.
4. Effective phonon frequency ω from the accepting mode (or from the phonon DOS
   of the DFPT run in ph.in).
5. Report S, λ, and the configuration-coordinate diagram (E_g(Q), E_e(Q)
   parabolas) — this is the direct evidence of STE formation the reviewers want.

## Tools
- CarrierCapture.jl / Nonrad — S, capture rates, CC diagrams from the two geometries.
- VASP/CP2K ΔSCF — the excited-state optimisation.
- ph.in (this folder) — the phonon frequencies feeding ω.
