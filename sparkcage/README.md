# SparkCage

Open-source feasibility models for a single-component, electrochemically read
de novo protein-switch biosensor: a LOCKR/LucCage-style cage-latch carrying a
methylene blue reporter on a carbon electrode.

Full analysis, including the software-substitution table and the risk
assessment, is in **[FEASIBILITY.md](FEASIBILITY.md)** (written in Chinese).

## Run it

```bash
pip install numpy scipy matplotlib
python3 sparkcage/run_demo.py     # figures + numeric tables, < 1 min on 4 CPUs
python3 sparkcage/validate.py     # 16 checks against published values / exact limits
```

No GPU, no licensed software, no Rosetta, no COMSOL.

## What each module is for

| Module | Replaces | What it does |
|---|---|---|
| `thermo/switch_model.py` | Rosetta ddG + hand algebra | Keyless three-state switch model with closed-form EC50 and dynamic range |
| `thermo/luccage_reference.py` | — | Exact reimplementation of the published LucCage ten-equation coupled-equilibrium model |
| `echem/swv.py` | COMSOL Electroanalysis | Square-wave voltammetry of a surface-confined redox reporter, integrated analytically per half-cycle |
| `echem/transport.py` | COMSOL Transport of Diluted Species | 1D diffusion to the electrode with Langmuir capture, plus the transport-regime diagnostic |
| `structure/analyze.py` | PyMOL / Rosetta analysis | Cage-latch contact topology and rigid-body hinge angles, numpy only |

## Three findings worth reading the report for

1. **COMSOL is the wrong tool here, not merely a replaceable one.** A
   surface-confined redox species has no spatial dimension, so its voltammetry
   has an exact closed-form solution. The paper this project would replicate
   (Dauphin-Ducharme et al., *Langmuir* 2017) had to approximate the monolayer
   with a 150 µm thin-layer cell because COMSOL cannot represent it directly.

2. **Dropping the key costs 42× in sensitivity.** LucCage is a two-component
   system whose key binding energy pays part of the cage-opening cost. A
   single-component electrochemical design loses that, and needs its own
   compensating interaction. The latch-electrode adsorption the quote already
   budgets molecular dynamics for is the natural candidate.

3. **Almost all of the conformational change is wasted.** The measured 4–30%
   signal change of real protein E-AB sensors implies the reporter moves only
   1–3 Å, against an 11–21 Å domain motion. Where the methylene blue is
   attached matters more than how large the protein's motion is.

## Honesty note

Trends, scalings and design rankings from these models are reliable. Absolute
numbers (currents, EC50, response times) depend on calibration parameters
(`k0_contact`, `β`, `Γ`, `σ`) and must be anchored to a measured voltammogram
before being quoted. Use the models to rank designs, not to predict an
absolute limit of detection.
