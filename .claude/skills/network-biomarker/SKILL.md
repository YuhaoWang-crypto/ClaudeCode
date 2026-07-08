---
name: network-biomarker
description: >-
  Apply the irreducibility/symmetry + dynamical-systems biomarker pipeline
  (graph automorphism → CRNT deficiency → elementary flux modes →
  critical-slowing / Lyapunov / bifurcation) to a gene-regulatory or metabolic
  network, to find irreducible core modules and detectable early-warning
  biomarkers. Use when analyzing a pathway or network for core modules,
  bistability / switch capacity, or tipping-point biomarkers; extending the
  grn_pipeline package with a new pathway; reproducing an exact literature
  model (BioModels); or producing the honesty-labeled summary report. Enforces
  rigorous ✅-rigorous vs ⚠️-hypothesis labeling on every claim.
---

# Network irreducibility & critical-transition biomarker analysis

A reusable methodology (and a working `grn_pipeline` package) that turns four
abstract mathematical tools into *computed* numbers on concrete, literature-
grounded biological networks — and always labels what is rigorous vs.
hypothetical.

## The core question this answers

Given a regulatory / metabolic network, find (a) its **irreducible core
modules** (via symmetry / fibration and flux irreducibility) and (b)
**model-free, measurable early-warning biomarkers** of a critical transition
(cell-fate flip, toxicity, resistance). The unifying object is the **leading
eigenvalue = largest Lyapunov exponent → 0 at a bifurcation**, with variance,
lag-1 autocorrelation and the DNB index rising alongside.

## The four tools + when each applies

| Tool | Question it answers | Module | Rigor |
|---|---|---|---|
| Graph automorphism → quotient; input-tree **fibration** | which nodes are redundant / collapse to one core node? | `m1_symmetry`, `m11_fibration` | ✅ deterministic |
| **CRNT deficiency** δ = n − ℓ − s | can this topology be a switch at all? (δ=0 weakly-reversible ⇒ never bistable) | `m2_crnt`, `m12_dualphos` | ✅ linear algebra + ODE |
| **Elementary flux modes** | what are the irreducible flux generators (the "prime factors")? | `m3_efm` | ✅ extreme rays + nnls |
| **Critical slowing / Lyapunov / DNB** | how far from tipping? what's the biomarker? | `m4_dnb_lyapunov`, `m16_erk_dnb` | ✅ signals; ⚠️ "λ-flip ⇒ toxicity" |

## Bifurcation taxonomy — read the RIGHT signature

The single most important operational rule. Three classes of dynamical
instability, three distinct early-warning signatures:

- **Saddle-node (switch / bistable):** variance ↑, autocorrelation ↑, recovery
  time τ ↑, **no oscillation**. → read variance + AR1. (`m19`, `m15/m16`)
- **Hopf (clock / oscillator):** variance ↑ **and a spectral peak sharpens at a
  fixed intrinsic frequency** (amplitude → 0). → read the **power spectrum**,
  not autocorrelation. (`m21`)
- **SNIC (excitable / cell-cycle, saddle-node ON a limit cycle):** finite-
  amplitude spikes whose **period → ∞ (frequency → 0)**; ISI mean + CV both
  grow. → read the **period / ISI trend**. (`m22`)

## Run the existing pipeline

```bash
pip install numpy scipy networkx matplotlib
# optional, for exact BioModels reproduction: pip install libroadrunner
python3 -m grn_pipeline.run_all          # full pipeline + all figures
python3 -m grn_pipeline.m15_markevich_mm  # or any single module
```

Every module exposes `report()` returning a dict of its computed results and
writing a figure to `figures/`. `run_all.py` aggregates them. See `README.md`
(module table) and `REPORT.md` (full Chinese write-up with all data tables).

## The non-negotiable discipline: honesty labeling

Every result carries one of:

- **✅ rigorous** — deterministic computation, decimal-exact literature
  reproduction, or a real measured statistic.
- **⚠️ hypothesis** — a plausible but unvalidated claim needing experimental
  calibration (e.g. "λ crosses zero ⇒ apoptosis/toxicity", illustrative
  couplings, approximate PK values).

Never blur the two. Negative and partial results (e.g. M17: real single-cell
ERK data was all supra-threshold, so the biomarker could not be decisively
confirmed) are reported as findings, not hidden. When a prior claim is
overturned by data (e.g. M10 overturned M7's ranking), say so explicitly.

## Common tasks → where to look

- **Add a new pathway** (bistable switch, oscillator, or mixed) →
  `reference/adding-a-pathway.md` + `assets/pathway_template.py`. Reuse the
  generic engines instead of re-deriving: `m19.bistable_window` / `m19.titration`
  (1-D switches), `m21.find_hopf` / `m21.power_spectrum` (oscillators),
  `m4.benettin_lle` (validated Lyapunov), `m17._ar1` (early-warning stat).
- **Reproduce an EXACT literature model** (decimal-exact, not just topology) →
  `reference/data-access.md`. The working recipe behind the agent proxy is the
  **biomodels GitHub mirror + libRoadRunner**; `biomodels.org`/EBI are blocked.
- **The math in depth** (definitions, theorems, gotchas that were real bugs) →
  `reference/methodology.md`.
- **Produce the summary report** → build a self-contained themed HTML report
  (see `figures/summary_report.html` as the reference implementation) that
  weaves the narrative arc + the computed data tables + a rigor ledger.

## Hard-won gotchas (these were real bugs)

- CRNT δ must be computed on the **effective** network: strip chemostatted
  species from the complexes, or δ is wrong (Schlögl looked like δ=0).
- FIM observable ranking must use **dimensionless relative** sensitivity, or
  concentrations win on units alone and ratio-biomarkers wrongly lose.
- For stationary titration data use **mean-subtraction**, not rolling detrend —
  rolling detrend removes the slow critical-slowing fluctuation you're measuring.
- On an OFF branch where the mean → 0, CV explodes spuriously; use **absolute
  residual SD**.
- Track the **least-stable / disappearing** branch (max λ) into a saddle-node,
  not the persistent one.
