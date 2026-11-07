---
name: glycan-steric-design
description: >-
  Program N-linked glycan composition through protein structure: score glycan
  burial with the two-pass DeltaRASA metric, rank designs on an ordinal
  underprocessing score, and gate candidates on sequon, occupancy and
  glycan-aware structure prediction. Use when designing or analysing a
  glycoprotein whose glycan state matters (oligomannose vs complex), scoring
  glycan shielding or accessibility on any structure, triaging de novo
  glycoprotein designs, planning Endo H / PNGase F / kifunensine validation,
  or extending the glycan_pipeline package. Enforces the same rigorous vs
  hypothesis labelling as network-biomarker.
---

# Steric programming of N-linked glycan composition

A reusable methodology (and a working `glycan_pipeline` package) adapting
Lunn-Halbert et al., *Programming N-linked glycan composition through de novo
protein design* (bioRxiv 2026.10.07.757460, Baker and King labs).

## The core question this answers

N-glycan composition is **not templated by the amino acid sequence**. After
OST transfers Glc3Man9GlcNAc2 onto an N-X-S/T sequon, glycosidases trim it,
and how far that trimming goes depends on whether processing enzymes can
physically reach the glycan. So: **can you encode a glycan's final composition
in the fold around it?** The paper's answer is yes - enclose the glycan and
processing arrests early, leaving oligomannose instead of complex-type.

The unifying object is **DeltaRASA**, a per-atom measure of how much of a
glycan's surface the protein takes away. Everything else ranks, gates or
validates against it.

## The steric ladder - the central empirical result

Progressively walling in one glycan on a modular repeat protein walked its
composition straight down the processing pathway:

| Enclosure | Construct | Glycan observed | Endo H |
|---|---|---|---|
| none, flat exposed surface | THRgly1 | complex, heterogeneous | resistant |
| low, one wall | THRgly2 | **Man8 / Man9** | cleaved |
| medium | THRgly3 | **Glc1Man9** | partial |
| high, fully enclosed | THRgly4 | **Glc1Man9** | partial |

Glc1Man9 is an ER intermediate essentially never seen on mature glycoproteins
(outside freaks like Asn917 of complement C3). Reaching it by geometry alone
is the strongest form of the claim. Note the cost: more enclosure also meant
**lower occupancy** - OST has to reach the site too.

## The four tools + when each applies

| Tool | Question it answers | Module | Rigor |
|---|---|---|---|
| **DeltaRASA** two-pass Shrake-Rupley | how much does the protein shield this glycan? | `g1_delta_rasa` | rigorous, verified vs Biopython to 1e-14 A^2 |
| **Confound control** (size, missing coords) | is my ranking measuring shielding or artefacts? | `g2_confounds` | rigorous |
| **Underprocessing score** (ordinal 1-20) | how far down the processing pathway is this site? | `g3_underprocessing` | engine rigorous; **ladder reconstructed** |
| **Design triage** (sequon / occupancy / confidence gate) | should this design be expressed at all? | `g4_triage` | rigorous, reproduces two paper claims |

## DeltaRASA - get the definition exactly right

```
pass 1  RASA of glycan atoms WITH protein present
pass 2  RASA of the same glycan with PROTEIN ATOMS REMOVED
        (the glycan's own atoms stay - this is the point)
DeltaRASA = |pass2 - pass1|
```

RASA normalises per atom by the isolated-atom maximum `4*pi*(r + r_probe)^2`,
probe 1.4 A. Rank designs on **mean per-atom DeltaRASA** and on the **number
of atoms with non-zero DeltaRASA**; the paper inspects roughly the top 5% by
either, then eyeballs them.

Pass 2 keeps the glycan's own atoms so that self-occlusion cancels and what
survives the subtraction is protein-dependent shielding only. Get that wrong
and the metric means nothing.

## Run the pipeline

```bash
pip install numpy biopython matplotlib
python3 -m glycan_pipeline.run_all            # all four modules + figures
python3 -m glycan_pipeline.g1_delta_rasa      # or any single module
```

Each module exposes `report()` returning a dict and writing a figure to
`figures/`. Structures are cached in `figures/_pdb_cache/`, so only the first
run needs network. Full run is ~45 s.

Score your own structure:

```python
from glycan_pipeline.glycans import load_model, find_glycan_trees
from glycan_pipeline.g1_delta_rasa import site_delta_rasa

model = load_model("my_design.pdb")
rows  = site_delta_rasa(model, find_glycan_trees(model), context="protein")
```

See `assets/score_design.py` for a copy-paste campaign scorer.

## The non-negotiable discipline: honesty labeling

Same two labels as `network-biomarker`:

- **rigorous** - deterministic computation, exact literature reproduction, or
  a real measured statistic.
- **hypothesis / caution** - plausible but unvalidated, reconstructed, or a
  lower bound.

In this package specifically:

- **rigorous**: all SASA/RASA/DeltaRASA geometry; the confound correlations;
  the Fisher reproduction (p = 2.0e-11 vs the paper's "p < 1e-10"); the
  AF2-gate reproduction (only O2 survives, as the paper states); the sequon
  audit (0 of 128 observed glycans sat off a canonical sequon).
- **caution**: the 1-20 underprocessing ladder is **reconstructed, not
  published**; `g3`'s glycan profiles are **illustrative**, built only from the
  paper's qualitative sentences; the neighbour-glycan shielding term is a
  **lower bound**.

Report negatives as findings. `g4` cannot test the NXT-over-NXS occupancy
preference from crystal structures (only 11 unoccupied sequons, p = 1.00) and
says so rather than implying support.

## Common tasks -> where to look

- **The metric in depth** (definition, normalisation, ceiling, what the two
  passes mean) -> `reference/methodology.md`.
- **Run an actual design campaign** (RFdiffusionAA -> LigandMPNN -> RF3/AF3 ->
  DeltaRASA filter, with the paper's exact parameters) ->
  `reference/design-pipeline.md` + `assets/score_design.py`.
- **Prove the glycan state at the bench** (Endo H / PNGase F / kifunensine /
  intact MS / LC-MS/MS decision tree) -> `reference/wet-validation.md`.
- **Score a new structure** -> `g1_delta_rasa.site_delta_rasa`.
- **Swap in real LC-MS/MS data** -> replace `g3.ILLUSTRATIVE_PROFILES`; the
  ladder is overridable via the `ladder=` argument everywhere.

## Hard-won gotchas (these were real bugs or real findings)

- **Mean per-atom DeltaRASA is size-confounded.** Spearman(atoms modelled,
  DeltaRASA) = **-0.38** on 5FYL; a one-sugar stub lies flat on the protein
  and scores high, a big tree dilutes its mean in bulk solvent. Standardising
  to a shared core drops it to **-0.07**. Graft **one identical glycan** (the
  paper uses GlycoShape Man9, or Man6 for AF3-server limits) onto every design
  before ranking. Never rank experimental trees against each other.
- **Use the biological assembly, not the asymmetric unit.** In 5FYL the ASU is
  one protomer, so every neighbouring-subunit glycan is simply absent and
  glycan-glycan shielding reads as zero. `pdbcache.fetch(id, assembly=1)`.
- **Crystallographic glycans are truncated** - median **55%** of a Man9GlcNAc2
  in 5FYL, 45% of sites below half. Any crowding term from deposited
  coordinates is a lower bound, so a crystal structure cannot demonstrate the
  HIV-style "intrinsic mannose patch" mechanism even where it is real.
- **The paper's DeltaRASA excludes neighbouring glycans on purpose.** Its
  designs carry one glycan, so DeltaRASA there means *protein-pocket burial*.
  A screen built on it will not find crowding-based oligomannose sites. Both
  contexts are available; pick deliberately.
- **DeltaRASA cannot reach 100%.** Bonded neighbours already occlude the
  isolated-atom reference. The paper's per-atom colour bar stops at 50%. It is
  a relative shielding score, not "fraction buried".
- **Do not filter glycoprotein designs with a glycan-blind predictor.** AF2
  without a glycan scored O1 and O3 at pLDDT 42.7 / 67.7 and RMSD 21.8 / 16.8 A
  - it reads a glycan-housing pocket as a huge empty cavity. Both were fine.
  Only O2 of the four would have survived the first-generation AF2 gate.
- **Gate on occupancy before composition.** A perfectly shielded site at 5%
  occupancy tells you nothing, and enclosure actively suppresses occupancy by
  blocking OST. Paper gates: >5% to analyse, >50% to select. Bias +2 positions
  to Thr (NXT) for occupancy.
- **Biopython's `protein_letters_3to1` is keyed UPPERCASE** (`"ASN"`, not
  `"Asn"`). Using `.capitalize()` returns zero sequons silently rather than
  raising. This cost a debugging round.
- **Sphere resolution matters more than usual.** DeltaRASA is a difference of
  two similar numbers, so Biopython's default `n_points=100` leaves ~1%
  relative error. The package defaults to 960; symmetry-copy replicates then
  agree to 0.012% (median), which is the real resolution of the metric -
  differences smaller than that are noise.
- **The 1500 underprocessing cutoff is not ladder-robust.** Orderings between
  stochastically-ordered profiles held in 17000/17000 random monotone ladder
  sweeps, but all four O-designs flip pass/fail depending on interior ladder
  spacing. Trust the ranking, not the absolute threshold.
