---
name: allosteric-biosensor
description: >-
  Design single-component allosteric protein switches (biosensors) by circularly
  permuting a small ML/de-novo ligand-binding domain and inserting it into a
  loop of a reporter enzyme (TEM-1 beta-lactamase, de-novo luciferase), then
  validate the constructs in silico with structure+binding prediction (Boltz).
  Reproduces the workflow of Guo, Baker & Alexandrov, "Artificial allosteric
  protein switches with machine-learning-designed receptors", Nat. Biotechnol.
  2026 (doi:10.1038/s41587-026-03081-9). Use when: building a biosensor / protein
  switch for a small molecule, peptide or protein analyte; circularly permuting a
  binder and grafting it into an enzyme reporter; screening a focused
  insertion/permutation library; scoring candidate chimeras for switchability;
  or reproducing/extending the biosensor_pipeline package with a new
  receptor / reporter / analyte. Enforces ✅-rigorous vs ⚠️-hypothesis labeling
  on every claim (a wet-lab kobs titration is the only ground-truth dynamic range).
---

# Artificial allosteric biosensors (ML-receptor protein switches)

A reusable methodology — plus a working `biosensor_pipeline` package — that turns
the paper's wet-lab design recipe into *deterministic, reproducible* sequence
engineering, and adds an *in-silico validation* layer that screens candidate
chimeras with structure prediction. Every result is labeled **✅ rigorous** or
**⚠️ hypothesis**; the two are never blurred.

## The core idea the paper established

You can build an efficient single-component biosensor from a **receptor** that
undergoes *no global conformational change* — as long as it is small and rigid
(exactly what modern ML/de-novo binders are). Ligand binding does not swing a
domain; it **reduces the conformational entropy** of the chimera, which restores
the **reporter enzyme's** catalytic activity (the ON state). This removes the
historical requirement for a ligand-induced macromolecular conformation change
and makes almost any ML binder a candidate receptor.

## The design recipe (paper Methods, "Design of chimeric proteins")

Four deterministic steps — implemented exactly in `biosensor_pipeline/design.py`:

1. **Pick a permutation site** — an arbitrary residue *inside a loop* of the
   receptor. (We pick loops rigorously from the deposited structure via
   `annotate_sse`; ✅.)
2. **Circularly permute** — remove that residue to create new N/C termini;
   join the *native* N and C termini with a flexible Gly/Ser linker "of
   sufficient length". → `cpReceptor = R[site+1:] + GSlinker + R[:site]`.
3. **Insert into the reporter** — graft the cpReceptor into a permissive surface
   loop of the reporter enzyme (TEM-1 β-lactamase **position 253**; or 196/197
   for logic gates), a single glycine on each side.
   → `chimera = reporter[:q] + G + cpReceptor + G + reporter[q:]`.
4. **Screen a tiny library** — because ML binders are small, the permutation
   library is "fewer than ten variants"; assay each for ligand-dependent
   activity and keep the switchable ones.

Extensions the paper demonstrates (supported as parameters / patterns):
YES gate (duplicate receptor at loops 41 & 197 → larger dynamic range),
AND gate (two orthogonal receptors at 41 & 197), auxiliary binding domain fused
by a linker to offset the circular-permutation affinity penalty, and alternative
reporters (de-novo luciferase LuxSit Pro, NanoLuc, glucose dehydrogenase).

## Reproduction status in this repo

| System | Receptor (PDB) | Analyte | Reporter | Role |
|---|---|---|---|---|
| `dig`   | DIG10.3 (4J9A) | digoxigenin (steroid glycoside) | TEM-1 / BLA-253 | reproduction |
| `dfhbi` | mFAP1 (6CZI)   | DFHBI (fluorogen) — a different analyte, unrelated fold | TEM-1 / BLA-253 | validation |
| `vitd`  | CDL2.2 (5IEN)  | vitamin D3 (secosteroid) — receptor **auto-mined** by `discover.py` | TEM-1 / BLA-253 | validation (mined) |

All use **real public sequences** and the **identical** recipe. In-silico binding
retention after circular permutation + insertion: DFHBI 1.01, digoxigenin 0.92,
vitamin D3 0.74 → **1.07 after site×linker tuning** (see reference/report).

## The full workflow (mine → graft → optimize → validate)

```bash
pip install biotite numpy matplotlib            # structure-derived steps + figures

# 0) MINE a receptor for an analyte you only have as a molecule (Tier-A PDB)
python3 -m biosensor_pipeline.discover "<SMILES>" --name X --max-len 200
python3 -m biosensor_pipeline.discover_batch     # a whole analyte panel -> receptor_catalog.json

# 1) DESIGN — deterministic, 100% reproducible  ✅
python3 -m biosensor_pipeline.run_repro          # build+verify libraries -> biosensor_out/
python3 -m biosensor_pipeline.structure          # recompute & verify loop annotations
python3 -m biosensor_pipeline.test_design        # correctness + reproducibility tests

# 2) VALIDATE in silico via Boltz MCP (paid, ~$0.05/prediction)
#    submit holo/apo/control (see reference/boltz-validation.md), then:
python3 -m biosensor_pipeline.analyze_boltz      # constellation + retention + proxy
python3 -m biosensor_pipeline.make_figure        # figures/biosensor_validation.png

# 3) OPTIMIZE — tune the receptor<->reporter linker and/or permutation site
python3 -m biosensor_pipeline.tune_linker vitd   # flank-linker series -> Boltz each

# 4) MECHANISM / beyond DR — Boltz can't give dynamic range; probe the coupling
python3 -m biosensor_pipeline.coupling           # apo vs holo active-site ordering (free)
python3 -m biosensor_pipeline.md_entropy         # OpenMM MD flexibility/entropy (CPU smoke)
# 4b) PRODUCTION MD on a cloud GPU (apo/holo ΔS + ligand parameterization):
modal run biosensor_pipeline/modal_app.py --apo-pdb ... --holo-pdb ... --ligand-sdf ... --smiles ...
```

The switch **dynamic range** is kinetic/entropic — no platform computes it
directly (see `reference/beyond-boltz.md`). `coupling.py` reads apo→holo
active-site pLDDT (allosteric-coupling hint, free); `md_entropy.py` runs OpenMM
MD for the physical conformational-entropy route (needs GPU + µs to converge).

Module map: `design.py` (circular permutation + insertion, ✅ exact) ·
`systems.py` (receptors/reporters/analytes) · `structure.py` (annotate_sse loop
sites; tag/terminal-mismatch tolerant) · `screen.py` (focused library) ·
`scoring.py` (✅ geometric metrics + ⚠️ switch proxy) · `boltz_io.py` (payloads) ·
`discover.py`/`discover_batch.py` (PDB receptor mining, contact-verified) ·
`tune_linker.py` (linker series) · `analyze_boltz.py` (results → scores) ·
`coupling.py` (apo/holo active-site ordering) · `md_entropy.py` (OpenMM MD, incl. production GPU apo/holo ΔS + ligand FF) · `modal_app.py` (cloud-GPU runner).

## What is rigorous vs. hypothesis (read this before quoting any number)

- **✅ rigorous**
  - The chimera **sequence construction** — deterministic string surgery,
    round-trip-verified (`design.verify_chimera`): reporter preserved, every
    receptor residue conserved except the one removed.
  - **Loop / permutation-site selection** — from `annotate_sse` on the deposited
    PDB (reproducible with `structure.py`).
  - **Catalytic-residue location** — TEM-1 S70/S130 etc. found by unambiguous
    sequence motifs (`STFK`, `SDN`), cross-checked vs Ambler numbering.
  - **Geometric measurements on a predicted model** — e.g. Cα–Cα distances of
    the catalytic constellation; a real number computed from atoms.
  - **Boltz confidence/affinity numbers** — as reported by the model.

- **⚠️ hypothesis**
  - The **"switch proxy"** score that ranks a library. Structure prediction does
    **not** compute an allosteric dynamic range. The proxy is an *illustrative,
    transparent* triage heuristic — never a predicted DR.
  - Interpreting an intact active-site constellation as "catalytically ON", or an
    apo→holo confidence shift as "allosteric coupling". Plausible, not proven.
  - Any claim that a construct "works" as a biosensor. **Only a wet-lab
    kobs(+ligand)/kobs(−ligand) titration** (paper, "Quantification of biosensor
    performance parameters") is ground truth.

Negative/ambiguous in-silico results are reported as findings, not hidden. A
low-confidence Boltz model or a disrupted active site is a *useful* result: it
deprioritizes that variant before any bench work — exactly the triage value the
paper gets from its small library.

## Extending to a new analyte / receptor / reporter

See `reference/adding-a-system.md`. In short: add a `Receptor` (sequence +
ligand SMILES + structure-derived `loop_sites`) and/or a `Reporter` (sequence +
catalytic residues + insertion loops) to `systems.py`, then
`build_library` + the Boltz validation work unchanged. This is the whole point:
the recipe is analyte-agnostic.

## References

- `reference/methodology.md` — the paper's mechanism, the recipe, logic gates,
  and the honesty-labeling contract in detail.
- `reference/boltz-validation.md` — how the in-silico screen maps design
  questions onto Boltz structure+binding predictions, and how to read them.
- `reference/adding-a-system.md` — plug in a new receptor/reporter/analyte.
- `reference/beyond-boltz.md` — computing dynamic range & the readout: what
  QM/MM, MD, FEP/MM-GBSA each give (and why none returns DR directly).
- `reference/md-gpu-protocol.md` — the production apo/holo MD + ligand
  parameterization protocol, env, and Modal cloud-GPU quickstart.
- `assets/system_template.py` — copy-paste skeleton for a new system.
