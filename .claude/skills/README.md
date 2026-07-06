# Skills — the pipeline as reusable Claude Code capabilities

Four skills package the technical platforms in this repo. In a Claude Code session
with this repo as the working directory they are auto-discovered; invoke by name
(the model picks them from the description, or type `/<name>`).

| skill | track | what it does | key entrypoint |
|---|---|---|---|
| **lipid-library-screen** | A | enumerate a combinatorial ionizable-lipid library → score with LiON → confidence-aware ranked shortlist (per-organ) | `analysis/enumerate_library*.py`, `lion_library.enumerate_michael_lipids` |
| **lion-modal** | A infra | train (lite / full-resilient) & screen the LiON D-MPNN on Modal GPUs | `modal_app/lion.py` |
| **target-evidence** | B | mine ChEMBL + DrugCLIP + humanPPI for a target, cross-compare + surface-enrichment | `analysis/run_target.py`, `lipidlib/targetpipe.py` |
| **lnp-delivery-kinetics** | mechanistic | ODE uptake→escape→translation model; expression dynamics & sensitivity | `lipidlib/kinetics.py`, `analysis/delivery_kinetics.py` |

Typical composition:
1. `target-evidence` — pick/validate a receptor + surface ligand (Track B).
2. `lion-modal` → `lipid-library-screen` — design the ionizable lipid (Track A),
   optionally per target organ.
3. `lnp-delivery-kinetics` — turn the chosen lipid's potency into expression
   dynamics and reason about cargo / dosing.

Each skill lists its prerequisites (deps, `MODAL_TOKEN_ID/SECRET` for GPU work),
exact commands, outputs, and caveats. Free/CPU: enumeration, mining, enrichment,
kinetics. GPU/Modal: LiON training + screening. Hosted MCP: Boltz, ChEMBL.
