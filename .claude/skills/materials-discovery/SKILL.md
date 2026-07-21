---
name: materials-discovery
description: >
  GNoME-style discovery of new stable crystalline materials — reproduce the workflow of
  Merchant et al., "Scaling deep learning for materials discovery" (Nature 2023,
  doi:10.1038/s41586-023-06735-9). Generate candidate crystals by ion substitution into
  known prototypes (structural pipeline) and by enumerating charge-neutral stoichiometries
  (compositional pipeline), prescreen them with a fast energy model, DFT-verify the top hits,
  and keep those on the convex hull (energy above hull ≈ 0). The convex-hull math and
  candidate generators are exact and unit-tested offline; energies come from a pluggable
  scorer — an offline ionicity surrogate for the demo, Meta FAIR UMA (fairchem) for the ML
  prescreen, or Quantum ESPRESSO on Modal for DFT verification. Use to screen for new stable
  compounds (Li-ion conductors, perovskite oxides, 2D layered, spinels), rank candidates by
  predicted stability, compute energy-above-hull / decomposition, or reproduce this paper.
  Triggers: discover new materials, crystal structure prediction, convex hull stability,
  energy above hull, ion substitution, GNoME, screen stable compounds, formation energy.
---

# Materials discovery (GNoME-style)

Reproduces the **discovery loop** of GNoME (Google DeepMind, Nature 2023): cheaply
predict which candidate crystals are stable, verify the best with DFT, keep the ones
on the convex hull. This skill provides the exact, tested **scaffold** (candidate
generation + convex-hull decision) and **pluggable scorers** for the energies.

## The loop (maps 1:1 to the paper)

```
generate candidates            structural (substitute ions into a prototype)
   │                           compositional (enumerate charge-neutral stoichiometries)
   ▼
cheap prescreen                score every candidate → energy above hull (E_hull)
   │                           [surrogate offline · UMA for real ML prescreen]
   ▼
verify top-k with DFT          Quantum ESPRESSO SCF on Modal (the active-learning label)
   │
   ▼
keep E_hull ≤ tol              convex-hull decision (exact LP) → new stable materials
   │
   └── retrain cheap model on new labels → repeat
```

## Quick start (offline, no cloud, no weights)

```bash
python3 -m materials_discovery.test_materials     # exact hull + generator tests
python3 -m materials_discovery.discover           # alkaline-earth perovskite demo
```

The demo generates ABO₃ perovskite candidates over {Ca,Sr,Ba}×{Ti,Zr}×O, scores them
with the offline ionicity **surrogate**, and ranks by E_hull vs the binary oxides —
recovering BaZrO₃/BaTiO₃/SrTiO₃/SrZrO₃ as the most stable. ⚠️ surrogate energies are a
heuristic, **not DFT** — they exercise the loop; real stability needs UMA/DFT.

## Going real (the two paid/heavy back-ends)

| stage | back-end | how |
|---|---|---|
| ML prescreen (GNoME's GNN) | **UMA** (Meta FAIR fairchem) | `score.uma_energy` → run in the **`uma-crystal-mof`** skill env; relax candidate + subtract elemental refs |
| DFT verification (the label) | **Quantum ESPRESSO** on Modal GPU | `score.qe_modal_energy` → **`qe-modal-bader-density`** / **`mlp-modal`** skills |
| triage / routing | — | **`materials-compute`** |

Swap the scorer: `discover(role_candidates, prototype, refs, scorer="uma")` for the
prescreen, then verify the surviving candidates with `qe_modal_energy`.

## Files
- `materials_discovery/composition.py` — formula algebra (no pymatgen dep)
- `materials_discovery/hull.py` — convex hull / E_hull / decomposition (exact LP)
- `materials_discovery/prototypes.py` — crystal prototypes, oxidation states, Shannon radii, tolerance factor, demo reference energies
- `materials_discovery/generate.py` — structural (substitution) + compositional pipelines with charge-neutrality & tolerance filters
- `materials_discovery/score.py` — surrogate / UMA / QE-on-Modal scorers (same signature)
- `materials_discovery/discover.py` — the active-learning loop + demo
- `materials_discovery/test_materials.py` — deterministic tests (hull signs, neutrality, tolerance)
- `reference/methodology.md` — GNoME mapping + honesty labels
- `reference/going-real.md` — wiring UMA + Modal DFT

## Honesty
- ✅ convex-hull math, decomposition, candidate generation, charge/tolerance filters — exact, unit-tested.
- ⚠️ the offline **surrogate** energy is an ionicity heuristic, not DFT — for workflow demonstration only.
- ✅ real discovery requires UMA (prescreen) + DFT (verify); the deciding number for any claim is a **DFT energy-above-hull**, and ultimately experimental synthesis — the same discipline the paper follows.
