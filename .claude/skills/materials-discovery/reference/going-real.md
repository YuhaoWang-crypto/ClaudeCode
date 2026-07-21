# Going real: wiring UMA (prescreen) + Quantum ESPRESSO on Modal (DFT verify)

The offline demo uses `score.surrogate_energy`. To make genuine discoveries, replace it
with the two real back-ends. Both keep the same signature
`formation_energy_per_atom(composition) -> eV/atom`, so `discover.screen()` is unchanged.

## 1. UMA (Meta FAIR fairchem) — the fast ML prescreen (GNoME's GNN)

Environment: the sibling **`uma-crystal-mof`** skill (fairchem-core + UMA weights, GPU).
UMA weights are gated (Hugging Face) — request access once.

```python
from fairchem.core import pretrained_mlip, FAIRChemCalculator
from ase.io import read
from ase.optimize import BFGS

def uma_formation_energy(cif_path, elemental_ev_per_atom):
    atoms = read(cif_path)
    atoms.calc = FAIRChemCalculator(pretrained_mlip.get_predict_unit("uma-s-1"),
                                    task_name="omat")
    BFGS(atoms).run(fmax=0.02, steps=200)            # relax the substituted structure
    e = atoms.get_potential_energy() / len(atoms)     # eV/atom (total)
    # formation energy = total/atom minus composition-weighted elemental references
    ef = e - sum(frac[el] * elemental_ev_per_atom[el] for el in frac)
    return ef
```

Use UMA to prescreen thousands of substituted candidates, keep the lowest-E_hull few
percent, then verify those with DFT.

## 2. Quantum ESPRESSO on Modal — the DFT verification label

Environment: **`qe-modal-bader-density`** / **`mlp-modal`** skills (build the QE image,
stage SSSP pseudopotentials, run `pw.x` SCF on a Modal GPU). For each survivor:

1. Build/relax the candidate cell (from the prototype lattice + substituted species).
2. `pw.x` SCF → total energy; repeat for the elemental references (or use a consistent
   reference table) → **DFT formation energy per atom**.
3. Feed the DFT number into `PhaseDiagram.e_above_hull` to get the **DFT E_hull**.
4. Keep E_hull ≲ 0.03 eV/atom. These are the defensible predicted-stable materials.

## 3. Reference energies for the hull

For a real hull, pull competing-phase formation energies from **Materials Project**
(same PBE/PBE+U settings as your candidates) instead of the demo's illustrative table,
or compute them with the same QE settings so everything is on one energy scale. Mixing
energy scales (e.g. MP GGA+U vs your GGA) silently corrupts E_hull — keep one scale.

## Cost discipline
- Prescreen is cheap (UMA is ~ms–s/candidate); DFT is the expense.
- Verify only the top-k by predicted E_hull, in batches; log how many candidates were
  dropped by the prescreen so coverage is honest.
- The final arbiter of a stability claim is DFT E_hull; the final arbiter of a *material*
  is experimental synthesis (as in the paper's companion A-Lab work).
