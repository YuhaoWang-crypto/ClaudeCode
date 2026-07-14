# Stage 1 — DFT (the "gold standard" reference data)

**Goal:** produce the first-principles numbers that everything downstream is
trained on and validated against.

| Quantity | What it gives you | QE tool |
|---|---|---|
| Total / interaction energy | Water & ion binding to the pore, charge transfer | `pw.x` (`scf_graphene.in`) |
| Bader / Löwdin charges | Charge redistribution at the pore edge | `pp.x` + Bader |
| **Pore-crossing barrier** `E_a` | The key selectivity number (water vs Na⁺ vs Cl⁻) | `neb.x` (`neb_na_pore.in`) |
| AIMD forces | Training labels for the MLP | `cp.x` / CP2K |

### How to run
- **Cloud (recommended to start):** `modal run modal_run_dft.py --infile qe/scf_graphene.in`
  — plane-wave DFT is **CPU-bound**, so this uses a fat CPU box, not a GPU.
- **HPC:** `sbatch ../slurm/qe_scf.slurm`

### Generating inputs from the built geometry
```python
from ase.io import read, write
atoms = read("../figures/graphene_pore.xyz")
write("scf_from_builder.in", atoms, format="espresso-in",
      input_data={...}, pseudopotentials={...})
```
The `__PLACEHOLDER__` fields in the templates (nat, cell, positions) are what
ASE fills in automatically.

### Functional choice
Water/ion physisorption needs dispersion: use **vdW-DF2** (set in the templates)
or PBE+D3. Plain PBE will underbind water and give the wrong barriers.
