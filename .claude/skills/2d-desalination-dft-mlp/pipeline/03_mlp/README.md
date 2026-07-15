# Stage 2 — Machine-Learning Potential (DFT accuracy, MD speed)

**Goal:** train a potential that reproduces the DFT energies/forces/barriers so
you can run **ns-scale NEMD** at ~classical cost. This is the recommended route
over ReaxFF for *physical* (non-reactive) RO desalination — see the top-level
README for the MLP-vs-ReaxFF rationale.

### Files
- `deepmd_input.json` — DeepPot-SE model/descriptor config (rcut 6 Å, force-weighted loss).
- `dpgen_param.json` — DP-GEN **active-learning** loop: explore → pick high-uncertainty frames → DFT-label → retrain. Minimises expensive DFT calls and, crucially, covers the **confined, high-pressure** states the NEMD piston will visit.
- `modal_train_deepmd.py` — train on a **Modal GPU** (training *is* GPU-bound).

### Workflow
1. Convert DFT outputs to DeePMD format with `dpdata`.
2. Train: `modal run modal_train_deepmd.py --input deepmd_input.json` (or `sbatch ../slurm/deepmd_train.slurm`).
3. **Validate the model against DFT** before trusting it:
   - energy/force RMSE via `dp test`;
   - re-compute the NEB barrier from Stage 1 with the MLP and check it matches.
4. Freeze → `graph-compress.pb`, hand to LAMMPS in Stage 3.

### DeePMD vs MACE vs NequIP bake-off
Three trainers on the **same DFT data**, so you can pick the winner on force RMSE
and NEB-barrier accuracy:

| Trainer | Config | Data format | GPU launcher |
|---|---|---|---|
| DeePMD-kit | `deepmd_input.json` | npy (native) | `modal_train_deepmd.py` |
| MACE | `mace_config.yaml` | extended-XYZ | `modal_train_mace.py` |
| NequIP/Allegro | `nequip_config.yaml` | extended-XYZ | `nequip-train …` |

`MACE`/`NequIP` are E(3)-equivariant and usually more **data-efficient** (fewer
DFT frames for a given force accuracy) — worth the bake-off, especially early
when DFT frames are scarce.

**Data conversion** (npy → extended-XYZ for MACE/NequIP):
```bash
python prepare_deepmd_data.py --to-extxyz ../data/train/pore_crossing \
       --out ../data/train/pore_crossing.xyz   # writes REF_energy + REF_forces
```
Keep `r_max`/`rcut = 6.0 Å` identical across all three for a fair comparison.
Validate every trained model the same way: energy/force RMSE **and** reproduce
the Stage-1 NEB barrier.
