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

### Alternatives to DeePMD
`MACE`, `NequIP`/`Allegro` (equivariant, often better data efficiency) — same
data, different trainer. Worth a bake-off if force RMSE stalls.
