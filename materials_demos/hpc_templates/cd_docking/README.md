# Cd(2+) metal-site docking / refinement (upgrades demo #5)

Demo #5 shortlists candidate Cd2+ sites geometrically. This folder takes the top
sites to a defensible binding model. Metal-ion docking is NOT AutoDock's default
strength (scoring functions are parameterised for organics), so use a metal-aware
path and confirm with QM/MM.

## Recommended pipeline

1. **Shortlist** (local demo #5): `cd_binding_site_predictor.py your.pdb`
   -> top clusters + residues.

2. **Prepare receptor** (AutoDockTools / Meeko):
   ```bash
   prepare_receptor -r receptor.pdb -o receptor.pdbqt
   # keep the predicted donor residues; assign protonation at your pH (propka)
   ```

3. **Place Cd2+** and dock with a metal-compatible setup:
   - AutoDock4 has a `Cd` atom type in `AD4_parameters.dat` (r_eqm/eps entries);
     add/verify the Cd row before `autogrid4`.
   - Center the grid box on the predicted site centroid (from demo #5 output),
     ~15-20 Å box.
   - Better for metals: **AutoDock-GPU** with a custom Cd map, or **GOLD** (has
     explicit metal coordination handling), or a restrained docking that biases
     Cd toward the shortlisted donor atoms.

4. **Refine the coordination sphere with QM/MM** (the honest, publishable step):
   - cut a cluster (Cd2+ + first-shell donors + backbone), optimise in ORCA/Gaussian
     (e.g. B3LYP-D3/def2-TZVP, Cd with def2 ECP), verify Cd–L distances 2.3-2.6 Å
     and coordination number 4-6.

5. **Design mutants for MST**: mutate the shortlisted donors (e.g. His->Ala,
   Glu->Gln) and compare predicted binding; these are the constructs to test.

## Files
- `dock_config_vina.txt` — AutoDock Vina box template (organic-scoring caveat noted).
- `qmmm_cluster.md` — how to cut and run the QM cluster refinement.
