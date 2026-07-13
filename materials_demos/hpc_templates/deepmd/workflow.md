# DeePMD-kit: from AIMD frames to macroscopic properties

## 1. Package AIMD into training data
```python
import dpdata
# CP2K:  fmt="cp2k/aimd_output"      VASP: fmt="vasp/outcar" or "vasp/xml"
ls = dpdata.LabeledSystem("cp2k_run_dir", fmt="cp2k/aimd_output")
ls.to("deepmd/npy", "deepmd_dataset/T700K")
```
(the local demo #4 `pack_deepmd()` does exactly this from ASE frames).

## 2. Train + freeze + compress
```bash
dp train  train.json
dp freeze -o graph.pb
dp compress -i graph.pb -o graph-compressed.pb
```
Watch `lcurve.out`: force RMSE should fall to ~few×10 meV/Å.

## 3. Validate the potential
```bash
dp test -m graph.pb -s ../deepmd_dataset/T900K -n 200
```
Also do **DP-GEN active learning** if the melt explores configs the initial
AIMD missed (recommended for transferable potentials across the 6 temperatures).

## 4. Long MD -> density & viscosity (the client's real deliverable)
Run LAMMPS with the trained potential:
```
pair_style deepmd graph-compressed.pb
pair_coeff * *
fix 1 all npt temp ${T} ${T} 0.1 iso 1.0 1.0 1.0   # NPT -> density(T)
# Green-Kubo (stress ACF) or Einstein (MSD) -> viscosity / diffusion
```
- **density(T)**: average box volume in NPT.
- **viscosity(T)**: Green-Kubo integral of the off-diagonal stress autocorrelation.
- **diffusion(T)**: MSD slope (the local demo #4 computes this on the proxy).
