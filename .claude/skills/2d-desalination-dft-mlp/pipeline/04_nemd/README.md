# Stage 3 — NEMD (water flux & salt rejection)

**Goal:** apply a pressure difference across the membrane and measure the two
numbers that define an RO membrane:

- **Water flux / permeance** — throughput of water across the pore.
- **Salt rejection** `R = 1 − c_permeate / c_feed`.

### Files
- `in.desalination.lammps` — LAMMPS NEMD: rigid membrane, feed **piston** applies
  the pressure, DeePMD potential from Stage 2, trajectory dumped for analysis.
- `analyze_flux_rejection.py` — engine-agnostic analysis. Importable functions:
  - `count_crossings` — net feed→permeate crossings (boolean-side counting, so a
    particle sitting exactly on the membrane plane is never double-counted);
  - `water_flux`, `salt_rejection`, `flux_to_permeance` (→ L·cm⁻²·day⁻¹·MPa⁻¹);
  - `load_lammps_dump` — parse a `custom id type x y z` dump with no heavy deps;
  - `summarize` — one call → all headline numbers.

### Run
```
lmp -in in.desalination.lammps          # or: sbatch ../slurm/lammps_nemd.slurm
python analyze_flux_rejection.py --demo  # self-check on synthetic data
```
Then analyse the real dump:
```python
import analyze_flux_rejection as a
zw, zi, t = a.load_lammps_dump("traj.lammpstrj", water_types=[3], ion_types=[4,5])
print(a.summarize(zw, zi, z_membrane=0.0, times_ns=t*5e-4,
                  area_A2=150., n_ions_feed_initial=20, pressure_MPa=100.))
```

### Calibrate the piston force
The per-atom `fz` in the LAMMPS script is a placeholder — compute it from your
target pressure: `f = P · A_box / N_piston` (watch units: metal = eV/Å).

### Validation
`../05_toy_validation/toy_nemd.py` exercises **exactly these analysis functions**
against a known ground truth (see the top-level README).
