# HPC input-file templates

Ready-to-adapt input files for the **production** calculations that the local
demos only proxy. Each folder maps to one client direction. Replace the marked
`<<...>>` placeholders (system, paths, pseudopotentials, resources) and submit.

> These are templates, not turnkey jobs. Always converge cutoff / k-points /
> timestep on YOUR system before trusting production numbers.

| Folder | Client direction | Produces | Local demo it upgrades |
|---|---|---|---|
| `moltensalt_cp2k/` | 熔盐 ML 势 (郭硕) | AIMD trajectory (energy+forces) | #4 |
| `moltensalt_vasp/` | 熔盐 ML 势 — VASP alternative | AIMD trajectory | #4 |
| `deepmd/` | 熔盐 ML 势 | trained DeePMD potential | #4 |
| `cd_docking/` | 蛋白-Cd²⁺ (bizlikery) | metal-site docking poses | #5 |
| `sers_au_slab/` | 手性 SERS on Au (旺仔秋秋唐) | adsorption E + Raman on Au(111) | #6 |
| `phonon_qe/` | 声子 / 电子-声子 / STE (审稿补算) | phonon spectrum + Huang-Rhys | #3 |

## Typical end-to-end (molten-salt ML potential example)

```
1. moltensalt_builder.py FLiNaK --natoms 100      # initial structure
2. moltensalt_cp2k/aimd.inp                        # AIMD at each T point
3. dpdata: parse AIMD -> deepmd/npy                # (pack_deepmd in demo #4)
4. deepmd/train.json -> dp train                   # train the ML potential
5. LAMMPS + dp pair_style                          # long MD -> density, viscosity
```

## Honesty note

Every template encodes a *reasonable default*, not a converged production
setting. Cutoffs, k-meshes, timesteps, and functionals must be validated per
system. Flag to the client which numbers are converged vs. placeholder.
