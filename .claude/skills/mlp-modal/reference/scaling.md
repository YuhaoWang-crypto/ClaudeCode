# Scaling from demo to production

The bundled demo (molten NaCl, 64 atoms, few steps) proves the pipeline. To reach
production quality, change parameters only — the code path is identical.

## Parameters
| Goal | Change |
|---|---|
| System size 50–100+ atoms | `--nx/--ny/--nz` (2×2×2 = 64; 3×3×3 rocksalt supercells etc.) |
| Multiple temperatures | `--temps "900,1000,1100,1200,1300,1400"` (production stage) |
| Proper equilibration + sampling | `--steps` from ~10 to several thousand (1 fs each ≈ ps-scale) |
| Longer/better training | `--dp-steps` 3000 → 10⁵–10⁶ |
| Faster training | in the decorator, `gpu="A100"` or `gpu="H100"` instead of `"T4"` |
| Faster DFT | raise `cpu=` on `run_cp2k`; or MPI-parallelize (see below) |

## Getting production accuracy (not just a working pipeline)
1. **Melt properly**: start from the crystal, run NVT well above T_melt to liquefy,
   then equilibrate at each target T for ~10 ps *before* collecting training frames.
   The demo skips this (heats an expanded lattice briefly) — fine for plumbing only.
2. **Sample decorrelated frames**: take frames every N steps (e.g. every 20–50 fs),
   not every step, to reduce correlation.
3. **Cover the phase space** the model will see: all 6 T points, and ideally a couple
   of densities/volumes if you want NPT properties.
4. **Validate the potential**: hold out frames; check energy/force RMSE; run an MD
   with the trained model (`dp` + LAMMPS/ASE) and compare RDF, MSD/diffusion,
   density against the DFT reference before trusting macroscopic numbers.

## Different chemistries
Change `build_nacl` (or drop in ASE/pymatgen to build any structure), update the
`&KIND` blocks (basis/pseudopotential per element), and the DeepMD `type_map`/`sel`.
For non-cubic or multi-species melts, generate coordinates however you like and feed
them into `make_cp2k_input`.

## Parallelizing the temperature sweep
The sweep runs temperatures sequentially for reliability. To parallelize, replace the
`for ... run_cp2k.remote(inp)` loop with `run_cp2k.map(inputs)` — but first raise
`memory=`, lower per-run `cpu=`, and confirm the account's concurrent-container limit,
or co-scheduled MPI runs abort (that's why the default is sequential).

## Other engines
Same Modal pattern works for Quantum ESPRESSO (conda-forge / NGC image, GPU via
OpenACC), VASP (user supplies licensed source, compile in image), phonopy (light —
it orchestrates the DFT engine's force calls). ORCA is CPU-only (no GPU benefit).
