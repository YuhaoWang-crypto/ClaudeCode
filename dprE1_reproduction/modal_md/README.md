# Part 3 — Protein–Ligand MD on Modal (making-it-rain port)

Reproduces the 500 ns MD + binding-energy stage of the paper (their CHARMm /
Discovery Studio protocol) with an **open-source OpenMM + AMBER** pipeline,
ported from the *making-it-rain* `Protein_ligand` notebook
(Arantes et al., https://github.com/pablo-arantes/making-it-rain) to **Modal**
serverless GPUs.

## What it does (faithful to making-it-rain protein-ligand)

| Stage | Tool / setting |
|---|---|
| Protein prep | PDBFixer (add missing atoms + H at pH 7) |
| Protein FF | AMBER **ff14SB** (or ff19SB) |
| Ligand FF | **GAFF2**, AM1-BCC charges (antechamber / parmchk2) |
| Topology | AmberTools **tleap** → `prmtop`/`inpcrd` |
| Water / ions | **TIP3P**, 12 Å box, Na⁺/Cl⁻ neutralise |
| Engine | **OpenMM 8**, LangevinMiddleIntegrator, 2 fs, HBonds, PME 1.0 nm |
| Protocol | minimise → NVT heat → NPT equilibrate → NPT production |
| Analysis | mdtraj: ligand / protein-backbone / cofactor RMSD, protein RMSF (paper Table 4) |

## Run it

```bash
pip install modal
modal token new                       # one-time, uses your Modal account

# 1) build docked starting complexes locally (uses Part 2 docking)
python modal_md/prepare_inputs.py                     # all candidates, DprE1_WT
#   -> modal_md/inputs/DprE1_WT__GTD_9.7/{protein.pdb, ligand.sdf}

# 2) launch MD on Modal GPUs
modal run modal_md/app.py --system DprE1_WT__GTD_9.7 --ns 5      # quick demo
modal run modal_md/app.py --all --ns 500                        # full paper protocol
```

Each job returns the Table-4 metrics and writes trajectory + `metrics.json` to a
Modal Volume (`dprE1-md-results`). Jobs fan out one GPU container per complex, so
all 11 candidates run in parallel.

## Cost / time note

500 ns of an explicit-solvent protein–ligand system is ~real compute: on an
A10G expect roughly a day per complex. Start with `--ns 5` to validate the setup,
then scale up. Swap `gpu="A10G"` for `"A100"` in `app.py` for ~2–3× throughput.

## Extending to match the paper exactly

- **FAD / HEM cofactor in the box.** `prepare_inputs.py` currently writes the
  protein without its cofactor (the making-it-rain protein-ligand flow
  parametrises a single GAFF2 ligand). To reproduce the paper's *cofactor RMSD*
  column, add FAD/HEM parameters in the `tleap` step — e.g. load a prebuilt FAD
  `lib`+`frcmod` (Manchester AMBER parameter database) before `combine`. The
  analysis already reports cofactor RMSD when those residues are present.
- **MM-GBSA ΔG (paper Table 4).** AmberTools `MMPBSA.py` is in the image; run it
  on the last 100 ns of the trajectory (as the paper does) with an `igb=5`/GBSW-
  like model to obtain ΔG_binding. A helper stub is in `app.py::analyse` comments.
- **Replicas / longer sampling.** Set `--ns 500` and launch multiple seeds by
  calling `simulate.spawn` repeatedly.
```
