# Host structures

- **`bcd.sdf`** — β-cyclodextrin, real 3D geometry (147 atoms, C₄₂H₇₀O₃₅) from
  the PDB Chemical Component Dictionary (ligand code `BCD`, ideal coordinates,
  RCSB `files.rcsb.org/ligands/download/BCD_ideal.sdf`). This is what
  `config/system.yaml → host.structure` points at; antechamber converts it to a
  charged mol2 during parameterization. **Shipped and ready — no action needed.**

To use a different host, drop an sdf/pdb/mol2 here and update `host.structure`.
- Modified hosts are generated on the fly by `src/fep_ti_ddg.py` from
  `config/modifications.yaml` (or supply a pre-built mol2 per modification).

The pipeline asserts the host net charge from config; make sure protonation is
correct (neutral β-CD; cationic for amino/ammonium modifications).
