# Host structures

Drop your host structure here as `bcd.mol2` (Tripos mol2, explicit hydrogens).

- `bcd.mol2` — β-cyclodextrin. A clean starting structure is available from the
  pAPRika β-CD tutorial, or build one from the PDB ligand code `BCD` / a GLYCAM
  build and export to mol2 with explicit H.
- Modified hosts are generated on the fly by `src/fep_ti_ddg.py` from
  `config/modifications.yaml` (or supply a pre-built mol2 per modification).

The pipeline asserts the host net charge from config; make sure protonation is
correct (neutral β-CD; cationic for amino/ammonium modifications).
