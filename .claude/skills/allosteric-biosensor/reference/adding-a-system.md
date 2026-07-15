# Adding a new receptor / reporter / analyte

The whole point of the recipe is that it is analyte-agnostic. Plugging in a new
detection target is a data edit in `systems.py`, not new logic.

## Add a new receptor (and its analyte)

1. Get a sequence. Prefer a deposited structure (RCSB FASTA:
   `https://www.rcsb.org/fasta/entry/<PDB>`); strip purification tags (His6 etc.).
2. Get the analyte SMILES (PubChem REST:
   `.../compound/name/<name>/property/CanonicalSMILES/JSON`).
3. Derive loop permutation sites from the structure:
   ```python
   from biosensor_pipeline.structure import annotate
   seq_modeled, offset, sse, loop_centers = annotate("<PDB>", "A", full_seq)
   ```
   `loop_centers` is the `<10`-member candidate list.
4. Add a `Receptor(...)` to `systems.py` with `seq`, `ligand_name`,
   `ligand_smiles`, `loop_sites=loop_centers`, `modeled_offset=offset`.

## Add a new reporter

Add a `Reporter(...)` with:
- `seq` — the mature enzyme sequence.
- `catalytic` — 0-indexed catalytic residues. Locate them by **sequence motif**
  (robust to numbering), not by literal residue number. Verify against the
  structure.
- `insertion_sites` — 0-indexed cut points at permissive **surface loops** (use
  `annotate` on the reporter structure; pick coil residues away from the active
  site). Map any literature residue number via the deposited structure's author
  numbering.

`scoring.map_reporter_index` assumes catalytic residues sit N-terminal to the
insertion (true for TEM-1/BLA-253). If you insert *upstream* of a catalytic
residue, that mapping already shifts it by the insert length — no change needed.

## Wire up a System

```python
SYSTEMS["mytag"] = System(
    key="mytag", role="validation",
    receptor=MY_RECEPTOR, reporter=TEM1,
    primary_insertion="253",
    description="...",
)
```

Then everything works unchanged:
```bash
python3 -m biosensor_pipeline.run_repro mytag     # build + verify library, emit payloads
python3 -m biosensor_pipeline.structure           # verify loop annotations reproduce
# submit Boltz holo/apo/control (reference/boltz-validation.md), then:
python3 -m biosensor_pipeline.analyze_boltz
```

## Logic gates (advanced)

- **YES gate** — build two chimeras inserting the *same* cpReceptor at two loops
  (e.g. TEM-1 `insertion_sites["41"]` and `["197"]`) into one construct. Extend
  `build_chimera` to take a list of (site, cpReceptor) inserts.
- **AND gate** — same, but two *orthogonal* receptors → activity requires both
  ligands. Validate each single-ligand and the double-ligand holo separately.
- **Auxiliary domain** — append a second binder to the chimera C-terminus with a
  flexible linker to offset the CP affinity penalty.

Keep the honesty labels: sequence construction is ✅; whether a gate actually
computes AND is a wet-lab measurement (⚠️ until assayed).
