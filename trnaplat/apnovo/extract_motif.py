"""Extract an AP Novo motif CIF from a deposited PDB entry, correctly.

Why this module replaces a hand-written `_atom_site` loop: a motif file that AP
Novo can load is a **complete mmCIF**, not an atom list. The working example
shipped in the repo (`examples/kemp_eliminase/*.cif`) carries `_entry`,
`_chem_comp`, `_entity`, `_entity_poly`, `_entity_poly_seq`,
`_pdbx_poly_seq_scheme`, `_pdbx_nonpoly_scheme`, `_struct_asym` and
`_atom_site`. ✅ Enumerated from that file, after three successive failures on
Modal taught the lesson one category at a time:

1. `ValueError: The CIF file does not start with the data_ field.`
   — the header was `data pylrs_pyl_motif`, with a space.
2. `KeyError: '_atom_site.pdbx_PDB_model_num'`
   — four of the fifteen columns AF3 reads were absent.
3. `KeyError: ('B', 1)`
   — the author→internal residue map is built from
   `struct.author_naming_scheme`, which only exists when the entity and
   `pdbx_*_scheme` categories are there to define it.

Hand-patching category by category is how you get a fourth failure. `gemmi`
(MIT) writes all of them from a real structure, so this reads the deposited
entry and subsets it.

## What it does

* keeps the chosen protein residues in chain A, under their **author**
  numbering, which is what the manifest's `A300,A302,…` refers to;
* moves the ligand into its own chain so the manifest's `/B1` resolves;
* drops waters and cryoprotectants (`HOH`, `EDO`), which are not motif;
* runs `setup_entities()` so the entity and scheme categories are generated
  rather than invented.

⚠️ 2Q7H holds everything in a single deposited chain A, including the ligand.
Splitting the ligand out is required by the manifest grammar, not a property of
the entry.

Run:  python3 -m trnaplat.apnovo.extract_motif --pdb 2q7h.cif --out motif.cif \\
          --residues 300,302,305,306,346,348,384,401,417 --ligand YLY
"""

from __future__ import annotations

import argparse
import pathlib

#: Not part of any motif: ordered water and the cryoprotectant from the drop.
DISCARD = frozenset({"HOH", "EDO", "DOD"})


def extract(pdb_path: pathlib.Path, residues: list[int], ligand: str | None,
            entry_id: str, protein_chain: str = "A",
            ligand_chain: str = "B") -> tuple[object, dict]:
    """Subset a deposited entry into a motif structure. Returns (structure, info)."""
    import gemmi

    source = gemmi.read_structure(str(pdb_path))
    # ⚠️ 2Q7H models MET344 in two alternate conformations. Keeping both doubles
    # that residue's atoms and silently moves any distance computed from it --
    # it flipped MET344 and GLU396 between subsites in the moiety split below.
    # Standard practice, and what the audited package did: keep one.
    source.remove_alternative_conformations()
    source.setup_entities()
    model = source[0]

    wanted = set(residues)
    found_protein: list = []
    found_ligand: list = []
    for chain in model:
        for res in chain:
            if res.name in DISCARD:
                continue
            info = gemmi.find_tabulated_residue(res.name)
            is_aa = bool(info and info.is_amino_acid())
            if is_aa and res.seqid.num in wanted:
                found_protein.append(res)
            elif ligand and res.name == ligand:
                found_ligand.append(res)

    missing = wanted - {r.seqid.num for r in found_protein}
    if missing:
        raise ValueError(
            f"{pdb_path.name} has no amino-acid residue at {sorted(missing)}; "
            "check the numbering against the deposited entry"
        )
    if ligand and not found_ligand:
        raise ValueError(f"{pdb_path.name} has no {ligand!r} residue")
    if len(found_ligand) > 1:
        raise ValueError(
            f"{pdb_path.name} has {len(found_ligand)} copies of {ligand!r}; "
            "the manifest's single /B1 reference would be ambiguous"
        )

    out = gemmi.Structure()
    out.name = entry_id
    out.spacegroup_hm = source.spacegroup_hm
    out.cell = source.cell
    out_model = gemmi.Model("1")

    chain_a = gemmi.Chain(protein_chain)
    for res in sorted(found_protein, key=lambda r: r.seqid.num):
        chain_a.add_residue(res.clone())
    out_model.add_chain(chain_a)

    if found_ligand:
        chain_b = gemmi.Chain(ligand_chain)
        lig = found_ligand[0].clone()
        lig.seqid = gemmi.SeqId(1, " ")      # the manifest refers to it as B1
        lig.het_flag = "H"
        chain_b.add_residue(lig)
        out_model.add_chain(chain_b)

    out.add_model(out_model)
    out.setup_entities()
    out.assign_label_seq_id()

    # `make_mmcif_document` only emits `_entity_poly_seq` and
    # `_pdbx_poly_seq_scheme` for a polymer entity that declares its sequence,
    # and those are two of the categories AF3 needs to build
    # `author_naming_scheme` -- the map whose absence produced KeyError ('B', 1).
    for entity in out.entities:
        if entity.entity_type == gemmi.EntityType.Polymer and not entity.full_sequence:
            names: list[str] = []
            for subchain in entity.subchains:
                for chain in out_model:
                    for res in chain:
                        if res.subchain == subchain:
                            names.append(res.name)
            entity.full_sequence = names or [r.name for r in chain_a]

    return out, {
        "protein_residues": [(r.name, r.seqid.num) for r in chain_a],
        "ligand": (found_ligand[0].name if found_ligand else None),
        "ligand_atoms": (len(found_ligand[0]) if found_ligand else 0),
        "protein_atoms": sum(len(r) for r in chain_a),
    }


def _add_scheme_loops(block, structure) -> None:
    """Append the `pdbx_*_scheme` loops gemmi does not write.

    ⚠️ gemmi's `MmcifOutputGroups` has no option for these two categories, so
    they are written here. They are the label↔author correspondence tables, and
    `_pdbx_poly_seq_scheme` is what lets AF3 build `author_naming_scheme`; its
    absence is the KeyError ('B', 1) in the module docstring. Column order
    follows the repo's working example.
    """
    import gemmi

    model = structure[0]
    entity_of = {}
    for entity in structure.entities:
        for subchain in entity.subchains:
            entity_of[subchain] = entity.name

    poly_rows, nonpoly_rows = [], []
    for chain in model:
        for res in chain:
            entity_id = entity_of.get(res.subchain, "1")
            info = gemmi.find_tabulated_residue(res.name)
            is_aa = bool(info and info.is_amino_acid())
            seq = str(res.label_seq) if res.label_seq is not None else "."
            if is_aa:
                poly_rows.append([
                    res.subchain or chain.name, str(res.seqid.num), entity_id,
                    "n", res.name, ".", str(res.seqid.num), chain.name, seq,
                ])
            else:
                nonpoly_rows.append([
                    res.subchain or chain.name, str(res.seqid.num), entity_id,
                    res.name, ".", str(res.seqid.num), chain.name,
                ])

    if poly_rows:
        loop = block.init_loop("_pdbx_poly_seq_scheme.", [
            "asym_id", "auth_seq_num", "entity_id", "hetero", "mon_id",
            "pdb_ins_code", "pdb_seq_num", "pdb_strand_id", "seq_id",
        ])
        for row in poly_rows:
            loop.add_row(row)
    if nonpoly_rows:
        loop = block.init_loop("_pdbx_nonpoly_scheme.", [
            "asym_id", "auth_seq_num", "entity_id", "mon_id", "pdb_ins_code",
            "pdb_seq_num", "pdb_strand_id",
        ])
        for row in nonpoly_rows:
            loop.add_row(row)


def write(structure, target: pathlib.Path) -> None:
    doc = structure.make_mmcif_document()
    _add_scheme_loops(doc.sole_block(), structure)
    target.write_text(doc.as_string())


def categories(path: pathlib.Path) -> list[str]:
    """mmCIF categories present, for comparison against a working example."""
    seen = []
    for line in path.read_text().splitlines():
        stripped = line.strip()
        if stripped.startswith("_") and "." in stripped:
            cat = stripped.split(".", 1)[0]
            if cat not in seen:
                seen.append(cat)
    return seen


#: Categories the repo's own working example carries. A motif missing one of
#: these is the shape AF3 refused three times above.
REQUIRED_CATEGORIES = [
    "_entry", "_chem_comp", "_entity", "_entity_poly", "_entity_poly_seq",
    "_pdbx_poly_seq_scheme", "_struct_asym", "_atom_site",
]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pdb", type=pathlib.Path, required=True,
                    help="deposited mmCIF, e.g. 2q7h.cif from files.rcsb.org")
    ap.add_argument("--out", type=pathlib.Path, required=True)
    ap.add_argument("--residues", required=True,
                    help="comma-separated author residue numbers")
    ap.add_argument("--ligand", default="YLY")
    ap.add_argument("--entry-id", default=None)
    args = ap.parse_args(argv)

    residues = [int(x) for x in args.residues.split(",")]
    entry_id = args.entry_id or args.out.stem
    structure, info = extract(args.pdb, residues, args.ligand, entry_id)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    write(structure, args.out)

    print(f"{args.out}")
    print(f"  protein : {len(info['protein_residues'])} residues,"
          f" {info['protein_atoms']} atoms")
    print(f"            {', '.join(f'{n}{i}' for n, i in info['protein_residues'])}")
    if info["ligand"]:
        print(f"  ligand  : {info['ligand']} as chain B residue 1,"
              f" {info['ligand_atoms']} atoms")
    present = categories(args.out)
    missing = [c for c in REQUIRED_CATEGORIES if c not in present]
    print(f"  mmCIF categories: {len(present)}")
    if missing:
        print(f"  ❌ missing {missing} -- AF3 will refuse this file")
        return 1
    print("  ✅ every category the repo's working example carries is present")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
