"""MOL2 and SDF/MOL reader for small molecule ligands."""
from __future__ import annotations
import numpy as np
from typing import List
from ..core import Atom, Bond, Molecule


def read_mol2(path: str) -> Molecule:
    """Parse a Tripos MOL2 file."""
    atoms: List[Atom] = []
    bonds: List[Bond] = []
    mol_name = path.split('/')[-1].replace('.mol2', '')

    with open(path) as f:
        lines = f.readlines()

    section = None
    for line in lines:
        line = line.rstrip()
        if line.startswith('@<TRIPOS>'):
            section = line[9:].strip()
            continue

        if section == 'MOLECULE' and not atoms and not mol_name.endswith('.mol2'):
            mol_name = line.strip() or mol_name
            section = None

        elif section == 'ATOM':
            parts = line.split()
            if len(parts) < 6:
                continue
            try:
                idx = int(parts[0]) - 1
                name = parts[1]
                x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
                atom_type = parts[5]
                charge = float(parts[-1]) if len(parts) >= 9 else 0.0
                # Element from atom_type (e.g. 'C.3' -> 'C', 'N.am' -> 'N')
                element = atom_type.split('.')[0].upper()
                # Residue info
                res_name = parts[7] if len(parts) > 7 else 'LIG'
                res_id = int(''.join(filter(str.isdigit, parts[6])) or '1') if len(parts) > 6 else 1
                atoms.append(Atom(
                    index=idx, name=name, element=element,
                    residue_name=res_name, residue_id=res_id, chain='L',
                    position=np.array([x, y, z], dtype=np.float64),
                    charge=charge, atom_type=atom_type, is_hetatm=True
                ))
            except (ValueError, IndexError):
                continue

        elif section == 'BOND':
            parts = line.split()
            if len(parts) < 4:
                continue
            try:
                a1, a2 = int(parts[1]) - 1, int(parts[2]) - 1
                btype = parts[3]
                order = {'1': 1, '2': 2, '3': 3, 'ar': 2, 'am': 1}.get(btype, 1)
                bonds.append(Bond(a1, a2, order))
            except (ValueError, IndexError):
                continue

    mol = Molecule(mol_name)
    mol.atoms = atoms
    mol.bonds = bonds
    return mol


def read_sdf(path: str) -> List[Molecule]:
    """Parse an SDF file, return list of molecules (one per record)."""
    molecules: List[Molecule] = []

    with open(path) as f:
        content = f.read()

    for record in content.split('$$$$'):
        record = record.strip()
        if not record:
            continue
        lines = record.split('\n')
        if len(lines) < 4:
            continue

        mol_name = lines[0].strip() or 'LIG'
        # counts line (line 4, index 3)
        try:
            counts = lines[3]
            n_atoms = int(counts[:3])
            n_bonds = int(counts[3:6])
        except (ValueError, IndexError):
            continue

        atoms: List[Atom] = []
        bonds: List[Bond] = []

        for i in range(4, 4 + n_atoms):
            try:
                parts = lines[i].split()
                x, y, z = float(parts[0]), float(parts[1]), float(parts[2])
                element = parts[3].upper()
                atoms.append(Atom(
                    index=i - 4, name=f"{element}{i-4}",
                    element=element, residue_name='LIG', residue_id=1,
                    chain='L',
                    position=np.array([x, y, z], dtype=np.float64),
                    is_hetatm=True
                ))
            except (ValueError, IndexError):
                continue

        for i in range(4 + n_atoms, 4 + n_atoms + n_bonds):
            try:
                parts = lines[i].split()
                a1, a2, order = int(parts[0]) - 1, int(parts[1]) - 1, int(parts[2])
                bonds.append(Bond(a1, a2, order))
            except (ValueError, IndexError):
                continue

        mol = Molecule(mol_name)
        mol.atoms = atoms
        mol.bonds = bonds
        molecules.append(mol)

    return molecules
