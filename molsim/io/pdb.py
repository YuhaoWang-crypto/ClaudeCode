"""PDB format reader and writer."""
from __future__ import annotations
import numpy as np
import re
from typing import List, Optional
from ..core import Atom, Bond, Molecule

_ELEM_RE = re.compile(r'[A-Z][a-z]?')


def _guess_element(atom_name: str, residue_name: str) -> str:
    """Guess element from PDB atom name."""
    name = atom_name.strip()
    # Remove leading digits (e.g. '1HB' -> 'HB')
    name = name.lstrip('0123456789')
    if not name:
        return 'C'
    # First letter after stripping is usually element for common atoms
    # PDB columns 13-14 are atom name; element is usually first non-digit char
    el = name[0]
    if len(name) >= 2 and name[1].islower():
        el = name[:2].upper()
    # Special cases
    if el in ('C', 'N', 'O', 'S', 'H', 'P', 'F'):
        return el
    if name[:2].upper() in ('CL', 'BR', 'FE', 'ZN', 'MG', 'CA', 'NA', 'MN', 'CU', 'CO', 'NI'):
        return name[:2].upper()
    return el


def read_pdb(path: str, model: int = 0) -> Molecule:
    """
    Parse a PDB file and return a Molecule.

    Parameters
    ----------
    path : str
        Path to .pdb file.
    model : int
        Which MODEL record to read (0 = first).
    """
    atoms: List[Atom] = []
    bonds: List[Bond] = []
    conect: dict = {}
    current_model = -1
    model_started = False
    index = 0

    with open(path) as f:
        for line in f:
            rec = line[:6].strip()

            if rec == 'MODEL':
                current_model += 1
                model_started = True
                if current_model > model:
                    break
                continue

            if rec == 'ENDMDL':
                if current_model == model:
                    break
                continue

            if rec not in ('ATOM', 'HETATM'):
                if rec == 'CONECT':
                    # Parse connectivity
                    fields = line.split()
                    if len(fields) >= 2:
                        src = int(fields[1]) - 1
                        for dst_str in fields[2:]:
                            dst = int(dst_str) - 1
                            if dst > src:
                                conect.setdefault(src, set()).add(dst)
                continue

            if model_started and current_model < model:
                continue

            try:
                is_hetatm = rec == 'HETATM'
                serial = int(line[6:11])
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                chain = line[21].strip() or 'A'
                res_id = int(line[22:26])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])

                # Element from columns 77-78 if present, else guess
                element_col = line[76:78].strip().upper() if len(line) > 76 else ''
                element = element_col if element_col else _guess_element(atom_name, res_name)

                atom = Atom(
                    index=index,
                    name=atom_name,
                    element=element,
                    residue_name=res_name,
                    residue_id=res_id,
                    chain=chain,
                    position=np.array([x, y, z], dtype=np.float64),
                    is_hetatm=is_hetatm,
                )
                atoms.append(atom)
                index += 1
            except (ValueError, IndexError):
                continue

    # Build bonds from CONECT records
    for src, dsts in conect.items():
        for dst in dsts:
            bonds.append(Bond(src, dst))

    mol = Molecule(name=path.split('/')[-1].replace('.pdb', ''))
    mol.atoms = atoms
    mol.bonds = bonds
    return mol


def write_pdb(mol: Molecule, path: str, remark: str = "") -> None:
    """Write molecule to PDB format."""
    with open(path, 'w') as f:
        if remark:
            f.write(f"REMARK {remark}\n")
        for atom in mol.atoms:
            rec = 'HETATM' if atom.is_hetatm else 'ATOM  '
            name = atom.name
            # PDB atom name formatting: 1-3 char names start at col 14
            if len(name) < 4:
                name_str = f" {name:<3s}"
            else:
                name_str = f"{name:<4s}"
            x, y, z = atom.position
            el = atom.element[:2] if len(atom.element) > 1 else f' {atom.element}'
            line = (
                f"{rec}{atom.index+1:5d} {name_str}{atom.residue_name:3s} "
                f"{atom.chain:1s}{atom.residue_id:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}"
                f"  1.00  0.00          {el:>2s}\n"
            )
            f.write(line)

        # CONECT records
        if mol.bonds:
            bond_map: dict = {}
            for b in mol.bonds:
                bond_map.setdefault(b.atom1, []).append(b.atom2)
                bond_map.setdefault(b.atom2, []).append(b.atom1)
            for src in sorted(bond_map):
                dsts = bond_map[src]
                for i in range(0, len(dsts), 4):
                    chunk = dsts[i:i+4]
                    ids = ''.join(f"{d+1:5d}" for d in chunk)
                    f.write(f"CONECT{src+1:5d}{ids}\n")

        f.write("END\n")
