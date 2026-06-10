"""
Simplified AMBER ff14SB / GAFF force field parameters.

VdW: Lennard-Jones epsilon (kcal/mol) and R_min/2 (Angstrom).
Partial charges: backbone from AMBER ff14SB, Gasteiger for unknowns.
"""
from __future__ import annotations
import numpy as np
from typing import Dict, List, Tuple
from ..core import Atom, Molecule

# fmt: off
# (epsilon kcal/mol, R_min/2 Angstrom) — AMBER GAFF / ff14SB values
VDW_PARAMS: Dict[str, Tuple[float, float]] = {
    # AMBER atom types
    'C':   (0.0860, 1.9080),   # sp2 carbon
    'CA':  (0.0860, 1.9080),   # aromatic C
    'CB':  (0.0860, 1.9080),
    'CC':  (0.0860, 1.9080),
    'CT':  (0.1094, 1.9080),   # sp3 carbon
    'CX':  (0.1094, 1.9080),
    'C8':  (0.1094, 1.9080),
    'N':   (0.1700, 1.8240),   # sp2 N
    'NA':  (0.1700, 1.8240),
    'N2':  (0.1700, 1.8240),
    'N3':  (0.1700, 1.8240),
    'NT':  (0.1700, 1.8240),
    'NY':  (0.1700, 1.8240),
    'O':   (0.2100, 1.6612),   # sp2 O
    'O2':  (0.2100, 1.6612),
    'OH':  (0.2104, 1.7210),   # hydroxyl O
    'OS':  (0.1700, 1.6837),   # ether O
    'S':   (0.2500, 2.0000),
    'SH':  (0.2500, 2.0000),
    'H':   (0.0157, 1.4870),
    'HC':  (0.0157, 1.4870),
    'HO':  (0.0000, 0.0001),
    'HN':  (0.0157, 1.4870),
    'HS':  (0.0157, 1.4870),
    'P':   (0.2000, 2.1000),
    'F':   (0.0610, 1.7500),
    'Cl':  (0.2650, 1.9480),
    'Br':  (0.3200, 2.2200),
    'I':   (0.4000, 2.3600),
    # Element fallbacks
    'CL':  (0.2650, 1.9480),
    'BR':  (0.3200, 2.2200),
}

# Element-based fallback VdW (for GAFF-assigned types)
_ELEM_VDW: Dict[str, Tuple[float, float]] = {
    'C':  (0.1094, 1.9080),
    'N':  (0.1700, 1.8240),
    'O':  (0.2100, 1.6612),
    'S':  (0.2500, 2.0000),
    'H':  (0.0157, 1.4870),
    'P':  (0.2000, 2.1000),
    'F':  (0.0610, 1.7500),
    'CL': (0.2650, 1.9480),
    'BR': (0.3200, 2.2200),
    'I':  (0.4000, 2.3600),
    'FE': (0.0130, 1.2000),
    'ZN': (0.0130, 1.0900),
    'MG': (0.8947, 0.7926),
    'CA': (0.4497, 1.7131),
}

# Gasteiger EN parameters: (a, b, c) for chi = a + b*q + c*q^2
# where q is the partial charge
_GASTEIGER: Dict[str, Tuple[float, float, float]] = {
    'C':  (7.98,  9.18,  1.88),
    'H':  (7.17,  6.24, -0.56),
    'O':  (8.50, 13.79,  1.00),
    'N':  (7.32, 12.36,  1.72),
    'S':  (6.60,  7.73,  0.50),
    'F':  (10.89, 17.37,  1.00),
    'CL': (8.29, 10.98,  1.00),
    'BR': (7.59,  9.81,  1.00),
    'I':  (6.98,  9.00,  1.00),
    'P':  (8.90,  8.24,  0.50),
}

# Backbone partial charges for standard amino acids (AMBER ff14SB)
# Format: residue_name -> {atom_name: charge}
AMBER_RESIDUE_CHARGES: Dict[str, Dict[str, float]] = {
    '_BB': {'N': -0.4157, 'H': 0.2719, 'CA': 0.0337, 'HA': 0.0823,
            'C': 0.5973, 'O': -0.5679},
    'ALA': {'CB': -0.1825, 'HB1': 0.0603, 'HB2': 0.0603, 'HB3': 0.0603},
    'GLY': {'HA2': 0.0698},
    'VAL': {'CB': 0.2985, 'HB': -0.0297, 'CG1': -0.3192, 'HG11': 0.0791,
            'HG12': 0.0791, 'HG13': 0.0791, 'CG2': -0.3192,
            'HG21': 0.0791, 'HG22': 0.0791, 'HG23': 0.0791},
    'LEU': {'CB': -0.2469, 'HB2': 0.0974, 'HB3': 0.0974, 'CG': 0.3502,
            'HG': -0.0361, 'CD1': -0.4121, 'HD11': 0.1000,
            'CD2': -0.4121, 'HD21': 0.1000},
    'ILE': {'CB': 0.1303, 'HB': 0.0187, 'CG1': -0.0430, 'CG2': -0.3204,
            'CD1': -0.0660},
    'PRO': {'CB': -0.2030, 'CG': 0.0189, 'CD': 0.0192},
    'PHE': {'CB': -0.1102, 'CG': 0.1108, 'CD1': -0.1801, 'CD2': -0.1801,
            'CE1': -0.2802, 'CE2': -0.2802, 'CZ': -0.1399, 'HZ': 0.1375},
    'TRP': {'CB': -0.0548, 'CG': -0.1415, 'CD1': -0.1638, 'CD2': 0.1243,
            'NE1': -0.3418, 'HE1': 0.3412, 'CE2': 0.1380, 'CE3': -0.2387,
            'CZ2': -0.2601, 'CZ3': -0.1972, 'CH2': -0.1134},
    'MET': {'CB': -0.0237, 'CG': 0.0018, 'SD': -0.2737, 'CE': -0.0536},
    'SER': {'CB': 0.2117, 'HB2': 0.0345, 'HB3': 0.0345,
            'OG': -0.6546, 'HG': 0.4275},
    'THR': {'CB': 0.3654, 'HB': 0.0043, 'OG1': -0.6761, 'HG1': 0.4102,
            'CG2': -0.2438},
    'CYS': {'CB': -0.1231, 'HB2': 0.1112, 'HB3': 0.1112,
            'SG': -0.3119, 'HG': 0.1933},
    'TYR': {'CB': -0.0152, 'CG': -0.0011, 'CD1': -0.1906, 'CD2': -0.1906,
            'CE1': -0.2341, 'CE2': -0.2341, 'CZ': 0.3226,
            'OH': -0.5579, 'HH': 0.3992},
    'ASN': {'CB': -0.2041, 'CG': 0.7130, 'OD1': -0.5931,
            'ND2': -0.9191, 'HD21': 0.4196, 'HD22': 0.4196},
    'GLN': {'CB': -0.0536, 'CG': -0.0645, 'CD': 0.6951,
            'OE1': -0.6086, 'NE2': -0.9407, 'HE21': 0.4251, 'HE22': 0.4251},
    'ASP': {'CB': -0.0302, 'CG': 0.7994, 'OD1': -0.8014, 'OD2': -0.8014},
    'GLU': {'CB': 0.0560, 'CG': 0.0136, 'CD': 0.8054,
            'OE1': -0.8188, 'OE2': -0.8188},
    'LYS': {'CB': -0.0094, 'CG': 0.0187, 'CD': -0.0479, 'CE': -0.0143,
            'NZ': -0.3854, 'HZ1': 0.3400, 'HZ2': 0.3400, 'HZ3': 0.3400},
    'ARG': {'CB': -0.0007, 'CG': 0.0390, 'CD': 0.0486, 'NE': -0.5295,
            'HE': 0.3456, 'CZ': 0.8076, 'NH1': -0.8627, 'NH2': -0.8627,
            'HH11': 0.4478, 'HH12': 0.4478, 'HH21': 0.4478, 'HH22': 0.4478},
    'HIS': {'CB': -0.0452, 'CG': -0.0266, 'ND1': -0.3811, 'HD1': 0.3649,
            'CD2': 0.1292, 'CE1': 0.2057, 'NE2': -0.5727},
}
# fmt: on


def get_vdw(atom_type: str, element: str) -> Tuple[float, float]:
    """Return (epsilon, R_min/2) for an atom."""
    if atom_type in VDW_PARAMS:
        return VDW_PARAMS[atom_type]
    el = element.upper()
    if el in _ELEM_VDW:
        return _ELEM_VDW[el]
    return (0.1094, 1.9080)   # generic carbon fallback


def assign_atom_types(mol: Molecule) -> None:
    """
    Assign simplified AMBER atom types to each atom in-place.
    Uses element + bond connectivity heuristics.
    """
    # Build connectivity count
    bond_count: Dict[int, int] = {a.index: 0 for a in mol.atoms}
    bond_partners: Dict[int, List[int]] = {a.index: [] for a in mol.atoms}
    for b in mol.bonds:
        bond_count[b.atom1] += 1
        bond_count[b.atom2] += 1
        bond_partners[b.atom1].append(b.atom2)
        bond_partners[b.atom2].append(b.atom1)

    atom_by_idx = {a.index: a for a in mol.atoms}

    for atom in mol.atoms:
        el = atom.element.upper()
        n_bonds = bond_count.get(atom.index, 0)
        partners = bond_partners.get(atom.index, [])
        partner_elements = [atom_by_idx[p].element.upper() for p in partners if p in atom_by_idx]

        if el == 'C':
            if n_bonds <= 2:
                atom.atom_type = 'C'    # sp (approximation)
            elif n_bonds == 3 or 'O' in partner_elements or 'N' in partner_elements:
                atom.atom_type = 'CA'   # aromatic / sp2
            else:
                atom.atom_type = 'CT'   # sp3
        elif el == 'N':
            if 'H' in partner_elements:
                atom.atom_type = 'N3' if n_bonds >= 3 else 'N'
            else:
                atom.atom_type = 'N2'
        elif el == 'O':
            if 'H' in partner_elements:
                atom.atom_type = 'OH'
            elif n_bonds == 1:
                atom.atom_type = 'O'
            else:
                atom.atom_type = 'OS'
        elif el == 'S':
            atom.atom_type = 'SH' if 'H' in partner_elements else 'S'
        elif el == 'H':
            if partner_elements:
                pe = partner_elements[0]
                if pe == 'C':
                    atom.atom_type = 'HC'
                elif pe == 'O':
                    atom.atom_type = 'HO'
                elif pe == 'N':
                    atom.atom_type = 'HN'
                elif pe == 'S':
                    atom.atom_type = 'HS'
                else:
                    atom.atom_type = 'H'
            else:
                atom.atom_type = 'H'
        elif el == 'P':
            atom.atom_type = 'P'
        elif el == 'F':
            atom.atom_type = 'F'
        elif el in ('CL', 'BR', 'I'):
            atom.atom_type = el.capitalize()
        else:
            atom.atom_type = el


def gasteiger_charges(mol: Molecule, n_iter: int = 6) -> None:
    """
    Compute Gasteiger-Marsili partial charges and assign to mol.atoms in-place.

    Works for small organic molecules. Proteins should use template charges.
    """
    atoms = mol.atoms
    n = len(atoms)
    q = np.zeros(n)

    # Build adjacency
    adj: List[List[int]] = [[] for _ in range(n)]
    idx_map = {a.index: i for i, a in enumerate(atoms)}
    for b in mol.bonds:
        i = idx_map.get(b.atom1)
        j = idx_map.get(b.atom2)
        if i is not None and j is not None:
            adj[i].append(j)
            adj[j].append(i)

    def chi(el: str, qi: float) -> float:
        a, b, c = _GASTEIGER.get(el.upper(), (7.0, 9.0, 1.5))
        return a + b * qi + c * qi * qi

    damp = 1.0
    for _ in range(n_iter):
        damp *= 0.5
        dq = np.zeros(n)
        for i, atom in enumerate(atoms):
            for j in adj[i]:
                el_i = atoms[i].element.upper()
                el_j = atoms[j].element.upper()
                chi_i = chi(el_i, q[i])
                chi_j = chi(el_j, q[j])
                denom = max(chi_i + chi_j, 0.001)
                transfer = (chi_j - chi_i) / denom * damp
                dq[i] += transfer
        q += dq

    # Assign back
    for i, atom in enumerate(atoms):
        atom.charge = float(q[i])


def apply_amber_charges(mol: Molecule) -> None:
    """Apply AMBER ff14SB backbone + sidechain charges to protein residues."""
    bb = AMBER_RESIDUE_CHARGES.get('_BB', {})
    for atom in mol.atoms:
        if not atom.is_hetatm:
            # Backbone first
            atom.charge = bb.get(atom.name, 0.0)
            # Sidechain override
            sc = AMBER_RESIDUE_CHARGES.get(atom.residue_name, {})
            if atom.name in sc:
                atom.charge = sc[atom.name]
