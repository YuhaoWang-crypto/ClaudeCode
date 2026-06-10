"""
Molecular contacts: hydrogen bonds, hydrophobic contacts, residue-residue pairs.
"""
from __future__ import annotations
import numpy as np
from typing import List, Dict, Tuple, Optional
from ..core import Molecule


def find_hbonds(
    donor_mol: Molecule,
    acceptor_mol: Optional[Molecule] = None,
    d_cutoff: float = 3.5,
    angle_cutoff: float = 120.0,
) -> List[Dict]:
    """
    Find hydrogen bonds between donor and acceptor molecules.
    If acceptor_mol is None, finds intramolecular H-bonds within donor_mol.

    Returns list of dicts: {donor_atom, acceptor_atom, distance, angle}.
    """
    if acceptor_mol is None:
        acceptor_mol = donor_mol

    DONORS = {'N', 'O', 'S'}
    ACCEPTORS = {'N', 'O', 'S', 'F'}
    H_ELEMENTS = {'H'}

    # Build donor H list: (H_atom, donor_atom)
    donor_Hs: List[Tuple] = []
    adj: Dict[int, List[int]] = {}
    idx_map = {a.index: i for i, a in enumerate(donor_mol.atoms)}
    for b in donor_mol.bonds:
        adj.setdefault(b.atom1, []).append(b.atom2)
        adj.setdefault(b.atom2, []).append(b.atom1)

    for i, atom in enumerate(donor_mol.atoms):
        if atom.element.upper() in H_ELEMENTS:
            for nb_idx in adj.get(atom.index, []):
                nb_pos = idx_map.get(nb_idx)
                if nb_pos is not None:
                    nb = donor_mol.atoms[nb_pos]
                    if nb.element.upper() in DONORS:
                        donor_Hs.append((i, nb_pos, atom, nb))

    # Acceptor atoms
    acceptors = [(i, a) for i, a in enumerate(acceptor_mol.atoms)
                 if a.element.upper() in ACCEPTORS]

    pos_acc = acceptor_mol.positions
    hbonds: List[Dict] = []

    for h_idx, d_idx, h_atom, d_atom in donor_Hs:
        h_pos = h_atom.position
        d_pos = d_atom.position
        for acc_idx, acc_atom in acceptors:
            # Skip self (intramolecular)
            if acceptor_mol is donor_mol and acc_idx in (h_idx, d_idx):
                continue
            acc_pos = pos_acc[acc_idx]
            dist_H_A = np.linalg.norm(h_pos - acc_pos)
            if dist_H_A > d_cutoff:
                continue
            # D-H...A angle
            vec_DH = h_pos - d_pos
            vec_HA = acc_pos - h_pos
            cos_angle = np.dot(vec_DH, vec_HA) / (
                np.linalg.norm(vec_DH) * np.linalg.norm(vec_HA) + 1e-12)
            angle = float(np.degrees(np.arccos(np.clip(cos_angle, -1, 1))))
            if angle < angle_cutoff:
                continue
            hbonds.append({
                'donor_atom': d_atom,
                'h_atom': h_atom,
                'acceptor_atom': acc_atom,
                'distance_H_A': float(dist_H_A),
                'distance_D_A': float(np.linalg.norm(d_pos - acc_pos)),
                'angle_DHA': angle,
            })

    return hbonds


def find_hydrophobic_contacts(
    mol1: Molecule,
    mol2: Optional[Molecule] = None,
    cutoff: float = 4.5,
) -> List[Dict]:
    """
    Find hydrophobic contacts (C-C pairs within cutoff).
    If mol2 is None, finds intramolecular contacts in mol1.
    """
    HYDROPHOBIC = {'C', 'S'}

    def hydro_atoms(mol):
        return [(i, a) for i, a in enumerate(mol.atoms)
                if a.element.upper() in HYDROPHOBIC]

    if mol2 is None:
        mol2 = mol1

    h1 = hydro_atoms(mol1)
    h2 = hydro_atoms(mol2)

    contacts = []
    pos2 = np.array([mol2.atoms[i].position for i, _ in h2])
    for i1, a1 in h1:
        dists = np.linalg.norm(pos2 - a1.position, axis=1)
        for k, (i2, a2) in enumerate(h2):
            if mol1 is mol2 and i1 >= i2:
                continue
            if dists[k] <= cutoff:
                contacts.append({
                    'atom1': a1,
                    'atom2': a2,
                    'distance': float(dists[k]),
                })
    return contacts


def residue_contacts(
    mol1: Molecule,
    mol2: Molecule,
    cutoff: float = 5.0,
    heavy_only: bool = True,
) -> List[Tuple[int, int, float]]:
    """
    Find all residue pairs between mol1 and mol2 with any atom within cutoff.

    Returns list of (residue_id_1, residue_id_2, min_distance).
    """
    atoms1 = [(a.residue_id, a.position) for a in mol1.atoms
              if not heavy_only or a.element.upper() != 'H']
    atoms2 = [(a.residue_id, a.position) for a in mol2.atoms
              if not heavy_only or a.element.upper() != 'H']

    if not atoms1 or not atoms2:
        return []

    rids1 = [x[0] for x in atoms1]
    pos1 = np.array([x[1] for x in atoms1])
    rids2 = [x[0] for x in atoms2]
    pos2 = np.array([x[1] for x in atoms2])

    # Min distance matrix per residue pair
    pair_min: Dict[Tuple[int, int], float] = {}
    for i, (rid1, p1) in enumerate(zip(rids1, pos1)):
        dists = np.linalg.norm(pos2 - p1, axis=1)
        for j, (rid2, d) in enumerate(zip(rids2, dists)):
            if d <= cutoff:
                key = (rid1, rid2)
                if key not in pair_min or d < pair_min[key]:
                    pair_min[key] = d

    return [(r1, r2, d) for (r1, r2), d in sorted(pair_min.items())]
