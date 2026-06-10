"""
Structural analysis: RMSD, RMSF, Rg, secondary structure, Ramachandran.
"""
from __future__ import annotations
import numpy as np
from typing import List, Optional, Dict, Tuple
from ..core import Molecule


# ─── RMSD ────────────────────────────────────────────────────────────────────

def rmsd(pos1: np.ndarray, pos2: np.ndarray) -> float:
    """Root-mean-square deviation between two (N,3) coordinate arrays."""
    assert pos1.shape == pos2.shape, "Coordinate arrays must have same shape"
    diff = pos1 - pos2
    return float(np.sqrt((diff ** 2).sum() / len(pos1)))


def kabsch_rotation(P: np.ndarray, Q: np.ndarray) -> np.ndarray:
    """
    Kabsch algorithm: find rotation matrix R that minimizes RMSD of P->Q.

    Both P, Q are (N, 3) centered at origin.
    Returns 3x3 rotation matrix.
    """
    H = P.T @ Q
    U, S, Vt = np.linalg.svd(H)
    d = np.linalg.det(Vt.T @ U.T)
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    return R


def aligned_rmsd(
    mol1: Molecule,
    mol2: Molecule,
    atom_names: Optional[List[str]] = None,
) -> Tuple[float, np.ndarray]:
    """
    Align mol2 onto mol1 using selected atoms, return (RMSD, aligned_positions).

    Parameters
    ----------
    mol1, mol2 : Molecule
        Must have the same number of atoms (or matching by selected names).
    atom_names : list of str, optional
        Atom names to use for alignment (e.g. ['CA']). Default: all heavy atoms.

    Returns
    -------
    (rmsd_value, aligned_full_coords)
    """
    if atom_names is None:
        # Use all non-hydrogen heavy atoms
        idx1 = [i for i, a in enumerate(mol1.atoms) if a.element.upper() != 'H']
        idx2 = [i for i, a in enumerate(mol2.atoms) if a.element.upper() != 'H']
    else:
        idx1 = [i for i, a in enumerate(mol1.atoms) if a.name in atom_names]
        idx2 = [i for i, a in enumerate(mol2.atoms) if a.name in atom_names]

    n = min(len(idx1), len(idx2))
    P = mol2.positions[idx2[:n]]
    Q = mol1.positions[idx1[:n]]

    P_c = P - P.mean(0)
    Q_c = Q - Q.mean(0)

    R = kabsch_rotation(P_c, Q_c)

    # Apply to full mol2 coords
    full_pos = mol2.positions.copy()
    com2 = mol2.positions[idx2[:n]].mean(0)
    com1 = mol1.positions[idx1[:n]].mean(0)
    full_pos = (full_pos - com2) @ R.T + com1

    aligned_ref = (P - P.mean(0)) @ R.T + Q.mean(0)
    rms = float(np.sqrt(((aligned_ref - Q) ** 2).sum() / n))
    return rms, full_pos


def rmsf(trajectory: List[np.ndarray], indices: Optional[List[int]] = None) -> np.ndarray:
    """
    Per-atom root-mean-square fluctuation from a trajectory.

    Parameters
    ----------
    trajectory : list of (N, 3) arrays
    indices : atom indices to compute (default: all)

    Returns
    -------
    (M,) array of RMSF values in Angstrom
    """
    if not trajectory:
        return np.array([])
    traj = np.array(trajectory)           # (T, N, 3)
    if indices is not None:
        traj = traj[:, indices, :]
    mean_pos = traj.mean(0)               # (M, 3)
    dev = traj - mean_pos                 # (T, M, 3)
    return np.sqrt((dev ** 2).sum(-1).mean(0))   # (M,)


def radius_of_gyration(mol: Molecule, heavy_only: bool = True) -> float:
    """Mass-weighted radius of gyration (Angstrom)."""
    atoms = mol.atoms
    if heavy_only:
        atoms = [a for a in atoms if a.element.upper() != 'H']
    if not atoms:
        return 0.0
    pos = np.array([a.position for a in atoms])
    mass = np.array([a.mass for a in atoms])
    com = (pos * mass[:, None]).sum(0) / mass.sum()
    r2 = ((pos - com) ** 2).sum(1)
    return float(np.sqrt((mass * r2).sum() / mass.sum()))


def distance_matrix(mol1: Molecule, mol2: Optional[Molecule] = None) -> np.ndarray:
    """
    Compute pairwise distance matrix.
    If mol2 is None, compute intra-mol1 distances.
    Returns (N1, N2) array.
    """
    pos1 = mol1.positions
    pos2 = mol2.positions if mol2 is not None else pos1
    # Efficient broadcasting
    diff = pos1[:, None, :] - pos2[None, :, :]   # (N1, N2, 3)
    return np.sqrt((diff ** 2).sum(-1))


# ─── Ramachandran angles ─────────────────────────────────────────────────────

def _dihedral(p0, p1, p2, p3) -> float:
    """Compute dihedral angle (radians) from 4 points."""
    b1 = p1 - p0
    b2 = p2 - p1
    b3 = p3 - p2
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)
    m1 = np.cross(n1, b2 / (np.linalg.norm(b2) + 1e-12))
    x = np.dot(n1, n2)
    y = np.dot(m1, n2)
    return np.arctan2(y, x)


def ramachandran_angles(mol: Molecule) -> Dict[int, Dict[str, float]]:
    """
    Compute phi/psi backbone dihedral angles for each residue.

    Returns dict: {residue_id: {'phi': deg, 'psi': deg, 'resname': str}}
    """
    # Index atoms by (residue_id, name)
    idx_map: Dict[Tuple[int, str], int] = {}
    for i, a in enumerate(mol.atoms):
        if a.element.upper() != 'H':
            idx_map[(a.residue_id, a.name)] = i

    res_ids = mol.residue_ids()
    pos = mol.positions
    result = {}

    for k, rid in enumerate(res_ids):
        entry: Dict[str, float] = {
            'resname': next((a.residue_name for a in mol.atoms
                            if a.residue_id == rid), '???'),
        }

        # phi: C(i-1) - N(i) - CA(i) - C(i)
        if k > 0:
            prev_rid = res_ids[k - 1]
            C_prev = idx_map.get((prev_rid, 'C'))
            N = idx_map.get((rid, 'N'))
            CA = idx_map.get((rid, 'CA'))
            C = idx_map.get((rid, 'C'))
            if all(x is not None for x in [C_prev, N, CA, C]):
                phi = _dihedral(pos[C_prev], pos[N], pos[CA], pos[C])
                entry['phi'] = float(np.degrees(phi))

        # psi: N(i) - CA(i) - C(i) - N(i+1)
        if k < len(res_ids) - 1:
            next_rid = res_ids[k + 1]
            N = idx_map.get((rid, 'N'))
            CA = idx_map.get((rid, 'CA'))
            C = idx_map.get((rid, 'C'))
            N_next = idx_map.get((next_rid, 'N'))
            if all(x is not None for x in [N, CA, C, N_next]):
                psi = _dihedral(pos[N], pos[CA], pos[C], pos[N_next])
                entry['psi'] = float(np.degrees(psi))

        if len(entry) > 1:     # has at least phi or psi
            result[rid] = entry

    return result


# ─── Secondary structure (DSSP-like) ─────────────────────────────────────────

def secondary_structure(mol: Molecule) -> Dict[int, str]:
    """
    Assign secondary structure to each residue using H-bond + dihedral criteria.

    Returns dict: {residue_id: 'H' (helix), 'E' (sheet), 'T' (turn), 'C' (coil)}
    """
    res_ids = sorted(set(a.residue_id for a in mol.atoms))
    n = len(res_ids)
    rid_to_idx = {rid: i for i, rid in enumerate(res_ids)}

    pos_map: Dict[Tuple[int, str], np.ndarray] = {}
    for a in mol.atoms:
        pos_map[(a.residue_id, a.name)] = a.position

    # Collect backbone atoms
    N_pos = {rid: pos_map.get((rid, 'N')) for rid in res_ids}
    H_pos = {rid: pos_map.get((rid, 'H')) or pos_map.get((rid, 'HN')) for rid in res_ids}
    CA_pos = {rid: pos_map.get((rid, 'CA')) for rid in res_ids}
    C_pos = {rid: pos_map.get((rid, 'C')) for rid in res_ids}
    O_pos = {rid: pos_map.get((rid, 'O')) for rid in res_ids}

    # H-bond detection: N-H(i) ... O=C(j)
    # donor: N-H, acceptor: C=O
    # DSSP criterion: E = q1*q2*(1/r_ON + 1/r_CH - 1/r_OH - 1/r_CN) < -0.5 kcal/mol
    # Simplified: just distance + angle check
    def is_hbond(i: int, j: int) -> bool:
        ri, rj = res_ids[i], res_ids[j]
        if abs(i - j) < 2:
            return False
        h = H_pos.get(ri) or N_pos.get(ri)
        n = N_pos.get(ri)
        o = O_pos.get(rj)
        c = C_pos.get(rj)
        if h is None or o is None:
            return False
        d_HO = np.linalg.norm(h - o)
        if d_HO > 2.5:
            return False
        if n is not None and c is not None:
            # Check N-H...O angle
            vec_NH = h - n
            vec_HO = o - h
            cos_a = np.dot(vec_NH, vec_HO) / (
                np.linalg.norm(vec_NH) * np.linalg.norm(vec_HO) + 1e-12)
            if cos_a < 0.0:    # angle > 90 deg
                return False
        return True

    # Detect alpha helices (i, i+4 H-bonds) and beta sheets (extended)
    ss = {rid: 'C' for rid in res_ids}

    # Alpha helix: 4 consecutive (i->i+4) H-bonds
    helix_set = set()
    for i in range(n - 4):
        if is_hbond(i + 4, i) and is_hbond(i + 3, i - 1 if i > 0 else i):
            for k in range(i, i + 5):
                helix_set.add(k)

    # Mark helices
    for k in helix_set:
        ss[res_ids[k]] = 'H'

    # Beta sheet (antiparallel): i-j H-bond pairs with i<j, |i-j|>4
    for i in range(n):
        for j in range(i + 4, min(i + 30, n)):
            if is_hbond(i, j) and is_hbond(j, i):
                if ss[res_ids[i]] != 'H':
                    ss[res_ids[i]] = 'E'
                if ss[res_ids[j]] != 'H':
                    ss[res_ids[j]] = 'E'

    # Turn: any 3-6 residue H-bond (i to i+2, i+3, i+4)
    for i in range(n):
        for delta in (2, 3, 4):
            j = i + delta
            if j < n and is_hbond(j, i):
                for k in range(i, j + 1):
                    if ss[res_ids[k]] == 'C':
                        ss[res_ids[k]] = 'T'

    return ss
