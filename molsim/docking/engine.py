"""
Molecular docking engine.

Algorithm: Monte Carlo random search + L-BFGS-B local minimization.
Scoring: simplified Vina-style (VdW + H-bond + hydrophobic + electrostatic).
"""
from __future__ import annotations
import numpy as np
from dataclasses import dataclass, field
from typing import List, Optional, Tuple
from scipy.optimize import minimize as scipy_minimize

from ..core import Molecule, Bond
from ..energy.calculator import EnergyCalculator, _vdw_arrays, COULOMB_CONST
from ..forcefield.params import assign_atom_types, gasteiger_charges


@dataclass
class DockingResult:
    """A single docked pose."""
    rank: int
    score: float              # kcal/mol (lower = better)
    vdw: float
    elec: float
    hbond: float
    positions: np.ndarray     # (N_lig, 3) Angstrom
    rmsd_from_top: float = 0.0

    def __repr__(self) -> str:
        return (f"DockingResult(rank={self.rank}, score={self.score:.2f}, "
                f"vdw={self.vdw:.2f}, elec={self.elec:.2f}, hbond={self.hbond:.2f})")


def rotation_matrix(axis: np.ndarray, angle: float) -> np.ndarray:
    """Rodrigues rotation formula: rotate angle (rad) around axis."""
    axis = axis / (np.linalg.norm(axis) + 1e-12)
    c, s = np.cos(angle), np.sin(angle)
    t = 1 - c
    x, y, z = axis
    return np.array([
        [t*x*x + c,   t*x*y - s*z, t*x*z + s*y],
        [t*x*y + s*z, t*y*y + c,   t*y*z - s*x],
        [t*x*z - s*y, t*y*z + s*x, t*z*z + c  ],
    ])


def random_rotation() -> np.ndarray:
    """Uniform random rotation matrix via quaternion."""
    q = np.random.randn(4)
    q /= np.linalg.norm(q)
    w, x, y, z = q
    return np.array([
        [1-2*(y*y+z*z), 2*(x*y-w*z),   2*(x*z+w*y)  ],
        [2*(x*y+w*z),   1-2*(x*x+z*z), 2*(y*z-w*x)  ],
        [2*(x*z-w*y),   2*(y*z+w*x),   1-2*(x*x+y*y)],
    ])


def find_rotatable_bonds(mol: Molecule) -> List[Tuple[int, int]]:
    """
    Identify rotatable bonds in a small molecule.
    A bond is rotatable if it is:
    - single bond
    - not in a ring (approximated as: both atoms have degree >= 2)
    - not terminal (both atoms have > 1 heavy-atom neighbor)
    """
    # Count heavy-atom neighbors
    heavy_degree: dict = {}
    for b in mol.bonds:
        e1 = mol.atoms[b.atom1].element.upper() if b.atom1 < len(mol.atoms) else 'H'
        e2 = mol.atoms[b.atom2].element.upper() if b.atom2 < len(mol.atoms) else 'H'
        if e2 != 'H':
            heavy_degree[b.atom1] = heavy_degree.get(b.atom1, 0) + 1
        if e1 != 'H':
            heavy_degree[b.atom2] = heavy_degree.get(b.atom2, 0) + 1

    rotatable = []
    for b in mol.bonds:
        if b.order != 1:
            continue
        e1 = mol.atoms[b.atom1].element.upper() if b.atom1 < len(mol.atoms) else 'H'
        e2 = mol.atoms[b.atom2].element.upper() if b.atom2 < len(mol.atoms) else 'H'
        if e1 == 'H' or e2 == 'H':
            continue
        if heavy_degree.get(b.atom1, 0) < 2 or heavy_degree.get(b.atom2, 0) < 2:
            continue
        rotatable.append((b.atom1, b.atom2))

    return rotatable


def _atoms_on_side(mol: Molecule, bond_atom1: int, bond_atom2: int) -> List[int]:
    """BFS to get atom indices on bond_atom2's side when bond is cut."""
    adj: dict = {}
    for b in mol.bonds:
        adj.setdefault(b.atom1, []).append(b.atom2)
        adj.setdefault(b.atom2, []).append(b.atom1)

    visited = {bond_atom1}
    queue = [bond_atom2]
    side = []
    while queue:
        node = queue.pop(0)
        if node in visited:
            continue
        visited.add(node)
        side.append(node)
        queue.extend(adj.get(node, []))
    return side


def _hbond_score(
    rec: Molecule,
    lig_coords: np.ndarray,
    lig: Molecule,
    cutoff: float = 3.5,
) -> float:
    """
    Simple H-bond scoring: -1.0 kcal/mol per donor-acceptor pair within cutoff.
    Donors: N-H, O-H in receptor. Acceptors: N, O in ligand (and vice versa).
    """
    HBOND_DONORS = {'N', 'O'}
    HBOND_ACCEPTORS = {'N', 'O', 'S', 'F'}

    score = 0.0
    pos_r = rec.positions

    # Receptor donors/acceptors
    rec_donor_idx = [i for i, a in enumerate(rec.atoms)
                     if a.element.upper() in HBOND_DONORS]
    rec_accept_idx = [i for i, a in enumerate(rec.atoms)
                      if a.element.upper() in HBOND_ACCEPTORS]

    # Ligand donors/acceptors
    lig_donor_idx = [i for i, a in enumerate(lig.atoms)
                     if a.element.upper() in HBOND_DONORS]
    lig_accept_idx = [i for i, a in enumerate(lig.atoms)
                      if a.element.upper() in HBOND_ACCEPTORS]

    def _count_pairs(rec_idx, lig_idx, pos_l):
        if not rec_idx or not lig_idx:
            return 0.0
        rpos = pos_r[rec_idx]   # (M, 3)
        lpos = pos_l[lig_idx]   # (K, 3)
        s = 0.0
        for lp in lpos:
            dists = np.linalg.norm(rpos - lp, axis=1)
            close = (dists < cutoff).sum()
            s -= close * 1.0
        return s

    score += _count_pairs(rec_donor_idx, lig_accept_idx, lig_coords)
    score += _count_pairs(rec_accept_idx, lig_donor_idx, lig_coords)
    return score


class DockingEngine:
    """
    Monte Carlo docking engine.

    Usage
    -----
    engine = DockingEngine(receptor, ligand)
    engine.set_binding_site(center, radius=12.0)
    results = engine.dock(n_poses=100, n_refine=10)
    """

    def __init__(
        self,
        receptor: Molecule,
        ligand: Molecule,
        prepare: bool = True,
    ):
        self.receptor = receptor
        self.ligand = ligand.clone()

        if prepare:
            assign_atom_types(self.receptor)
            assign_atom_types(self.ligand)
            # Use AMBER charges for protein, Gasteiger for ligand
            from ..forcefield.params import apply_amber_charges
            apply_amber_charges(self.receptor)
            gasteiger_charges(self.ligand)

        self._calc = EnergyCalculator(receptor=self.receptor, ligand=self.ligand)
        self._site_center: Optional[np.ndarray] = None
        self._site_radius: float = 15.0
        self._rot_bonds: List[Tuple[int, int]] = find_rotatable_bonds(self.ligand)

    def set_binding_site(self, center: np.ndarray, radius: float = 12.0) -> None:
        """Define binding site by center and radius (Angstrom)."""
        self._site_center = np.asarray(center, dtype=float)
        self._site_radius = radius

    def set_binding_site_from_residues(self, residue_ids: List[int]) -> None:
        """Define binding site from a list of receptor residue IDs."""
        atoms = [a for a in self.receptor.atoms if a.residue_id in residue_ids]
        if not atoms:
            raise ValueError("No atoms found for given residue IDs.")
        pos = np.array([a.position for a in atoms])
        self._site_center = pos.mean(0)
        self._site_radius = float(np.linalg.norm(pos - self._site_center, axis=1).max()) + 5.0

    def _random_pose(self) -> np.ndarray:
        """Generate a random ligand pose within the binding site."""
        lig = self.ligand
        coords = lig.positions.copy()
        com = coords.mean(0)
        coords -= com

        # Random rotation
        R = random_rotation()
        coords = coords @ R.T

        # Random translation within sphere
        r = self._site_radius * np.cbrt(np.random.rand())
        direction = np.random.randn(3)
        direction /= np.linalg.norm(direction) + 1e-12
        t = self._site_center + r * direction
        coords += t

        return coords

    def _apply_torsions(self, coords: np.ndarray, torsions: np.ndarray) -> np.ndarray:
        """Apply torsion angles to rotatable bonds."""
        coords = coords.copy()
        for k, (a1, a2) in enumerate(self._rot_bonds[:len(torsions)]):
            if a1 >= len(coords) or a2 >= len(coords):
                continue
            angle = torsions[k]
            axis = coords[a2] - coords[a1]
            R = rotation_matrix(axis, angle)
            moving = _atoms_on_side(self.ligand, a1, a2)
            pivot = coords[a2]
            for m in moving:
                if m < len(coords):
                    coords[m] = R @ (coords[m] - pivot) + pivot
        return coords

    def _score_pose(self, coords: np.ndarray) -> Tuple[float, float, float, float]:
        """Return (total, vdw, elec, hbond)."""
        e = self._calc.score(coords)
        # Decompose
        pos_r = self._calc._pos_r
        eps_r = self._calc._eps_r_arr
        rmin_r = self._calc._rmin_r_arr
        q_r = self._calc._q_r
        cutoff2 = self._calc.cutoff ** 2
        eps_rr = self._calc.epsilon_r

        eps_l, rmin_l = _vdw_arrays(self.ligand)
        q_l = self.ligand.charges

        e_vdw = 0.0
        e_elec = 0.0
        for j in range(len(coords)):
            diff = pos_r - coords[j]
            r2 = (diff * diff).sum(1)
            mask = r2 < cutoff2
            rv = np.sqrt(r2[mask])
            rv = np.maximum(rv, 0.5)
            eps_ij = np.sqrt(eps_r[mask] * eps_l[j])
            rmin_ij = rmin_r[mask] + rmin_l[j]
            x = rmin_ij / rv
            e_vdw += (eps_ij * (x**12 - 2.0 * x**6)).sum()
            e_elec += (COULOMB_CONST * q_r[mask] * q_l[j] / (eps_rr * rv)).sum()

        hb = _hbond_score(self.receptor, coords, self.ligand)
        return float(e_vdw + e_elec + hb), float(e_vdw), float(e_elec), float(hb)

    def _optimize_pose(self, coords: np.ndarray, max_iter: int = 200) -> np.ndarray:
        """Local L-BFGS-B minimization of a pose."""
        x0 = coords.flatten()

        def func_grad(x):
            c = x.reshape(-1, 3)
            e, g = self._calc.score_and_grad(c)
            return e, g.flatten()

        result = scipy_minimize(
            func_grad, x0, method='L-BFGS-B', jac=True,
            options={'maxiter': max_iter, 'ftol': 1e-5, 'gtol': 1e-4}
        )
        return result.x.reshape(-1, 3)

    def dock(
        self,
        n_random: int = 2000,
        n_refine: int = 20,
        seed: Optional[int] = None,
        verbose: bool = False,
    ) -> List[DockingResult]:
        """
        Run docking search.

        Parameters
        ----------
        n_random : int
            Number of random poses to generate and score.
        n_refine : int
            Top-N poses to locally minimize.
        seed : int, optional
            Random seed for reproducibility.
        verbose : bool
            Print progress.

        Returns
        -------
        List[DockingResult]
            Sorted by score (best first).
        """
        if seed is not None:
            np.random.seed(seed)

        if self._site_center is None:
            # Default: center of receptor
            self._site_center = self.receptor.positions.mean(0)

        if verbose:
            print(f"Sampling {n_random} random poses...")

        # Phase 1: Random sampling
        scored: List[Tuple[float, np.ndarray]] = []
        for i in range(n_random):
            coords = self._random_pose()
            e, *_ = self._score_pose(coords)
            # Reject severe clashes
            if e < 5000:
                scored.append((e, coords))

        if not scored:
            return []

        scored.sort(key=lambda x: x[0])
        top = scored[:n_refine]

        if verbose:
            print(f"Refining top {len(top)} poses...")

        # Phase 2: Local minimization
        results: List[DockingResult] = []
        for rank, (e_raw, coords) in enumerate(top):
            opt_coords = self._optimize_pose(coords)
            total, vdw, elec, hb = self._score_pose(opt_coords)
            results.append(DockingResult(
                rank=rank + 1,
                score=total,
                vdw=vdw,
                elec=elec,
                hbond=hb,
                positions=opt_coords,
            ))

        results.sort(key=lambda r: r.score)
        for i, r in enumerate(results):
            r.rank = i + 1

        # Compute RMSD from top pose
        top_pos = results[0].positions if results else None
        if top_pos is not None:
            for r in results[1:]:
                diff = r.positions - top_pos
                r.rmsd_from_top = float(np.sqrt((diff ** 2).sum(1).mean()))

        if verbose:
            print(f"\nTop docking results:")
            for r in results[:5]:
                print(f"  Rank {r.rank}: score={r.score:.2f}, "
                      f"vdw={r.vdw:.2f}, elec={r.elec:.2f}, hbond={r.hbond:.2f}")

        return results

    def apply_pose(self, result: DockingResult) -> Molecule:
        """Return ligand molecule with coordinates from a docking result."""
        lig = self.ligand.clone()
        lig.positions = result.positions
        return lig
