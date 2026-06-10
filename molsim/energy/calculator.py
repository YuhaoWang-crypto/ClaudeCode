"""
Non-bonded energy calculation (LJ + Coulomb).

All energies in kcal/mol, distances in Angstrom.
"""
from __future__ import annotations
import numpy as np
from typing import Optional, Tuple, List
from ..core import Molecule
from ..forcefield.params import get_vdw

COULOMB_CONST = 332.0636   # e^2/(Angstrom * kcal/mol)
DEFAULT_CUTOFF = 12.0       # Angstrom
DEFAULT_EPSILON_R = 4.0     # distance-dependent dielectric factor


def _vdw_arrays(mol: Molecule) -> Tuple[np.ndarray, np.ndarray]:
    """Return epsilon and R_min/2 arrays for all atoms."""
    eps = np.empty(mol.n_atoms)
    rmin = np.empty(mol.n_atoms)
    for i, atom in enumerate(mol.atoms):
        e, r = get_vdw(atom.atom_type, atom.element)
        eps[i] = e
        rmin[i] = r
    return eps, rmin


def nonbonded_energy(
    mol: Molecule,
    cutoff: float = DEFAULT_CUTOFF,
    epsilon_r: float = DEFAULT_EPSILON_R,
    exclude_12: bool = True,
    exclude_13: bool = True,
) -> Tuple[float, float]:
    """
    Compute intramolecular VdW and Coulomb energy.

    Returns (E_vdw, E_elec) in kcal/mol.
    """
    pos = mol.positions
    n = mol.n_atoms
    eps_arr, rmin_arr = _vdw_arrays(mol)
    charges = mol.charges

    # Build exclusion pairs (1-2 and 1-3 bonds)
    excluded: set = set()
    if exclude_12 or exclude_13:
        adj: dict = {i: set() for i in range(n)}
        idx_map = {a.index: i for i, a in enumerate(mol.atoms)}
        for b in mol.bonds:
            i, j = idx_map.get(b.atom1, -1), idx_map.get(b.atom2, -1)
            if i >= 0 and j >= 0:
                adj[i].add(j)
                adj[j].add(i)
        if exclude_12:
            for i, nbrs in adj.items():
                for j in nbrs:
                    if i < j:
                        excluded.add((i, j))
        if exclude_13:
            for i, nbrs in adj.items():
                for j in nbrs:
                    for k in adj[j]:
                        if k != i and min(i, k) < max(i, k):
                            excluded.add((min(i, k), max(i, k)))

    e_vdw = 0.0
    e_elec = 0.0

    for i in range(n - 1):
        diff = pos[i+1:] - pos[i]
        r2 = (diff * diff).sum(1)
        mask = r2 < cutoff * cutoff

        indices_j = np.where(mask)[0] + i + 1

        for jj, j in enumerate(indices_j):
            if (i, j) in excluded:
                continue
            r = np.sqrt(r2[j - i - 1])
            if r < 0.5:     # hard clash guard
                r = 0.5

            # Lorentz-Berthelot combining rules
            eps_ij = np.sqrt(eps_arr[i] * eps_arr[j])
            rmin_ij = rmin_arr[i] + rmin_arr[j]
            x = rmin_ij / r
            e_vdw += eps_ij * (x**12 - 2.0 * x**6)

            # Coulomb with distance-dependent dielectric eps_r * r
            if charges[i] != 0.0 or charges[j] != 0.0:
                e_elec += COULOMB_CONST * charges[i] * charges[j] / (epsilon_r * r)

    return e_vdw, e_elec


def interaction_energy(
    receptor: Molecule,
    ligand: Molecule,
    cutoff: float = DEFAULT_CUTOFF,
    epsilon_r: float = DEFAULT_EPSILON_R,
) -> dict:
    """
    Compute receptor-ligand interaction energy.

    Returns dict with keys: vdw, elec, total (all kcal/mol).
    """
    pos_r = receptor.positions
    pos_l = ligand.positions
    eps_r_arr, rmin_r_arr = _vdw_arrays(receptor)
    eps_l_arr, rmin_l_arr = _vdw_arrays(ligand)
    q_r = receptor.charges
    q_l = ligand.charges

    e_vdw = 0.0
    e_elec = 0.0

    for j in range(ligand.n_atoms):
        diff = pos_r - pos_l[j]
        r2 = (diff * diff).sum(1)
        mask = r2 < cutoff * cutoff
        r_vals = np.sqrt(r2[mask])
        r_vals = np.maximum(r_vals, 0.5)

        eps_ij = np.sqrt(eps_r_arr[mask] * eps_l_arr[j])
        rmin_ij = rmin_r_arr[mask] + rmin_l_arr[j]
        x = rmin_ij / r_vals
        e_vdw += (eps_ij * (x**12 - 2.0 * x**6)).sum()

        e_elec += (COULOMB_CONST * q_r[mask] * q_l[j]
                   / (epsilon_r * r_vals)).sum()

    return {'vdw': float(e_vdw), 'elec': float(e_elec), 'total': float(e_vdw + e_elec)}


class EnergyCalculator:
    """
    Convenient wrapper for energy + gradient computation.

    Used by the minimizer and docking engine.
    """

    def __init__(
        self,
        receptor: Optional[Molecule] = None,
        ligand: Optional[Molecule] = None,
        cutoff: float = DEFAULT_CUTOFF,
        epsilon_r: float = DEFAULT_EPSILON_R,
    ):
        self.receptor = receptor
        self.ligand = ligand
        self.cutoff = cutoff
        self.epsilon_r = epsilon_r

        if receptor is not None:
            self._pos_r = receptor.positions
            self._eps_r_arr, self._rmin_r_arr = _vdw_arrays(receptor)
            self._q_r = receptor.charges

    def score(self, coords: np.ndarray) -> float:
        """
        Score ligand at given (N_lig, 3) coordinates against pre-loaded receptor.
        Fast path for docking inner loop.
        """
        pos_r = self._pos_r
        eps_r = self._eps_r_arr
        rmin_r = self._rmin_r_arr
        q_r = self._q_r
        cutoff2 = self.cutoff * self.cutoff
        eps_rr = self.epsilon_r

        # Pre-compute ligand VdW and charges once
        if self.ligand is not None:
            eps_l, rmin_l = _vdw_arrays(self.ligand)
            q_l = self.ligand.charges
        else:
            eps_l = np.ones(len(coords)) * 0.1
            rmin_l = np.ones(len(coords)) * 1.9
            q_l = np.zeros(len(coords))

        e = 0.0
        for j in range(len(coords)):
            diff = pos_r - coords[j]
            r2 = (diff * diff).sum(1)
            mask = r2 < cutoff2
            if not mask.any():
                continue
            rv = np.sqrt(r2[mask])
            rv = np.maximum(rv, 0.5)
            eps_ij = np.sqrt(eps_r[mask] * eps_l[j])
            rmin_ij = rmin_r[mask] + rmin_l[j]
            x = rmin_ij / rv
            vdw = (eps_ij * (x**12 - 2.0 * x**6)).sum()
            elec = (COULOMB_CONST * q_r[mask] * q_l[j] / (eps_rr * rv)).sum()
            e += vdw + elec

        return float(e)

    def score_and_grad(self, coords: np.ndarray) -> Tuple[float, np.ndarray]:
        """Return (energy, gradient wrt ligand coords)."""
        pos_r = self._pos_r
        eps_r = self._eps_r_arr
        rmin_r = self._rmin_r_arr
        q_r = self._q_r
        cutoff2 = self.cutoff * self.cutoff
        eps_rr = self.epsilon_r

        if self.ligand is not None:
            eps_l, rmin_l = _vdw_arrays(self.ligand)
            q_l = self.ligand.charges
        else:
            eps_l = np.ones(len(coords)) * 0.1
            rmin_l = np.ones(len(coords)) * 1.9
            q_l = np.zeros(len(coords))

        e = 0.0
        grad = np.zeros_like(coords)

        for j in range(len(coords)):
            diff = pos_r - coords[j]   # (N_rec, 3)
            r2 = (diff * diff).sum(1)
            mask = r2 < cutoff2
            if not mask.any():
                continue
            rv = np.sqrt(r2[mask])
            rv = np.maximum(rv, 0.5)
            eps_ij = np.sqrt(eps_r[mask] * eps_l[j])
            rmin_ij = rmin_r[mask] + rmin_l[j]
            x = rmin_ij / rv

            # Energy
            vdw_e = eps_ij * (x**12 - 2.0 * x**6)
            elec_e = COULOMB_CONST * q_r[mask] * q_l[j] / (eps_rr * rv)
            e += vdw_e.sum() + elec_e.sum()

            # Gradient wrt coords[j]
            dvdw_dr = eps_ij * (-12.0 * x**12 + 12.0 * x**6) / rv
            delec_dr = -COULOMB_CONST * q_r[mask] * q_l[j] / (eps_rr * rv**2)
            # force on j from receptor atoms: -dE/d(r_j) = +sum over r * unit_vec
            unit = diff[mask] / rv[:, None]   # direction from j to r
            grad[j] += ((dvdw_dr + delec_dr)[:, None] * unit).sum(0)

        return e, grad
