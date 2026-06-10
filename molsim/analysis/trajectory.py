"""Simple trajectory container for MD / conformational sampling results."""
from __future__ import annotations
import numpy as np
from typing import List, Optional
from ..core import Molecule
from .structure import rmsd, rmsf, radius_of_gyration


class Trajectory:
    """
    Stores a sequence of coordinate frames for a molecule.

    Parameters
    ----------
    mol : Molecule
        Template molecule (topology is shared, coordinates vary per frame).
    """

    def __init__(self, mol: Molecule):
        self.mol = mol
        self._frames: List[np.ndarray] = []    # list of (N, 3) arrays
        self._energies: List[float] = []

    def add_frame(self, positions: np.ndarray, energy: Optional[float] = None) -> None:
        self._frames.append(positions.copy())
        self._energies.append(energy if energy is not None else float('nan'))

    @property
    def n_frames(self) -> int:
        return len(self._frames)

    @property
    def positions(self) -> np.ndarray:
        """(T, N, 3) array of all frames."""
        return np.array(self._frames)

    @property
    def energies(self) -> np.ndarray:
        return np.array(self._energies)

    def rmsd_series(self, ref_frame: int = 0, atom_names: Optional[List[str]] = None) -> np.ndarray:
        """Per-frame RMSD relative to reference frame."""
        ref = self._frames[ref_frame]
        if atom_names:
            idx = [i for i, a in enumerate(self.mol.atoms) if a.name in atom_names]
            ref = ref[idx]
            return np.array([rmsd(ref, self._frames[t][idx]) for t in range(self.n_frames)])
        return np.array([rmsd(ref, self._frames[t]) for t in range(self.n_frames)])

    def rmsf_per_atom(self, indices: Optional[List[int]] = None) -> np.ndarray:
        return rmsf(self._frames, indices)

    def rmsf_per_residue(self) -> dict:
        """RMSF averaged per residue (CA atoms if available, else all heavy)."""
        ca_idx = [i for i, a in enumerate(self.mol.atoms) if a.name == 'CA']
        if not ca_idx:
            ca_idx = [i for i, a in enumerate(self.mol.atoms)
                      if a.element.upper() != 'H']
        vals = self.rmsf_per_atom(ca_idx)
        res_ids = [self.mol.atoms[i].residue_id for i in ca_idx]
        return dict(zip(res_ids, vals.tolist()))

    def rg_series(self) -> np.ndarray:
        """Radius of gyration per frame."""
        result = []
        orig = self.mol.positions.copy()
        for frame in self._frames:
            self.mol.positions = frame
            result.append(radius_of_gyration(self.mol))
        self.mol.positions = orig
        return np.array(result)

    def get_frame(self, i: int) -> Molecule:
        """Return molecule with coordinates from frame i."""
        m = self.mol.clone()
        m.positions = self._frames[i]
        return m

    def __len__(self) -> int:
        return self.n_frames

    def __repr__(self) -> str:
        return f"Trajectory({self.n_frames} frames, {self.mol.n_atoms} atoms)"
