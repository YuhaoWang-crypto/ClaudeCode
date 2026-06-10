"""Core data structures for molecular simulation."""
from __future__ import annotations
import numpy as np
from dataclasses import dataclass, field
from typing import List, Optional, Dict, Tuple


# fmt: off
ELEMENT_MASS: Dict[str, float] = {
    'H': 1.008,  'C': 12.011, 'N': 14.007, 'O': 15.999, 'S': 32.06,
    'P': 30.974, 'F': 18.998, 'CL': 35.45,  'BR': 79.904, 'I': 126.90,
    'FE': 55.845,'ZN': 65.38, 'MG': 24.305,'CA': 40.078, 'NA': 22.990,
    'K': 39.098, 'MN': 54.938,'CU': 63.546,'CO': 58.933, 'NI': 58.693,
}
ELEMENT_RADIUS: Dict[str, float] = {   # VdW radii in Angstrom
    'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52, 'S': 1.80,
    'P': 1.80, 'F': 1.47, 'CL': 1.75, 'BR': 1.85, 'I': 1.98,
    'FE': 1.80,'ZN': 1.39,'MG': 1.73, 'CA': 1.74, 'NA': 2.27,
    'K': 2.75,
}
# fmt: on


@dataclass
class Atom:
    index: int
    name: str
    element: str
    residue_name: str
    residue_id: int
    chain: str
    position: np.ndarray        # (3,) Angstrom
    charge: float = 0.0
    atom_type: str = ""
    is_hetatm: bool = False

    @property
    def mass(self) -> float:
        return ELEMENT_MASS.get(self.element.upper(), 12.0)

    @property
    def vdw_radius(self) -> float:
        return ELEMENT_RADIUS.get(self.element.upper(), 1.70)

    def __repr__(self) -> str:
        return (f"Atom({self.index} {self.name} {self.residue_name}"
                f"{self.residue_id}:{self.chain})")


@dataclass
class Bond:
    atom1: int          # atom indices
    atom2: int
    order: int = 1      # 1=single, 2=double, 3=triple, 1.5=aromatic


@dataclass
class Molecule:
    """Container for atoms + connectivity."""
    name: str
    atoms: List[Atom] = field(default_factory=list)
    bonds: List[Bond] = field(default_factory=list)

    @property
    def n_atoms(self) -> int:
        return len(self.atoms)

    @property
    def positions(self) -> np.ndarray:
        """Return (N,3) coordinate array."""
        return np.array([a.position for a in self.atoms])

    @positions.setter
    def positions(self, coords: np.ndarray) -> None:
        for i, atom in enumerate(self.atoms):
            atom.position = coords[i].copy()

    @property
    def charges(self) -> np.ndarray:
        return np.array([a.charge for a in self.atoms])

    @property
    def masses(self) -> np.ndarray:
        return np.array([a.mass for a in self.atoms])

    @property
    def center_of_mass(self) -> np.ndarray:
        m = self.masses
        return (self.positions * m[:, None]).sum(0) / m.sum()

    @property
    def elements(self) -> List[str]:
        return [a.element for a in self.atoms]

    def residue_ids(self) -> List[int]:
        seen, result = set(), []
        for a in self.atoms:
            if a.residue_id not in seen:
                seen.add(a.residue_id)
                result.append(a.residue_id)
        return result

    def select(self, **kwargs) -> 'Molecule':
        """Select atoms by field value (e.g. element='C', residue_name='ALA')."""
        atoms = self.atoms
        for k, v in kwargs.items():
            if isinstance(v, (list, tuple, set)):
                atoms = [a for a in atoms if getattr(a, k) in v]
            else:
                atoms = [a for a in atoms if getattr(a, k) == v]
        mol = Molecule(self.name + "_sel")
        mol.atoms = atoms
        return mol

    def get_ca_atoms(self) -> List[Atom]:
        return [a for a in self.atoms if a.name == 'CA']

    def get_backbone_atoms(self) -> List[Atom]:
        return [a for a in self.atoms if a.name in ('N', 'CA', 'C', 'O')]

    def translate(self, vec: np.ndarray) -> None:
        for a in self.atoms:
            a.position += vec

    def rotate(self, R: np.ndarray, center: Optional[np.ndarray] = None) -> None:
        """Rotate by 3x3 rotation matrix around center (default: origin)."""
        pos = self.positions
        if center is not None:
            pos -= center
        pos = pos @ R.T
        if center is not None:
            pos += center
        self.positions = pos

    def clone(self) -> 'Molecule':
        from copy import deepcopy
        return deepcopy(self)

    def __repr__(self) -> str:
        return f"Molecule('{self.name}', {self.n_atoms} atoms)"
