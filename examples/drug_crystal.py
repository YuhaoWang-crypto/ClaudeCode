"""Organic / drug molecular-crystal relaxation & polymorph ranking (UMA `omc`).

Given one or more candidate crystal structures (CIF) of the same molecule --
e.g. different polymorphs of a drug -- this relaxes each cell with the `omc`
(Open Molecular Crystals) head and ranks them by lattice energy per molecule.
The lowest-energy structure is the model's predicted most stable polymorph.

This is the per-structure evaluation step of a crystal-structure-prediction
(CSP) workflow; UMA is the potential behind FastCSP (arXiv:2508.02641).

Usage:
    python examples/drug_crystal.py --cif form_I.cif form_II.cif form_III.cif
    python examples/drug_crystal.py --cif *.cif --relax-cell
"""

from __future__ import annotations

import argparse
import sys

from ase.filters import FrechetCellFilter
from ase.io import read
from ase.optimize import FIRE

sys.path.insert(0, "src")
from fairchem_predict.loader import get_calculator  # noqa: E402


def count_molecules(atoms) -> int:
    """Number of formula units, assuming one connected molecule repeats.
    Falls back to 1 so energies stay comparable if detection is unreliable."""
    try:
        from ase.build.tools import niggli_reduce  # noqa: F401
        from scipy.sparse.csgraph import connected_components

        from ase.neighborlist import NeighborList, natural_cutoffs

        nl = NeighborList(natural_cutoffs(atoms), self_interaction=False, bothways=True)
        nl.update(atoms)
        n = len(atoms)
        import numpy as np

        adj = np.zeros((n, n))
        for i in range(n):
            for j in nl.get_neighbors(i)[0]:
                adj[i, j] = 1
        ncomp, _ = connected_components(adj, directed=False)
        return max(ncomp, 1)
    except Exception:  # noqa: BLE001
        return 1


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cif", nargs="+", required=True, help="candidate crystal CIFs")
    ap.add_argument("--relax-cell", action="store_true",
                    help="relax lattice + atoms (else atoms only)")
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--fmax", type=float, default=0.03)
    args = ap.parse_args()

    calc = get_calculator("omc", device=args.device)

    rows = []
    for path in args.cif:
        atoms = read(path)
        atoms.calc = calc
        target = FrechetCellFilter(atoms) if args.relax_cell else atoms
        FIRE(target, logfile="-").run(fmax=args.fmax, steps=300)
        z = count_molecules(atoms)
        e_total = atoms.get_potential_energy()
        rows.append((path, e_total, z, e_total / z))

    ref = min(r[3] for r in rows)
    print("\n=== Polymorph ranking (per molecule) ===")
    print(f"{'structure':30s} {'E/mol (eV)':>12s} {'ΔE (kJ/mol)':>12s} {'Z':>3s}")
    for path, _e, z, epm in sorted(rows, key=lambda r: r[3]):
        dkj = (epm - ref) * 96.485
        print(f"{path[:30]:30s} {epm:12.4f} {dkj:12.3f} {z:3d}")
    print("lowest ΔE = model's predicted most stable polymorph")


if __name__ == "__main__":
    main()
