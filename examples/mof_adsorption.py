"""MOF water / CO2 adsorption energy with the UMA `odac` head.

Computes the adsorption (binding) energy of a guest molecule in a MOF:

    E_ads = E(MOF + guest) - E(MOF) - E(guest)

A negative value means binding is favourable. This is the core quantity behind
MOF water uptake / CO2 capture screening (the ODAC = Open Direct Air Capture
task the model was trained on).

Usage:
    python examples/mof_adsorption.py --mof path/to/mof.cif --guest H2O
    python examples/mof_adsorption.py --mof path/to/mof.cif --guest CO2 --relax

Provide your own MOF CIF (e.g. from the CoRE MOF / QMOF databases). The guest
is inserted at the largest pore cavity as a simple starting guess; for
production work you should sample several insertion sites.
"""

from __future__ import annotations

import argparse
import sys

import numpy as np
from ase import Atoms
from ase.build import molecule
from ase.io import read
from ase.optimize import FIRE

sys.path.insert(0, "src")
from fairchem_predict.loader import get_calculator  # noqa: E402


def largest_void(mof: Atoms) -> np.ndarray:
    """Crude estimate of an empty pore point: the grid point in the cell
    farthest from any framework atom. Good enough as an insertion seed."""
    cell = mof.cell.array
    n = 12
    frac = np.stack(
        np.meshgrid(*(np.linspace(0, 1, n, endpoint=False),) * 3, indexing="ij"),
        axis=-1,
    ).reshape(-1, 3)
    pts = frac @ cell
    fpos = mof.get_positions()
    # min distance from each grid point to any atom (ignores PBC images: ok as seed)
    d = np.min(np.linalg.norm(pts[:, None, :] - fpos[None, :, :], axis=-1), axis=1)
    return pts[int(np.argmax(d))]


def energy(atoms: Atoms, calc, relax: bool, fmax: float = 0.05, steps: int = 100):
    atoms = atoms.copy()
    atoms.calc = calc
    if relax:
        FIRE(atoms, logfile="-").run(fmax=fmax, steps=steps)
    return atoms.get_potential_energy()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mof", required=True, help="MOF structure file (CIF, ...)")
    ap.add_argument("--guest", default="H2O", help="ASE molecule name, e.g. H2O, CO2")
    ap.add_argument("--relax", action="store_true", help="relax geometries first")
    ap.add_argument("--device", default="cpu")
    args = ap.parse_args()

    calc = get_calculator("odac", device=args.device)

    mof = read(args.mof)
    guest = molecule(args.guest)

    # Build MOF + guest.
    loaded = mof.copy()
    guest_in = guest.copy()
    guest_in.translate(largest_void(mof) - guest_in.get_center_of_mass())
    loaded += guest_in

    e_mof = energy(mof, calc, args.relax)
    e_guest = energy(guest, calc, args.relax)     # gas-phase reference
    e_loaded = energy(loaded, calc, args.relax)
    e_ads = e_loaded - e_mof - e_guest

    print("\n=== MOF adsorption ===")
    print(f"guest              : {args.guest}")
    print(f"E(MOF)             : {e_mof:12.4f} eV")
    print(f"E(guest)           : {e_guest:12.4f} eV")
    print(f"E(MOF+guest)       : {e_loaded:12.4f} eV")
    print(f"E_ads              : {e_ads:12.4f} eV  ({e_ads * 96.485:8.2f} kJ/mol)")
    print("negative => binding favourable" if e_ads < 0 else "positive => unbound")


if __name__ == "__main__":
    main()
