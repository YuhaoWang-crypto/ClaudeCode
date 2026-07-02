"""Inorganic material relaxation + energy with the UMA `omat` head (OMat24).

Relaxes a bulk crystal to its local minimum and reports total / per-atom energy.
This is the OMat24 use case: fast DFT-quality energies & forces for inorganic
materials screening, stability, and geometry optimisation.

Usage:
    python examples/material_relax.py --cif structure.cif --relax-cell
"""

from __future__ import annotations

import argparse
import sys

from ase.filters import FrechetCellFilter
from ase.io import read
from ase.optimize import FIRE

sys.path.insert(0, "src")
from fairchem_predict.loader import get_calculator  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cif", required=True, help="crystal structure file")
    ap.add_argument("--relax-cell", action="store_true")
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--fmax", type=float, default=0.02)
    # OMAT24 EquiformerV2 checkpoints also exist; UMA `omat` is the current default.
    ap.add_argument("--model", default="uma-s-1p1")
    args = ap.parse_args()

    calc = get_calculator("omat", model=args.model, device=args.device)
    atoms = read(args.cif)
    atoms.calc = calc

    e0 = atoms.get_potential_energy()
    target = FrechetCellFilter(atoms) if args.relax_cell else atoms
    FIRE(target, logfile="-").run(fmax=args.fmax, steps=300)
    e1 = atoms.get_potential_energy()

    print("\n=== Inorganic material relaxation ===")
    print(f"formula        : {atoms.get_chemical_formula()}")
    print(f"E initial      : {e0:12.4f} eV  ({e0/len(atoms):8.4f} eV/atom)")
    print(f"E relaxed      : {e1:12.4f} eV  ({e1/len(atoms):8.4f} eV/atom)")
    print(f"E gained       : {e1 - e0:12.4f} eV")


if __name__ == "__main__":
    main()
