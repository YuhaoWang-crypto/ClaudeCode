#!/usr/bin/env python3
"""
Build a nanoporous graphene membrane for RO desalination studies.

Creates a graphene sheet supercell, carves a circular pore of a target
radius, and (optionally) passivates the pore edge with hydrogen. The
output is written both as an extended-XYZ (for inspection / ASE) and as a
LAMMPS data file (for classical / DeePMD MD).

This runs with a plain ASE install -- no DFT/MD engine required. It is the
STAGE-1 geometry generator: the same sheet is later used for
  * QE/CP2K DFT single points and NEB pore-crossing barriers (02_dft)
  * DeePMD training-structure generation (03_mlp)
  * LAMMPS NEMD water flux / salt rejection runs (04_nemd)

Usage
-----
    python build_graphene_pore.py --nx 6 --ny 4 --pore-radius 3.5 \
        --passivate H --out ../figures/graphene_pore

Notes
-----
Pore radius is the physically meaningful knob: sub-nm pores (~0.5-1.0 nm
diameter) are the regime where graphene rejects hydrated Na+/Cl- while
passing water. Start around r = 3.0-4.5 Angstrom and scan.
"""
import argparse
import numpy as np
from ase.build import graphene
from ase import Atom
from ase.io import write


def build_pore(nx, ny, pore_radius, passivate=None, vacuum=20.0):
    """Return an ASE Atoms graphene sheet with a circular pore removed.

    Parameters
    ----------
    nx, ny : int
        Supercell repetitions of the 2-atom graphene primitive cell.
    pore_radius : float
        Radius (Angstrom) of atoms removed around the sheet centre.
    passivate : str or None
        Element symbol to cap dangling edge carbons (e.g. 'H'), or None.
    vacuum : float
        Vacuum padding (Angstrom) added along z on each side.
    """
    sheet = graphene(size=(nx, ny, 1), vacuum=vacuum / 2.0)
    # graphene() centres the sheet in z; get the in-plane centre.
    cx = sheet.cell[0, 0] / 2.0 + sheet.cell[1, 0] / 2.0
    cy = sheet.cell[1, 1] / 2.0
    z0 = sheet.positions[:, 2].mean()

    # Identify atoms inside the pore radius (in-plane distance).
    dx = sheet.positions[:, 0] - cx
    dy = sheet.positions[:, 1] - cy
    r = np.sqrt(dx**2 + dy**2)
    keep = r > pore_radius
    removed = int((~keep).sum())
    sheet = sheet[keep]

    edge_atoms = []
    if passivate:
        edge_atoms = _find_edge_carbons(sheet, cx, cy, pore_radius, z0)
        for pos in edge_atoms:
            sheet.append(Atom(passivate, position=pos))

    info = {
        "n_carbon": int((np.array(sheet.get_chemical_symbols()) == "C").sum()),
        "n_removed": removed,
        "n_passivant": len(edge_atoms),
        "pore_radius": pore_radius,
        "cell": sheet.cell.lengths().tolist(),
    }
    return sheet, info


def _find_edge_carbons(sheet, cx, cy, pore_radius, z0, cc_cut=1.9, nb_full=3):
    """Place passivant atoms pointing radially inward from under-coordinated
    PORE-edge carbons.

    Coordination is counted with the minimum-image convention so that the
    outer carbons sitting on the periodic cell boundary are NOT mistaken for
    dangling bonds -- only carbons lining the pore should be passivated. We
    additionally restrict candidates to a ring just outside the pore radius.
    """
    pos = sheet.positions
    sym = np.array(sheet.get_chemical_symbols())
    c_idx = np.where(sym == "C")[0]
    cpos = pos[c_idx]
    cell = np.array(sheet.cell)
    inv = np.linalg.inv(cell)

    r_from_centre = np.sqrt((cpos[:, 0] - cx) ** 2 + (cpos[:, 1] - cy) ** 2)
    ring = r_from_centre < pore_radius + 1.6  # only carbons lining the pore

    new_positions = []
    for local_i, i in enumerate(c_idx):
        if not ring[local_i]:
            continue
        # minimum-image displacement of every carbon relative to atom i
        d = cpos - cpos[local_i]
        frac = d @ inv.T
        frac -= np.round(frac)
        dmin = frac @ cell
        dist = np.linalg.norm(dmin, axis=1)
        n_neighbours = int(((dist > 0.1) & (dist < cc_cut)).sum())
        if n_neighbours < nb_full:
            # point the passivant back toward the pore centre (into the void)
            rvec = np.array([cx - cpos[local_i, 0], cy - cpos[local_i, 1], 0.0])
            norm = np.linalg.norm(rvec)
            rhat = rvec / norm if norm > 1e-6 else np.array([1.0, 0.0, 0.0])
            new_positions.append(cpos[local_i] + 1.09 * rhat)  # C-H ~1.09 A
    return new_positions


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--nx", type=int, default=6)
    ap.add_argument("--ny", type=int, default=4)
    ap.add_argument("--pore-radius", type=float, default=3.5,
                    help="Pore radius in Angstrom")
    ap.add_argument("--passivate", default="H",
                    help="Edge passivant element, or 'none'")
    ap.add_argument("--vacuum", type=float, default=20.0)
    ap.add_argument("--out", default="graphene_pore",
                    help="Output basename (writes .xyz and .lammps-data)")
    args = ap.parse_args()

    passivate = None if args.passivate.lower() == "none" else args.passivate
    sheet, info = build_pore(args.nx, args.ny, args.pore_radius,
                             passivate=passivate, vacuum=args.vacuum)

    write(f"{args.out}.xyz", sheet)
    try:
        write(f"{args.out}.lammps-data", sheet, format="lammps-data")
    except Exception as exc:  # pragma: no cover - depends on ASE version
        print(f"[warn] LAMMPS data write skipped: {exc}")

    print("Built nanoporous graphene:")
    for k, v in info.items():
        print(f"  {k:12s}: {v}")
    print(f"Wrote {args.out}.xyz and {args.out}.lammps-data")


if __name__ == "__main__":
    main()
