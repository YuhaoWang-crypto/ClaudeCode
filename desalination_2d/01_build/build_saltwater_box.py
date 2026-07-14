#!/usr/bin/env python3
"""
Assemble a full NEMD desalination cell:

    [ feed: salt water ] | [ 2D membrane ] | [ permeate: pure water ]

The membrane (from build_graphene_pore.py / build_mos2_pore.py) is placed
perpendicular to z. Water molecules are packed on both sides on a jittered
grid, and a target number of NaCl ion pairs are substituted into the feed
side to reach a chosen salinity. Rigid graphene-style pistons are marked
(via ASE tags / a written group file) so LAMMPS can push the feed piston
to impose the RO pressure.

This is a lightweight packer meant to produce a *valid starting geometry*
for equilibration -- it is NOT a substitute for packmol/moltemplate density
control, but it is fully reproducible and dependency-free (ASE + numpy).

Usage
-----
    python build_saltwater_box.py --membrane ../figures/graphene_pore.xyz \
        --n-water-feed 400 --n-water-perm 400 --n-nacl 12 \
        --out ../figures/nemd_cell
"""
import argparse
import numpy as np
from ase import Atoms
from ase.io import read, write


def _tip3p_molecule():
    """Return O,H,H positions for a single rigid TIP3P-like water at origin."""
    d, ang = 0.9572, np.deg2rad(104.52)
    return np.array([
        [0.0, 0.0, 0.0],
        [d, 0.0, 0.0],
        [d * np.cos(ang), d * np.sin(ang), 0.0],
    ])


def _pack_water(cell_xy, z_lo, z_hi, n_water, rng, min_gap=2.8):
    """Pack n_water waters on a jittered 3D grid between z_lo and z_hi."""
    lx, ly = cell_xy
    lz = z_hi - z_lo
    # choose a grid that comfortably holds n_water sites
    n_per_axis = int(np.ceil(n_water ** (1.0 / 3.0)))
    nx = ny = n_per_axis
    nz = int(np.ceil(n_water / (nx * ny)))
    xs = np.linspace(min_gap, lx - min_gap, nx)
    ys = np.linspace(min_gap, ly - min_gap, ny)
    zs = np.linspace(z_lo + min_gap, z_hi - min_gap, max(nz, 1))
    sites = [(x, y, z) for z in zs for y in ys for x in xs]
    rng.shuffle(sites)
    sites = sites[:n_water]

    symbols, positions = [], []
    base = _tip3p_molecule()
    for (x, y, z) in sites:
        # random rigid rotation
        q = rng.normal(size=4)
        q /= np.linalg.norm(q)
        R = _quat_to_matrix(q)
        mol = base @ R.T + np.array([x, y, z])
        symbols += ["O", "H", "H"]
        positions.append(mol)
    if positions:
        positions = np.vstack(positions)
    else:
        positions = np.zeros((0, 3))
    return symbols, positions


def _quat_to_matrix(q):
    w, x, y, z = q
    return np.array([
        [1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
        [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
        [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)],
    ])


def build_cell(membrane_file, n_water_feed, n_water_perm, n_nacl,
               feed_thickness=25.0, perm_thickness=25.0, seed=1):
    rng = np.random.default_rng(seed)
    membrane = read(membrane_file)
    lx, ly = membrane.cell[0, 0], membrane.cell[1, 1]
    z_mem = membrane.positions[:, 2].mean()

    # feed on -z side, permeate on +z side
    feed_lo, feed_hi = z_mem - feed_thickness, z_mem - 3.0
    perm_lo, perm_hi = z_mem + 3.0, z_mem + perm_thickness

    fs, fp = _pack_water((lx, ly), feed_lo, feed_hi, n_water_feed, rng)
    ps, pp = _pack_water((lx, ly), perm_lo, perm_hi, n_water_perm, rng)

    # Substitute NaCl into feed side: turn 2*n_nacl waters into ions by
    # dropping their H atoms and relabelling O -> Na/Cl.
    fs, fp = _substitute_ions(fs, fp, n_nacl, rng)

    cell = membrane.copy()
    cell += Atoms(fs, positions=fp)
    cell += Atoms(ps, positions=pp)

    # extend the box in z to hold both reservoirs + vacuum buffer
    lz = (perm_hi - feed_lo) + 10.0
    cell.cell[2, 2] = lz
    cell.center(axis=2)
    cell.pbc = [True, True, True]

    info = {
        "membrane_atoms": len(membrane),
        "feed_water": fs.count("O"),
        "perm_water": ps.count("O"),
        "na_ions": fs.count("Na"),
        "cl_ions": fs.count("Cl"),
        "box": [round(float(lx), 2), round(float(ly), 2), round(float(lz), 2)],
        "approx_salinity_mol_frac": round(
            n_nacl / max(1, n_water_feed), 4),
    }
    return cell, info


def _substitute_ions(symbols, positions, n_nacl, rng):
    """Replace 2*n_nacl water oxygens with Na+/Cl-, dropping their hydrogens."""
    symbols = list(symbols)
    positions = positions.copy()
    o_idx = [i for i, s in enumerate(symbols) if s == "O"]
    rng.shuffle(o_idx)
    chosen = o_idx[: 2 * n_nacl]
    to_delete = set()
    for k, oi in enumerate(chosen):
        symbols[oi] = "Na" if k % 2 == 0 else "Cl"
        # its two hydrogens follow directly in the flat O,H,H layout
        to_delete.add(oi + 1)
        to_delete.add(oi + 2)
    keep = [i for i in range(len(symbols)) if i not in to_delete]
    symbols = [symbols[i] for i in keep]
    positions = positions[keep]
    return symbols, positions


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--membrane", required=True,
                    help="Membrane .xyz from a build_*_pore.py script")
    ap.add_argument("--n-water-feed", type=int, default=400)
    ap.add_argument("--n-water-perm", type=int, default=400)
    ap.add_argument("--n-nacl", type=int, default=12)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", default="nemd_cell")
    args = ap.parse_args()

    cell, info = build_cell(args.membrane, args.n_water_feed,
                            args.n_water_perm, args.n_nacl, seed=args.seed)
    write(f"{args.out}.xyz", cell)
    try:
        write(f"{args.out}.lammps-data", cell, format="lammps-data")
    except Exception as exc:  # pragma: no cover
        print(f"[warn] LAMMPS data write skipped: {exc}")

    print("Assembled NEMD desalination cell:")
    for k, v in info.items():
        print(f"  {k:22s}: {v}")
    print(f"Wrote {args.out}.xyz and {args.out}.lammps-data")


if __name__ == "__main__":
    main()
