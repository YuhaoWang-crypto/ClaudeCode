#!/usr/bin/env python3
"""
Milestone 4 helper -- convert a target RO pressure into the per-atom piston
force (LAMMPS 'metal' units: eV/Angstrom) needed in in.desalination.lammps.

Physics
-------
A rigid piston of N atoms pushing on cross-sectional area A produces pressure
    P = F_total / A = N * f_atom / A     =>     f_atom = P * A / N

Unit bookkeeping (metal units): energy in eV, length in A, so force in eV/A.
    1 Pa = 1 N/m^2 = 1 J/m^2/m ... we convert via:
    f[eV/A] = P[Pa] * A[A^2] / N  * (1 J = 6.2415e18 eV) * (1 m = 1e10 A)^-2 ...
We do it cleanly through SI: force_SI[N] = P[Pa]*A[m^2]/N, then N -> eV/A.

    python piston_force.py --data ../figures/nemd_cell.lammps-data \
        --pressure-mpa 100 --n-piston 60
"""
import argparse
import numpy as np

# unit conversions
EV = 1.602176634e-19          # J per eV
ANG = 1.0e-10                 # m per Angstrom
N_PER_EV_PER_ANG = EV / ANG   # 1 eV/A = 1.602e-9 N  => N to eV/A: divide


def force_per_atom(pressure_pa, area_ang2, n_piston):
    """Return per-atom piston force in eV/Angstrom."""
    area_m2 = area_ang2 * ANG**2
    f_total_N = pressure_pa * area_m2         # N
    f_atom_N = f_total_N / n_piston           # N per atom
    f_atom_eV_A = f_atom_N / N_PER_EV_PER_ANG  # -> eV/A
    return f_atom_eV_A


def read_box_area(data_file):
    """Read xlo/xhi/ylo/yhi from a LAMMPS data file -> in-plane area (A^2)."""
    xlo = xhi = ylo = yhi = None
    with open(data_file) as fh:
        for line in fh:
            if "xlo xhi" in line:
                xlo, xhi = map(float, line.split()[:2])
            elif "ylo yhi" in line:
                ylo, yhi = map(float, line.split()[:2])
            if xlo is not None and ylo is not None:
                break
    if None in (xlo, xhi, ylo, yhi):
        raise SystemExit("Could not read box bounds from data file.")
    return (xhi - xlo) * (yhi - ylo)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", help="LAMMPS data file (reads box area)")
    ap.add_argument("--area-ang2", type=float,
                    help="Override in-plane area (A^2) instead of --data")
    ap.add_argument("--pressure-mpa", type=float, default=100.0)
    ap.add_argument("--n-piston", type=int, required=True,
                    help="Number of atoms in the feed piston group")
    args = ap.parse_args()

    if args.area_ang2:
        area = args.area_ang2
    elif args.data:
        area = read_box_area(args.data)
    else:
        raise SystemExit("Provide --data or --area-ang2")

    p_pa = args.pressure_mpa * 1e6
    f = force_per_atom(p_pa, area, args.n_piston)
    print(f"in-plane area      : {area:.2f} A^2")
    print(f"target pressure    : {args.pressure_mpa} MPa")
    print(f"piston atoms (N)   : {args.n_piston}")
    print(f"total force        : {p_pa*area*ANG**2:.3e} N  "
          f"= {p_pa*area*ANG**2/N_PER_EV_PER_ANG:.4f} eV/A")
    print(f"\n  -> set in LAMMPS:  variable fz equal {f:.6f}   # eV/A per piston atom")


if __name__ == "__main__":
    main()
