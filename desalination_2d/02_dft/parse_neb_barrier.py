#!/usr/bin/env python3
"""
Parse a Quantum ESPRESSO NEB result and report the energy barrier E_a.

neb.x writes a `<prefix>.dat` file with columns:
    reaction_coordinate   energy(eV, rel. to image 1)   error
This reads that file (or `neb.out`) and returns:
    * forward barrier  E_a = max(E) - E(initial)
    * reverse barrier          = max(E) - E(final)
    * reaction energy  dE      = E(final) - E(initial)

    python parse_neb_barrier.py --dat neb_na.dat

The forward barrier is THE selectivity number: compare E_a(Na+) vs E_a(Cl-)
vs E_a(H2O). A good desalination pore has a low water barrier and high ion
barriers.
"""
import argparse
import numpy as np


def parse_dat(path):
    """Read a neb .dat file -> (reaction_coord, energy_eV)."""
    coords, energies = [], []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            try:
                coords.append(float(parts[0]))
                energies.append(float(parts[1]))
            except (ValueError, IndexError):
                continue
    return np.array(coords), np.array(energies)


def barrier(coords, energies):
    e = np.asarray(energies)
    imax = int(np.argmax(e))
    e_forward = float(e[imax] - e[0])
    e_reverse = float(e[imax] - e[-1])
    de = float(e[-1] - e[0])
    return {
        "forward_barrier_eV": round(e_forward, 4),
        "reverse_barrier_eV": round(e_reverse, 4),
        "reaction_energy_eV": round(de, 4),
        "saddle_image": imax + 1,
        "n_images": len(e),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dat", required=True, help="neb.x <prefix>.dat file")
    args = ap.parse_args()
    coords, energies = parse_dat(args.dat)
    if len(energies) < 2:
        raise SystemExit(f"Could not parse >=2 images from {args.dat}")
    res = barrier(coords, energies)
    print(f"NEB barrier from {args.dat}:")
    for k, v in res.items():
        print(f"  {k:20s}: {v}")
    print(f"\n  -> forward barrier E_a = {res['forward_barrier_eV']} eV "
          f"({res['forward_barrier_eV']*96.485:.1f} kJ/mol)")


if __name__ == "__main__":
    main()
