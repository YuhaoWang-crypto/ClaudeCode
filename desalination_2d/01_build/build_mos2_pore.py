#!/usr/bin/env python3
"""
Build a single-layer MoS2 (2H) nanoporous membrane for desalination.

MoS2 is the STAGE-2 material. Unlike graphene it is a three-atom-thick
S-Mo-S sandwich, so a "pore" removes atoms from all three sub-planes. The
chemistry of the pore rim (Mo-terminated vs S-terminated) strongly changes
water flux, so we expose that choice.

Runs with a plain ASE install. Output feeds the same DFT -> MLP -> NEMD
pipeline as the graphene builder.

Usage
-----
    python build_mos2_pore.py --nx 6 --ny 6 --pore-radius 3.5 \
        --rim auto --out ../figures/mos2_pore
"""
import argparse
import numpy as np
from ase.build import mx2
from ase.io import write


def build_pore(nx, ny, pore_radius, rim="auto", vacuum=20.0, a=3.18, thickness=3.19):
    """Return a single-layer 2H-MoS2 sheet with a circular pore removed.

    rim : 'auto' removes every atom within pore_radius (mixed rim);
          'Mo' additionally trims S atoms so the rim is Mo-rich;
          'S'  additionally trims Mo so the rim is S-rich.
    """
    sheet = mx2(formula="MoS2", kind="2H", a=a, thickness=thickness,
                size=(nx, ny, 1), vacuum=vacuum / 2.0)
    cx = sheet.cell[0, 0] / 2.0 + sheet.cell[1, 0] / 2.0
    cy = sheet.cell[1, 1] / 2.0

    dx = sheet.positions[:, 0] - cx
    dy = sheet.positions[:, 1] - cy
    r = np.sqrt(dx**2 + dy**2)
    sym = np.array(sheet.get_chemical_symbols())

    remove = r <= pore_radius
    if rim == "Mo":
        # widen the pore for S so the innermost surviving ring is Mo
        remove |= (sym == "S") & (r <= pore_radius + 0.8)
    elif rim == "S":
        remove |= (sym == "Mo") & (r <= pore_radius + 0.8)

    keep = ~remove
    removed = int(remove.sum())
    n_mo_before = int((sym == "Mo").sum())
    sheet = sheet[keep]
    sym2 = np.array(sheet.get_chemical_symbols())

    info = {
        "n_Mo": int((sym2 == "Mo").sum()),
        "n_S": int((sym2 == "S").sum()),
        "n_removed": removed,
        "S_to_Mo_ratio": round(float((sym2 == "S").sum()) /
                               max(1, (sym2 == "Mo").sum()), 3),
        "rim": rim,
        "pore_radius": pore_radius,
        "cell": sheet.cell.lengths().tolist(),
    }
    return sheet, info


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--nx", type=int, default=6)
    ap.add_argument("--ny", type=int, default=6)
    ap.add_argument("--pore-radius", type=float, default=3.5)
    ap.add_argument("--rim", choices=["auto", "Mo", "S"], default="auto")
    ap.add_argument("--vacuum", type=float, default=20.0)
    ap.add_argument("--out", default="mos2_pore")
    args = ap.parse_args()

    sheet, info = build_pore(args.nx, args.ny, args.pore_radius,
                             rim=args.rim, vacuum=args.vacuum)
    write(f"{args.out}.xyz", sheet)
    try:
        write(f"{args.out}.lammps-data", sheet, format="lammps-data",
              specorder=["Mo", "S"])
    except Exception as exc:  # pragma: no cover
        print(f"[warn] LAMMPS data write skipped: {exc}")

    print("Built nanoporous MoS2:")
    for k, v in info.items():
        print(f"  {k:14s}: {v}")
    print(f"Wrote {args.out}.xyz and {args.out}.lammps-data")


if __name__ == "__main__":
    main()
