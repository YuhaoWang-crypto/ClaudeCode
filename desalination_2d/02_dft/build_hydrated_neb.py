#!/usr/bin/env python3
"""
Step 2 (real-DFT side) -- build CI-NEB endpoints for a HYDRATED ion crossing the
pore, the physically-correct desalination barrier that bare-ion NEB misses
(see 06_dft_surrogate_study/REPORT.md section 3.5).

    python build_hydrated_neb.py --membrane ../figures/graphene_pore_r25.xyz \
        --ion Na --n-water 6 --out qe/neb_na_hydrated.in

Endpoints:
  initial : the ion + its full first solvation shell sitting in the feed, a few
            A above the pore
  final   : the ion in the pore plane -- during optimisation the shell partially
            strips (dehydration), which is exactly the barrier we want

The resulting neb.x run gives the DFT dehydration barrier; feed it to
parse_neb_barrier.py. Use a charged supercell (tot_charge is set from the ion).
For a quantitative free energy, follow up with umbrella sampling / metadynamics
in MD using the trained MLP.
"""
import argparse
import numpy as np
from ase.io import read, write
from ase import Atoms

# reuse the hydrated-cluster geometry + the NEB writer already in the repo
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__),
                                "..", "06_dft_surrogate_study"))
import run_hydrated_pmf as hp        # noqa: E402
import build_neb_endpoints as bne    # noqa: E402

PSEUDOS = {
    "C": "C.pbe-n-kjpaw_psl.1.0.0.UPF",
    "H": "H.pbe-kjpaw_psl.1.0.0.UPF",
    "O": "O.pbe-n-kjpaw_psl.1.0.0.UPF",
    "Na": "Na.pbe-spn-kjpaw_psl.1.0.0.UPF",
    "Cl": "Cl.pbe-n-kjpaw_psl.1.0.0.UPF",
    "Mo": "Mo.pbe-spn-kjpaw_psl.1.0.0.UPF",
    "S": "S.pbe-n-kjpaw_psl.1.0.0.UPF",
}


def make_endpoints(membrane_file, ion, n_water, height=4.5):
    membrane = read(membrane_file)
    cx, cy, z0 = bne.pore_centre(membrane) if hasattr(bne, "pore_centre") \
        else _centre(membrane)
    cluster, anchor = hp.build_hydrated_ion(ion, n_water)

    initial = membrane.copy()
    c0 = cluster.copy()
    c0.translate([cx, cy, z0 + height] - c0.positions[anchor])
    initial += c0

    final = membrane.copy()
    cf = cluster.copy()
    cf.translate([cx, cy, z0] - cf.positions[anchor])  # ion in pore plane
    final += cf
    return initial, final


def _centre(membrane):
    cx = membrane.cell[0, 0] / 2 + membrane.cell[1, 0] / 2
    cy = membrane.cell[1, 1] / 2
    z0 = membrane.positions[:, 2].mean()
    return cx, cy, z0


def write_neb(initial, final, out, ion, ecutwfc=60.0, num_images=9,
              kpts=(3, 3, 1)):
    """Reuse the NEB writer from build_neb_endpoints, but with hydrated
    endpoints and the O/H pseudos included."""
    # patch the module's pseudo table to include O for the water shell
    bne.PSEUDOS = PSEUDOS
    bne.CHARGE.setdefault(ion, +1.0 if ion in ("Na", "K", "Li") else -1.0)
    return bne.write_neb_input(initial, final, out, ion,
                               ecutwfc=ecutwfc, num_images=num_images,
                               kpts=kpts)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--membrane", required=True)
    ap.add_argument("--ion", default="Na", choices=["Na", "Cl", "K", "Li"])
    ap.add_argument("--n-water", type=int, default=6)
    ap.add_argument("--height", type=float, default=4.5)
    ap.add_argument("--num-images", type=int, default=9)
    ap.add_argument("--ecutwfc", type=float, default=60.0)
    ap.add_argument("--out", default="qe/neb_hydrated.in")
    args = ap.parse_args()

    initial, final = make_endpoints(args.membrane, args.ion, args.n_water,
                                    args.height)
    write("neb_hyd_initial.xyz", initial)
    write("neb_hyd_final.xyz", final)
    write_neb(initial, final, args.out, args.ion,
              ecutwfc=args.ecutwfc, num_images=args.num_images)
    print(f"Wrote {args.out}")
    print(f"  {args.ion}.(H2O){args.n_water} hydrated NEB, nat={len(initial)}, "
          f"images={args.num_images}")
    print("  endpoints: neb_hyd_initial.xyz (feed) -> neb_hyd_final.xyz (pore)")
    print("  run: neb.x -inp", args.out, "> neb.out ; "
          "parse with parse_neb_barrier.py")


if __name__ == "__main__":
    main()
