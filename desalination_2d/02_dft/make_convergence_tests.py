#!/usr/bin/env python3
"""
Step 1 helper -- generate a grid of QE SCF inputs to converge ecutwfc and the
k-point mesh before trusting any energy/barrier.

    python make_convergence_tests.py --xyz ../figures/graphene_pore_r25.xyz \
        --outdir qe/conv_graphene

Produces one input per (ecutwfc, kpts) pair. Run them all, then plot total
energy vs ecutwfc and vs k-points; pick the smallest values where the energy is
converged to your tolerance (e.g. 1 meV/atom). Only then run production SCF/NEB.
"""
import argparse
import os
import make_qe_input as mq


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--xyz", required=True)
    ap.add_argument("--outdir", default="qe/conv")
    ap.add_argument("--ecutwfc", type=float, nargs="+",
                    default=[40, 50, 60, 70, 80])
    ap.add_argument("--kmeshes", nargs="+",
                    default=["2x2x1", "3x3x1", "4x4x1"])
    ap.add_argument("--functional", default="vdw-df2")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    made = []
    # ecutwfc sweep at a fixed reference k-mesh (the middle one)
    ref_k = args.kmeshes[len(args.kmeshes) // 2]
    kx, ky, kz = (int(v) for v in ref_k.split("x"))
    for ec in args.ecutwfc:
        out = os.path.join(args.outdir, f"scf_ecut{int(ec)}_k{ref_k}.in")
        mq.make(args.xyz, out, ecutwfc=ec, kpts=(kx, ky, kz),
                functional=args.functional)
        made.append(out)

    # k-mesh sweep at the largest ecutwfc (safest)
    ec_ref = max(args.ecutwfc)
    for km in args.kmeshes:
        kx, ky, kz = (int(v) for v in km.split("x"))
        out = os.path.join(args.outdir, f"scf_ecut{int(ec_ref)}_k{km}.in")
        mq.make(args.xyz, out, ecutwfc=ec_ref, kpts=(kx, ky, kz),
                functional=args.functional)
        made.append(out)

    print(f"\nGenerated {len(made)} convergence inputs in {args.outdir}")
    print("Run each (Modal or HPC), then extract energies with:")
    print("  grep '!    total energy' scf_*.out")
    print("and plot E/atom vs ecutwfc and vs k-mesh; pick the converged knee.")


if __name__ == "__main__":
    main()
