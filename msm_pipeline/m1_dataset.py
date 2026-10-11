"""
Module 1 - fetch the dataset and report what is actually in it.

No modelling here; this establishes the provenance and the sampling budget
that every later claim has to live within.
"""
import numpy as np

from .common import (APO_STARTS, DT_NS, LIGAND_NAMES, TRAJ_NS, VARIANTS,
                     apo_indices, banner, holo_indices, load)


def main():
    banner("M1  dataset: TEM beta-lactamase, 48 x 500 ns (Zenodo 6565237)")
    X, names = load()
    n_traj, n_fr = X.shape[0], X.shape[1]
    print(f"  trajectories      {n_traj}")
    print(f"  frames each       {n_fr}  ({DT_NS:.4f} ns apart, {TRAJ_NS:.0f} ns per run)")
    print(f"  C-alpha atoms     {X.shape[2]}  (Ambler 26-290, 239/253 absent)")
    print(f"  aggregate         {n_traj * TRAJ_NS / 1000:.0f} us")
    apo, holo = apo_indices(names), holo_indices(names)
    print(f"  apo / holo        {len(apo)} / {len(holo)}")
    print(f"  apo design        {len(VARIANTS)} variants x {len(APO_STARTS)} crystal starts"
          f"  -> paired by starting structure")
    print(f"  holo design       {len(VARIANTS)} variants x 4 antibiotics "
          f"({', '.join(LIGAND_NAMES.values())})")

    print("\n  sampling budget per variant (apo):")
    for v in VARIANTS:
        ks = apo_indices(names, v)
        print(f"    {v:7s} {len(ks)} runs  {len(ks)*TRAJ_NS/1000:.0f} us")

    print("\n  NOTE  these are 32 *independent* 500 ns runs, not one long run and")
    print("        not adaptively reseeded.  M3 shows what that costs an MSM.")
    return {"n_traj": n_traj, "aggregate_us": n_traj * TRAJ_NS / 1000}


if __name__ == "__main__":
    main()
