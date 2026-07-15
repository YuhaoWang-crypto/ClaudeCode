#!/usr/bin/env python3
"""
Steps 3-4 (local demonstration) -- generate a small labelled training set and
prove the DFT -> data -> MLP pipeline end to end.

Real workflow: label AIMD/DFT frames (Modal/HPC) with QE, convert with dpdata.
Here we stand in for that expensive step with the CHGNet universal potential:
we build diverse graphene+water configurations (random displacements + a few
MD-like kicks), label energy+forces with CHGNet, and write extended-XYZ that
MACE / NequIP read directly (and DeePMD via prepare_deepmd_data.py). Swapping
CHGNet for QE labels is a one-line change; the data plumbing is identical.

    python generate_training_frames.py --n-frames 80 \
        --out-train ../data/train/demo_frames.xyz \
        --out-valid ../data/valid/demo_frames.xyz
"""
import argparse
import os
import warnings
warnings.filterwarnings("ignore")
import numpy as np

from ase.build import graphene
from ase import Atoms
from ase.io import write


def _add_waters(atoms, n_water, rng, zgap=3.2):
    """Drop a few rigid waters above the sheet at random in-plane positions."""
    lx, ly = atoms.cell[0, 0], atoms.cell[1, 1]
    z0 = atoms.positions[:, 2].max()
    d, ang = 0.9572, np.deg2rad(104.52)
    base = np.array([[0, 0, 0], [d, 0, 0],
                     [d * np.cos(ang), d * np.sin(ang), 0]])
    for _ in range(n_water):
        shift = np.array([rng.uniform(2, lx - 2), rng.uniform(2, ly - 2),
                          z0 + zgap + rng.uniform(0, 2)])
        q = rng.normal(size=4); q /= np.linalg.norm(q)
        from ase.build import minimize_rotation_and_translation  # noqa
        w, x, y, z = q
        R = np.array([
            [1 - 2*(y*y+z*z), 2*(x*y-z*w), 2*(x*z+y*w)],
            [2*(x*y+z*w), 1 - 2*(x*x+z*z), 2*(y*z-x*w)],
            [2*(x*z-y*w), 2*(y*z+x*w), 1 - 2*(x*x+y*y)]])
        pos = base @ R.T + shift
        atoms += Atoms("OHH", positions=pos)
    return atoms


def generate(n_frames, n_water=2, seed=1):
    from chgnet.model.dynamics import CHGNetCalculator
    calc = CHGNetCalculator()
    rng = np.random.default_rng(seed)
    sheet = graphene(size=(4, 4, 1), vacuum=8.0)

    frames = []
    for i in range(n_frames):
        atoms = sheet.copy()
        atoms = _add_waters(atoms, n_water, rng)
        # thermal-like disorder: random displacements
        amp = rng.uniform(0.03, 0.12)
        atoms.positions += rng.normal(0, amp, size=atoms.positions.shape)
        atoms.calc = calc
        atoms.info["REF_energy"] = float(atoms.get_potential_energy())
        atoms.arrays["REF_forces"] = atoms.get_forces()
        frames.append(atoms)
        if (i + 1) % 20 == 0:
            print(f"  labelled {i+1}/{n_frames} frames", flush=True)
    return frames


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n-frames", type=int, default=80)
    ap.add_argument("--n-water", type=int, default=2)
    ap.add_argument("--valid-fraction", type=float, default=0.2)
    ap.add_argument("--out-train", default="../data/train/demo_frames.xyz")
    ap.add_argument("--out-valid", default="../data/valid/demo_frames.xyz")
    args = ap.parse_args()

    frames = generate(args.n_frames, args.n_water)
    n_val = max(1, int(len(frames) * args.valid_fraction))
    valid, train = frames[:n_val], frames[n_val:]
    for path, fr in [(args.out_train, train), (args.out_valid, valid)]:
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
        write(path, fr, format="extxyz")
    e = np.array([f.info["REF_energy"] for f in frames])
    print(f"\nWrote {len(train)} train + {len(valid)} valid frames "
          f"({len(frames[0])} atoms each)")
    print(f"  energy range: {e.min():.2f} .. {e.max():.2f} eV "
          f"(spread {e.max()-e.min():.2f} eV)")
    print(f"  train: {args.out_train}\n  valid: {args.out_valid}")


if __name__ == "__main__":
    main()
