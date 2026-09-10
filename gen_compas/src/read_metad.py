"""Reconstruct the metadynamics free-energy surface from the shared bias files.

Lets the reference FES be inspected while the walkers are still running.
"""

import glob
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common  # noqa: E402

BIAS_FACTOR = 8.0


def load_fes(bias_dir=None):
    bias_dir = bias_dir or os.path.join(common.RESULTS, "metad_bias")
    files = glob.glob(os.path.join(bias_dir, "bias_*.npy"))
    if not files:
        return None, None, None
    # OpenMM writes bias_<walker-id>_<index>.npy; keep the highest index
    # for each walker and sum the walkers' biases.
    latest = {}
    for f in files:
        base = os.path.basename(f)[:-4]
        _, walker, index = base.split("_")
        n = int(index)
        if walker not in latest or n > latest[walker][0]:
            latest[walker] = (n, f)
    total = None
    for walker, (n, f) in sorted(latest.items()):
        b = np.load(f)
        total = b if total is None else total + b
    # well-tempered: F = -(gamma/(gamma-1)) * V_bias
    fes = -BIAS_FACTOR / (BIAS_FACTOR - 1.0) * total
    fes = fes * 0.239005736  # kJ/mol -> kcal/mol
    fes -= fes.min()
    edges = np.linspace(-180.0, 180.0, fes.shape[0] + 1)
    centers = 0.5 * (edges[1:] + edges[:-1])
    # Metadynamics grids are indexed [psi, phi] (last variable varies fastest)
    return fes, centers, sorted(latest.items())


def main():
    fes, centers, files = load_fes()
    if fes is None:
        print("no bias files yet")
        return
    print(f"grid {fes.shape}, walkers/files: {[f for _, (n, f) in files]}")
    print(f"span = {fes.max() - fes.min():.2f} kcal/mol")
    # report the minima
    from scipy import ndimage
    sm = ndimage.gaussian_filter(fes, 1.0, mode="wrap")
    mins = (sm == ndimage.minimum_filter(sm, size=7, mode="wrap")) & (sm < 6.0)
    ij = np.argwhere(mins)
    print("\n minima (phi, psi, F kcal/mol):")
    for i, j in ij:
        print(f"   phi={centers[j]:7.1f}  psi={centers[i]:7.1f}   F={fes[i, j]:5.2f}")


if __name__ == "__main__":
    main()
