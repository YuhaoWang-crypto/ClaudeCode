"""Render a 2D slice of the differential charge density Δρ from a cube file,
analogous to the report's Figures 1-3 (electron accumulation vs depletion)."""
import os, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
RDIR = os.path.join(HERE, "..", "results")
FDIR = os.path.join(HERE, "..", "figures")
os.makedirs(FDIR, exist_ok=True)
BOHR = 0.52917721067

def read_cube(fn):
    with open(fn) as f:
        f.readline(); f.readline()
        natm, ox, oy, oz = f.readline().split()[:4]
        natm = int(natm); origin = np.array([float(ox), float(oy), float(oz)])
        n = []; vecs = []
        for _ in range(3):
            parts = f.readline().split()
            n.append(int(parts[0])); vecs.append([float(x) for x in parts[1:4]])
        n = np.array(n); vecs = np.array(vecs)
        atoms = []
        for _ in range(natm):
            p = f.readline().split()
            atoms.append((int(p[0]), np.array([float(p[2]), float(p[3]), float(p[4])])))
        data = np.fromstring(" ".join(f.read().split()), sep=" ").reshape(n)
    return data, n, origin, vecs, atoms

def main():
    fn = os.path.join(RDIR, "demo_delta_rho.cube")
    data, n, origin, vecs, atoms = read_cube(fn)
    # integrate out y -> project onto x-z plane (adsorbate binds along z)
    # take the y-slice through the adsorbate binding oxygen for a clean cut
    zc = [p for Z, p in atoms]
    # sum over y near molecular plane (a few central slices) for signal
    j0 = n[1] // 2
    sl = data[:, j0-3:j0+4, :].sum(axis=1)
    x = origin[0] + np.arange(n[0]) * vecs[0][0]
    z = origin[2] + np.arange(n[2]) * vecs[2][2]
    X, Zg = np.meshgrid(x * BOHR, z * BOHR, indexing="ij")

    vmax = np.percentile(np.abs(sl), 99.5)
    fig, ax = plt.subplots(figsize=(6, 6))
    pcm = ax.pcolormesh(X, Zg, sl, cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                        shading="auto")
    # overlay atoms projected on x-z
    col = {22: "#9aa0a6", 8: "#d93025", 6: "#202124", 1: "#e8eaed"}
    for Z, p in atoms:
        ax.scatter(p[0]*BOHR, p[2]*BOHR, s=120 if Z > 2 else 40,
                   c=col.get(Z, "#888"), edgecolors="k", linewidths=0.6, zorder=3)
    ax.set_xlabel("x (Å)"); ax.set_ylabel("z (Å)")
    ax.set_title("Differential charge density Δρ (HCOOH / Ti₂O₄)\n"
                 "blue = electron accumulation, red = depletion")
    cb = fig.colorbar(pcm, ax=ax, shrink=0.8); cb.set_label("Δρ (a.u., y-projected)")
    ax.set_aspect("equal")
    out = os.path.join(FDIR, "demo_delta_rho_slice.png")
    fig.tight_layout(); fig.savefig(out, dpi=150)
    print("wrote", out)

if __name__ == "__main__":
    main()
