"""Does the generative model actually earn its place?

Takes the conformations a finished Gen-COMPAS run collected, and compares the
intermediates proposed by the trained diffusion model (DDIM inversion +
latent interpolation) against the obvious baseline: straight-line
interpolation between the same pairs of end-state structures in the same
feature space.

Reports, for each, the fraction of proposals that are physically realisable
and where they land in (phi, psi), and draws the comparison.
"""

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common  # noqa: E402
import mdops  # noqa: E402
import models  # noqa: E402
from analysis import load_store  # noqa: E402

plt.rcParams.update({"figure.dpi": 150, "savefig.dpi": 150, "font.size": 8})


def phi_of(gen, template):
    full = common.defeaturize_into(np.tile(template, (len(gen), 1, 1)), gen)
    return common.phi_psi(full)


def main():
    tag = sys.argv[1] if len(sys.argv) > 1 else "main"
    n_pairs = 60
    lambdas = np.linspace(0.15, 0.85, 7)

    store = load_store(tag)
    rng = np.random.default_rng(0)
    phi_all, _ = common.phi_psi(store["coords"])
    core = common.which_basin(phi_all)
    a = store["coords"][core == 0]
    b = store["coords"][core == 1]
    print(f"[compare] state A frames {len(a)}, state B frames {len(b)}")

    feats = common.featurize(store["coords"])
    dm = models.DiffusionModel(feats.shape[1], n_steps=400, hidden=384, seed=0)
    dm.fit(feats[rng.choice(len(feats), min(6000, len(feats)), replace=False)],
           steps=8000)

    fa = common.featurize(a[rng.choice(len(a), n_pairs)])
    fb = common.featurize(b[rng.choice(len(b), n_pairs)])

    gen_d = dm.interpolate(fa, fb, lambdas, n_infer=40)
    gen_l = np.concatenate([(1 - lam) * fa + lam * fb for lam in lambdas])

    out = {}
    for name, g in (("diffusion", gen_d), ("linear", gen_l)):
        ok = np.array([mdops.target_is_physical(x.reshape(-1, 3)) for x in g])
        ph, ps = phi_of(g, store["coords"][0])
        out[name] = (g, ok, ph, ps)
        near = np.abs(((ph + 180) % 360) - 180) < 40  # near the phi ~ 0 saddle
        print(f"[compare] {name:9s}: {100 * ok.mean():5.1f}% physical, "
              f"{100 * (ok & near).mean():5.1f}% physical *and* in the "
              f"barrier region (|phi| < 40)")

    _, psi_all = common.phi_psi(store["coords"])
    fig, axes = plt.subplots(1, 2, figsize=(7.6, 3.4), sharex=True,
                             sharey=True)
    for ax, name in zip(axes, ("diffusion", "linear")):
        g, ok, ph, ps = out[name]
        ax.scatter(phi_all[::5], psi_all[::5], s=1.5, c="#d8dee9",
                   label="sampled conformations")
        ax.scatter(ph[~ok], ps[~ok], s=9, c="#c53030", marker="x", lw=0.7,
                   label="proposal, not realisable")
        ax.scatter(ph[ok], ps[ok], s=9, c="#2b6cb0",
                   label="proposal, physical")
        ax.set_xlim(-180, 180)
        ax.set_ylim(-180, 180)
        ax.set_aspect("equal")
        ax.set_xlabel(r"$\phi$ (deg)")
        ax.set_title(f"{name} interpolation\n"
                     f"{100 * ok.mean():.0f}% physically realisable")
    axes[0].set_ylabel(r"$\psi$ (deg)")
    axes[0].legend(fontsize=6, loc="lower left", framealpha=0.9)
    fig.tight_layout()
    p = os.path.join(common.FIGURES, f"generators_{tag}.png")
    fig.savefig(p, bbox_inches="tight")
    print("wrote", p)


if __name__ == "__main__":
    main()
