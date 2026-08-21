"""Render a binder-target complex as a publication-style figure (no PyMOL needed).

Draws a depth-shaded backbone tube for the target and the binder from CA traces,
highlights the reference epitope, and adds an interface contact map.

Usage:
    python render_complex.py complex.cif --target-chain A --binder-chain X --out fig.png
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from Bio.PDB import MMCIFParser, PDBParser

BG = "#1a1a19"
TARGET = "#8a8a84"
EPITOPE = "#3987e5"
BINDER = "#d95926"
CONTACT = "#c98500"


def load(path: Path):
    p = MMCIFParser(QUIET=True) if path.suffix.lower() in (".cif", ".mmcif") else PDBParser(QUIET=True)
    return p.get_structure("m", str(path))[0]


def ca_trace(chain):
    nums, xyz = [], []
    for r in chain:
        if r.id[0] == " " and "CA" in r:
            nums.append(r.id[1])
            xyz.append(r["CA"].coord)
    return np.array(nums), np.array(xyz)


def smooth(points, n=8):
    """Catmull-Rom-ish subdivision for a smooth tube."""
    if len(points) < 4:
        return points
    out = []
    p = np.vstack([points[0], points, points[-1]])
    for i in range(len(p) - 3):
        p0, p1, p2, p3 = p[i:i + 4]
        for t in np.linspace(0, 1, n, endpoint=False):
            t2, t3 = t * t, t * t * t
            out.append(0.5 * ((2 * p1) + (-p0 + p2) * t +
                              (2 * p0 - 5 * p1 + 4 * p2 - p3) * t2 +
                              (-p0 + 3 * p1 - 3 * p2 + p3) * t3))
    out.append(points[-1])
    return np.array(out)


def principal_view(coords):
    """Rotate so the two largest-variance axes are in-plane; depth along the third."""
    c = coords - coords.mean(0)
    _, _, vt = np.linalg.svd(c, full_matrices=False)
    return vt


def draw_tube(ax, pts, colors, base_lw, depth, zorder):
    segs = np.stack([pts[:-1, :2], pts[1:, :2]], axis=1)
    d = (depth - depth.min()) / max(1e-6, float(depth.max() - depth.min()))
    lw = base_lw * (0.55 + 0.75 * d[:-1])
    rgba = np.array([matplotlib.colors.to_rgba(c) for c in colors[:-1]])
    rgba[:, :3] *= (0.45 + 0.55 * d[:-1])[:, None]
    # dark outline for separation
    ax.add_collection(LineCollection(segs, colors="#05070a", linewidths=lw + 2.6,
                                     capstyle="round", zorder=zorder - 0.1))
    ax.add_collection(LineCollection(segs, colors=rgba, linewidths=lw,
                                     capstyle="round", zorder=zorder))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("complex", type=Path)
    ap.add_argument("--target-chain", default="A")
    ap.add_argument("--binder-chain", default="X")
    ap.add_argument("--offset", type=int, default=0)
    ap.add_argument("--epitope", choices=["siteI", "siteII"], default="siteI")
    ap.add_argument("--title", default="")
    ap.add_argument("--subtitle", default="")
    ap.add_argument("--out", type=Path, default=Path("complex.png"))
    args = ap.parse_args()

    model = load(args.complex)
    tnum, txyz = ca_trace(model[args.target_chain])
    bnum, bxyz = ca_trace(model[args.binder_chain])

    EPITOPES = {
        "siteI": ({30, 33, 54, 61, 66, 69, 73, 74, 75, 78, 172, 175, 178, 179, 180, 182, 183},
                  "site I epitope (IL-6Rα footprint)"),
        "siteII": ({19, 24, 27, 28, 30, 31, 34, 110, 111, 113, 114, 117, 118, 121, 124, 125, 128},
                   "site II epitope (gp130 footprint)"),
    }
    site1, epitope_label = EPITOPES[args.epitope]

    # contacts (CA-CA within 10 A as a display proxy; real contacts come from the analysis script)
    dmat = np.linalg.norm(txyz[:, None] - bxyz[None], axis=-1)
    tgt_contact = dmat.min(1) < 9.0
    bnd_contact = dmat.min(0) < 9.0

    rot = principal_view(np.vstack([txyz, bxyz]))
    center = np.vstack([txyz, bxyz]).mean(0)
    T = (txyz - center) @ rot.T
    B = (bxyz - center) @ rot.T

    fig = plt.figure(figsize=(13, 6.2), facecolor=BG)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.55, 1], wspace=0.18,
                          left=0.03, right=0.97, top=0.86, bottom=0.11)

    ax = fig.add_subplot(gs[0, 0], facecolor=BG)
    ax.set_aspect("equal")
    ax.axis("off")

    tcol = [EPITOPE if (n + args.offset) in site1 else TARGET for n in tnum]
    ts = smooth(T)
    tcol_s = np.repeat(tcol, len(ts) // max(1, len(T)))
    tcol_s = list(tcol_s) + [tcol[-1]] * (len(ts) - len(tcol_s))
    draw_tube(ax, ts, np.array(tcol_s, dtype=object), 4.2, ts[:, 2], zorder=2)

    bs = smooth(B)
    bcol_s = [BINDER] * len(bs)
    draw_tube(ax, bs, np.array(bcol_s, dtype=object), 6.0, bs[:, 2], zorder=3)

    # mark contacting target residues
    ax.scatter(T[tgt_contact, 0], T[tgt_contact, 1], s=18, c=EPITOPE, alpha=0.55,
               edgecolors="none", zorder=1.5)
    ax.scatter(B[bnd_contact, 0], B[bnd_contact, 1], s=18, c=CONTACT, alpha=0.6,
               edgecolors="none", zorder=3.5)

    pad = 6
    allp = np.vstack([ts, bs])
    ax.set_xlim(allp[:, 0].min() - pad, allp[:, 0].max() + pad)
    ax.set_ylim(allp[:, 1].min() - pad, allp[:, 1].max() + pad)

    handles = [
        plt.Line2D([], [], color=TARGET, lw=4, label="IL-6 (target)"),
        plt.Line2D([], [], color=EPITOPE, lw=4, label=epitope_label),
        plt.Line2D([], [], color=BINDER, lw=5, label="designed miniprotein"),
    ]
    leg = ax.legend(handles=handles, loc="lower left", frameon=False, fontsize=9.5,
                    labelcolor="#c3c2b7")
    for t in leg.get_texts():
        t.set_color("#c3c2b7")

    # contact map
    ax2 = fig.add_subplot(gs[0, 1], facecolor="#1a1a19")
    im = ax2.imshow(dmat.T, cmap="magma_r", vmin=3, vmax=25, aspect="auto", origin="lower")
    ax2.set_xlabel("IL-6 residue (1P9M numbering)", color="#c3c2b7", fontsize=9.5)
    ax2.set_ylabel("binder residue", color="#c3c2b7", fontsize=9.5)
    xt = np.linspace(0, len(tnum) - 1, 8).astype(int)
    ax2.set_xticks(xt)
    ax2.set_xticklabels([str(tnum[i] + args.offset) for i in xt])
    ax2.tick_params(colors="#8f8e86", labelsize=8.5)
    for s in ax2.spines.values():
        s.set_color("#30363d")
    for n in sorted(site1):
        idx = np.where(tnum + args.offset == n)[0]
        if len(idx):
            ax2.axvline(idx[0], color=EPITOPE, alpha=0.28, lw=1.0)
    cb = fig.colorbar(im, ax=ax2, fraction=0.046, pad=0.02)
    cb.set_label("Cα–Cα distance (Å)", color="#c3c2b7", fontsize=9)
    cb.ax.tick_params(colors="#8f8e86", labelsize=8)
    cb.outline.set_edgecolor("#30363d")
    ax2.set_title("interface contact map", color="#c3c2b7", fontsize=10.5, pad=8)

    fig.suptitle(args.title, color="#ffffff", fontsize=15, x=0.03, ha="left", y=0.965)
    if args.subtitle:
        fig.text(0.03, 0.905, args.subtitle, color="#8f8e86", fontsize=10, ha="left")

    fig.savefig(args.out, dpi=190, facecolor=BG)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
