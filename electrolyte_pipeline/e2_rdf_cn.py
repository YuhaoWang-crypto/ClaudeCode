"""
E2 — RDF and coordination number: what a peak can and cannot tell   [digest p5]

For each cation X (Li+ or Na+) and partner atom set Y:
  1. g_XY(r) with the partner number density rho_Y of the *whole box*
     (uniform-fluid normalisation, as stated in the digest footnote),
  2. r_min = first minimum after the first peak  -> shell criterion, reported,
  3. N(r) = 4*pi*rho_Y * int_0^r g s^2 ds   (density-weighted integral),
  4. direct-count check: mean number of Y within r_min (must equal N(r_min)),
  5. P(n): distribution of per-ion counts, at atom level AND molecule level
     (a DME can donate two O; "O-atom CN" != "DME-molecule CN"),
  6. block statistics over the trajectory (stability of the mean),
  7. the most populated (n_O_DME, n_O_TFSI) shell as a representative
     configuration written to a PDB.
"""
from __future__ import annotations

import json
import os

import numpy as np

from electrolyte_pipeline.e0_systems import WORK, FIG, SERIES
from electrolyte_pipeline.traj import Traj

PAIRS = [("DME", "O", "O(DME)"), ("TFSI", "O", "O(TFSI)"), ("TFSI", "N", "N(TFSI)"),
         ("TFSI", "F", "F(TFSI)")]


def rdf_and_counts(tr: Traj, cat_idx, part_idx, part_mol, r_max=8.0, dr=0.02, stride=1, nblocks=5):
    from MDAnalysis.lib.distances import distance_array
    edges = np.arange(0, r_max + dr, dr)
    hist = np.zeros(len(edges) - 1)
    vols = []
    nfr = 0
    per_frame_d = []            # list of (n_cat, n_part) distance arrays -> too big; keep counts only
    frame_counts_atom = []      # per frame: per-cation counts within a provisional cut (filled after r_min known)
    dists = []
    for ts in tr.frames(stride):
        box = tr.box(ts)
        d = distance_array(tr.u.atoms.positions[cat_idx], tr.u.atoms.positions[part_idx],
                           box=ts.dimensions)
        hist += np.histogram(d, bins=edges)[0]
        vols.append(np.prod(box))
        dists.append(d.astype(np.float32))
        nfr += 1
    vols = np.array(vols)
    rho = len(part_idx) / vols.mean()                        # partner number density (Å^-3)
    r = 0.5 * (edges[1:] + edges[:-1])
    shell = 4 * np.pi * r**2 * dr
    g = hist / (nfr * len(cat_idx) * rho * shell)
    N_r = np.cumsum(4 * np.pi * rho * g * r**2 * dr)         # eq. on digest p5
    # first peak / first minimum (on a lightly smoothed g so that noise in a
    # short trajectory does not create spurious minima right after the peak)
    ipk = int(np.argmax(g * (r < 4.5)))
    gs = np.convolve(g, np.ones(7) / 7, mode="same")
    w = int(0.25 / dr)                                        # local-min window ±0.25 Å
    imin = None
    for i in range(ipk + int(0.3 / dr), min(len(gs) - w, ipk + int(3.0 / dr))):
        if gs[i] <= gs[i - w:i + w + 1].min() and gs[i] < 0.5 * gs[ipk]:
            imin = i
            break
    if imin is None:                                          # fallback: global min in window
        lo, hi = ipk + int(0.3 / dr), ipk + int(3.0 / dr)
        imin = lo + int(np.argmin(gs[lo:hi]))
    r_min, r_peak = float(r[imin]), float(r[ipk])
    # direct counts per cation per frame within r_min (atom level + molecule level)
    D = np.stack(dists)                                       # (frames, ncat, npart)
    within = D < r_min
    n_atom = within.sum(axis=2)                               # (frames, ncat)
    mol_of = np.asarray(part_mol)
    n_mol = np.zeros_like(n_atom)
    for f in range(D.shape[0]):
        for c in range(D.shape[1]):
            n_mol[f, c] = len(np.unique(mol_of[within[f, c]]))
    blocks = np.array_split(np.arange(D.shape[0]), nblocks)
    bm = np.array([n_atom[b].mean() for b in blocks])
    pn_atom = np.bincount(n_atom.ravel(), minlength=12)[:12] / n_atom.size
    pn_mol = np.bincount(n_mol.ravel(), minlength=12)[:12] / n_mol.size
    return {"r": r, "g": g, "N_r": N_r, "r_peak": r_peak, "r_min": r_min,
            "g_peak": float(g[ipk]), "N_at_rmin": float(N_r[imin]),
            "cn_direct": float(n_atom.mean()), "cn_block_means": bm.round(3).tolist(),
            "cn_sem": float(bm.std(ddof=1) / np.sqrt(nblocks)),
            "cn_mol": float(n_mol.mean()), "P_n_atom": pn_atom.tolist(), "P_n_mol": pn_mol.tolist(),
            "n_frames": int(nfr), "rho_partner_A3": float(rho),
            "_n_atom": n_atom, "_D": D}


def analyse(label: str, workdir: str = WORK, stride: int = 1, write_rep=True) -> dict:
    tr = Traj(label, workdir)
    cat = tr.sel(species=tr.cation)
    res = {"label": label, "cation": tr.cation, "frame_ps": tr.frame_ps, "pairs": {}}
    shells = {}
    for sp, el, name in PAIRS:
        if sp not in tr.counts():
            continue
        part = tr.sel(species=sp, element=el)
        out = rdf_and_counts(tr, cat, part, tr.molid[part], stride=stride)
        shells[name] = out
        res["pairs"][name] = {k: v for k, v in out.items() if not k.startswith("_") and k not in ("r", "g", "N_r")}
        res["pairs"][name]["curve"] = {"r": out["r"].round(3).tolist(), "g": out["g"].round(4).tolist(),
                                       "N_r": out["N_r"].round(4).tolist()}
    # joint shell composition (O_DME, O_TFSI) per cation per frame
    if "O(TFSI)" in shells:
        nd, nt = shells["O(DME)"]["_n_atom"], shells["O(TFSI)"]["_n_atom"]
        combos = {}
        for a, b in zip(nd.ravel(), nt.ravel()):
            combos[(int(a), int(b))] = combos.get((int(a), int(b)), 0) + 1
        top = sorted(combos.items(), key=lambda kv: -kv[1])[:6]
        res["joint_shell_top"] = [{"n_O_DME": k[0], "n_O_TFSI": k[1], "frac": v / nd.size} for k, v in top]
        # total O coordination
        tot = nd + nt
        res["total_O_CN"] = float(tot.mean())
        res["P_total_O"] = (np.bincount(tot.ravel(), minlength=10)[:10] / tot.size).tolist()
        if write_rep:
            _write_representative(tr, cat, shells, top[0][0], label)
    with open(os.path.join(tr.dir, "e2_rdf.json"), "w") as fh:
        json.dump(res, fh)
    _plot(res, label)
    return res


def _write_representative(tr, cat, shells, combo, label):
    """Write the first Li whose shell matches the most populated (n_O_DME, n_O_TFSI)
    together with all DME / TFSI molecules that have an O within r_min."""
    nd, nt = shells["O(DME)"]["_n_atom"], shells["O(TFSI)"]["_n_atom"]
    hits = np.argwhere((nd == combo[0]) & (nt == combo[1]))
    if len(hits) == 0:
        return
    f, c = hits[len(hits) // 2]
    ts = tr.u.trajectory[tr.start + f]
    ion = cat[c]
    mols = set()
    for name in ("O(DME)", "O(TFSI)"):
        sp = "DME" if "DME" in name else "TFSI"
        part = tr.sel(species=sp, element="O")
        d = shells[name]["_D"][f, c]
        mols |= set(tr.molid[part[d < shells[name]["r_min"]]].tolist())
    idx = np.concatenate([[ion], np.where(np.isin(tr.molid, list(mols)))[0]])
    ag = tr.u.atoms[idx]
    # unwrap around the ion
    pos = ag.positions - tr.u.atoms.positions[ion]
    box = tr.box(ts)
    pos -= box * np.round(pos / box)
    ag.positions = pos + tr.u.atoms.positions[ion]
    ag.write(os.path.join(FIG, f"e2_shell_{label}.pdb"))


def _plot(res, label):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
    for name, p in res["pairs"].items():
        r, g, N = map(np.array, (p["curve"]["r"], p["curve"]["g"], p["curve"]["N_r"]))
        ax[0].plot(r, g, label=f"{name}  r_min={p['r_min']:.2f} Å")
        ax[1].plot(r, N, label=f"{name}  N(r_min)={p['N_at_rmin']:.2f}")
        ax[1].axvline(p["r_min"], ls=":", lw=0.8, color="gray")
    ax[0].set_xlabel("r (Å)"); ax[0].set_ylabel(f"g_{res['cation']}-Y(r)"); ax[0].set_xlim(1, 8); ax[0].legend(fontsize=7)
    ax[1].set_xlabel("r (Å)"); ax[1].set_ylabel("N(r) = 4πρ∫g s² ds"); ax[1].set_xlim(1, 8); ax[1].set_ylim(0, 8); ax[1].legend(fontsize=7)
    for name in ("O(DME)", "O(TFSI)"):
        if name in res["pairs"]:
            p = res["pairs"][name]
            ax[2].bar(np.arange(12) - 0.2 + 0.4 * (name == "O(TFSI)"), p["P_n_atom"], width=0.4,
                      label=f"{name} atoms  <n>={p['cn_direct']:.2f}±{p['cn_sem']:.2f}")
    if "O(DME)" in res["pairs"]:
        ax[2].plot(np.arange(12), res["pairs"]["O(DME)"]["P_n_mol"], "k.--",
                   label=f"DME molecules  <n>={res['pairs']['O(DME)']['cn_mol']:.2f}")
    ax[2].set_xlabel("n within r_min"); ax[2].set_ylabel("P(n)"); ax[2].legend(fontsize=7); ax[2].set_xlim(-0.5, 9.5)
    fig.suptitle(f"E2  {label}: RDF (uniform normalisation) · integrated N(r) · CN distribution", fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, f"e2_rdf_{label}.png"), dpi=140)
    plt.close(fig)


def report(labels=None, workdir: str = WORK) -> dict:
    labels = labels or [c.label for c in SERIES if os.path.exists(os.path.join(workdir, c.label, "prod.dcd"))]
    out = {}
    print("E2  RDF / CN   (r_min = first minimum of g(r); CN by density-weighted integral AND direct count)")
    print(f"{'system':12s} {'pair':9s} {'r_peak':>6s} {'r_min':>6s} {'N(r_min)':>8s} {'direct':>7s} "
          f"{'±sem':>5s} {'mol-CN':>6s}  blocks")
    for lab in labels:
        r = analyse(lab, workdir)
        out[lab] = r
        for name, p in r["pairs"].items():
            print(f"{lab:12s} {name:9s} {p['r_peak']:6.2f} {p['r_min']:6.2f} {p['N_at_rmin']:8.2f} "
                  f"{p['cn_direct']:7.2f} {p['cn_sem']:5.2f} {p['cn_mol']:6.2f}  {p['cn_block_means']}")
        if "joint_shell_top" in r:
            top = ", ".join(f"({t['n_O_DME']},{t['n_O_TFSI']}):{t['frac']:.2f}" for t in r["joint_shell_top"][:4])
            print(f"{'':12s} joint (n_O_DME, n_O_TFSI) top: {top};  total O CN = {r['total_O_CN']:.2f}")
    _plot_series(out)
    return out


def _plot_series(out):
    """Li vs Na and concentration trend of P(n) — digest p5 panel C."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    if not out:
        return
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
    for lab, r in out.items():
        if "P_total_O" in r:
            ax[0].plot(np.arange(10), r["P_total_O"], "o-", ms=3, label=f"{lab} <n>={r['total_O_CN']:.2f}")
        if "O(DME)" in r["pairs"]:
            p = r["pairs"]["O(DME)"]
            ax[1].plot(p["curve"]["r"], p["curve"]["g"], label=lab)
    ax[0].set_xlabel("total O (DME+TFSI) within r_min"); ax[0].set_ylabel("P(n)"); ax[0].legend(fontsize=7)
    ax[1].set_xlabel("r (Å)"); ax[1].set_ylabel("g cation–O(DME)"); ax[1].set_xlim(1, 8); ax[1].legend(fontsize=7)
    fig.suptitle("E2  coordination-number distribution and cation–O(DME) RDF across the series", fontsize=10)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "e2_series.png"), dpi=140); plt.close(fig)


if __name__ == "__main__":
    import sys
    report(sys.argv[1:] or None)
