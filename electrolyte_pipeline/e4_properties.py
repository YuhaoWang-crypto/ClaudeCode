"""
E4 — bulk properties that can be compared with experiment   [digest p8]

  density        from the NPT production log (block SEM)              vs. experiment
  D_i            MSD of molecular centres of mass, Einstein relation,
                 linear fit on the log-log-slope≈1 window, block errors
  σ_NE           Nernst–Einstein from D_+ and D_-   (ideal, uncorrelated ions)
  σ_EH           Einstein–Helfand collective charge displacement  (includes
                 ion correlations; ratio σ_EH/σ_NE = "ionicity")
  ε_solvent      Kirkwood–Fröhlich from the fluctuation of the DME dipole sum
                 (ionic contribution EXCLUDED — labelled as such)
  viscosity      from e4_nemd_viscosity (periodic-perturbation NEMD), read
                 here if present, so that everything for one system is in one
                 table with the experimental reference from E0.

Charges used for σ are the *simulation* charges (scaled by 0.8): this is the
conductivity of the model, and it is stated as such.
"""
from __future__ import annotations

import json
import os

import numpy as np

from electrolyte_pipeline.e0_systems import WORK, FIG, SERIES, EXPERIMENT, SERIES_BY_LABEL
from electrolyte_pipeline.traj import Traj

KB = 1.380649e-23
E = 1.602176634e-19


def _unwrapped_com(tr: Traj, species: str, stride=1):
    """Centre-of-mass trajectories (Å) of all molecules of `species`, unwrapped
    with the minimum-image convention frame to frame (box may change under NPT)."""
    idx = tr.sel(species=species)
    mols = tr.molid[idx]
    umol, inv = np.unique(mols, return_inverse=True)
    m = tr.mass[idx]
    M = np.zeros(len(umol)); np.add.at(M, inv, m)
    coms, boxes = [], []
    prev_pos = None
    for ts in tr.frames(stride):
        pos = tr.u.atoms.positions[idx].astype(float)
        box = tr.box(ts)
        if prev_pos is not None:                       # unwrap atoms first
            d = pos - prev_pos
            pos -= box * np.round(d / box)
        prev_pos = pos.copy()
        c = np.zeros((len(umol), 3)); np.add.at(c, inv, pos * m[:, None]); c /= M[:, None]
        coms.append(c); boxes.append(box)
    return np.array(coms), np.array(boxes)          # (frames, nmol, 3), (frames, 3)


def _msd(X, max_lag_frac=0.5, nlags=60):
    """Einstein MSD averaged over time origins (FFT-free, lag subset)."""
    nfr = X.shape[0]
    lags = np.unique(np.linspace(1, int(nfr * max_lag_frac), nlags).astype(int))
    msd = np.array([np.mean(np.sum((X[l:] - X[:-l]) ** 2, axis=-1)) for l in lags])
    return lags, msd


def _fit_D(t_ps, msd_A2, lo=0.3, hi=0.9):
    """D from the slope of MSD(t) in the window [lo,hi] of the lag range, plus
    the log-log slope on the same window (should be ≈1 for diffusive)."""
    i0, i1 = int(len(t_ps) * lo), int(len(t_ps) * hi)
    p = np.polyfit(t_ps[i0:i1], msd_A2[i0:i1], 1)
    ll = np.polyfit(np.log(t_ps[i0:i1]), np.log(msd_A2[i0:i1]), 1)[0]
    D = p[0] / 6 * 1e-20 / 1e-12                   # Å²/ps -> m²/s
    return D, float(ll)


def analyse(label: str, workdir: str = WORK, stride: int = 1, nblocks: int = 4) -> dict:
    import pandas as pd
    tr = Traj(label, workdir)
    T = tr.record["dynamics"]["temperature_K"]
    dt = tr.frame_ps * stride
    res = {"label": label, "T_K": T, "frame_ps": dt}

    # --- density -----------------------------------------------------------
    df = pd.read_csv(os.path.join(tr.dir, "prod.csv")); df.columns = [c.strip('#" ') for c in df.columns]
    rho = df["Density (g/mL)"].values
    bm = np.array([b.mean() for b in np.array_split(rho, 5)])
    res["density"] = {"mean": float(rho.mean()), "sem": float(bm.std(ddof=1) / np.sqrt(5)),
                      "label": "✅ block SEM"}
    V = df["Box Volume (nm^3)"].values.mean() * 1e-27          # m³
    counts = tr.counts()
    if tr.cation in counts:
        res["molarity_from_box"] = counts[tr.cation] / (6.02214076e23 * V * 1e3)

    # --- MSD / D per species -------------------------------------------------
    species = [s for s in counts]
    D, msd_curves = {}, {}
    coms = {}
    if tr.n_frames < 50:
        res["diffusion"] = {}
        res["dielectric_solvent"] = {"eps": float("nan"), "label": "❌ too few frames"}
        res["note"] = f"❌ only {tr.n_frames} frames: transport not evaluated (smoke run)"
        with open(os.path.join(tr.dir, "e4_properties.json"), "w") as fh:
            json.dump(res, fh, indent=1, default=float)
        return res
    for sp in species:
        X, boxes = _unwrapped_com(tr, sp, stride)
        coms[sp] = X
        lags, msd = _msd(X)
        t = lags * dt
        Dv, slope = _fit_D(t, msd)
        # block errors: split trajectory into nblocks and refit
        Db = []
        for b in np.array_split(np.arange(X.shape[0]), nblocks):
            lb, mb = _msd(X[b]); Db.append(_fit_D(lb * dt, mb)[0])
        Db = np.array(Db)
        D[sp] = {"D_m2_s": Dv, "sem": float(Db.std(ddof=1) / np.sqrt(nblocks)), "loglog_slope": slope,
                 "fit_window_ps": [float(t[int(len(t) * 0.3)]), float(t[int(len(t) * 0.9)])],
                 "label": "✅ diffusive" if 0.9 < slope < 1.1 else f"⚠️ log-log slope {slope:.2f} ≠ 1 (sub-diffusive window; extend trajectory)"}
        msd_curves[sp] = (t, msd)
    res["diffusion"] = D

    # --- conductivity: Nernst–Einstein and Einstein–Helfand --------------------
    if tr.cation in counts and "TFSI" in counts:
        qc = float(tr.charge[tr.sel(species=tr.cation)][0])
        qa = float(np.sum(tr.charge[tr.sel(species="TFSI")]) / counts["TFSI"])
        n = counts[tr.cation] / V
        sig_NE = n * E**2 / (KB * T) * (qc**2 * D[tr.cation]["D_m2_s"] + qa**2 * D["TFSI"]["D_m2_s"])
        # collective charge displacement  P(t) = sum_i q_i r_i(t)
        P = qc * coms[tr.cation].sum(axis=1) + qa * coms["TFSI"].sum(axis=1)        # (frames, 3) Å·e
        lags, mP = _msd(P[:, None, :])
        t = lags * dt
        i0, i1 = int(len(t) * 0.3), int(len(t) * 0.9)
        slope = np.polyfit(t[i0:i1], mP[i0:i1], 1)[0]                            # Å² e² / ps
        sig_EH = slope * 1e-20 / 1e-12 * E**2 / (6 * V * KB * T)
        Vb = V
        res["conductivity"] = {"sigma_NE_S_m": sig_NE, "sigma_EH_S_m": sig_EH,
                               "ionicity_EH_over_NE": sig_EH / sig_NE if sig_NE > 0 else None,
                               "charges_used": {"cation": qc, "anion": qa},
                               "label": "⚠️ σ_EH from a single collective observable (no block error possible "
                                        "without independent trajectories); σ_NE inherits D errors"}
        msd_curves["charge"] = (t, mP)

    # --- dielectric (solvent dipole fluctuation, ions excluded) -----------------
    idx = tr.sel(species="DME")
    q = tr.charge[idx]; mols = tr.molid[idx]
    umol, inv = np.unique(mols, return_inverse=True)
    Mtot = []
    for ts in tr.frames(stride * 2):
        pos = tr.u.atoms.positions[idx].astype(float)
        box = tr.box(ts)
        # make each molecule whole (relative to its first atom) before dipole
        first = np.zeros((len(umol), 3)); np.maximum.at(first, inv, 0)  # placeholder
        ref = pos[np.searchsorted(umol, umol)]  # not used
        # simpler: per molecule, unwrap around first atom
        d = np.zeros_like(pos)
        starts = np.r_[0, np.where(np.diff(mols))[0] + 1]
        for s0 in starts:
            s1 = s0 + np.searchsorted(mols[s0:], mols[s0], side="right")
            p = pos[s0:s1] - pos[s0]
            p -= box * np.round(p / box)
            d[s0:s1] = p
        mu = np.zeros((len(umol), 3)); np.add.at(mu, inv, q[:, None] * d)
        Mtot.append(mu.sum(axis=0))
    Mtot = np.array(Mtot) * E * 1e-10                                             # C·m
    fluct = np.mean(np.sum(Mtot**2, axis=1)) - np.sum(np.mean(Mtot, axis=0) ** 2)
    eps0 = 8.8541878128e-12
    eps = 1 + fluct / (3 * eps0 * V * KB * T)
    res["dielectric_solvent"] = {"eps": float(eps),
                                 "label": "⚠️ DME dipole fluctuation only (ions excluded; conducting-system "
                                          "correction not applied); non-polarisable FF, ε∞=1 assumed"}

    # --- viscosity from NEMD if present -----------------------------------------
    nemd = os.path.join(tr.dir, "e4_nemd.json")
    if os.path.exists(nemd):
        res["viscosity"] = json.load(open(nemd))
    res["experiment"] = EXPERIMENT.get(label, {})
    with open(os.path.join(tr.dir, "e4_properties.json"), "w") as fh:
        json.dump(res, fh, indent=1, default=float)
    _plot(res, msd_curves, label)
    return res


def _plot(res, curves, label):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
    for sp, (t, m) in curves.items():
        if sp == "charge":
            continue
        d = res["diffusion"][sp]
        ax[0].loglog(t, m, label=f"{sp}: D={d['D_m2_s']*1e10:.2f}e-10 m²/s, slope {d['loglog_slope']:.2f}")
    ax[0].set_xlabel("t (ps)"); ax[0].set_ylabel("MSD of COM (Å²)"); ax[0].legend(fontsize=7)
    if "charge" in curves:
        t, m = curves["charge"]
        ax[1].plot(t, m, "k-")
        c = res["conductivity"]
        ax[1].set_title(f"σ_NE={c['sigma_NE_S_m']:.3f}  σ_EH={c['sigma_EH_S_m']:.3f} S/m  "
                        f"(ratio {c['ionicity_EH_over_NE']:.2f})", fontsize=9)
    ax[1].set_xlabel("t (ps)"); ax[1].set_ylabel("<|ΣqΔr|²> (Å² e²)")
    fig.suptitle(f"E4  {label}: self-diffusion (Einstein) and collective charge displacement", fontsize=10)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, f"e4_msd_{label}.png"), dpi=140); plt.close(fig)


def report(labels=None, workdir: str = WORK) -> dict:
    labels = labels or [c.label for c in SERIES if os.path.exists(os.path.join(workdir, c.label, "prod.dcd"))]
    out = {}
    print("E4  bulk properties (model charges; see labels in e4_properties.json)")
    print(f"{'system':12s} {'M':>5s} {'rho sim':>8s} {'rho exp':>7s} {'D_cat':>7s} {'D_TFSI':>7s} {'D_DME':>7s} "
          f"{'sNE':>6s} {'sEH':>6s} {'ratio':>5s} {'eps_solv':>8s} {'eta sim':>8s} {'eta exp':>7s}")
    for lab in labels:
        r = analyse(lab, workdir); out[lab] = r
        d = r["diffusion"]; ex = r["experiment"]; c = r.get("conductivity", {})
        cat = [s for s in d if s in ("Li", "Na")]
        f = lambda s: f"{d[s]['D_m2_s']*1e10:7.2f}" if s in d else f"{'-':>7s}"
        eta = r.get("viscosity", {}).get("eta_mPa_s")
        print(f"{lab:12s} {r.get('molarity_from_box', 0):5.2f} {r['density']['mean']:8.4f} "
              f"{ex.get('density_g_cm3', float('nan')):7.4f} {f(cat[0]) if cat else '      -'} {f('TFSI')} {f('DME')} "
              f"{c.get('sigma_NE_S_m', float('nan')):6.3f} {c.get('sigma_EH_S_m', float('nan')):6.3f} "
              f"{(c.get('ionicity_EH_over_NE') or float('nan')):5.2f} {r['dielectric_solvent']['eps']:8.2f} "
              f"{(eta if eta is not None else float('nan')):8.3f} {ex.get('viscosity_mPa_s', float('nan')):7.2f}")
    print("   D in 1e-10 m²/s; σ in S/m; η in mPa·s. Experimental viscosities are eye-read from Chem.Rev. Fig.25.")
    _plot_series(out)
    return out


def _plot_series(out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    labs = [l for l in out if l in SERIES_BY_LABEL]
    if not labs:
        return
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
    M = [out[l].get("molarity_from_box", 0) for l in labs]
    li = [i for i, l in enumerate(labs) if not l.startswith("Na")]
    ax[0].errorbar([M[i] for i in li], [out[labs[i]]["density"]["mean"] for i in li],
                   yerr=[out[labs[i]]["density"]["sem"] for i in li], fmt="o-", label="MD (Sage, q×0.8)")
    ex = [(M[i], EXPERIMENT[labs[i]]["density_g_cm3"]) for i in li if "density_g_cm3" in EXPERIMENT.get(labs[i], {})]
    if ex:
        ax[0].plot(*zip(*ex), "s", mfc="none", color="k", label="experiment")
    ax[0].set_xlabel("M (from box)"); ax[0].set_ylabel("density (g/mL)"); ax[0].legend(fontsize=8)
    for sp, mk in (("Li", "o"), ("TFSI", "^"), ("DME", "s")):
        xs = [M[i] for i in li if sp in out[labs[i]]["diffusion"]]
        ys = [out[labs[i]]["diffusion"][sp]["D_m2_s"] * 1e10 for i in li if sp in out[labs[i]]["diffusion"]]
        ax[1].semilogy(xs, ys, mk + "-", label=sp)
    ax[1].set_xlabel("M"); ax[1].set_ylabel("D (1e-10 m²/s)"); ax[1].legend(fontsize=8)
    ys = [out[labs[i]].get("viscosity", {}).get("eta_mPa_s") for i in li]
    xs = [M[i] for i, y in zip(li, ys) if y is not None]
    ys = [y for y in ys if y is not None]
    if ys:
        ax[2].semilogy(xs, ys, "o-", label="MD NEMD")
    ex = [(M[i], EXPERIMENT[labs[i]]["viscosity_mPa_s"]) for i in li if "viscosity_mPa_s" in EXPERIMENT.get(labs[i], {})]
    if ex:
        ax[2].semilogy(*zip(*ex), "s", mfc="none", color="k", label="exp (Fig.25, eye-read)")
    ax[2].set_xlabel("M"); ax[2].set_ylabel("viscosity (mPa·s)"); ax[2].legend(fontsize=8)
    fig.suptitle("E4  MD vs experiment: density · self-diffusion · viscosity (LiTFSI/DME, 298 K)", fontsize=10)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "e4_series.png"), dpi=140); plt.close(fig)


if __name__ == "__main__":
    import sys
    report(sys.argv[1:] or None)
