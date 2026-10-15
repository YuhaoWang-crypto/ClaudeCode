"""Observables, statistical errors and the time-step selection rule for the CG trajectories."""
import glob
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
MD_DIR = os.environ.get("MB_MD_DIR", os.path.join(HERE, "md_runs"))
Z_MIN = 0.35            # MB radius = contact height (nm)


def load(name):
    d = np.load(os.path.join(MD_DIR, name + ".npz"))
    return dict(t_ps=d["t_ps"], z=d["z_mb"].astype(float), rg=d["rg_dna"].astype(float),
                r_end=d["r_end"].astype(float), nc=d["n_contact"].astype(float), meta=json.loads(str(d["meta"])))


def tau_int(x, dt, c=5.0):
    """Integrated autocorrelation time (Sokal automatic window), same units as dt."""
    x = np.asarray(x, float) - np.mean(x)
    n = len(x)
    f = np.fft.rfft(x, 2 * n)
    acf = np.fft.irfft(f * np.conj(f))[:n]
    if acf[0] <= 0:
        return dt * 0.5
    acf /= acf[0]
    tau = 0.5
    for m in range(1, n):
        tau += acf[m]
        if m >= c * tau:
            break
    return max(tau, 0.5) * dt


def replica_observables(run, beta=5.0):
    """Per-replica time averages (arrays are (frames, replicas)); each replica is one independent sample."""
    z = run["z"]
    w = np.exp(-beta * (z - Z_MIN))
    dt = run["t_ps"][1] - run["t_ps"][0]
    return dict(z_mean=z.mean(0), rg_mean=run["rg"].mean(0), p_contact=(z < Z_MIN + 0.2).mean(0),
                w_mean=w.mean(0), nc_mean=run["nc"].mean(0),
                tau_z_ps=np.array([tau_int(z[:, r], dt) for r in range(z.shape[1])]))


def pooled(names, beta=5.0):
    """Mean over all replicas of the given runs, with the across-replica standard error."""
    per = [replica_observables(load(n), beta) for n in names]
    keys = per[0].keys()
    allv = {k: np.concatenate([p[k] for p in per]) for k in keys}
    m = {k: float(v.mean()) for k, v in allv.items()}
    se = {k: float(v.std(ddof=1) / np.sqrt(len(v))) for k, v in allv.items()}
    return m, se, per


def dt_selection(candidates=(1.0, 0.5, 0.25, 0.125), tol=0.05):
    """Fastest dt whose halving changes every observable by <= tol (relative)."""
    table = {}
    for dt in candidates:
        names = sorted(os.path.basename(p)[:-4] for p in glob.glob(os.path.join(MD_DIR, f"dt{dt}_s*.npz")))
        if names:
            table[dt] = pooled(names)
    decision, rows = None, []
    for dt in candidates:
        half = dt / 2
        if dt not in table or half not in table:
            continue
        m1, se1, _ = table[dt]
        m2, se2, _ = table[half]
        rel = {k: (m1[k] - m2[k]) / m2[k] for k in m1}
        zsc = {k: (m1[k] - m2[k]) / np.hypot(se1[k], se2[k]) for k in m1}
        ok = all(abs(v) <= tol for v in rel.values())
        rows.append(dict(dt=dt, half=half, rel=rel, z=zsc, pass_=ok))
        if ok and decision is None:
            decision = dt
    return decision, rows, {dt: (t[0], t[1]) for dt, t in table.items()}


if __name__ == "__main__":
    dec, rows, tab = dt_selection()
    for dt, (m, se) in tab.items():
        print(f"dt={dt}: " + ", ".join(f"{k}={m[k]:.4g}±{se[k]:.2g}" for k in m))
    for r in rows:
        print(f"{r['dt']} vs {r['half']}: pass={r['pass_']} " +
              ", ".join(f"{k}:{100*v:+.1f}% (z={r['z'][k]:+.1f})" for k, v in r["rel"].items()))
    print("selected dt:", dec)
