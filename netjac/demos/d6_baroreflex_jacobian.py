"""
D6 — A data-driven physiological Jacobian: the baroreflex loop (heart
period RR <-> systolic pressure SBP) in REAL beat-to-beat recordings,
across healthy ageing (PhysioNet Autonomic Aging, Schumann & Bär 2021).

This is the physiological twin of the user's 2-variable examples: a
negative feedback loop (SBP up -> baroreflex lengthens RR; RR up -> more
diastolic run-off -> pressure falls) that, like Sturis/Goodwin in D5, can
oscillate only because the loop contains DELAYS (sympathetic vascular
response ~ seconds). Its signature is the Mayer wave (~0.1 Hz).

Per subject (resting, ~8-45 min):
  * beat series -> 2 Hz, band-limited (0.02 Hz high-pass);
  * VAR(1) on (RR, SBP): the direct "pairwise" 2x2 Jacobian
        J = logm(A) / dt
    (eigenvalues: can a 2x2 loop alone explain the 0.1 Hz wave?);
  * VAR(p), p by BIC (p <= 20 samples = 10 s): the delay-embedded system,
    eigenvalues -> continuous lambda = ln(mu)/dt; the dominant mode in the
    LF band (0.04-0.15 Hz) = Mayer-wave mode: frequency and damping
    (|Re lambda|, damping ratio zeta);
  * baroreflex sensitivity (BRS, ms/mmHg) by the transfer-function method
    in LF where coherence > 0.5 — the POSITIVE CONTROL (BRS falls with age
    is textbook; if we did not see it the pipeline would be suspect).
Age: Spearman + rank regression adjusted for sex and recording device.
"""
import os

import numpy as np

from netjac import FIG

FS = 2.0
LF = (0.04, 0.15)


def prepare(npz):
    from scipy.interpolate import interp1d
    from scipy.signal import butter, filtfilt
    d = np.load(npz)
    t, rr, sbp, good = d["t"], d["rr"], d["sbp"], d["good"].astype(bool)
    if good.mean() < 0.85 or good.sum() < 300:
        return None
    t, rr, sbp = t[good], rr[good], sbp[good]
    tg = np.arange(t[0], t[-1], 1 / FS)
    if len(tg) < 600:                                  # < 5 min usable
        return None
    R = interp1d(t, rr, kind="cubic")(tg)
    P = interp1d(t, sbp, kind="cubic")(tg)
    b, a = butter(2, 0.02 / (FS / 2), "high")
    Z = np.column_stack([filtfilt(b, a, R), filtfilt(b, a, P)])
    return {"t": tg, "Z": Z, "rr_mean": rr.mean(), "sbp_mean": sbp.mean()}


def fit_var(Z, p):
    """OLS VAR(p); returns coefficient list A_1..A_p, residual cov, BIC."""
    T, k = Z.shape
    Y = Z[p:]
    X = np.hstack([Z[p - i:T - i] for i in range(1, p + 1)])
    B, *_ = np.linalg.lstsq(X, Y, rcond=None)
    E = Y - X @ B
    S = E.T @ E / len(Y)
    bic = np.log(np.linalg.det(S)) + p * k * k * np.log(len(Y)) / len(Y)
    A = [B[(i * k):(i + 1) * k].T for i in range(p)]
    return A, S, bic


def companion(A):
    k, p = A[0].shape[0], len(A)
    C = np.zeros((k * p, k * p))
    C[:k] = np.hstack(A)
    C[k:, :-k] = np.eye(k * (p - 1))
    return C


def brs_tf(Z):
    from scipy.signal import csd, welch, coherence
    f, Pxx = welch(Z[:, 1], FS, nperseg=256)
    _, Pxy = csd(Z[:, 1], Z[:, 0], FS, nperseg=256)
    _, Coh = coherence(Z[:, 1], Z[:, 0], FS, nperseg=256)
    band = (f >= LF[0]) & (f <= LF[1])
    ok = band & (Coh > 0.5)
    if ok.sum() == 0:
        return np.nan, Coh[band].mean()
    return np.mean(np.abs(Pxy[ok]) / Pxx[ok]), Coh[band].mean()


def var_spectrum(A, S, freqs, dt):
    """Parametric cross-spectral matrix of a VAR at the given frequencies."""
    k = A[0].shape[0]
    out = []
    for f in freqs:
        z = np.exp(-2j * np.pi * f * dt)
        H = np.eye(k, dtype=complex)
        for i, Ai in enumerate(A, start=1):
            H = H - Ai * z ** i
        Hi = np.linalg.inv(H)
        out.append(Hi @ S @ Hi.conj().T)
    return np.array(out)


def analyse_subject(npz):
    from scipy.linalg import logm
    pr = prepare(npz)
    if pr is None:
        return None
    Z = pr["Z"]
    Zs = Z / Z.std(0)                       # standardise for the VARs
    dt = 1 / FS
    # (i) the pairwise 2x2 Jacobian from VAR(1)
    A1, _, _ = fit_var(Zs, 1)
    J = np.real(logm(A1[0])) / dt
    lamJ = np.linalg.eigvals(J)
    # (ii) delay-embedded VAR(p)
    fits = [(p,) + fit_var(Zs, p) for p in range(1, 21)]
    p, A, S, _ = min(fits, key=lambda r: r[3])
    mu = np.linalg.eigvals(companion(A))
    lam = np.log(mu.astype(complex)) / dt
    f = np.abs(lam.imag) / (2 * np.pi)
    lf = (f >= LF[0]) & (f <= LF[1]) & (lam.imag > 0)
    # Mayer wave = peak of the SBP spectrum in LF; take the pole nearest it
    fr = np.linspace(LF[0], LF[1], 111)
    Ssbp = var_spectrum(A, S, fr, dt)[:, 1, 1].real
    f_peak = fr[np.argmax(Ssbp)]
    peak_inside = 2 < np.argmax(Ssbp) < len(fr) - 3
    out = {"p": p, "rho": np.abs(mu).max(), "J": J, "lamJ": lamJ,
           "f_peak": f_peak, "peak_inside": bool(peak_inside),
           "J_complex": bool(np.abs(lamJ.imag).max() > 1e-9),
           "rr": pr["rr_mean"], "sbp": pr["sbp_mean"]}
    if lf.any():
        j = np.argmin(np.where(lf, np.abs(f - f_peak), np.inf))
        L = lam[j]
        out.update(lf_freq=abs(L.imag) / (2 * np.pi), lf_decay=-L.real,
                   lf_zeta=-L.real / abs(L))
    else:
        out.update(lf_freq=np.nan, lf_decay=np.nan, lf_zeta=np.nan)
    out["brs"], out["coh"] = brs_tf(Z)
    return out


def _adjusted_rank_assoc(y, age, covs):
    """Spearman-type partial association of y with age given covariates."""
    from scipy.stats import rankdata, pearsonr
    ry, ra = rankdata(y), rankdata(age)
    X = np.column_stack([np.ones(len(y))] + [np.asarray(c, float)
                                             for c in covs])
    ey = ry - X @ np.linalg.lstsq(X, ry, rcond=None)[0]
    ea = ra - X @ np.linalg.lstsq(X, ra, rcond=None)[0]
    return pearsonr(ey, ea)


def report(make_fig=True, ids=None, workers=6):
    import pandas as pd
    from scipy.stats import spearmanr
    from netjac.data import autonomic_aging as aa
    info = aa.subject_info()
    ids = ids or aa.default_selection()
    aa.fetch_cohort(ids, workers=workers)        # cached after first run
    rows = []
    for sid in ids:
        p = os.path.join(aa.DIR, f"{sid}.npz")
        if not os.path.exists(p):
            continue
        try:
            r = analyse_subject(p)
        except Exception:                           # noqa: BLE001
            r = None
        if r is not None:
            r["ID"] = sid
            rows.append(r)
    R = pd.DataFrame(rows).merge(info, on="ID").dropna(subset=["age"])
    R["group"] = pd.cut(R.age, [0, 30, 50, 100], labels=["18-29", "30-49",
                                                          "50+"])
    feats = {"brs": "BRS (ms/mmHg) [positive control]",
             "f_peak": "SBP spectral LF peak (Hz)",
             "lf_freq": "Mayer-mode frequency (Hz)",
             "lf_decay": "Mayer-mode decay |Re λ| (1/s)",
             "lf_zeta": "Mayer-mode damping ratio ζ",
             "rho": "spectral radius of VAR (stability margin)",
             "rr": "mean RR (ms)", "sbp": "mean SBP (mmHg)"}
    assoc = {}
    for c in feats:
        d = R.dropna(subset=[c, "Sex", "Device"])
        s = spearmanr(d[c], d.age)
        pa = _adjusted_rank_assoc(d[c].values, d.age.values,
                                  [d.Sex.values, d.Device.values])
        assoc[c] = (s.statistic, s.pvalue, pa.statistic, pa.pvalue, len(d))
    med = R.groupby("group", observed=True)[list(feats)].median()
    n_by = R.groupby("group", observed=True).size()
    jc = R.J_complex.mean()
    Jmed = np.median(np.stack(R.J.values), axis=0)
    lf_found = R.lf_freq.notna().mean()
    print("=" * 68)
    print("D6  BAROREFLEX JACOBIAN FROM REAL BEAT-TO-BEAT DATA vs AGE")
    print("=" * 68)
    print(f"subjects analysed: {len(R)} (" + ", ".join(
        f"{g}: {n}" for g, n in n_by.items()) + ")")
    print(f"[i] pairwise 2x2 Jacobian (VAR(1), standardised units, 1/s), "
          f"median:\n    {np.round(Jmed, 3).tolist()}")
    print(f"    its eigenvalues are complex (can oscillate) in "
          f"{jc*100:.0f} % of subjects")
    print(f"    SBP LF spectral peak strictly inside 0.04-0.15 Hz in "
          f"{R.peak_inside.mean()*100:.0f} % of subjects")
    print(f"[ii] delay-embedded VAR(p): median p = {R.p.median():.0f} "
          f"({R.p.median()/FS:.1f} s of memory); an LF (Mayer) mode found in "
          f"{lf_found*100:.0f} %; median frequency {R.lf_freq.median():.3f} Hz")
    print("[iii] age associations (Spearman; partial | sex + device):")
    for c, lab in feats.items():
        s, p, pa, ppa, n = assoc[c]
        print(f"    {lab:44s} rho {s:+.2f} (p={p:.1g}); adj {pa:+.2f} "
              f"(p={ppa:.1g}); n={n}")
    print("    group medians:")
    print(med.round(4).to_string())
    if make_fig:
        _figure(R, feats)
    return {"R": R, "assoc": assoc, "median": med, "Jmed": Jmed,
            "J_complex": jc, "lf_found": lf_found}


def _figure(R, feats):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(FIG, exist_ok=True)
    fig, ax = plt.subplots(2, 3, figsize=(16, 9))
    for a, c in zip(ax.ravel()[:5], ["brs", "lf_freq", "lf_decay", "lf_zeta",
                                     "rho"]):
        d = R.dropna(subset=[c])
        a.scatter(d.age + np.random.default_rng(0).uniform(-1.5, 1.5, len(d)),
                  d[c], s=10, alpha=.6, c=d.Device, cmap="coolwarm")
        for g, sub in d.groupby("group", observed=True):
            a.plot(sub.age.median(), sub[c].median(), "ks", ms=8)
        a.set_xlabel("age (group midpoint, jittered)")
        a.set_ylabel(feats[c])
        if c == "brs":
            a.set_yscale("log")
    a = ax[1, 2]
    lam = np.concatenate(R.lamJ.values)
    a.scatter(lam.real, lam.imag, s=6, alpha=.4, label="2x2 Jacobian (VAR(1))")
    d = R.dropna(subset=["lf_freq"])
    a.scatter(-d.lf_decay, 2 * np.pi * d.lf_freq, s=8, alpha=.6, c="C3",
              label="Mayer mode (VAR(p), delays)")
    a.axvline(0, color="k", lw=.5)
    a.set_xlabel("Re λ (1/s)"); a.set_ylabel("Im λ (rad/s)")
    a.set_title("pairwise 2x2 vs delay-embedded eigenvalues")
    a.legend(fontsize=8)
    fig.suptitle("D6  real baroreflex loop (PhysioNet Autonomic Aging): "
                 "colour = recording device", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "d6_baroreflex_jacobian.png"), dpi=120)
    plt.close(fig)


if __name__ == "__main__":
    report()
