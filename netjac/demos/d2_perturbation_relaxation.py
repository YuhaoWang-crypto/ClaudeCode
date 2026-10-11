"""
D2 — Perturbation-relaxation diagnostics: push the system, read |Re λ| off the
recovery.

  A  Exact model, stable node (Teusink 2000, BIOMD0000000064): a G6P bolus,
     only 2 noisy observables recorded -> Hankel-DMD recovers the slowest
     eigenvalue of the 16-species Jacobian.
  B  Exact model with a Hopf bifurcation (Wolf 2000 yeast glycolytic
     oscillations, BIOMD0000000206; glucose influx k0 as the "ageing" knob):
     approaching the Hopf the recovery time τ = 1/|Re λ| diverges while the
     return becomes ringing at the intrinsic frequency; past it the steady
     state is unstable and the system sits on a limit cycle. Eigenvalues
     classify stable node / stable focus / oscillatory; the limit cycle is
     confirmed by simulation (eigenvalues alone cannot certify it).
     Single observable = NADH (the autofluorescence actually measured in
     yeast oscillation experiments).
  C  REAL DATA (CGMacros, PhysioNet; 45 adults, healthy / prediabetes / T2D by
     HbA1c): every logged meal is a natural perturbation. Per meal, fit the
     post-peak glucose decay rate k (= the relaxation eigenvalue of the
     glucose-insulin loop, cf. Bergman 1979 minimal model); per person take the
     median. Also the passive critical-slowing indicators (night-time CGM
     lag-1 autocorrelation and SD; Scheffer 2018 PNAS, Olde Rikkert 2016 CCM).
"""
import os
import numpy as np

from netjac import FIG
from netjac.core import KineticModel
from netjac import spectrum as sp, relax as rx, discriminate as di


# ------------------------------------------------------------------- part A
def part_a(noise=0.01, reps=40, seed=0):
    """Bolus on G6P, observe G6P + F16P with `noise` (fraction of the peak
    deviation), repeat `reps` noise realisations; global 2-mode fit."""
    m = KineticModel("BIOMD0000000064", "Teusink2000")
    x = m.steady_state()
    lam_true = np.sort(np.linalg.eigvals(m.reduced_jacobian(x)).real)[::-1]
    obs = ["G6P", "F16P"]
    t, Y = rx.pulse_response(m, x, "G6P", 0.2 * x[m.state_ids.index("G6P")],
                             t_end=3.0, n=601, observe=obs)
    rng = np.random.default_rng(seed)
    sel = t >= 0.2                       # after the fast (|λ|>50) transients
    dmd_clean = rx.hankel_dmd(Y[sel], t[1] - t[0], rank=2, delays=20)[0].real
    est = []
    for _ in range(reps):
        Yn = Y + noise * np.abs(Y).max(0) * rng.standard_normal(Y.shape)
        est.append(rx.fit_modes(t[sel], Yn[sel], [-3.0, -10.0])["rates"][0])
    est = np.array(est)
    Yshow = Y + noise * np.abs(Y).max(0) * rng.standard_normal(Y.shape)
    return {"lam_true": lam_true[0], "lam2_true": lam_true[1],
            "dmd_clean": dmd_clean, "lam_est": np.median(est),
            "iqr": tuple(np.percentile(est, [25, 75])), "t": t, "Y": Yshow,
            "obs": obs, "rel_err": abs(np.median(est) - lam_true[0])
            / abs(lam_true[0]), "noise": noise}


# ------------------------------------------------------------------- part B
def _lead(m, x):
    return sp.classify(m.reduced_jacobian(x))


def find_hopf(m, x_ref, lo, hi, p="k0", tol=1e-4):
    def re(k):
        m.set_params({p: k})
        return _lead(m, m.fixed_point(x_ref))["lambda_lead"].real
    a, b = lo, hi
    assert re(a) < 0 < re(b)
    while b - a > tol * abs(b):
        c = 0.5 * (a + b)
        a, b = (c, b) if re(c) < 0 else (a, c)
    return 0.5 * (a + b)


def part_b(seed=1):
    m = KineticModel("BIOMD0000000206", "Wolf2000")
    k0_def = m.get("k0")
    m.rr.simulate(0, 50, 10)
    x_ref = m.fixed_point(m.x())
    k_hopf = find_hopf(m, x_ref, 0.5 * k0_def, k0_def)
    rng = np.random.default_rng(seed)
    rows = []
    for frac in [0.6, 0.8, 0.85, 0.9, 0.93, 0.95, 0.97, 0.98, 0.99, 1.02, 1.1,
                 1.3]:
        k0 = frac * k_hopf
        m.set_params({"k0": k0})
        x = m.fixed_point(x_ref)
        cl = _lead(m, x)
        row = {"frac": frac, "k0": k0, "lam": cl["lambda_lead"],
               "kind": cl["kind"], "tau": cl["tau_relax"]}
        ev = cl["eigenvalues"]
        cp = ev[ev.imag > 1e-6]
        crit = cp[np.argmax(cp.real)]            # the pair that crosses
        row["crit"] = crit
        row["tau_crit"] = 1 / abs(crit.real)
        if cl["lambda_lead"].real < 0:
            # pulse on glucose, observe NADH only, 2 % noise; DMD to
            # initialise, then the unbiased global mode fit
            t, Y = rx.pulse_response(m, x, "s1", 0.1 * x[0], 1.5, 1500, ["na"])
            Yn = Y + 0.02 * np.abs(Y).max() * rng.standard_normal(Y.shape)
            lam = rx.hankel_dmd(Yn, t[1] - t[0], rank=6, delays=60)
            oc = lam[lam.imag > 1]
            s0, w0 = ((min(oc[np.argmax(oc.real)].real, -0.5),
                       oc[np.argmax(oc.real)].imag) if len(oc) else (-20, 55))
            fit = rx.fit_modes(t, Yn, [-5.0, -10.0], [(s0, w0)])
            sfit, wfit = fit["osc"][0]
            ok = abs(wfit - crit.imag) < 0.2 * crit.imag and sfit > -100
            row.update(osc_fit=complex(sfit, wfit) if ok else None,
                       t=t, y=Yn[:, 0])
        else:
            # beyond Hopf: confirm limit cycle by simulation
            m.set_x(x * 1.01)
            r = m.rr.simulate(0, 60, 6001, selections=["time", "na"])
            tail = r[3000:, 1]
            row.update(amp=tail.max() - tail.min(), t=r[:, 0], y=r[:, 1])
        rows.append(row)
    m.set_params({"k0": k0_def})
    return {"k_hopf": k_hopf, "k0_default": k0_def, "rows": rows}


# ------------------------------------------------------------------- part C
def cgm_metrics(min_carbs=20, carb_band=None):
    import pandas as pd
    from netjac.data import cgmacros
    bio, traces = cgmacros.load()
    bio = bio.copy()
    bio["group"] = pd.cut(bio["A1c PDL (Lab)"], [0, 5.7, 6.5, 99], right=False,
                          labels=["healthy", "prediabetes", "T2D"])
    rows = []
    for pid, df in traces.items():
        df = df.copy()
        df["t"] = pd.to_datetime(df["Timestamp"])
        df = df.set_index("t")
        g5 = df["Dexcom GL"].astype(float).resample("5min").mean() \
            .interpolate(limit=3)
        carbs = df["Carbs"].fillna(0)
        ok = df["Meal Type"].notna() & (carbs >= min_carbs)
        if carb_band:
            ok &= (carbs >= carb_band[0]) & (carbs <= carb_band[1])
        meal_t = df.index[ok].sort_values()
        all_meals = df.index[df["Meal Type"].notna()].sort_values()
        ks, halfs, peaks = [], [], []
        for t0 in meal_t:
            later = all_meals[all_meals > t0]
            if len(later) and later[0] - t0 < pd.Timedelta("3h"):
                continue                       # isolated meals only
            seg = g5[t0 - pd.Timedelta("30min"): t0 + pd.Timedelta("3h")]
            if len(seg) < 36 or seg.isna().mean() > 0.1:
                continue
            base = seg[:t0].mean()
            post = seg[t0:]
            ip = int(np.nanargmax(post.values))
            tpk = (post.index[ip] - t0).total_seconds() / 60
            if tpk > 120 or post.values[ip] - base < 15:
                continue
            dec = post.iloc[ip:].dropna()
            tt = (dec.index - dec.index[0]).total_seconds().values / 60
            if len(tt) < 12:
                continue
            fit = rx.fit_decay(tt, dec.values)
            if fit["r2"] < 0.7:
                continue
            ks.append(fit["k"])
            peaks.append(post.values[ip] - base)
            half = base + (post.values[ip] - base) / 2
            below = np.where(dec.values <= half)[0]
            halfs.append(tt[below[0]] if len(below) else np.nan)
        night = g5[(g5.index.hour >= 1) & (g5.index.hour < 5)]
        ar, sd = [], []
        for _, seg in night.groupby(night.index.date):
            s = seg.dropna().values
            if len(s) < 40:
                continue
            i = np.arange(len(s))
            r = s - np.polyval(np.polyfit(i, s, 1), i)
            ar.append(np.corrcoef(r[:-1], r[1:])[0, 1])
            sd.append(r.std())
        rows.append({"subject": pid, "n_meals": len(ks),
                     "k": np.median(ks) if ks else np.nan,
                     "t_half": (np.nanmedian(halfs) if np.isfinite(halfs).any()
                                else np.nan) if halfs else np.nan,
                     "peak": np.median(peaks) if peaks else np.nan,
                     "ar1_night": np.median(ar) if ar else np.nan,
                     "sd_night": np.median(sd) if sd else np.nan})
    R = pd.DataFrame(rows).merge(
        bio[["subject", "group", "A1c PDL (Lab)", "Fasting GLU - PDL (Lab)",
             "BMI", "Age"]], on="subject")
    return R.rename(columns={"A1c PDL (Lab)": "hba1c",
                             "Fasting GLU - PDL (Lab)": "fasting_glucose"})


def _partial_spearman(x, y, z):
    from scipy.stats import rankdata, pearsonr
    rx_, ry, rz = rankdata(x), rankdata(y), rankdata(z)
    ex = rx_ - np.polyval(np.polyfit(rz, rx_, 1), rz)
    ey = ry - np.polyval(np.polyfit(rz, ry, 1), rz)
    return pearsonr(ex, ey)


def part_c():
    from scipy.stats import spearmanr
    R = cgm_metrics()
    feats = ["k", "t_half", "peak", "ar1_night", "sd_night",
             "fasting_glucose"]
    comps = {}
    for a, b in [("T2D", "healthy"), ("prediabetes", "healthy"),
                 ("T2D", "prediabetes")]:
        A, B = R[R.group == a], R[R.group == b]
        comps[(a, b)] = {f: di.auc(A[f].dropna(), B[f].dropna(), 2000)
                         for f in feats}
    corr = {}
    for f in feats[:-1]:
        s = R[[f, "hba1c", "fasting_glucose"]].dropna()
        corr[f] = (spearmanr(s[f], s["hba1c"]),
                   _partial_spearman(s[f].values, s["hba1c"].values,
                                     s["fasting_glucose"].values))
    med = R.groupby("group", observed=True)[["n_meals", "k", "t_half",
                                             "ar1_night", "sd_night",
                                             "fasting_glucose"]].median()
    Rb = cgm_metrics(carb_band=(30, 80))
    sb = Rb[["k", "hba1c"]].dropna()
    band = spearmanr(sb["k"], sb["hba1c"])
    return {"table": R, "comps": comps, "corr": corr, "median": med,
            "carb_band_corr": band, "n_band": len(sb)}


def report(make_fig=True, with_real=True):
    A = part_a()
    B = part_b()
    C = None
    if with_real:
        try:
            C = part_c()
        except Exception as e:                    # network / data access
            print(f"D2-C real data skipped ({type(e).__name__}: {e})")
    print("=" * 68)
    print("D2  PERTURBATION-RELAXATION")
    print("=" * 68)
    print(f"[A] Teusink 16-species Jacobian: slowest λ1 = {A['lam_true']:.3f}, "
          f"λ2 = {A['lam2_true']:.3f} /min")
    print(f"    G6P bolus, observe only G6P+F16P: noise-free Hankel-DMD λ1 = "
          f"{A['dmd_clean']:.3f}; with {A['noise']*100:.0f} % noise, global "
          f"2-mode fit λ1 = {A['lam_est']:.3f} (IQR {A['iqr'][0]:.2f}.."
          f"{A['iqr'][1]:.2f}; error {A['rel_err']*100:.1f} %)")
    print(f"[B] Wolf 2000: Hopf at k0* = {B['k_hopf']:.3f} (SBML default "
          f"{B['k0_default']:.0f} lies inside the oscillatory regime)")
    for r in B["rows"]:
        lam = r["lam"]
        if "amp" in r:
            extra = f"limit cycle, NADH amplitude {r['amp']:.3f}"
        elif r.get("osc_fit") is not None:
            extra = (f"NADH-pulse fit of ringing: {r['osc_fit'].real:+.2f}"
                     f"{r['osc_fit'].imag:+.2f}i")
        else:
            extra = "ringing too damped to resolve from NADH"
        c = r["crit"]
        print(f"    k0/k0*={r['frac']:.2f}: lead λ={lam.real:+.2f}"
              f"{lam.imag:+.2f}i [{r['kind'][:13]}]; critical pair "
              f"{c.real:+.2f}{c.imag:+.2f}i τc={r['tau_crit']:.3f}; {extra}")
    if C:
        print(f"[C] CGMacros REAL data, n = {len(C['table'])} "
              "(group by HbA1c, ADA cut-offs)")
        med = C["median"]
        for g, row in med.iterrows():
            print(f"    {g:11s} median k = {row['k']:.4f}/min (τ = "
                  f"{1/row['k']:.0f} min), t½ = {row['t_half']:.0f} min, "
                  f"night AR1 = {row['ar1_night']:.3f}, meals/person = "
                  f"{row['n_meals']:.0f}")
        for (a, b), d in C["comps"].items():
            print(f"    {a} vs {b}: " + ", ".join(
                f"{f} {v['auc']:.2f}[{v['ci'][0]:.2f},{v['ci'][1]:.2f}]"
                for f, v in d.items()))
        for f, (s, ps) in C["corr"].items():
            print(f"    {f:9s} Spearman vs HbA1c {s.statistic:+.2f} "
                  f"(p={s.pvalue:.1g}); partial | fasting glucose "
                  f"{ps.statistic:+.2f} (p={ps.pvalue:.1g})")
        cb = C["carb_band_corr"]
        print(f"    sensitivity (meals 30-80 g carbs only, n={C['n_band']}): "
              f"k vs HbA1c Spearman {cb.statistic:+.2f} (p={cb.pvalue:.1g})")
    if make_fig:
        _figure(A, B, C)
    return {"A": A, "B": B, "C": C}


def _figure(A, B, C):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(FIG, exist_ok=True)
    fig, ax = plt.subplots(2, 3, figsize=(16, 9))
    a = ax[0, 0]
    for j, o in enumerate(A["obs"]):
        a.plot(A["t"], A["Y"][:, j], label=o)
    a.set_xlabel("min"); a.set_ylabel("deviation (mM)")
    a.set_title(f"A  Teusink bolus: fitted λ1={A['lam_est']:.2f} vs true "
                f"{A['lam_true']:.2f}")
    a.legend(fontsize=8)
    a = ax[0, 1]
    rows = [r for r in B["rows"] if r["lam"].real < 0]
    a.semilogy([r["frac"] for r in rows], [r["tau_crit"] for r in rows],
               "o-", label="τc = 1/|Re λc| critical pair (exact)")
    fr = [r for r in rows if r.get("osc_fit") is not None]
    a.semilogy([r["frac"] for r in fr],
               [1 / abs(r["osc_fit"].real) for r in fr], "x", ms=9,
               label="from NADH pulse (fit)")
    a.semilogy([r["frac"] for r in rows], [r["tau"] for r in rows], ":",
               color="gray", label="1/|Re λ_lead| (non-critical real mode)")
    a.axvline(1, color="r", ls=":")
    a.set_xlabel("k0 / k0*(Hopf)"); a.set_ylabel("damping time (min)")
    a.set_title("B  critical slowing of the ringing mode"); a.legend(fontsize=7)
    a = ax[0, 2]
    off = 0
    for r in B["rows"]:
        if r["frac"] in (0.8, 0.97, 0.99, 1.1):
            sel = r["t"] <= 1.0
            y = r["y"][sel] - np.median(r["y"][sel])
            a.plot(r["t"][sel], y / np.abs(y).max() + off,
                   label=f"k0/k0*={r['frac']}")
            off -= 2.2
    a.set_yticks([])
    a.set_xlabel("min"); a.set_ylabel("NADH (normalised, offset)")
    a.set_title("B  node → ringing focus → limit cycle"); a.legend(fontsize=8)
    if C:
        T = C["table"]
        cols = {"healthy": "C2", "prediabetes": "C1", "T2D": "C3"}
        for k, (f, lab) in enumerate([("k", "meal recovery rate k (1/min)"),
                                      ("ar1_night", "night-time CGM AR1")]):
            a = ax[1, k]
            for i, g in enumerate(["healthy", "prediabetes", "T2D"]):
                v = T[T.group == g][f].dropna()
                a.scatter(np.full(len(v), i) + np.random.default_rng(i)
                          .uniform(-.12, .12, len(v)), v, color=cols[g], s=18)
                a.plot([i - .25, i + .25], [v.median()] * 2, color="k")
            a.set_xticks(range(3))
            a.set_xticklabels(["healthy", "prediabetes", "T2D"])
            a.set_ylabel(lab)
            a.set_title(f"C  REAL CGMacros: {lab.split(' (')[0]}")
        a = ax[1, 2]
        for g in ["healthy", "prediabetes", "T2D"]:
            s = T[T.group == g]
            a.scatter(s["hba1c"], s["k"], color=cols[g], label=g, s=18)
        a.set_xlabel("HbA1c (%)"); a.set_ylabel("k (1/min)")
        s = C["corr"]["k"][0]
        a.set_title(f"C  k vs HbA1c: Spearman {s.statistic:+.2f} "
                    f"(p={s.pvalue:.0e})")
        a.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "d2_perturbation_relaxation.png"), dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    report()
