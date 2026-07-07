"""
Module 19 — Migrate the critical-slowing biomarker engine to MANY pathways.

M12/M15/M16/M18 built the dynamical pipeline (bistable switch -> saddle-node ->
critical slowing -> titration biomarker) for MAPK/ERK. This module abstracts it
into a GENERIC engine and applies it to a library of canonical bistable-switch
models spanning diverse pathways. The universal dynamical fact -- every
saddle-node bifurcation produces critical slowing (leading eigenvalue -> 0,
recovery time and fluctuation variance/autocorrelation diverge) -- is
demonstrated concretely, pathway by pathway.

For each 1-D switch dx/dt = f(x; p) with a physiological control parameter p:
  * steady states = roots of f; stability from f'(x) (= the Lyapunov exponent);
  * bistable window in p located by bisection of the saddle-node boundaries;
  * a Langevin titration across the boundary -> residual SD and lag-1
    autocorrelation, which should PEAK at the near-threshold dose.

Model forms are canonical mechanistic switches (positive feedback / covalent
cycle / substrate-depletion), parameterised per pathway; the point is the
SHARED biomarker geometry, not bespoke parameter fits.
"""
import numpy as np
from scipy.optimize import brentq


# ----------------------------- model library -----------------------------
# Each: f(x, p), xmax, control-param range [p_lo, p_hi] chosen to cross a
# saddle-node, param name, pathway, indication, mechanism class.
def _hill(x, K, n):
    xn = x ** n
    return xn / (K ** n + xn)


MODELS = [
    dict(id="MAPK_ERK", pathway="EGFR/RAS/RAF/MEK/ERK",
         indication="RAS-driven cancer", mech="positive feedback",
         pname="drug-induced degradation",
         f=lambda x, p: 0.2 + 3.0 * _hill(x, 1.0, 4) - (1.0 + p) * x,
         xmax=4.0, prange=(0.0, 0.95)),
    dict(id="TF_master", pathway="master-TF fate switch (e.g. PU.1/MyoD)",
         indication="differentiation / stemness", mech="self-activation",
         pname="repressor / inhibitor",
         f=lambda x, p: 0.1 + 2.2 * _hill(x, 0.8, 3) - (1.0 + p) * x,
         xmax=3.0, prange=(0.0, 1.1)),
    dict(id="Rb_E2F", pathway="Rb-E2F restriction point",
         indication="cell-cycle entry / cancer", mech="E2F self-activation",
         pname="growth-factor withdrawal",
         f=lambda x, p: 0.05 + 3.0 * _hill(x, 1.2, 4) - (1.0 + p) * x,
         xmax=4.0, prange=(0.0, 1.2)),
    dict(id="CaMKII", pathway="CaMKII autophosphorylation memory switch",
         indication="synaptic memory / LTP", mech="ultrasensitive pos. fb",
         pname="phosphatase / Ca2+ withdrawal",
         f=lambda x, p: 0.02 + 2.0 * _hill(x, 1.0, 8) - (1.0 + p) * x,
         xmax=3.0, prange=(0.0, 0.95)),
    dict(id="Cdc42_polarity", pathway="Cdc42 GTPase polarity",
         indication="cell polarity / migration", mech="pos. fb (n=2)",
         pname="GAP / inhibitor",
         f=lambda x, p: 0.1 + 2.0 * _hill(x, 1.0, 2) - (1.0 + p) * x,
         xmax=3.0, prange=(0.0, 0.4)),
    dict(id="Apoptosis", pathway="caspase / MOMP apoptotic switch",
         indication="cell death / chemotherapy", mech="irreversible pos. fb",
         pname="death stimulus (rising)",
         f=lambda x, p: p + 2.5 * _hill(x, 1.0, 4) - 1.0 * x,
         xmax=4.0, prange=(0.35, 0.0)),           # rising stimulus -> ON
    dict(id="Cdc2_mitosis", pathway="Cdc2/CyclinB mitotic trigger",
         indication="mitotic entry", mech="double-neg. fb (Cdc25/Wee1)",
         pname="cyclin synthesis (rising)",
         f=lambda x, p: p + 2.2 * _hill(x, 0.7, 3) - 1.0 * x,
         xmax=3.0, prange=(0.5, 0.0)),
    dict(id="Schlogl", pathway="Schlogl chemical switch",
         indication="generic mass-action bistability", mech="substrate deple.",
         pname="reactant supply c",
         f=lambda x, c: -1.0 * x ** 3 + 6.0 * x ** 2 - 11.0 * x + c,
         xmax=4.0, prange=(5.2, 6.8)),
    dict(id="Lac_operon", pathway="lac operon induction",
         indication="metabolic switch", mech="pos. fb (n=2)",
         pname="external inducer withdrawal",
         f=lambda x, p: 0.05 + 2.0 * _hill(x, 1.0, 2) - (1.0 + p) * x,
         xmax=3.0, prange=(0.0, 0.45)),
    dict(id="Wnt_bcat", pathway="Wnt / beta-catenin",
         indication="colorectal cancer / stemness", mech="pos. fb",
         pname="destruction-complex activity",
         f=lambda x, p: 0.15 + 2.5 * _hill(x, 1.0, 3) - (1.0 + p) * x,
         xmax=3.5, prange=(0.0, 0.7)),
]


# ----------------------------- generic engine ----------------------------
def _roots(f, p, xmax, n=3000):
    xs = np.linspace(1e-4, xmax, n)
    v = np.array([f(x, p) for x in xs])
    out = []
    for i in range(len(xs) - 1):
        if v[i] == 0 or v[i] * v[i + 1] < 0:
            xr = brentq(lambda x: f(x, p), xs[i], xs[i + 1])
            h = 1e-4 * max(1.0, xr)
            lam = (f(xr + h, p) - f(xr - h, p)) / (2 * h)   # = f'(x*) = LLE
            out.append({"x": xr, "lam": lam, "stable": lam < 0})
    return out


def n_stable(f, p, xmax):
    return sum(r["stable"] for r in _roots(f, p, xmax))


def bistable_window(f, prange, xmax, n=60):
    ps = np.linspace(prange[0], prange[1], n)
    inside = [p for p in ps if n_stable(f, p, xmax) >= 2]
    if not inside:
        return None
    lo, hi = min(inside), max(inside)
    step = (prange[1] - prange[0]) / (n - 1)

    def edge(a, b):                                # a outside, b inside
        for _ in range(28):
            m = 0.5 * (a + b)
            if n_stable(f, m, xmax) >= 2:
                b = m
            else:
                a = m
        return 0.5 * (a + b)
    return (edge(lo - step, lo), edge(hi + step, hi))


def titration(model, ncells=18, T=1500.0, dt=0.05, sigma=0.06, seed=0):
    """Follow the disappearing stable branch by continuation across the
    upper saddle-node; return residual SD and lag-1 autocorrelation vs p."""
    f, xmax, pr = model["f"], model["xmax"], model["prange"]
    win = bistable_window(f, pr, xmax)
    if win is None:
        return None
    # sweep the disappearing branch from the deep-stable end up to (just
    # inside) its saddle-node; cluster points near the fold where lambda -> 0.
    p_star, p_floor = max(win), min(win)
    p_deep = pr[0] if abs(pr[0] - p_floor) < abs(pr[1] - p_floor) else pr[1]
    u = np.linspace(0, 1, 12)
    frac = 1 - (1 - u) ** 2                          # dense near the fold
    ps = p_deep + frac * 0.99 * (p_star - p_deep)
    rng = np.random.default_rng(seed)
    rows = []
    x_prev = None
    for p in ps:
        st = [r for r in _roots(f, p, xmax) if r["stable"]]
        if not st:
            break
        # follow the branch approaching the saddle-node = the LEAST stable
        # stable root (lambda closest to 0); this is the one that will vanish.
        r = max(st, key=lambda r: r["lam"])
        x_prev = r["x"]
        # Langevin around this branch
        nstep = int(T / dt); sq = np.sqrt(dt); rec = []
        x = r["x"]
        for i in range(nstep):
            x = x + dt * f(x, p) + sigma * sq * rng.normal()
            x = min(max(x, 0.0), xmax)
            if i % 15 == 0:
                rec.append(x)
        rec = np.array(rec)
        rec = rec + rng.normal(0, 0.4 * sigma, size=len(rec))   # meas. noise
        res = rec - rec.mean()
        ar1 = (np.corrcoef(res[:-1], res[1:])[0, 1]
               if res.std() > 0 else np.nan)
        rows.append({"p": p, "x": r["x"], "lam": r["lam"],
                     "tau": (-1.0 / r["lam"] if r["lam"] < 0 else np.inf),
                     "sd": res.std(), "ar1": ar1})
    return rows


def analyze(model):
    f, xmax, pr = model["f"], model["xmax"], model["prange"]
    win = bistable_window(f, pr, xmax)
    tit = titration(model)
    res = {"id": model["id"], "window": win, "bistable": win is not None}
    if tit:
        sds = np.array([r["sd"] for r in tit])
        ars = np.array([r["ar1"] for r in tit])
        taus = [r["tau"] for r in tit if np.isfinite(r["tau"])]
        # critical-slowing signature = rise from the deep-stable end (index 0)
        # toward the saddle-node (last points). Robust to interior-vs-edge.
        near = np.nanmedian(sds[-3:]); deep = np.nanmedian(sds[:3])
        near_a = np.nanmedian(ars[-3:]); deep_a = np.nanmedian(ars[:3])
        res.update({
            "sd_peak_ratio": float(near / (deep + 1e-9)),
            "ar1_peak_ratio": float(near_a / (abs(deep_a) + 1e-9)),
            "max_tau": max(taus) if taus else np.inf,
            "sd_rise": float(near > 1.3 * deep),
            "tit": tit})
    return res


def report(fig_path="figures/m19_switch_library.png"):
    print("=" * 72)
    print("MODULE 19 — Critical-slowing biomarker migrated across pathways")
    print("=" * 72)
    results = [analyze(m) for m in MODELS]
    mmap = {m["id"]: m for m in MODELS}
    print(f"{'pathway/switch':26s} {'bistable':>8} {'window(p)':>16} "
          f"{'maxTau':>8} {'SDpeak':>7} {'AR1pk':>6} {'detect':>7}")
    ndet = 0
    for r in results:
        m = mmap[r["id"]]
        if not r["bistable"]:
            print(f"{r['id']:26s} {'no':>8}")
            continue
        w = r["window"]
        det = r.get("sd_peak_ratio", 0) > 1.3 and r.get("ar1_peak_ratio", 0) > 1.05
        ndet += int(det)
        print(f"{r['id']:26s} {'yes':>8} "
              f"[{min(w):.2f},{max(w):.2f}]".rjust(17)
              + f" {r['max_tau']:8.0f} {r.get('sd_peak_ratio',0):7.1f} "
              f"{r.get('ar1_peak_ratio',0):6.2f} "
              f"{'YES' if det else 'weak':>7}")
    print(f"\n{ndet}/{len(results)} pathway switches show the critical-slowing")
    print("biomarker (variance/autocorr peak at the saddle-node). The SAME")
    print("engine (M15/M16/M18) transfers to every bistable pathway; the")
    print("biomarker geometry is universal to the saddle-node bifurcation.")
    _figure(results, mmap, fig_path)
    print(f"\nfigure written to: {fig_path}")
    return {"results": [{k: v for k, v in r.items() if k != "tit"}
                        for r in results], "n_detected": ndet}


def _figure(results, mmap, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    good = [r for r in results if r.get("tit")]
    ncol = 5
    nrow = int(np.ceil(len(good) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.0 * ncol, 2.7 * nrow))
    axes = np.array(axes).reshape(-1)
    for ax, r in zip(axes, good):
        tit = r["tit"]
        ps = [t["p"] for t in tit]
        sd = [t["sd"] for t in tit]
        ax.plot(ps, sd, "o-", color="#16a085", ms=3)
        axb = ax.twinx()
        axb.plot(ps, [t["ar1"] for t in tit], "s-", color="#8e44ad", ms=3)
        axb.set_yticks([])
        ax.set_title(f"{r['id']}\nSD x{r.get('sd_peak_ratio',0):.1f}",
                     fontsize=8)
        ax.set_xticks([])
        ax.set_yticks([])
    for ax in axes[len(good):]:
        ax.axis("off")
    fig.suptitle("Critical-slowing biomarker across pathway switches — "
                 "residual SD (green) & lag-1 autocorr (purple) peak at each "
                 "saddle-node", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(path, dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    report()
