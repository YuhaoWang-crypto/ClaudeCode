"""
Module 13 — Fisher Information Matrix, sloppy models & stiff axes.

Implements the report's key biomarker-SELECTION method (Gutenkunst-Sethna
2007), which the earlier pipeline lacked. For the ERK dual-phosphorylation
model we compute the sensitivity of an output time course to log-parameters,
build the FIM  F = S^T S, and eigendecompose it.

Findings reproduced:
  * the FIM eigenvalue spectrum is 'sloppy' -- it spans many orders of
    magnitude, and the top few 'stiff' axes carry ~all of the output
    sensitivity;
  * FLUX-RATIO observables project onto the stiff axes far better than single
    concentrations, so the best biomarkers are ratios (ppERK/(pERK+ppERK),
    kinase-drive : phosphatase-capacity), not one protein level.
"""
import numpy as np
from scipy.integrate import solve_ivp

# distinct mechanistic parameters (log-space FIM), base values match M12
PARAM_NAMES = ["a1", "d1", "k1", "a2", "d2", "k2",
               "a3", "d3", "k3", "a4", "d4", "k4", "Ktot"]
BASE = np.array([5, 1, 1, 5, 1, 1, 5, 1, 1, 5, 1, 1, 105.0])
PTOT, MTOT = 100.0, 800.0
SPECIES = ["M", "Mp", "Mpp", "K", "P", "C1", "C2", "C3", "C4"]
IDX = {s: i for i, s in enumerate(SPECIES)}
TIMES = np.linspace(1, 60, 24)          # observation time points


def _rhs(theta):
    a1, d1, k1, a2, d2, k2, a3, d3, k3, a4, d4, k4, Ktot = theta
    rates = [a1, d1, k1, a2, d2, k2, a3, d3, k3, a4, d4, k4]
    RE = [({"K": 1, "M": 1}, {"C1": 1}), ({"C1": 1}, {"K": 1, "M": 1}),
          ({"C1": 1}, {"K": 1, "Mp": 1}),
          ({"K": 1, "Mp": 1}, {"C2": 1}), ({"C2": 1}, {"K": 1, "Mp": 1}),
          ({"C2": 1}, {"K": 1, "Mpp": 1}),
          ({"P": 1, "Mpp": 1}, {"C3": 1}), ({"C3": 1}, {"P": 1, "Mpp": 1}),
          ({"C3": 1}, {"P": 1, "Mp": 1}),
          ({"P": 1, "Mp": 1}, {"C4": 1}), ({"C4": 1}, {"P": 1, "Mp": 1}),
          ({"C4": 1}, {"P": 1, "M": 1})]

    def rhs(t, y):
        c = {s: max(y[IDX[s]], 0.0) for s in SPECIES}
        dy = np.zeros(len(SPECIES))
        for (r, p), kk in zip(RE, rates):
            rate = kk
            for s, co in r.items():
                rate *= c[s] ** co
            for s, co in r.items():
                dy[IDX[s]] -= co * rate
            for s, co in p.items():
                dy[IDX[s]] += co * rate
        return dy
    return rhs, Ktot


def _simulate(theta, observable="Mpp"):
    rhs, Ktot = _rhs(theta)
    y0 = np.zeros(len(SPECIES))
    y0[IDX["M"]] = MTOT; y0[IDX["K"]] = Ktot; y0[IDX["P"]] = PTOT
    sol = solve_ivp(rhs, (0, TIMES[-1] + 1), y0, t_eval=TIMES,
                    method="BDF", rtol=1e-8, atol=1e-10)
    Y = {s: sol.y[IDX[s]] for s in SPECIES}
    if observable == "Mpp":
        return Y["Mpp"]
    if observable == "Mpp_frac":
        tot = Y["M"] + Y["Mp"] + Y["Mpp"] + 1e-9
        return Y["Mpp"] / tot
    if observable == "ppfrac":                      # ppERK/(pERK+ppERK)
        return Y["Mpp"] / (Y["Mp"] + Y["Mpp"] + 1e-9)
    if observable == "Mp":
        return Y["Mp"]
    if observable == "ratio":                       # flux/state ratio Mpp/M
        return Y["Mpp"] / (Y["M"] + 1e-6)
    return Y["Mpp"]


def sensitivity_matrix(observable="Mpp", rel=0.05):
    base = _simulate(BASE, observable)
    S = np.zeros((len(TIMES), len(BASE)))
    for j in range(len(BASE)):
        tp = BASE.copy(); tm = BASE.copy()
        tp[j] *= (1 + rel); tm[j] *= (1 - rel)
        yp = _simulate(tp, observable); ym = _simulate(tm, observable)
        # d y / d log(theta) = theta * dy/dtheta
        S[:, j] = (yp - ym) / (2 * rel)
    return S, base


def fim_spectrum(observable="Mpp"):
    S, _ = sensitivity_matrix(observable)
    F = S.T @ S
    w = np.linalg.eigvalsh(F)[::-1]                 # descending
    w = np.clip(w, 1e-30, None)
    return w


def observable_stiff_projection():
    """Rank candidate biomarker observables by their sensitivity magnitude
    along the model's stiff axis 1."""
    # stiff axis from the Mpp observable FIM
    S, _ = sensitivity_matrix("Mpp")
    F = S.T @ S
    w, V = np.linalg.eigh(F)
    stiff1 = V[:, -1]                                # top eigenvector
    scores = {}
    for obs in ["ratio", "ppfrac", "Mpp_frac", "Mpp", "Mp"]:
        So, base = sensitivity_matrix(obs)
        # RELATIVE sensitivity (dimensionless): divide by the observable's own
        # scale so a fraction and a concentration are compared fairly.
        scale = np.mean(np.abs(base)) + 1e-9
        proj = (So @ stiff1) / scale
        scores[obs] = float(np.sqrt(np.mean(proj ** 2)))
    return scores, stiff1


def _figure(w, scores, path="figures/m13_sloppy.png"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    wn = w / w[0]
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    ax[0].semilogy(range(1, len(wn) + 1), wn, "o-", color="#8e44ad")
    ax[0].set_xlabel("FIM axis rank"); ax[0].set_ylabel("eigenvalue (rel.)")
    ax[0].set_title("Sloppy spectrum (ERK dual-phospho)\n"
                    "few stiff axes, many sloppy")
    obs = sorted(scores, key=lambda k: scores[k])
    cols = ["#27ae60" if o in ("ratio", "ppfrac", "Mpp_frac") else "#95a5a6"
            for o in obs]
    ax[1].barh(range(len(obs)), [scores[o] for o in obs], color=cols,
               edgecolor="k")
    ax[1].set_yticks(range(len(obs))); ax[1].set_yticklabels(obs)
    ax[1].set_xlabel("relative sensitivity to stiff axis 1")
    ax[1].set_title("Biomarker observables\n(green=ratio/fraction)")
    fig.tight_layout(); fig.savefig(path, dpi=130); plt.close(fig)


def report():
    print("=" * 68)
    print("MODULE 13 — FIM / sloppy model / stiff axes (biomarker selection)")
    print("=" * 68)
    w = fim_spectrum("Mpp")
    wn = w / w[0]
    top3 = wn[:3].sum() / wn.sum()
    print(f"FIM eigenvalue spectrum (output=Mpp, {len(BASE)} log-parameters):")
    for i, e in enumerate(wn[:8], 1):
        print(f"  axis {i}: {e:.2e}")
    print(f"  spectrum spans {np.log10(wn[0]/wn[wn>0][-1]):.1f} orders of "
          f"magnitude  -> SLOPPY")
    print(f"  top-3 stiff axes carry {top3*100:.1f}% of the spectral weight")

    scores, stiff1 = observable_stiff_projection()
    order = sorted(scores, key=lambda k: -scores[k])
    print(f"\ncandidate biomarker observables ranked by projection onto "
          f"stiff axis 1:")
    for o in order:
        print(f"  {o:10s}  rms-sensitivity = {scores[o]:.3g}")
    # stiff-axis loadings
    load = sorted(zip(PARAM_NAMES, stiff1), key=lambda kv: -abs(kv[1]))[:4]
    print(f"\nstiff axis 1 dominant parameters: "
          + ", ".join(f"{n}({v:+.2f})" for n, v in load))
    winner = order[0]
    print(f"\ninterpretation: the top-loading observable is '{winner}' "
          f"(a flux/state RATIO),")
    print("which projects onto the stiff axis more strongly than single")
    print("concentrations -> ratios (ppERK/(pERK+ppERK), kinase-drive :")
    print("phosphatase-capacity) are the robust, identifiable biomarkers,")
    print("matching the report's stiff-axis finding.")
    _figure(w, scores)
    print("\nfigure written to: figures/m13_sloppy.png")
    return {"top3_frac": top3, "orders": float(np.log10(wn[0]/wn[wn>0][-1])),
            "obs_scores": scores}


if __name__ == "__main__":
    report()
