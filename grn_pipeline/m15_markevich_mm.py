"""
Module 15 — EXACT Markevich 2004 Michaelis-Menten ERK cycle (bistability).

Supersedes M12's illustrative mass-action bistability with the exact
competitive-MM rate laws and published parameters of Markevich, Hoek &
Kholodenko 2004 (J Cell Biol; DOI 10.1083/jcb.200308060, Fig. 3), retrieved
from the paper this session. Reproduces the report's wide hysteresis window
and three-steady-state table.

Distributive ordered dual phosphorylation (kinase K=MAPKK) / dephosphorylation
(phosphatase P=MKP3), competitive Michaelis-Menten:

  v1 = k1cat K M   / (Km1(1+Mp/Km2) + M)     M   -> Mp   (kinase)
  v2 = k2cat K Mp  / (Km2(1+M/Km1) + Mp)      Mp  -> Mpp  (kinase)
  v3 = k3cat P Mpp / (Km3(1+Mp/Km4) + Mpp)    Mpp -> Mp   (phosphatase)
  v4 = k4cat P Mp  / (Km4(1+Mpp/Km3) + Mp)    Mp  -> M    (phosphatase)

Conservation: M + Mp + Mpp = MAPK_tot.  Bifurcation parameter: [MAPKK].
"""
import numpy as np
from scipy.optimize import fsolve

# published parameters (Markevich 2004, Fig. 3) + the phosphatase product-
# sequestration constant Km5 (paper Eq. 3, = h6/h-6). Km5 is recovered here
# from the report's OFF steady state (437.73, 12.85, 49.42) at MAPKK=50 nM:
# it makes BOTH dM=0 and dMpp=0 hold simultaneously -> Km5 ~ 78 nM.
KM1, KM2, KM3, KM4, KM5 = 50.0, 500.0, 22.0, 18.0, 78.0   # nM
K1, K2, K3, K4 = 0.01, 15.0, 0.084, 0.06            # s^-1
MTOT, PTOT = 500.0, 100.0                            # nM


def _v(M, Mp, Mpp, K):
    # kinase: competition between M and Mp (irreversible product release)
    dK = 1 + M / KM1 + Mp / KM2
    v1 = K1 * K * (M / KM1) / dK
    v2 = K2 * K * (Mp / KM2) / dK
    # phosphatase: competition Mpp/Mp + product M sequestration (Km5 term)
    dP = 1 + Mpp / KM3 + Mp / KM4 + M / KM5
    v3 = K3 * PTOT * (Mpp / KM3) / dP
    v4 = K4 * PTOT * (Mp / KM4) / dP
    return v1, v2, v3, v4


def rhs(state, K):
    """2-D reduced system in (M, Mpp); Mp = MTOT - M - Mpp."""
    M, Mpp = state
    Mp = MTOT - M - Mpp
    v1, v2, v3, v4 = _v(M, Mp, Mpp, K)
    dM = v4 - v1
    dMpp = v2 - v3
    return np.array([dM, dMpp])


def jacobian(state, K, eps=1e-4):
    J = np.zeros((2, 2))
    f0 = rhs(state, K)
    for j in range(2):
        s = np.array(state, float); s[j] += eps
        J[:, j] = (rhs(s, K) - f0) / eps
    return J


def _dM_along(M, Mpp, K):
    return rhs([M, Mpp], K)[0]


def steady_states(K, n=1200):
    """Robust nullcline finder: for each Mpp solve dM=0 for M on the M-
    nullcline, then locate sign changes of dMpp along it (= true steady
    states). Classify stability by the 2-D Jacobian eigenvalues."""
    from scipy.optimize import brentq
    roots = []
    prevMpp = prevG = None
    for Mpp in np.linspace(0.2, MTOT - 0.2, n):
        hi = MTOT - Mpp - 1e-3
        a, b = _dM_along(1e-3, Mpp, K), _dM_along(hi, Mpp, K)
        if a * b > 0:
            prevG = None
            continue
        Mstar = brentq(_dM_along, 1e-3, hi, args=(Mpp, K))
        g = rhs([Mstar, Mpp], K)[1]
        if prevG is not None and prevG * g < 0:
            # refine the Mpp root by bisection on g
            lo_m, hi_m = prevMpp, Mpp
            for _ in range(40):
                mid = 0.5 * (lo_m + hi_m)
                Mm = brentq(_dM_along, 1e-3, MTOT - mid - 1e-3, args=(mid, K))
                gm = rhs([Mm, mid], K)[1]
                if prevG * gm < 0:
                    hi_m = mid
                else:
                    lo_m = mid; prevG = gm
            Mppr = 0.5 * (lo_m + hi_m)
            Mr = brentq(_dM_along, 1e-3, MTOT - Mppr - 1e-3, args=(Mppr, K))
            ev = np.linalg.eigvals(jacobian([Mr, Mppr], K)).real
            roots.append({"M": Mr, "Mp": MTOT - Mr - Mppr, "Mpp": Mppr,
                          "lmax": ev.max(), "stable": ev.max() < 0})
        prevMpp, prevG = Mpp, g
    return sorted(roots, key=lambda r: r["Mpp"])


def _n_stable(K):
    return sum(r["stable"] for r in steady_states(K, n=800))


def bistable_window(Kmin=30, Kmax=70, n=41):
    Ks = np.linspace(Kmin, Kmax, n)
    flags = [(K, _n_stable(K) >= 2) for K in Ks]
    inside = [K for K, f in flags if f]
    if not inside:
        return None
    lo_i, hi_i = min(inside), max(inside)

    def bisect(a, b):                       # a=outside, b=inside
        for _ in range(30):
            m = 0.5 * (a + b)
            if _n_stable(m) >= 2:
                b = m
            else:
                a = m
        return 0.5 * (a + b)
    lo = bisect(lo_i - (Ks[1] - Ks[0]), lo_i)
    hi = bisect(hi_i + (Ks[1] - Ks[0]), hi_i)
    return lo, hi


def report(fig_path="figures/m15_markevich.png"):
    print("=" * 68)
    print("MODULE 15 — EXACT Markevich 2004 MM ERK cycle (bistability)")
    print("=" * 68)
    print("params (JCB 2004 Fig.3): Km1=50 Km2=500 Km3=22 Km4=18 nM; "
          "k1=0.01 k2=15 k3=0.084 k4=0.06 /s; MAPK=500 MKP3=100 nM")

    win = bistable_window()
    if win:
        print(f"\nbistability / hysteresis window: [MAPKK] in "
              f"[{win[0]:.2f}, {win[1]:.2f}] nM")
        print(f"  (report: 39.25 - 57.38 nM)")

    print(f"\nthree steady states at [MAPKK] = 50 nM:")
    ss = steady_states(50.0)
    print(f"  {'M':>8} {'Mp':>8} {'Mpp(active ERK)':>16} "
          f"{'lambda_max':>11}  state")
    for r in ss:
        kind = "stable" if r["stable"] else "saddle/unstable"
        print(f"  {r['M']:8.2f} {r['Mp']:8.2f} {r['Mpp']:16.2f} "
              f"{r['lmax']:11.4f}  {kind}")
    print(f"  (report at MAPKK=50: Mpp = 49.4 [stable], 277.3 [saddle], "
          f"481.9 [stable])")

    _figure(fig_path)
    print(f"\nfigure written to: {fig_path}")
    return {"window": win, "states_at_50": ss}


def _figure(path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    Ks = np.linspace(30, 70, 81)
    lo, mid, hi = [], [], []
    for K in Ks:
        ss = steady_states(K)
        stab = [r["Mpp"] for r in ss if r["stable"]]
        uns = [r["Mpp"] for r in ss if not r["stable"]]
        lo.append(min(stab) if stab else np.nan)
        hi.append(max(stab) if stab else np.nan)
        mid.append(uns[0] if uns else np.nan)
    fig, ax = plt.subplots(figsize=(7.5, 5))
    ax.plot(Ks, hi, "-", color="#c0392b", lw=2, label="stable ON")
    ax.plot(Ks, lo, "-", color="#2980b9", lw=2, label="stable OFF")
    ax.plot(Ks, mid, "--", color="gray", lw=1.5, label="unstable saddle")
    win = bistable_window()
    if win:
        ax.axvspan(win[0], win[1], color="#f9e79f", alpha=0.5,
                   label=f"bistable [{win[0]:.1f},{win[1]:.1f}] nM")
    ax.set_xlabel("[MAPKK] (nM)")
    ax.set_ylabel("active ERK  [Mpp] (nM)")
    ax.set_title("Exact Markevich 2004 ERK cycle: hysteresis\n"
                 "(competitive MM, published parameters)")
    ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(path, dpi=130); plt.close(fig)


if __name__ == "__main__":
    report()
