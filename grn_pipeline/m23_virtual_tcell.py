"""M23 - a mechanistic "virtual T cell": kinetic proofreading -> ERK switch -> IL-2.

This module exists to answer one question honestly: *what can a computational
T cell actually predict, and where does it stop?*

Three layers, each with an explicit rigour label:

  L1 (RIGOROUS)   McKeithan kinetic proofreading (PNAS 1995, 92:5042-5046).
                  Steady-state occupancy of the N-th modified TCR-pMHC complex.
                  Analytic result: d log(C_N)/d log(tau) = N*k_off/(k_off+k_p) + 1,
                  which -> N+1 when k_off >> k_p and collapses to 1 when it does
                  not. The chain amplifies a lifetime difference into an
                  (N+1)-th power output difference, but only in that regime.
                  Verified here by finite difference against the closed form.

  L2 (LITERATURE-GROUNDED TOPOLOGY, ILLUSTRATIVE PARAMETERS)
                  Altan-Bonnet & Germain (PLoS Biol 2005, 3:e356) added a
                  SHP-1 negative feedback and an ERK positive feedback on top
                  of proofreading. The topology - not the parameter set - is
                  what is reproduced here: a digital (bistable) ERK response
                  with a sharp threshold in ligand QUALITY that is only weakly
                  moved by ligand QUANTITY. Parameters are illustrative.

  L3 (RELATIVE ONLY)
                  IL-2 output is read out as a monotone function of the ERK-ON
                  fraction. This gives a RANKING of ligands and a threshold
                  position. It does NOT give pg/mL, and it is not calibrated to
                  any donor, hybridoma clone or assay.

What this module deliberately does not do: predict absolute IL-2 concentration,
predict a human donor's response, substitute for an MLR, an ELISPOT, a tetramer
precursor-frequency count, or HLA typing. See REPORT.md M23 for why.

Run:  python3 -m grn_pipeline.m23_virtual_tcell
"""

from __future__ import annotations

import os

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

FIGDIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "figures")


# ---------------------------------------------------------------- L1: proofreading
def kpr_output(tau_s: float, L: float, n_steps: int, k_p: float = 1.0,
               k_on: float = 1.0, R: float = 1.0) -> float:
    """Steady-state occupancy of the final (N-th) proofreading state.

    McKeithan's chain:  R + L  <-> C_0 -> C_1 -> ... -> C_N,
    every C_i unbinding at k_off = 1/tau, each forward step at k_p.

    Steady state gives  C_0 = k_on*L*R/(k_off + k_p)  and
    C_i = C_{i-1} * k_p/(k_off + k_p) for 0 < i < N, with the terminal state
    only losing ligand:  C_N = C_{N-1} * k_p/k_off.
    """
    k_off = 1.0 / tau_s
    phi = k_p / (k_off + k_p)
    c0 = k_on * L * R / (k_off + k_p)
    if n_steps == 0:
        return c0
    return c0 * phi ** (n_steps - 1) * (k_p / k_off)


def discrimination_exponent(tau_s: float, n_steps: int, k_p: float = 1.0,
                            eps: float = 1e-4) -> float:
    """d log C_N / d log tau, by central difference. Should approach N+1 for
    k_off >> k_p and 1 for k_off << k_p (the chain saturates)."""
    lo = kpr_output(tau_s * (1 - eps), 1.0, n_steps, k_p)
    hi = kpr_output(tau_s * (1 + eps), 1.0, n_steps, k_p)
    return (np.log(hi) - np.log(lo)) / (np.log(tau_s * (1 + eps)) - np.log(tau_s * (1 - eps)))


# ------------------------------------------------- L2: proofreading + feedback (ERK)
class VirtualTCell:
    """Proofreading output drives an ERK module with positive feedback and a
    SHP-1-like negative feedback. Topology after Altan-Bonnet & Germain 2005;
    parameters illustrative, chosen only so the qualitative regimes appear.

    States:  E = active ERK fraction [0,1],  S = active SHP-1 [0,1].
    """

    def __init__(self, n_steps: int = 4, k_p: float = 1.0,
                 alpha: float = 12.0,      # drive from proofreading output
                 beta: float = 12.0,       # ERK positive feedback strength
                 K_e: float = 0.45,        # feedback half-point
                 hill: float = 4.0,        # feedback cooperativity
                 gamma: float = 1.0,       # ERK deactivation
                 s_gain: float = 8.0,      # SHP-1 production from EARLY complexes
                 s_decay: float = 1.0,
                 s_inhib: float = 6.0):    # SHP-1 inhibition of ERK activation
        self.n_steps, self.k_p = n_steps, k_p
        self.alpha, self.beta, self.K_e, self.hill = alpha, beta, K_e, hill
        self.gamma = gamma
        self.s_gain, self.s_decay, self.s_inhib = s_gain, s_decay, s_inhib

    def drives(self, tau_s: float, L: float):
        """(specific drive, non-specific drive). The specific drive is the
        proofread output; the non-specific one is total early occupancy, which
        is what feeds the negative feedback - this is the ingredient that makes
        the response depend on ligand QUALITY rather than just occupancy."""
        k_off = 1.0 / tau_s
        c_total_early = L / (k_off + self.k_p)          # ~ C_0, quality-blind
        return kpr_output(tau_s, L, self.n_steps, self.k_p), c_total_early

    def rhs(self, _t, y, tau_s, L):
        E, S = y
        spec, nonspec = self.drives(tau_s, L)
        fb = self.beta * E ** self.hill / (self.K_e ** self.hill + E ** self.hill)
        act = (self.alpha * spec + fb) / (1.0 + self.s_inhib * S)
        dE = act * (1.0 - E) - self.gamma * E
        dS = self.s_gain * nonspec * (1.0 - S) - self.s_decay * S
        return [dE, dS]

    def steady_state(self, tau_s: float, L: float, e0: float = 0.0, t_end: float = 400.0):
        sol = solve_ivp(self.rhs, (0, t_end), [e0, 0.0], args=(tau_s, L),
                        method="LSODA", rtol=1e-8, atol=1e-10, dense_output=False)
        return float(sol.y[0, -1]), float(sol.y[1, -1])

    def bistable(self, tau_s: float, L: float, tol: float = 0.05) -> bool:
        """Two stable ERK states at the same input -> a digital cell."""
        lo, _ = self.steady_state(tau_s, L, e0=0.0)
        hi, _ = self.steady_state(tau_s, L, e0=1.0)
        return (hi - lo) > tol

    def tau_threshold(self, L: float, lo: float = 0.05, hi: float = 100.0):
        """Ligand-lifetime threshold for ERK-ON, from the low-E branch
        (i.e. a resting cell being triggered). None if no crossing in range."""
        f = lambda t: self.steady_state(t, L, e0=0.0)[0] - 0.5
        if f(lo) * f(hi) > 0:
            return None
        return brentq(f, lo, hi, xtol=1e-3)


def il2_readout(E_on_fraction: float, e_half: float = 0.5, n: float = 3.0) -> float:
    """RELATIVE IL-2 (arbitrary units, 0-1). A monotone saturating function of
    the ERK-ON signal. There is no honest way to put pg/mL on this axis without
    a calibration curve from the assay it is meant to replace."""
    x = max(E_on_fraction, 0.0)
    return x ** n / (e_half ** n + x ** n)


# -------------------------------------------------------------------------- report
def main() -> dict:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out: dict = {}
    print("=" * 74)
    print("M23  virtual T cell: kinetic proofreading -> ERK switch -> relative IL-2")
    print("=" * 74)

    # ---- L1 discrimination exponent -------------------------------------------
    print("\n[L1 RIGOROUS] McKeithan proofreading: discrimination exponent")
    print("  exact: d log C_N/d log tau = N*k_off/(k_off+k_p) + 1  -> N+1 as k_off >> k_p")
    rows = []
    for n in (0, 1, 2, 4, 6):
        d_lim = discrimination_exponent(1e-3, n)   # k_off = 1000 >> k_p = 1
        d_mid = discrimination_exponent(0.1, n)    # k_off = 10
        d_slow = discrimination_exponent(50.0, n)  # k_off = 0.02 << k_p
        rows.append((n, d_lim, d_mid, d_slow))
        print(f"    N={n}:  {d_lim:6.3f} (k_off=1000)   {d_mid:6.3f} (k_off=10)   "
              f"{d_slow:6.3f} (k_off=0.02)")
    out["exponents"] = rows
    for n, d_lim, _, d_slow in rows:
        assert abs(d_lim - (n + 1)) < 0.05, f"N={n}: exponent should approach N+1"
        if n >= 1:  # N=0 has no proofreading step, so it saturates at 0 instead
            assert abs(d_slow - 1.0) < 0.15, f"N={n}: exponent should saturate at 1"
    print("  -> in the fast-off-rate regime the chain converts a 2x lifetime")
    print(f"     difference into 2^(N+1): N=4 gives {2**5:.0f}x output difference.")
    print("  -> in the slow-off-rate regime the exponent collapses to 1: proofreading")
    print("     stops discriminating once every complex survives the chain.")

    # ---- L2 digital response ---------------------------------------------------
    print("\n[L2 ILLUSTRATIVE] ERK response: quality threshold vs quantity")
    cell = VirtualTCell(n_steps=4)
    taus = np.logspace(np.log10(0.2), np.log10(60), 60)
    doses = [0.3, 1.0, 3.0, 10.0]
    curves, thresholds = {}, {}
    for L in doses:
        curves[L] = np.array([cell.steady_state(t, L, e0=0.0)[0] for t in taus])
        thresholds[L] = cell.tau_threshold(L)
        th = thresholds[L]
        print(f"    ligand dose L={L:5.1f}:  tau_threshold = "
              + (f"{th:6.2f} s" if th else "  none in range"))
    valid = [v for v in thresholds.values() if v]
    if len(valid) > 1:
        spread = max(valid) / min(valid)
        dose_span = max(doses) / min(doses)
        print(f"  -> a {dose_span:.0f}x change in ligand QUANTITY moves the quality "
              f"threshold {spread:.2f}x.")
        out["threshold_spread"] = spread
        out["dose_span"] = dose_span
    out["thresholds"] = thresholds

    # hysteresis: scan for the tau window where a resting and a pre-activated
    # cell settle on different ERK states at the SAME ligand (memory of prior
    # stimulation - the reason "resting vs restimulated" is not one function).
    scan = np.logspace(np.log10(0.02), np.log10(60), 70)
    bist_taus = [t for t in scan if cell.bistable(t, 1.0)]
    if bist_taus:
        lo_b, hi_b = min(bist_taus), max(bist_taus)
        print(f"    bistable window at L=1: tau in [{lo_b:.2f}, {hi_b:.2f}] s "
              f"({len(bist_taus)}/{len(scan)} grid points)")
        print("    -> inside it the SAME ligand gives ERK-off or ERK-on depending on")
        print("       whether the cell was previously activated (hysteresis).")
        out["bistable_window"] = (lo_b, hi_b)
    else:
        print("    no bistable window found on this grid at L=1")
        out["bistable_window"] = None

    # ---- L3 ligand ranking ------------------------------------------------------
    print("\n[L3 RELATIVE ONLY] altered-peptide-ligand series -> relative IL-2")
    print("  (tau values are a SWEEP, not measured constants for any named peptide)")
    apl = [("agonist",        30.0),
           ("weak agonist",   10.0),
           ("partial agonist", 4.0),
           ("weak partial",    2.0),
           ("antagonist",      0.8),
           ("null",            0.3)]
    il2 = []
    for name, tau in apl:
        E, S = cell.steady_state(tau, 1.0, e0=0.0)
        il2.append((name, tau, E, il2_readout(E)))
        print(f"    {name:16s} tau={tau:5.1f}s  ERK={E:5.3f}  IL-2(rel)={il2[-1][3]:5.3f}")
    out["apl"] = il2
    ranks = [r[3] for r in il2]
    assert all(ranks[i] >= ranks[i + 1] - 1e-9 for i in range(len(ranks) - 1)), \
        "IL-2 ranking must be monotone in ligand lifetime"
    print("  -> ranking is monotone in tau. The ORDER is a prediction;")
    print("     the y-axis is arbitrary units and cannot be read as pg/mL.")

    # ---- figure -----------------------------------------------------------------
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))

    ax = axes[0]
    for n in (0, 1, 2, 4):
        y = [kpr_output(t, 1.0, n) for t in taus]
        ax.loglog(taus, y, label=f"N={n}")
    ax.set_xlabel("pMHC-TCR lifetime tau (s)")
    ax.set_ylabel("proofread output $C_N$ (a.u.)")
    ax.set_title("L1 rigorous: proofreading\nsharpens lifetime discrimination")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")

    ax = axes[1]
    for L in doses:
        ax.semilogx(taus, curves[L], label=f"L={L}")
        if thresholds[L]:
            ax.axvline(thresholds[L], ls=":", lw=0.8, color="grey")
    ax.set_xlabel("pMHC-TCR lifetime tau (s)")
    ax.set_ylabel("steady-state active ERK")
    ax.set_title("L2 illustrative: digital response,\nthreshold set by quality not quantity")
    ax.legend(fontsize=8, title="ligand dose")
    ax.grid(alpha=0.3)

    ax = axes[2]
    names = [r[0] for r in il2]
    vals = [r[3] for r in il2]
    ax.barh(range(len(names))[::-1], vals, color="#4c72b0")
    ax.set_yticks(range(len(names))[::-1])
    ax.set_yticklabels(names, fontsize=9)
    ax.set_xlabel("relative IL-2 (arbitrary units)")
    ax.set_title("L3 relative only: ligand RANKING\n(not pg/mL, not donor-specific)")
    ax.grid(alpha=0.3, axis="x")

    fig.suptitle("M23 virtual T cell - what a mechanistic model can and cannot predict",
                 fontsize=12)
    fig.tight_layout()
    os.makedirs(FIGDIR, exist_ok=True)
    path = os.path.join(FIGDIR, "m23_virtual_tcell.png")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"\n  figure -> {path}")

    print("\n" + "-" * 74)
    print("HONESTY LINE")
    print("  computed here : lifetime-discrimination exponent; a quality threshold;")
    print("                  a monotone ligand ranking; a hysteresis window.")
    print("  NOT computed  : pg/mL IL-2, a named peptide's tau, a donor's response,")
    print("                  MLR alloreactivity, ELISPOT spot counts, precursor")
    print("                  frequency, or anyone's HLA type.")
    print("-" * 74)
    return out


# the rest of the pipeline calls modules through report()
report = main


if __name__ == "__main__":
    main()
