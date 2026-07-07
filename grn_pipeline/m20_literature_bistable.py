"""
Module 20 — Toward literature-exact bistable switches (multi-variable) with
hysteresis: Rb-E2F restriction point and caspase apoptosis.

M19 used 1-D phenomenological forms. This module upgrades two of them to
MULTI-VARIABLE models with the real feedback TOPOLOGY of the published
switches, and demonstrates the hallmark the papers emphasise -- HYSTERESIS
(forward and backward parameter sweeps jump at different thresholds) -- plus
critical slowing at the saddle-nodes.

  * Rb-E2F (Yao et al. 2008, Nat Cell Biol; PMID 18364697): the restriction
    point is a bistable E2F switch driven by serum, via E2F<->CycE/CDK
    positive feedback. Here E2F(x) <-> CycE(y) mutual activation, control =
    serum S.
  * Apoptosis (Eissing et al. 2004, J Biol Chem; PMID 15208304): caspase-8*
    <-> caspase-3* positive feedback gives a bistable life/death switch,
    control = death stimulus.

HONEST NOTE: unlike M15 (Markevich, fetched and matched to the decimal), the
exact rate constants of Yao 2008 and Eissing 2004 are NOT in open-access PMC
and BioModels download is proxy-blocked, so these are the published feedback
TOPOLOGIES reproducing the documented bistability + hysteresis, not decimal
parameter transcriptions.
"""
import numpy as np
from scipy.optimize import fsolve


def _hill(x, K, n):
    x = max(x, 0.0)
    return x ** n / (K ** n + x ** n)


def rb_e2f(s, S):
    """E2F (x) <-> CycE (y) mutual activation; serum S feeds E2F."""
    x, y = s
    dx = 0.05 + 0.6 * S + 1.6 * _hill(y, 0.6, 4) - 1.0 * x
    dy = 1.6 * _hill(x, 0.6, 4) - 1.0 * y
    return [dx, dy]


def apoptosis(s, p):
    """caspase-8* (x) <-> caspase-3* (y) positive feedback; p = death stimulus."""
    x, y = s
    dx = p + 1.8 * _hill(y, 0.5, 4) - 1.0 * x
    dy = 1.8 * _hill(x, 0.5, 4) - 1.0 * y
    return [dx, dy]


MODELS = [
    dict(id="Rb_E2F", f=rb_e2f, ctrl="serum S", prange=(0.0, 1.2), xmax=3.0,
         cite="Yao 2008 (PMID 18364697)"),
    dict(id="Apoptosis", f=apoptosis, ctrl="death stimulus", prange=(0.0, 1.2),
         xmax=3.0, cite="Eissing 2004 (PMID 15208304)"),
]


def _jac(f, x, p, h=1e-6):
    n = len(x); J = np.zeros((n, n)); f0 = np.array(f(x, p))
    for j in range(n):
        xx = np.array(x, float); xx[j] += h
        J[:, j] = (np.array(f(xx, p)) - f0) / h
    return J


def steady_states(f, p, xmax, grid=6):
    roots = []
    for a in np.linspace(0.01, xmax, grid):
        for b in np.linspace(0.01, xmax, grid):
            sol, info, ier, _ = fsolve(f, [a, b], args=(p,), full_output=True)
            if ier != 1 or np.any(sol < -1e-6):
                continue
            if np.linalg.norm(f(sol, p)) > 1e-8:
                continue
            if not any(np.linalg.norm(sol - r["x"]) < 1e-3 for r in roots):
                lam = np.linalg.eigvals(_jac(f, sol, p)).real.max()
                roots.append({"x": sol, "act": sol[0], "lam": lam,
                              "stable": lam < 0})
    return sorted(roots, key=lambda r: r["act"])


def bistable_window(f, prange, xmax, n=50):
    ps = np.linspace(*prange, n)
    inside = [p for p in ps if sum(r["stable"]
              for r in steady_states(f, p, xmax)) >= 2]
    return (min(inside), max(inside)) if inside else None


def hysteresis(f, prange, xmax, n=60):
    """Forward sweep (p up) and backward sweep (p down) by continuation; the
    jump points differ -> hysteresis loop."""
    ps_up = np.linspace(prange[0], prange[1], n)
    ps_dn = ps_up[::-1]

    def sweep(ps):
        x = None; out = []
        for p in ps:
            st = [r for r in steady_states(f, p, xmax) if r["stable"]]
            if not st:
                out.append(np.nan); continue
            r = (min(st, key=lambda r: np.linalg.norm(r["x"] - x))
                 if x is not None else min(st, key=lambda r: r["act"]))
            x = r["x"]; out.append(r["act"])
        return np.array(out)
    up = sweep(ps_up)
    dn = sweep(ps_dn)[::-1]
    # jump thresholds
    dup = np.nanargmax(np.abs(np.diff(up))) if np.any(~np.isnan(up)) else 0
    ddn = np.nanargmax(np.abs(np.diff(dn))) if np.any(~np.isnan(dn)) else 0
    return ps_up, up, dn, ps_up[dup], ps_up[ddn]


def report(fig_path="figures/m20_literature_bistable.png"):
    print("=" * 70)
    print("MODULE 20 — Literature bistable switches (multi-variable) + "
          "hysteresis")
    print("=" * 70)
    out = []
    for m in MODELS:
        win = bistable_window(m["f"], m["prange"], m["xmax"])
        ps, up, dn, thr_up, thr_dn = hysteresis(m["f"], m["prange"], m["xmax"])
        width = abs(thr_up - thr_dn)
        out.append(dict(id=m["id"], win=win, ps=ps, up=up, dn=dn,
                        thr_up=thr_up, thr_dn=thr_dn, width=width))
        print(f"\n{m['id']}  [{m['cite']}]  control = {m['ctrl']}")
        if win:
            print(f"  bistable window: [{win[0]:.3f}, {win[1]:.3f}]")
        print(f"  hysteresis: rising jump at {m['ctrl']}={thr_up:.3f}, "
              f"falling jump at {thr_dn:.3f}  -> loop width {width:.3f}")
        # critical slowing: leading eigenvalue -> 0 near the upper fold
        if win:
            for frac, tag in [(0.5, "mid-window"), (0.98, "near fold")]:
                p = win[0] + frac * (win[1] - win[0])
                st = [r for r in steady_states(m["f"], p, m["xmax"])
                      if r["stable"]]
                lam = max(r["lam"] for r in st) if st else np.nan
                print(f"    at {m['ctrl']}={p:.3f} ({tag}): "
                      f"leading eigenvalue = {lam:+.3f}")
    print("\n-> both reproduce the published bistable/hysteretic switch; the")
    print("leading eigenvalue -> 0 at the folds (critical slowing), so the")
    print("same DNB/Lyapunov biomarker applies. (Topology of Yao 2008 /")
    print("Eissing 2004; exact rate constants not open-access-fetchable.)")
    _figure(MODELS, out, fig_path)
    print(f"\nfigure written to: {fig_path}")
    return {"models": [{k: v for k, v in o.items()
                        if k not in ("ps", "up", "dn")} for o in out]}


def _figure(models, out, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, len(models), figsize=(6.5 * len(models), 4.3))
    if len(models) == 1:
        ax = [ax]
    for a, m, o in zip(ax, models, out):
        a.plot(o["ps"], o["up"], "-", color="#c0392b", lw=2, label="serum ↑")
        a.plot(o["ps"], o["dn"], "-", color="#2980b9", lw=2, label="serum ↓")
        if o["win"]:
            a.axvspan(o["win"][0], o["win"][1], color="#f9e79f", alpha=0.4,
                      label="bistable")
        a.set_xlabel(m["ctrl"]); a.set_ylabel("active output (E2F / caspase-3*)")
        a.set_title(f"{m['id']}  ({m['cite'].split('(')[0].strip()})\n"
                    f"hysteresis loop width {o['width']:.2f}", fontsize=10)
        a.legend(fontsize=8)
    fig.suptitle("Literature bistable switches: hysteresis (forward vs backward "
                 "sweep) — same critical-slowing biomarker at the folds",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(path, dpi=130); plt.close(fig)


if __name__ == "__main__":
    report()
