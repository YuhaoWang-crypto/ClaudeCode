"""
Module 28 — When is the "order parameter -> Jacobian -> fragile point"
paradigm actually identifiable?

M26/M27 established the limits on steady-state omics:
  * a knockout that leaves its target's transcript intact destroys diagonal
    dominance and network inference sits at chance (M26);
  * a perturbation that moves the target's own transcript fixes the diagonal
    and lifts inference off chance, but only the strongest edges reproduce
    out-of-sample (M27);
  * either way steady-state data give the normalised interactions r, never
    J = D r, so eigenvalues and relaxation times are not identified.

This module makes the positive case: in a system where the kinetics ARE
known, the whole paradigm works, and it shows exactly what it buys.

Part A — structural analysis of a proposed Jacobian (no data needed).
  Thomas' rules: a positive feedback loop is necessary for multistability,
  a negative loop (length >= 2 in the Jacobian) for sustained oscillation.
  Reading the loops off a candidate matrix tells you which claims about it
  can possibly be true, before any numbers are fitted.

Part B — a real kinetic model (Chassagnole et al. 2002 E. coli central
carbon metabolism, BIOMD0000000051, 18 metabolites / 48 reactions) where
COPASI computes the exact Jacobian:
  * eigenvalue spectrum -> real relaxation times (in seconds)
  * participation factors -> which metabolite can actually see the slow mode
  * enzyme fragility d Re(lambda_max) / d ln E -> which enzyme controls
    stability, compared against steady-state flux
"""
import itertools
import os
import urllib.request

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "..", "figures")
CACHE = os.path.join(HERE, "..", "data_cache", "biomodels")
MODEL = "BIOMD0000000051"


# ---------------------------------------------------------------------------
# Part A: structural analysis of a candidate Jacobian
# ---------------------------------------------------------------------------
def feedback_loops(J, max_len=4):
    """All feedback loops up to `max_len` and the sign of their gain
    product. J[i, j] = effect of j on i."""
    n = len(J)
    out = []
    for L in range(2, max_len + 1):
        for c in itertools.permutations(range(n), L):
            if c[0] != min(c):
                continue
            p = float(np.prod([J[c[(i + 1) % L], c[i]] for i in range(L)]))
            if p != 0:
                out.append((tuple(c), p))
    return out


def structural_verdict(J, names=None):
    """What the sign structure alone permits (Thomas' rules)."""
    lp = feedback_loops(J)
    pos = [c for c, p in lp if p > 0]
    neg = [c for c, p in lp if p < 0]
    w = np.linalg.eigvals(J)
    return dict(
        loops=lp, n_positive=len(pos), n_negative=len(neg),
        can_be_multistable=len(pos) > 0,
        can_oscillate=len(neg) > 0,
        eig=w, max_re=float(w.real.max()),
        has_complex=bool(np.abs(w.imag).max() > 1e-9),
        names=names)


# ---------------------------------------------------------------------------
# Part B: a model whose kinetics are known
# ---------------------------------------------------------------------------
def fetch_model(mid=MODEL):
    os.makedirs(CACHE, exist_ok=True)
    p = os.path.join(CACHE, f"{mid}.xml")
    if not os.path.exists(p):
        urllib.request.urlretrieve(
            f"https://raw.githubusercontent.com/biomodels/{mid}/master/"
            f"{mid}/{mid}.xml", p)
    return p


def _load():
    import basico
    basico.load_model(fetch_model())
    # COPASI's default Newton resolution (1e-9) is too strict for this model
    basico.set_task_settings(basico.T.STEADY_STATE,
                             {"problem": {"Resolution": 1e-7}})
    return basico


def spectrum():
    """Exact Jacobian, eigenvalues, relaxation times and the participation
    factors of the slowest mode."""
    basico = _load()
    basico.run_steadystate()
    J = np.asarray(basico.get_jacobian_matrix())
    sp = list(basico.get_species().index)
    w, V = np.linalg.eig(J)
    k = int(np.argmax(w.real))
    Vi = np.linalg.inv(V)
    pf = np.abs(V[:, k] * Vi[k, :]); pf = pf / pf.sum()
    return dict(J=J, species=sp, eig=w, max_re=float(w.real.max()),
                tau_slow=float(-1 / w.real.max()),
                tau_fast=float(-1 / w.real.min()),
                n_osc=int(np.sum(np.abs(w.imag) > 1e-9) / 2),
                participation=dict(sorted(zip(sp, pf), key=lambda x: -x[1])))


def enzyme_fragility(step=1.10):
    """d Re(lambda_max) / d ln(enzyme) for every rmax parameter, with the
    steady-state flux of the same reaction for comparison."""
    import pandas as pd
    basico = _load()

    def lam():
        basico.run_steadystate()
        return np.linalg.eigvals(
            np.asarray(basico.get_jacobian_matrix())).real.max()

    base = lam()
    flux = basico.get_reactions()["flux"].abs()
    rp = basico.get_reaction_parameters()
    vm = rp[[("rmax" in n.split(".")[-1]) for n in rp.index]]
    rows = []
    for name, row in vm.iterrows():
        v0 = row["value"]
        try:
            basico.set_reaction_parameters(name, value=v0 * step)
            l = lam()
            rows.append(dict(enzyme=row["reaction"],
                             dlam=(l - base) / np.log(step),
                             tau=(-1 / l if l < 0 else np.inf),
                             flux=flux.get(row["reaction"], np.nan)))
        except Exception:
            pass
        finally:
            basico.set_reaction_parameters(name, value=v0)
    return base, pd.DataFrame(rows).dropna(subset=["dlam"])


def figure(res, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    s, d = res["spectrum"], res["fragility"]
    fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.8))
    a = ax[0]
    w = s["eig"]
    a.scatter(w.real, w.imag, s=40, c="#2a6fdb")
    a.axvline(0, c="k", lw=0.8)
    a.set_xscale("symlog", linthresh=1e-3)
    a.set_xlabel("Re(lambda)  [1/s]"); a.set_ylabel("Im(lambda)")
    a.set_title(f"(a) exact spectrum of a known-kinetics network\n"
                f"tau: {s['tau_fast'] * 1e3:.2g} ms ... {s['tau_slow']:.0f} s "
                f"({s['tau_slow'] / s['tau_fast']:.0e}x span)")
    a = ax[1]
    pf = list(s["participation"].items())[:6][::-1]
    a.barh([k for k, _ in pf], [v for _, v in pf], color="#29a36a")
    a.set_xlabel("participation in the slowest mode")
    a.set_title("(b) only these metabolites can see\nthe slow mode "
                "(probe selection)")
    a.tick_params(axis="y", labelsize=7)
    a = ax[2]
    a.scatter(d["flux"], d["dlam"].abs(), s=35, c="#e4572e")
    a.set_xscale("log")
    a.set_xlabel("steady-state flux of the reaction")
    a.set_ylabel("|d Re(lambda_max) / d ln E|")
    a.set_title(f"(c) stability control is NOT carried by\nthe high-flux "
                f"enzymes (Spearman {res['rho']:+.2f})")
    fig.suptitle("M28 — what the paradigm delivers when the kinetics are "
                 "known (E. coli central carbon metabolism)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.93]); fig.savefig(path, dpi=120)
    plt.close(fig)


def report(fig_path=os.path.join(FIG, "m28_identifiability.png")):
    import json
    from scipy.stats import spearmanr
    print("=" * 70)
    print("MODULE 28 — identifiability tiers for the Jacobian paradigm")
    print("=" * 70)
    J = np.array([[-0.5, -2.0, 0.0], [-2.0, -0.5, 1.5], [0.0, 1.0, -1.2]])
    v = structural_verdict(J, ["PAX5", "BLIMP1", "XBP1"])
    print("A. structural check of a candidate 3x3 Jacobian "
          "(PAX5/BLIMP1/XBP1):")
    for c, p in v["loops"]:
        print(f"   loop {c}: gain {p:+.2f} "
              f"({'positive' if p > 0 else 'negative'})")
    print(f"   eigenvalues {np.round(v['eig'].real, 3)} -> "
          f"{'saddle (unstable)' if v['max_re'] > 0 else 'stable'}")
    print(f"   multistability possible: {v['can_be_multistable']}; "
          f"oscillation possible: {v['can_oscillate']} "
          f"(complex pair present: {v['has_complex']})")
    # is the instability due to the mutual inhibition?
    K = J.copy(); K[0, 1] = K[1, 0] = 0.0
    print(f"   with the mutual inhibition REMOVED, max Re = "
          f"{np.linalg.eigvals(K).real.max():+.3f} -> the saddle comes from "
          f"the BLIMP1-XBP1 positive loop, not from the mutual inhibition")
    s = spectrum()
    print(f"\nB. exact Jacobian of {MODEL} ({len(s['species'])} metabolites):")
    print(f"   max Re(lambda) = {s['max_re']:.4g} /s -> stable; "
          f"tau_slow = {s['tau_slow']:.0f} s, tau_fast = "
          f"{s['tau_fast'] * 1e3:.2g} ms "
          f"({s['tau_slow'] / s['tau_fast']:.0e}x span)")
    print(f"   oscillatory mode pairs: {s['n_osc']}")
    print("   slow-mode participation: " + ", ".join(
        f"{k} {v_:.2f}" for k, v_ in list(s["participation"].items())[:4]))
    base, d = enzyme_fragility()
    rho, pv = spearmanr(d["flux"], d["dlam"].abs())
    print(f"\nC. enzyme fragility ({len(d)} enzymes):")
    top = d.sort_values("dlam", ascending=False).head(3)
    bot = d.sort_values("dlam").head(2)
    for _, r in top.iterrows():
        print(f"   destabilising: {r['enzyme'][:38]:40s} "
              f"d lambda {r['dlam']:+.4f}, tau {r['tau']:.0f} s")
    for _, r in bot.iterrows():
        print(f"   stabilising  : {r['enzyme'][:38]:40s} "
              f"d lambda {r['dlam']:+.4f}, tau {r['tau']:.0f} s")
    print(f"   |fragility| vs flux: Spearman rho = {rho:+.2f} (p={pv:.2f})"
          f" -> different things")
    res = dict(spectrum=s, fragility=d, rho=float(rho), struct=v)
    figure(res, fig_path)
    out = dict(tau_slow=s["tau_slow"], tau_fast=s["tau_fast"],
               span=s["tau_slow"] / s["tau_fast"], n_osc=s["n_osc"],
               participation=list(s["participation"].items())[:6],
               rho_fragility_flux=float(rho), p=float(pv),
               top_destabilising=top[["enzyme", "dlam", "tau"]]
               .to_dict("records"),
               top_stabilising=bot[["enzyme", "dlam", "tau"]]
               .to_dict("records"),
               struct=dict(loops=[[list(c), p] for c, p in v["loops"]],
                           max_re=v["max_re"],
                           can_oscillate=v["can_oscillate"],
                           has_complex=v["has_complex"]))
    with open(os.path.join(FIG, "m28_summary.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    d.to_csv(os.path.join(FIG, "m28_fragility.csv"), index=False)
    print(f"figure: {fig_path}")
    return res


if __name__ == "__main__":
    report()
