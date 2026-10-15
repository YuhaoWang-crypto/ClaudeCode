"""Tables and figures from results/cv_results_dt{dt}.json and results/r2_gates_dt{dt}.json."""
import json
import os
import sys

import numpy as np

import fit_cv as fc

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")
FIG = os.path.join(HERE, "figures")
ROUTE_ORDER = ("R0", "R1", "R2", "2pop")


def curve_table(res, fold_kind, beta=5.0, kind="B"):
    """Held-out per-curve metrics: {route: {cond: metrics}}."""
    tab = {}
    if kind == "B":
        for r in res["fitB"]:
            if r["fold_kind"] == fold_kind and r["beta"] == beta:
                tab.setdefault(r["route"], {}).update(r["test_metrics"])
    else:
        for r in res["fitA"]:
            if r["fold_kind"] == fold_kind:
                for rn, v in r["routes"].items():
                    tab.setdefault(rn, {}).update(v["test_metrics"])
    return tab


def agg(tab, key="nrmse_net"):
    return {rn: float(np.mean([m[key] for m in v.values()])) for rn, v in tab.items()}


def md_table(tab, title, conds):
    routes = [r for r in ROUTE_ORDER if r in tab]
    lines = [f"**{title}**", "", "| held-out curve | " + " | ".join(f"{r} NRMSE_net (peak err)" for r in routes) + " |",
             "|---|" + "---|" * len(routes)]
    for c in conds:
        lines.append(f"| {c} | " + " | ".join(f"{100*tab[r][c]['nrmse_net']:.1f}% ({100*tab[r][c]['peak_rel_err_net']:+.0f}%)"
                                              for r in routes) + " |")
    a = agg(tab)
    af = agg(tab, "nrmse_fwd")
    ab = agg(tab, "nrmse_bwd")
    lines.append("| **mean (equal weight)** | " + " | ".join(f"**{100*a[r]:.1f}%**" for r in routes) + " |")
    lines.append("| mean NRMSE fwd / bwd | " + " | ".join(f"{100*af[r]:.1f}% / {100*ab[r]:.1f}%" for r in routes) + " |")
    lines.append("| curves meeting NRMSE≤30 % and |peak err|≤30 % | " + " | ".join(
        f"{sum(1 for m in tab[r].values() if m['nrmse_net'] <= 0.30 and abs(m['peak_rel_err_net']) <= 0.30)}/{len(tab[r])}"
        for r in routes) + " |")
    if "R0" in a:
        lines.append("| relative change vs R0 | " + " | ".join(f"{100*(a[r]/a['R0']-1):+.1f}%" for r in routes) + " |")
    return "\n".join(lines)


CONDS = [f"f{f}_A{a}" for f in (20, 100, 500) for a in (25, 50, 75)]


def params_table(res, beta=5.0):
    lines = ["| route | fold | params (best) | Q (C) | train cost | starts within 1 % | at bound |", "|---|---|---|---|---|---|---|"]
    for r in sorted((r for r in res["fitB"] if r["beta"] == beta), key=lambda r: (ROUTE_ORDER.index(r["route"]), r["fold"])):
        p = ", ".join(f"{n}={v:.3g}" for n, v in zip(r["param_names"], r["fit"]["x"]))
        lines.append(f"| {r['route']} | {r['fold']} | {p} | {r['fit']['Q']:.3g} | {r['fit']['cost']:.4g} | "
                     f"{r['fit']['n_near_optimal']}/{r['fit']['n_starts']} | {', '.join(r['fit']['at_bound']) or '—'} |")
    return "\n".join(lines)


def figures(res, dt):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(FIG, exist_ok=True)
    obs = fc.data()
    col = {"R0": "#1f77b4", "R1": "#d62728", "R2": "#2ca02c", "2pop": "#9467bd"}
    # held-out predictions (frequency folds, fit B, beta 5)
    fig, axs = plt.subplots(3, 3, figsize=(13, 10), sharex=True)
    for r in res["fitB"]:
        if r["fold_kind"] != "freq" or r["beta"] != 5.0:
            continue
        for (f, a) in map(tuple, r["test"]):
            i, j = fc.FREQS.index(f), fc.AMPS.index(round(a, 3))
            ax = axs[i, j]
            key = f"f{int(f)}_A{int(round(a*1000))}"
            E = obs[(f, round(a, 3))]["E"]
            ls = {"R0": "-", "R1": "--", "R2": ":", "2pop": "-."}[r["route"]]
            for c in fc.COMPS:
                ax.plot(E, r["curves"][key][c], ls, color=col[r["route"]], lw=1.3,
                        label=r["route"] if c == "net" else None)
    for (f, a), d in obs.items():
        ax = axs[fc.FREQS.index(f), fc.AMPS.index(round(a, 3))]
        for c in fc.COMPS:
            ax.plot(d["E"], d[c], "o", ms=2.2, color="k", label="data (digitized)" if c == "net" else None)
        ax.set_title(f"{int(f)} Hz, {int(round(a*1000))} mV (held out)", fontsize=10)
        ax.invert_xaxis()
        ax.axhline(0, color="0.8", lw=0.5)
    axs[0, 0].legend(fontsize=8)
    for ax in axs[:, 0]:
        ax.set_ylabel("current (µA, cathodic +)")
    for ax in axs[-1]:
        ax.set_xlabel("E (V vs Ag|AgCl)")
    fig.suptitle("Leave-one-frequency-out predictions (fit B, β = 5 nm⁻¹): forward (+), backward (−), net", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, f"heldout_freq_dt{dt}.png"), dpi=120)
    plt.close(fig)
    # NRMSE bars
    fig, axs = plt.subplots(1, 2, figsize=(13, 4))
    for ax, fk in zip(axs, ("freq", "amp")):
        tab = curve_table(res, fk)
        routes = [r for r in ROUTE_ORDER if r in tab]
        x = np.arange(len(CONDS))
        for k, rn in enumerate(routes):
            ax.bar(x + (k - 1.5) * 0.2, [100 * tab[rn][c]["nrmse_net"] for c in CONDS], 0.2, color=col[rn], label=rn)
        ax.axhline(30, color="k", ls="--", lw=0.8)
        ax.axhline(3, color="0.5", ls=":", lw=0.8)
        ax.set_xticks(x)
        ax.set_xticklabels(CONDS, rotation=45, fontsize=8)
        ax.set_ylabel("held-out NRMSE of net (%)")
        ax.set_title({"freq": "leave-one-frequency-out (primary)", "amp": "leave-one-amplitude-out"}[fk])
    axs[0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, f"nrmse_dt{dt}.png"), dpi=120)
    plt.close(fig)


def main(dt):
    res = json.load(open(os.path.join(RES, f"cv_results_dt{dt}.json")))
    parts = []
    for fk, title in (("freq", "Primary: leave-one-frequency-out, fit B (per-route calibration), β = 5 nm⁻¹"),
                      ("amp", "Secondary: leave-one-amplitude-out, fit B, β = 5 nm⁻¹")):
        parts.append(md_table(curve_table(res, fk), title, CONDS))
    for fk in ("freq", "amp"):
        parts.append(md_table(curve_table(res, fk, kind="A"), f"Fit A (shared θ, s = 1), {fk} folds", CONDS))
    for b in (2.0, 10.0):
        parts.append(md_table(curve_table(res, "freq", beta=b), f"β = {b:g} nm⁻¹ sensitivity, frequency folds, fit B", CONDS))
    parts.append(params_table(res))
    g = [r for r in res["fitB"] if r["route"] == "R1"]
    parts.append("R1 40 vs 80 cells (same parameters): max |Δcurve|/range = " +
                 ", ".join(f"{r['fold']}(β{r['beta']:g}) {100*r['grid80_max_rel_diff']:.3f}%" for r in g))
    if res.get("profiles"):
        lines = ["| route | parameter | value | cost |", "|---|---|---|---|"]
        for p in sorted(res["profiles"], key=lambda p: (p["route"], p["value"])):
            lines.append(f"| {p['route']} | {p['param']} | {p['value']:.2f} | {p['cost']:.5g} |")
        parts.append("\n".join(lines))
    parts.append("Timescales: " + json.dumps(res["timescales"]["5.0"] if "5.0" in res["timescales"] else res["timescales"], default=float))
    md = "\n\n".join(parts)
    open(os.path.join(RES, f"tables_dt{dt}.md"), "w").write(md)
    figures(res, dt)
    print(md)


if __name__ == "__main__":
    main(float(sys.argv[1]) if len(sys.argv) > 1 else 0.5)
