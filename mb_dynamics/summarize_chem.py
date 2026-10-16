"""Round-2 tables and figures: routes within each chemistry, and chemistries compared (with parameter counts)."""
import json
import os

import numpy as np

import fit_cv as fc
from summarize import CONDS, md_table

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")
FIG = os.path.join(HERE, "figures")
ROUTES = ("R0", "R1", "R2", "2pop")
NPAR = {"C0": {"R0": 4, "R1": 5, "R2": 5, "2pop": 5}}   # incl. Q


def load():
    r1 = json.load(open(os.path.join(RES, "cv_results_dt0.5.json")))
    r2 = json.load(open(os.path.join(RES, "chem_results_dt0.5.json")))
    rows = []
    for r in r1["fitB"]:
        if r["beta"] == 5.0:
            rows.append(dict(r, chem="C0", K_p=None, fold=r["fold"], n_params=NPAR["C0"][r["route"]]))
    rows += r2["fitB"]
    return rows, r1, r2


def tab(rows, chem, kind, K_p=10.0):
    t = {}
    for r in rows:
        if r["chem"] != chem or (chem == "C2" and r["K_p"] != K_p):
            continue
        if (kind == "freq" and r["fold"].startswith("test_f")) or (kind == "amp" and r["fold"].startswith("test_A")):
            t.setdefault(r["route"], {}).update(r["test_metrics"])
    return t


def mean(t, rn, key="nrmse_net"):
    return float(np.mean([m[key] for m in t[rn].values()]))


def n_ok(t, rn):
    return sum(1 for m in t[rn].values() if m["nrmse_net"] <= 0.30 and abs(m["peak_rel_err_net"]) <= 0.30)


def chem_table(rows):
    lines = ["| chemistry | route | parameters (incl. Q) | full-data train cost | held-out NRMSE_net, freq folds | curves OK (freq) | held-out NRMSE_net, amp folds | curves OK (amp) |",
             "|---|---|---|---|---|---|---|---|"]
    for chem, kp in (("C0", None), ("C1", None), ("C2", 10.0), ("C2", 1.0), ("C2", 100.0)):
        for rn in ROUTES:
            full = [r for r in rows if r["chem"] == chem and r["route"] == rn and r["fold"] == "full"
                    and (chem != "C2" or r["K_p"] == kp)]
            if not full:
                continue
            tf, ta = tab(rows, chem, "freq", kp or 10.0), tab(rows, chem, "amp", kp or 10.0)
            label = chem + (f" (K_p={kp:g})" if chem == "C2" else "")
            fa = f"{100*mean(ta, rn):.1f}%" if rn in ta else "—"
            oa = f"{n_ok(ta, rn)}/9" if rn in ta else "—"
            lines.append(f"| {label} | {rn} | {full[0]['n_params']} | {full[0]['fit']['cost']:.4f} | {100*mean(tf, rn):.1f}% | "
                         f"{n_ok(tf, rn)}/9 | {fa} | {oa} |")
    return "\n".join(lines)


def params_table(rows):
    lines = ["| chem | route | fold | parameters | Q (C) | starts within 1 % | at bound | R1 grid 40/80 Δ |", "|---|---|---|---|---|---|---|---|"]
    for r in sorted((r for r in rows if r["chem"] != "C0" and (r["chem"] != "C2" or r["K_p"] == 10.0)),
                    key=lambda r: (r["chem"], ROUTES.index(r["route"]), r["fold"])):
        p = ", ".join(f"{n}={v:.3g}" for n, v in zip(r["param_names"], r["fit"]["x"]))
        g = (f"{100*r['grid40_max_rel_diff']:.2f}% / {100*r['grid80_max_rel_diff']:.2f}%" if "grid40_max_rel_diff" in r else "—")
        lines.append(f"| {r['chem']} | {r['route']} | {r['fold']} | {p} | {r['fit']['Q']:.3g} | "
                     f"{r['fit']['n_near_optimal']}/{r['fit']['n_starts']} | {', '.join(r['fit']['at_bound']) or '—'} | {g} |")
    return "\n".join(lines)


def fitA_table(r2, chem):
    t = {}
    for r in r2["fitA"]:
        if r["chem"] == chem and r["fold"].startswith("test_f"):
            for rn, v in r["routes"].items():
                t.setdefault(rn, {}).update(v["test_metrics"])
    return md_table(t, f"Fit A (shared θ, s = 1), {chem}, frequency folds", CONDS)


def figure(rows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    obs = fc.data()
    sel = [("C0", "R0", "#888888", "-", "C0 concerted 2e⁻ (R0)"), ("C1", "R0", "#1f77b4", "-", "C1 EE (R0)"),
           ("C2", "R0", "#d62728", "-", "C2 E–C–E (R0)"), ("C2", "R1", "#2ca02c", "--", "C2 E–C–E (R1)")]
    fig, axs = plt.subplots(3, 3, figsize=(13, 10), sharex=True)
    for chem, rn, col, ls, lab in sel:
        for r in rows:
            if r["chem"] != chem or r["route"] != rn or not r["fold"].startswith("test_f") or (chem == "C2" and r["K_p"] != 10.0):
                continue
            for (f, a) in map(tuple, r["test"]):
                ax = axs[fc.FREQS.index(f), fc.AMPS.index(round(a, 3))]
                key = f"f{int(f)}_A{int(round(a*1000))}"
                E = obs[(f, round(a, 3))]["E"]
                for c in fc.COMPS:
                    ax.plot(E, r["curves"][key][c], ls, color=col, lw=1.3, label=lab if c == "net" else None)
    for (f, a), d in obs.items():
        ax = axs[fc.FREQS.index(f), fc.AMPS.index(round(a, 3))]
        for c in fc.COMPS:
            ax.plot(d["E"], d[c], "o", ms=2.2, color="k", label="data (digitized)" if c == "net" else None)
        ax.set_title(f"{int(f)} Hz, {int(round(a*1000))} mV (held out)", fontsize=10)
        ax.invert_xaxis()
        ax.axhline(0, color="0.85", lw=0.5)
    axs[0, 0].legend(fontsize=7)
    for ax in axs[:, 0]:
        ax.set_ylabel("current (µA, cathodic +)")
    for ax in axs[-1]:
        ax.set_xlabel("E (V vs Ag|AgCl)")
    fig.suptitle("Leave-one-frequency-out predictions by chemistry (fit B, β = 5 nm⁻¹)", fontsize=11)
    fig.tight_layout()
    os.makedirs(FIG, exist_ok=True)
    fig.savefig(os.path.join(FIG, "heldout_freq_chem.png"), dpi=120)
    plt.close(fig)


def main():
    rows, r1, r2 = load()
    parts = ["## Chemistry × route summary (held-out, fit B, β = 5 nm⁻¹)", chem_table(rows)]
    for chem in ("C1", "C2"):
        for kind, title in (("freq", "leave-one-frequency-out"), ("amp", "leave-one-amplitude-out")):
            parts.append(md_table(tab(rows, chem, kind), f"{chem}: {title}, fit B", CONDS))
        parts.append(fitA_table(r2, chem))
    parts.append("## Parameters (round 2)")
    parts.append(params_table(rows))
    md = "\n\n".join(parts)
    open(os.path.join(RES, "tables_chem.md"), "w").write(md)
    figure(rows)
    print(md)


if __name__ == "__main__":
    main()
