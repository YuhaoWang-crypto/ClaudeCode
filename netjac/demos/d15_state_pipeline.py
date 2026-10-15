"""
D15  The state-restoration pipeline (netjac/statepipe.py) run end to end on
     three real intervention datasets, with pre-specified checks:

  AD   GSE130588 dupilumab vs placebo   (positive control: the fitted direct
       input b must point at the Th2 module, the drug's mechanism; the
       placebo arm must stay at the disease steady state)
  PSO  GSE85034 adalimumab / methotrexate (out-of-sample check: data contain
       no IL-17 drug, yet IL-17 blockade is the most effective psoriasis
       therapy; is the IL-17 keratinocyte module the best predicted lever?)
  IBD  GSE16879 infliximab, colon       (two time points only; reference =
       control mucosa, no personal reference)

Each prints a 'state card': steady state, order parameters, core network,
dynamics, minimal adjustment. AD is also run with data-driven modules
(no curated list) to test transfer to an indication without prior modules.

Run:  python3 -m netjac.demos.d15_state_pipeline
"""
import os
import pickle

import numpy as np

from netjac import DATA_CACHE, FIG
from netjac import statepipe as sp
from netjac.data import geo
from netjac.demos import d14_treatment_restoration as d14

RES = os.path.join(DATA_CACHE, "d15_results.pkl")
SKIN = {k: v for k, v in d14.SKIN.items() if not k.startswith("housekeeping")}
GUT = {k: v for k, v in d14.GUT.items() if not k.startswith("housekeeping")}


def data_ad():
    G, M = d14.load_ad()
    M = M[(M.tissue == "LS") | ((M.tissue == "NL") & (M.time == "Week0"))]
    tmap = {"Week0": 0.0, "Week4": 4.0, "Week16": 16.0}
    easi = {}
    for s, m in M[M.tissue == "LS"].groupby("subject id"):
        e = dict(zip(m.time, m.easi))
        if "Week0" in e and "Week16" in e and e["Week16"] not in ("NA", ""):
            easi[s] = float(e["Week16"]) <= 0.25 * float(e["Week0"])
    D = sp.StateData(X=G[M.index].T.values, features=np.array(G.index), subject=M["subject id"].values,
                     kind=np.where(M.tissue == "NL", "ref", "state"), time=M.time.map(tmap).values.astype(float),
                     arm=dict(zip(M["subject id"], M.treatment)), outcome=easi)
    return D


def data_pso():
    P = d14.load_psoriasis()
    X, subj, kind, time, arm, out, grp = [], [], [], [], {}, {}, {}
    tmap = {"NL": None, "LS": 0.0, "WK 1": 1.0, "WK 2": 2.0, "WK 4": 4.0, "WK 16": 16.0}
    for p in P["pts"]:
        arm[p["id"]] = p["drug"]
        grp[p["id"]] = p["drug"]
        if "WK 16" in p["pasi"] and "LS" in p["pasi"]:
            out[p["id"]] = p["pasi"]["WK 16"] <= 0.25 * p["pasi"]["LS"]
        for k, x in p["x"].items():
            if k not in tmap:
                continue
            X.append(x)
            subj.append(p["id"])
            kind.append("ref" if k == "NL" else "state")
            time.append(np.nan if k == "NL" else tmap[k])
    return sp.StateData(X=np.array(X), features=P["genes"], subject=np.array(subj), kind=np.array(kind),
                        time=np.array(time), arm=arm, outcome=out, group=grp)


def data_ibd():
    X, M = geo.series_matrix("GSE16879")
    G = d14._genes(X, "GPL570", log=True)
    M = M.assign(pid=M.title.str.replace(r"_(before|after)T$", "", regex=True),
                 when=M.title.str.extract(r"_(before|after)T$")[0])
    M = M[(M.tissue == "Colon")]
    ctrl = M.title.str.match(r"^CO\d+$")
    kind = np.where(ctrl, "ref", "state")
    time = np.where(M.when == "after", 5.0, 0.0)
    subj = np.where(ctrl, M.title, M.pid)
    keep = ctrl | M.pid.duplicated(keep=False)
    resp = {p: v == "Yes" for p, v in zip(M.pid, M["response to infliximab"]) if v in ("Yes", "No")}
    grp = {p: p[:2] if p.startswith("UC") else "CD" for p in resp}
    return sp.StateData(X=G[M.index[keep]].T.values, features=np.array(G.index), subject=subj[keep],
                        kind=kind[keep], time=time[keep], arm={p: "infliximab" for p in resp},
                        outcome=resp, group=grp)


def run():
    R = {}
    D = data_ad()
    R["AD"] = sp.run(D, SKIN)
    S = sp.steady_state(D)
    auto = sp.auto_modules(D, S, k=6)
    R["AD auto"] = sp.run(D, auto)
    R["AD auto modules"] = {k: v[:12] for k, v in auto.items()}
    # where do the curated modules fall among the auto modules
    R["AD auto overlap"] = {a: {m: len(set(g) & set(SKIN[m])) for m in SKIN if set(g) & set(SKIN[m])}
                            for a, g in auto.items()}
    R["PSO"] = sp.run(data_pso(), SKIN)
    R["IBD"] = sp.run(data_ibd(), GUT)
    pickle.dump(R, open(RES, "wb"))
    return R


def report(R):
    for k in ("AD", "AD auto", "PSO", "IBD"):
        sp.print_card(R[k], k)
        if k == "AD auto":
            print("  auto modules (first genes):")
            for a, g in R["AD auto modules"].items():
                print(f"    {a}: {g} | curated overlap {R['AD auto overlap'][a]}")


def figures(R):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(FIG, exist_ok=True)
    fig, ax = plt.subplots(2, 3, figsize=(17, 9.5))
    for row, key in enumerate(("AD", "PSO")):
        c = R[key]
        d, a = c["dyn"], c["adjust"]
        n = d["names"]
        short = [x.split(" (")[0].replace("IL-17 keratinocyte response", "IL-17 keratinocyte") for x in n]
        ax0 = ax[row, 0]
        im = ax0.imshow(d["J"], cmap="RdBu_r", vmin=-abs(d["J"]).max(), vmax=abs(d["J"]).max())
        ax0.set_xticks(range(len(n)), short, rotation=45, ha="right", fontsize=7)
        ax0.set_yticks(range(len(n)), short, fontsize=7)
        ax0.set_title(f"{key}: module Jacobian J (per week)\nLOSO RMSE {d['cv_rmse']:.2f} vs no-change {d['rmse_no_change']:.2f}",
                      fontsize=10)
        fig.colorbar(im, ax=ax0, fraction=0.046)
        ax1 = ax[row, 1]
        o = np.argsort(d["b"])
        ax1.barh(range(len(n)), d["b"][o], color=["C3" if i == o[0] else "0.6" for i in o])
        ax1.set_yticks(range(len(n)), [short[i] for i in o], fontsize=8)
        ax1.invert_yaxis()
        ax1.axvline(0, color="k", lw=0.6)
        ax1.set_xlabel("direct input b of the treatment (per week)")
        ax1.set_title(f"{key}: where the intervention acts directly", fontsize=10)
        ax2 = ax[row, 2]
        sing = sorted(a["singles"].items(), key=lambda kv: -kv[1]["restoration"])
        ax2.barh(range(len(sing)), [v["restoration"] for _, v in sing], color="C2")
        ax2.set_yticks(range(len(sing)), [k.split(" (")[0] for k, _ in sing], fontsize=8)
        ax2.invert_yaxis()
        ax2.set_xlim(0, 1)
        nset = len(a["path"][-1]["set"])
        ax2.set_xlabel("predicted restoration of the whole state")
        flag = "" if a["reliable"] else " - UNRELIABLE (steady state not identifiable)"
        ax2.set_title(f"{key}: single-module levers = bottlenecks, not causal drivers{flag}\n"
                      f"minimal set to 90%: {nset} of {len(a['singles'])} modules", fontsize=9)
    fig.tight_layout()
    p = os.path.join(FIG, "d15_state_pipeline.png")
    fig.savefig(p, dpi=110)
    plt.close(fig)
    return p


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--report", action="store_true")
    a = ap.parse_args()
    R = pickle.load(open(RES, "rb")) if a.report and os.path.exists(RES) else run()
    report(R)
    print(figures(R))
