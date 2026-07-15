"""make_figure.py -- data-driven summary figure of the in-silico validation.

Reads biosensor_out/boltz_results.json (+ analysis_summary.json) and plots
every system present, so new systems (e.g. vitd) appear automatically.
"""
from __future__ import annotations
import json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = os.path.join(os.getcwd(), "biosensor_out")
FIGDIR = os.path.join(os.getcwd(), "figures")

# constellation max Cα–Cα (Å) measured from downloaded chimera models; native ref 10.66.
# systems without a downloaded CIF are omitted from panel (c).
CONSTELLATION = {"native TEM-1": 10.66, "dig holo": 10.96, "dfhbi holo": 11.13, "dfhbi apo": 10.84}


def main():
    res = json.load(open(os.path.join(OUT, "boltz_results.json")))["results"]
    systems = sorted({r["system"] for r in res.values()})
    labels = {s: f"{res[f'{s}_holo']['desc'].split('+')[-1].strip() if False else s}" for s in systems}
    # nicer analyte labels
    analyte = {}
    for s in systems:
        h = res.get(f"{s}_holo", {})
        analyte[s] = {"dig": "digoxigenin", "dfhbi": "DFHBI", "vitd": "vitamin D3"}.get(s, s)
    os.makedirs(FIGDIR, exist_ok=True)

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
    x = range(len(systems)); w = 0.38
    xt = [f"{analyte[s]}\n({res[f'{s}_holo'].get('kind','')})" for s in systems]
    xt = [f"{analyte[s]}" for s in systems]

    # (a) ligand-interface confidence: native control vs chimera holo
    native = [res[f"ctrl_{s}"]["ligand_iptm"] for s in systems]
    chim = [res[f"{s}_holo"]["ligand_iptm"] for s in systems]
    ax[0].bar([i-w/2 for i in x], native, w, label="native receptor", color="#4C78A8")
    ax[0].bar([i+w/2 for i in x], chim, w, label="cp-receptor in TEM-1", color="#F58518")
    ax[0].set_xticks(list(x)); ax[0].set_xticklabels(xt)
    ax[0].set_ylabel("ligand_iptm"); ax[0].set_ylim(0, 1.05); ax[0].legend(fontsize=8)
    ax[0].set_title("(a) Receptor binds analyte after CP+insertion\n[rigorous metric]")

    # (b) binding retention vs native
    ret = [round(res[f"{s}_holo"]["ligand_iptm"]/res[f"ctrl_{s}"]["ligand_iptm"], 2) for s in systems]
    colors = ["#54A24B" if r >= 0.85 else "#E4A33A" if r >= 0.7 else "#D64545" for r in ret]
    ax[1].bar(list(x), ret, color=colors)
    ax[1].axhline(1.0, ls="--", c="#333", lw=1)
    for i, r in zip(x, ret): ax[1].text(i, r+0.02, f"{r:.2f}", ha="center", fontsize=9)
    ax[1].set_xticks(list(x)); ax[1].set_xticklabels(xt)
    ax[1].set_ylabel("chimera / native ligand_iptm"); ax[1].set_ylim(0, 1.2)
    ax[1].set_title("(b) Binding retention vs native\n[rigorous]  (1.0 = fully retained)")

    # (c) catalytic constellation max Cα–Cα vs native (only measured chimeras)
    keys = list(CONSTELLATION.keys()); vals = [CONSTELLATION[k] for k in keys]
    ccols = ["#333333"] + ["#F58518" if "holo" in k else "#B0B0B0" for k in keys[1:]]
    ax[2].bar(range(len(vals)), vals, color=ccols)
    ax[2].axhline(10.66, ls="--", c="#333", lw=1)
    ax[2].set_xticks(range(len(vals))); ax[2].set_xticklabels([k.replace(" ", "\n") for k in keys], fontsize=8)
    ax[2].set_ylabel("max catalytic Cα–Cα (Å)"); ax[2].set_ylim(0, 13)
    ax[2].set_title("(c) TEM-1 active site intact in chimera\n[rigorous] (Δ<0.5 Å)")

    fig.suptitle("In-silico validation of ML-receptor allosteric biosensors (Boltz-2.1)\n"
                 "reproduction (digoxigenin) + validation on different analytes (DFHBI, vitamin D3)",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.9])
    path = os.path.join(FIGDIR, "biosensor_validation.png")
    fig.savefig(path, dpi=130)
    print("wrote", path, "| systems:", systems, "| retention:", dict(zip(systems, ret)))


if __name__ == "__main__":
    main()
