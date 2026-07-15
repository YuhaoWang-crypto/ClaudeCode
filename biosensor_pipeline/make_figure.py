"""make_figure.py -- summary figure of the in-silico validation."""
from __future__ import annotations
import json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = os.path.join(os.getcwd(), "biosensor_out")
FIGDIR = os.path.join(os.getcwd(), "figures")


def main():
    with open(os.path.join(OUT, "boltz_results.json")) as fh:
        res = json.load(fh)["results"]
    os.makedirs(FIGDIR, exist_ok=True)

    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))

    # (a) ligand-interface confidence: native control vs chimera holo
    systems = [("dig", "digoxigenin\n(reproduction)"), ("dfhbi", "DFHBI\n(validation)")]
    x = range(len(systems))
    native = [res[f"ctrl_{k}"]["ligand_iptm"] for k, _ in systems]
    chim = [res[f"{k}_holo"]["ligand_iptm"] for k, _ in systems]
    w = 0.35
    ax[0].bar([i - w/2 for i in x], native, w, label="native receptor", color="#4C78A8")
    ax[0].bar([i + w/2 for i in x], chim, w, label="cp-receptor in TEM-1", color="#F58518")
    ax[0].set_xticks(list(x)); ax[0].set_xticklabels([s for _, s in systems])
    ax[0].set_ylabel("ligand_iptm (interface confidence)")
    ax[0].set_title("(a) Receptor still binds analyte\nafter CP + insertion  ✅")
    ax[0].set_ylim(0, 1.05); ax[0].legend(fontsize=8)

    # (b) apo vs holo complex pLDDT
    apo = [res[f"{k}_apo"]["complex_plddt"] for k, _ in systems]
    holo = [res[f"{k}_holo"]["complex_plddt"] for k, _ in systems]
    ax[1].bar([i - w/2 for i in x], apo, w, label="apo (OFF scaffold)", color="#B0B0B0")
    ax[1].bar([i + w/2 for i in x], holo, w, label="holo (+ analyte)", color="#54A24B")
    ax[1].set_xticks(list(x)); ax[1].set_xticklabels([s for _, s in systems])
    ax[1].set_ylabel("complex pLDDT")
    ax[1].set_title("(b) Fold confidence apo vs holo\n[interpretive]")
    ax[1].set_ylim(0, 1.05); ax[1].legend(fontsize=8)

    # (c) catalytic constellation max Cα–Cα vs native
    labels = ["native\nTEM-1", "dig\nholo", "dfhbi\nholo", "dfhbi\napo"]
    vals = [10.66, 10.96, 11.13, 10.84]
    colors = ["#333333", "#F58518", "#F58518", "#B0B0B0"]
    ax[2].bar(range(len(vals)), vals, color=colors)
    ax[2].axhline(10.66, ls="--", c="#333333", lw=1)
    ax[2].set_xticks(range(len(vals))); ax[2].set_xticklabels(labels)
    ax[2].set_ylabel("max catalytic Cα–Cα (Å)")
    ax[2].set_title("(c) TEM-1 active site intact\nin chimera [rigorous] (<0.5A)")
    ax[2].set_ylim(0, 13)

    fig.suptitle("In-silico validation of ML-receptor allosteric biosensors "
                 "(Boltz-2.1) — reproduction (digoxigenin) + different analyte (DFHBI)",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    path = os.path.join(FIGDIR, "biosensor_validation.png")
    fig.savefig(path, dpi=130)
    print("wrote", path)


if __name__ == "__main__":
    main()
