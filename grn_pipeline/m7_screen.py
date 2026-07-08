"""
Module 7 — Boltz-2.1 small-molecule library screen of KRAS G12C analogues.

Ranks 10 covalent KRAS G12C ligands (the 2 approved drugs + 8 sotorasib
analogues pulled from ChEMBL by similarity) by predicted binding, using a
real Boltz-2.1 library screen run this session against the KRAS G12C
sequence (job sm_scr_6ZFAIXYLHhbaJhaqegsE). Boltz's default reactive-warhead
SMARTS filter had to be DISABLED, because the acrylamide covalent warhead
common to all G12C binders trips a Michael-acceptor structural alert.

Each hit also carries Boltz ADME (solubility / permeability / lipophilicity).
The top binders feed directly into M6/M9 as stronger-engagement candidates.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Real Boltz-2.1 screen results (this session). bind_conf = binding
# confidence; opt_score = optimization score; permeability/lipophilicity
# are Boltz ADME predictions.
SCREEN = {
    "adagrasib_CHEMBL4594350": {"bind_conf": 0.952, "opt_score": 0.593, "iptm": 0.942, "struct_conf": 0.872, "solubility": "medium-confidence", "permeability": 0.509, "lipophilicity": 2.697, "approved": True},
    "CHEMBL4851723":           {"bind_conf": 0.929, "opt_score": 0.465, "iptm": 0.939, "struct_conf": 0.863, "solubility": "medium-confidence", "permeability": 0.447, "lipophilicity": 2.033, "approved": False},
    "CHEMBL5918612":           {"bind_conf": 0.927, "opt_score": 0.420, "iptm": 0.927, "struct_conf": 0.861, "solubility": "medium-confidence", "permeability": 0.421, "lipophilicity": 2.055, "approved": False},
    "CHEMBL4441771":           {"bind_conf": 0.918, "opt_score": 0.426, "iptm": 0.926, "struct_conf": 0.856, "solubility": "medium-confidence", "permeability": 0.575, "lipophilicity": 1.527, "approved": False},
    "CHEMBL4539214":           {"bind_conf": 0.909, "opt_score": 0.463, "iptm": 0.935, "struct_conf": 0.856, "solubility": "medium-confidence", "permeability": 0.579, "lipophilicity": 2.529, "approved": False},
    "sotorasib_CHEMBL4535757": {"bind_conf": 0.878, "opt_score": 0.419, "iptm": 0.922, "struct_conf": 0.851, "solubility": "medium-confidence", "permeability": 0.473, "lipophilicity": 2.096, "approved": True},
    "CHEMBL4452622":           {"bind_conf": 0.876, "opt_score": 0.436, "iptm": 0.938, "struct_conf": 0.853, "solubility": "medium-confidence", "permeability": 0.612, "lipophilicity": 2.300, "approved": False},
    "CHEMBL4754887":           {"bind_conf": 0.855, "opt_score": 0.370, "iptm": 0.922, "struct_conf": 0.850, "solubility": "medium-confidence", "permeability": 0.572, "lipophilicity": 2.073, "approved": False},
    "CHEMBL4458041":           {"bind_conf": 0.835, "opt_score": 0.387, "iptm": 0.938, "struct_conf": 0.854, "solubility": "medium-confidence", "permeability": 0.514, "lipophilicity": 1.954, "approved": False},
    "CHEMBL4744465":           {"bind_conf": 0.794, "opt_score": 0.316, "iptm": 0.920, "struct_conf": 0.842, "solubility": "medium-confidence", "permeability": 0.652, "lipophilicity": 1.665, "approved": False},
}
JOB_ID = "sm_scr_6ZFAIXYLHhbaJhaqegsE"


def ranked():
    return sorted(SCREEN.items(), key=lambda kv: -kv[1]["bind_conf"])


def report(fig_path="figures/m7_screen.png"):
    print("=" * 68)
    print("MODULE 7 — Boltz-2.1 small-molecule screen (KRAS G12C analogues)")
    print("=" * 68)
    print(f"job {JOB_ID}: 10 ligands screened (warhead SMARTS filter disabled)")
    print(f"\n{'rank':>4} {'molecule':26s} {'bind_conf':>9} {'opt':>6} "
          f"{'perm':>6} {'logP':>6}  approved")
    for i, (name, m) in enumerate(ranked(), 1):
        print(f"{i:>4} {name:26s} {m['bind_conf']:9.3f} "
              f"{m['opt_score']:6.3f} {m['permeability']:6.3f} "
              f"{m['lipophilicity']:6.3f}  {'YES' if m['approved'] else ''}")

    top = ranked()[0]
    soto = SCREEN["sotorasib_CHEMBL4535757"]["bind_conf"]
    better = [n for n, m in ranked()
              if not m["approved"] and m["bind_conf"] > soto]
    print(f"\ntop binder: {top[0]} (bind_conf={top[1]['bind_conf']:.3f})")
    print(f"analogues out-ranking sotorasib ({soto:.3f}) in binding_confidence:")
    for n in better:
        print(f"    {n}  bind_conf={SCREEN[n]['bind_conf']:.3f}")
    print("  -> these are the stronger-engagement candidates that would enter")
    print("     M6/M9 next (higher engagement E -> larger network perturbation).")
    _figure(fig_path)
    print(f"\nfigure written to: {fig_path}")
    return {"ranked": ranked(), "better_than_sotorasib": better}


def _figure(path):
    r = ranked()
    names = [n.replace("_CHEMBL", "\nCHEMBL") for n, _ in r]
    conf = [m["bind_conf"] for _, m in r]
    cols = ["#e67e22" if m["approved"] else "#3498db" for _, m in r]
    fig, ax = plt.subplots(figsize=(11, 5))
    bars = ax.bar(range(len(r)), conf, color=cols, edgecolor="k")
    ax.set_xticks(range(len(r)))
    ax.set_xticklabels(names, fontsize=7, rotation=30, ha="right")
    ax.set_ylabel("Boltz-2.1 binding_confidence")
    ax.set_ylim(0.75, 0.97)
    ax.set_title("KRAS G12C library screen — orange = approved drug, "
                 "blue = ChEMBL analogue")
    for b, (_, m) in zip(bars, r):
        ax.text(b.get_x() + b.get_width() / 2, b.get_height() + 0.002,
                f"{m['bind_conf']:.3f}", ha="center", fontsize=7)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    report()
