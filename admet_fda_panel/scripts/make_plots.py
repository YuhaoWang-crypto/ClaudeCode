"""Figures for the ADMET / safety-panel report."""

from __future__ import annotations
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

OUT = "/home/user/results/admet_fda_panel"
FIG = f"{OUT}/figures"
CSV = f"{OUT}/all_properties.csv"

from matplotlib import font_manager as _fm
for _p in ("/usr/share/fonts/truetype/wqy/wqy-zenhei.ttc",):
    try:
        _fm.fontManager.addfont(_p)
    except Exception:
        pass

plt.rcParams.update({
    "font.sans-serif": ["WenQuanYi Zen Hei", "DejaVu Sans"],
    "axes.unicode_minus": False,
    "figure.dpi": 160, "savefig.dpi": 160, "font.size": 9,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.6,
    "axes.axisbelow": True, "figure.facecolor": "white",
    "axes.labelsize": 9.5, "axes.titlesize": 10.5, "legend.frameon": False,
})

INK = "#1f2933"
BLUE = "#2f6f9f"
ORANGE = "#d9822b"
RED = "#c0392b"
GREEN = "#3d8b62"
GREY = "#98a2ad"


def save(fig, stem):
    os.makedirs(FIG, exist_ok=True)
    for ext in ("png", "svg"):
        fig.savefig(f"{FIG}/{stem}.{ext}", bbox_inches="tight",
                    facecolor="white")
    plt.close(fig)
    print(f"  wrote {stem}.png / .svg")


# --------------------------------------------------------------------------- #
def fig1_physchem(d):
    fig, axes = plt.subplots(2, 3, figsize=(10.5, 6))
    specs = [
        ("MW", "Molecular weight (Da)", [500], BLUE),
        ("cLogP", "cLogP (Crippen)", [5], BLUE),
        ("TPSA", "TPSA (Å²)", [140], BLUE),
        ("HBD", "H-bond donors", [5], BLUE),
        ("RotB", "Rotatable bonds", [10], BLUE),
        ("QED", "QED drug-likeness", [0.5], GREEN),
    ]
    for ax, (col, label, cuts, colour) in zip(axes.ravel(), specs):
        v = d[col].values
        bins = min(14, max(8, int(np.sqrt(len(v)) * 2.2)))
        ax.hist(v, bins=bins, color=colour, alpha=0.78, edgecolor="white",
                linewidth=0.8)
        for c in cuts:
            ax.axvline(c, color=RED, ls="--", lw=1.1)
            ax.text(c, ax.get_ylim()[1] * 0.93, f" {c:g}", color=RED,
                    fontsize=7.5, va="top")
        ax.axvline(np.median(v), color=INK, lw=1.3)
        ax.set_xlabel(label)
        ax.set_ylabel("Compounds")
        ax.set_title(f"median {np.median(v):.1f}", fontsize=8.5, color=GREY,
                     loc="right")
    fig.suptitle("Physicochemical property distribution — 30 FDA-approved drugs",
                 fontsize=11.5, y=0.99)
    fig.text(0.5, -0.015, "Dashed red = rule-of-five / Veber threshold;  "
             "solid dark = panel median", ha="center", fontsize=7.5, color=GREY)
    fig.tight_layout()
    save(fig, "fig1_physicochemical_overview")


def fig2_lipinski(d):
    fig, ax = plt.subplots(figsize=(8.2, 6))
    ax.add_patch(Rectangle((0, -2), 500, 7, facecolor=GREEN, alpha=0.07,
                           zorder=0, lw=0))
    ax.axvline(500, color=RED, ls="--", lw=1.1)
    ax.axhline(5, color=RED, ls="--", lw=1.1)

    viol = d["Lipinski_violations"].values
    colours = [GREEN if v == 0 else ORANGE if v == 1 else RED for v in viol]
    sizes = 40 + d["QED"].values * 260
    ax.scatter(d["MW"], d["cLogP"], s=sizes, c=colours, alpha=0.8,
               edgecolor="white", linewidth=1.0, zorder=3)

    label_these = set(d.nlargest(4, "MW")["name"]) | set(d.nlargest(3, "cLogP")["name"]) \
        | set(d.nsmallest(2, "MW")["name"]) | set(d.nsmallest(2, "cLogP")["name"])
    for _, r in d.iterrows():
        if r["name"] in label_these:
            ax.annotate(r["name"], (r["MW"], r["cLogP"]),
                        textcoords="offset points", xytext=(6, 4),
                        fontsize=7.2, color=INK)
    ax.set_xlabel("Molecular weight (Da)")
    ax.set_ylabel("cLogP")
    ax.set_title("Lipinski chemical space — bubble area ∝ QED, colour = Ro5 violations")
    handles = [plt.Line2D([], [], marker="o", ls="", color=c, markersize=8,
                          markeredgecolor="white", label=l)
               for c, l in [(GREEN, "0 violations"), (ORANGE, "1 violation"),
                            (RED, "≥2 violations")]]
    ax.legend(handles=handles, loc="upper left", fontsize=8)
    ax.text(0.99, 0.02, "shaded = rule-of-five compliant quadrant",
            transform=ax.transAxes, ha="right", fontsize=7.5, color=GREY)
    fig.tight_layout()
    save(fig, "fig2_lipinski_space")


def fig3_qed(d):
    s = d.sort_values("QED")
    fig, ax = plt.subplots(figsize=(7.6, 8.4))
    colours = [GREEN if q >= 0.67 else BLUE if q >= 0.5 else
               ORANGE if q >= 0.34 else RED for q in s["QED"]]
    ax.barh(s["name"], s["QED"], color=colours, alpha=0.88, height=0.72)
    for y, (q, a) in enumerate(zip(s["QED"], s["alert_total"])):
        ax.text(q + 0.012, y, f"{q:.2f}" + (f"  ⚑{a}" if a else ""),
                va="center", fontsize=7.3, color=INK)
    for x, lab in [(0.34, "poor"), (0.5, "moderate"), (0.67, "attractive")]:
        ax.axvline(x, color=GREY, ls=":", lw=0.9)
        ax.text(x, len(s) - 0.2, f" {lab}", fontsize=7, color=GREY, rotation=90,
                va="top")
    ax.set_xlim(0, 1.03)
    ax.set_xlabel("QED (quantitative estimate of drug-likeness)")
    ax.set_title("Developability ranking by QED\n⚑ = number of structural alerts",
                 fontsize=10.5)
    ax.grid(axis="y", visible=False)
    fig.tight_layout()
    save(fig, "fig3_qed_developability")


HEAT_ENDPOINTS = [
    ("hERG", "hERG blockade"),
    ("DILI", "Drug-induced liver injury"),
    ("AMES", "Ames mutagenicity"),
    ("ClinTox", "Clinical-trial toxicity"),
    ("Carcinogens_Lagunin", "Carcinogenicity"),
    ("Skin_Reaction", "Skin reaction"),
    ("CYP1A2_Veith", "CYP1A2 inhibition"),
    ("CYP2C9_Veith", "CYP2C9 inhibition"),
    ("CYP2C19_Veith", "CYP2C19 inhibition"),
    ("CYP2D6_Veith", "CYP2D6 inhibition"),
    ("CYP3A4_Veith", "CYP3A4 inhibition"),
    ("Pgp_Broccatelli", "P-gp inhibition"),
    ("BBB_Martins", "BBB penetration"),
    ("HIA_Hou", "Human intestinal absorption"),
    ("Bioavailability_Ma", "Oral bioavailability"),
    ("NR-AR", "NR: androgen receptor"),
    ("NR-ER", "NR: estrogen receptor"),
    ("NR-Aromatase", "NR: aromatase"),
    ("NR-PPAR-gamma", "NR: PPAR-γ"),
    ("NR-AhR", "NR: aryl-hydrocarbon receptor"),
    ("SR-ARE", "SR: oxidative stress (ARE)"),
    ("SR-MMP", "SR: mitochondrial toxicity"),
    ("SR-p53", "SR: p53 DNA damage"),
]


def fig4_heatmap(d):
    cols = [c for c, _ in HEAT_ENDPOINTS]
    labels = [l for _, l in HEAT_ENDPOINTS]
    m = d.set_index("name")[cols]
    m = m.loc[m.mean(axis=1).sort_values(ascending=False).index]

    fig, ax = plt.subplots(figsize=(11.5, 8.2))
    im = ax.imshow(m.values, aspect="auto", cmap="RdYlBu_r", vmin=0, vmax=1)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=55, ha="right", fontsize=7.6)
    ax.set_yticks(range(len(m)))
    ax.set_yticklabels(m.index, fontsize=7.6)
    ax.grid(False)
    for i in range(len(m)):
        for j in range(len(cols)):
            v = m.values[i, j]
            if v >= 0.85 or v <= 0.06:
                ax.text(j, i, f"{v:.2f}".lstrip("0"), ha="center", va="center",
                        fontsize=5.4, color="white" if v >= 0.85 else INK)
    # separate the toxicity / DMPK / Tox21 blocks
    for x in (5.5, 14.5):
        ax.axvline(x, color=INK, lw=1.4)
    cb = fig.colorbar(im, ax=ax, fraction=0.018, pad=0.012)
    cb.set_label("Predicted probability", fontsize=8.5)
    ax.set_title("ADMET-AI endpoint heatmap — rows sorted by mean predicted risk\n"
                 "blocks: toxicity | DMPK & absorption | Tox21 nuclear-receptor "
                 "& stress-response", fontsize=10.5)
    fig.tight_layout()
    save(fig, "fig4_admet_heatmap")


def fig5_validation():
    v = pd.read_csv(f"{OUT}/demo_validation.csv")
    pos = v[v.role == "positive control"]
    neg = v[v.role == "negative control"]
    order = ["hERG", "CYP3A4_Veith", "CYP2D6_Veith", "CYP1A2_Veith",
             "Pgp_Broccatelli", "BBB_Martins", "NR-Aromatase", "NR-ER",
             "NR-PPAR-gamma", "NR-AR"]
    fig, ax = plt.subplots(figsize=(9.6, 6.2))
    for i, ep in enumerate(order):
        p = pos[pos.endpoint == ep]["score"].values
        n = neg[neg.endpoint == ep]["score"].values
        ax.scatter(n, [i] * len(n), s=52, c=GREY, marker="o", alpha=0.75,
                   edgecolor="white", zorder=3, label="negative control" if i == 0 else None)
        ok = p > 0.5
        ax.scatter(p[ok], [i] * ok.sum(), s=95, c=GREEN, marker="D",
                   edgecolor="white", zorder=4,
                   label="positive control recovered" if i == 0 else None)
        ax.scatter(p[~ok], [i] * (~ok).sum(), s=95, c=RED, marker="X",
                   edgecolor="white", zorder=4,
                   label="positive control missed" if i == 0 else None)
        ax.axhline(i, color=GREY, lw=0.5, alpha=0.4, zorder=1)
    ax.axvline(0.5, color=INK, ls="--", lw=1.1)
    ax.text(0.5, len(order) - 0.35, " decision threshold 0.5", fontsize=7.5,
            color=INK)
    ax.axhspan(5.5, len(order) - 0.5, color=RED, alpha=0.05, zorder=0)
    ax.text(0.98, 8.0, "Tox21 nuclear-receptor endpoints\n"
            "fail to recover approved pharmacology", ha="right", fontsize=8,
            color=RED)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(order, fontsize=8.5)
    ax.set_xlim(-0.03, 1.03)
    ax.set_xlabel("Predicted probability")
    ax.set_title("Demo validation — known positive and negative control drugs\n"
                 "per ADMET-AI endpoint", fontsize=10.5)
    ax.legend(loc="lower right", fontsize=8)
    ax.grid(axis="y", visible=False)
    fig.tight_layout()
    save(fig, "fig5_demo_validation")


def fig6_coverage():
    """Safety-panel computational coverage, by target family."""
    from safety_panel import PANEL_FAMILIES
    fam = [f["family"] for f in PANEL_FAMILIES]
    b44 = [f["bowes44"] for f in PANEL_FAMILIES]
    cov = [f["covered"] for f in PANEL_FAMILIES]
    val = [f["validated"] for f in PANEL_FAMILIES]
    y = np.arange(len(fam))
    fig, ax = plt.subplots(figsize=(9.0, 5.0))
    ax.barh(y, b44, color=GREY, alpha=0.45, height=0.64,
            label="面板靶点数")
    ax.barh(y, cov, color=ORANGE, height=0.64, label="存在名义对应的计算终点")
    ax.barh(y, val, color=GREEN, height=0.64, label="且通过阳性对照验证")
    for i, (a, c, v) in enumerate(zip(b44, cov, val)):
        note = f"{c}/{a} 可映射"
        if c and not v:
            note += "（验证失败）"
        elif v:
            note += "，已验证"
        ax.text(a + 0.35, i, note, va="center", fontsize=8, color=INK)
    ax.set_yticks(y)
    ax.set_yticklabels(fam, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlim(0, 30)
    ax.set_xlabel("靶点数")
    ax.set_title("Bowes-44 安全药理面板的计算覆盖度（按靶点家族）\n"
                 f"44 个靶点中仅 {sum(cov)} 个存在对应计算终点，"
                 f"其中仅 {sum(val)} 个（hERG）通过验证", fontsize=10.5)
    ax.legend(fontsize=8.5, loc="lower right")
    ax.grid(axis="y", visible=False)
    fig.tight_layout()
    save(fig, "fig6_panel_coverage")


def generate_all_plots():
    d = pd.read_csv(CSV)
    print("Generating figures:")
    fig1_physchem(d)
    fig2_lipinski(d)
    fig3_qed(d)
    fig4_heatmap(d)
    fig5_validation()
    fig6_coverage()
    print("✓ Figures complete")


if __name__ == "__main__":
    import sys
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    generate_all_plots()
