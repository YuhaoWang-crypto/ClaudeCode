#!/usr/bin/env python3
"""Figures for the report. Every panel is generated from results/*.tsv."""
import csv
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import load_config, read_fasta, data_path, results_path, figures_path  # noqa: E402

INK = "#16202b"
MUTED = "#6b7a8c"
GRID = "#dde3ea"
ACCENT = "#2f6f9f"
WARN = "#c9752b"
BAD = "#b23b3b"
GOOD = "#3d8a6b"
ROLE_COLOR = {
    "test_article": BAD,
    "benchmark_ligand": ACCENT,
    "clinical_anchor": "#7a5ea8",
    "class_comparator": MUTED,
    "positive_control": WARN,
    "negative_control_self": GOOD,
}

plt.rcParams.update({
    "figure.dpi": 160, "savefig.dpi": 160, "savefig.bbox": "tight",
    "font.family": "DejaVu Sans", "font.size": 8.5,
    "axes.edgecolor": GRID, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.titlesize": 10, "axes.titleweight": "bold", "axes.titlecolor": INK,
    "figure.facecolor": "white", "axes.facecolor": "white",
})


def tsv(name):
    with open(results_path(name)) as f:
        return list(csv.DictReader(f, delimiter="\t"))


# --------------------------------------------------------------------------
def fig_panel_coverage():
    with open(results_path("m2_panel.json")) as f:
        d = json.load(f)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.5, 3.6),
                                   gridspec_kw={"width_ratios": [1.15, 1]})

    n = [c["n"] for c in d["curve"]]
    cov = [c["weighted_coverage"] * 100 for c in d["curve"]]
    ax1.plot(n, cov, "-o", color=ACCENT, ms=3.5, lw=1.6, zorder=3)
    ax1.axhline(d["target"] * 100, color=GOOD, ls="--", lw=1.2, zorder=2)
    ax1.axhline(d["legacy_weighted_coverage"] * 100, color=BAD, ls=":", lw=1.4, zorder=2)
    ax1.text(1.2, d["target"] * 100 + 1.4, f"target {d['target']*100:.0f}%",
             color=GOOD, fontsize=7.5)
    ax1.text(1.2, d["legacy_weighted_coverage"] * 100 - 4.2,
             f"legacy 15-molecule DR panel  {d['legacy_weighted_coverage']*100:.1f}%",
             color=BAD, fontsize=7.5)
    ax1.scatter([n[-1]], [cov[-1]], s=60, facecolor="white", edgecolor=ACCENT,
                lw=1.8, zorder=4)
    ax1.annotate(f"{d['panel_size_drb1']} DRB1 → {cov[-1]:.1f}%",
                 (n[-1], cov[-1]), textcoords="offset points", xytext=(-8, -16),
                 ha="right", fontsize=8, color=ACCENT, fontweight="bold")
    ax1.set_xlabel("DRB1 molecules in panel")
    ax1.set_ylabel("weighted US/EU phenotypic coverage (%)")
    ax1.set_title("Panel size is chosen, not assumed")
    ax1.set_ylim(0, 104)
    ax1.set_xticks(range(0, len(n) + 1, 2))
    ax1.grid(axis="y", color=GRID, lw=0.6)
    ax1.set_axisbelow(True)

    pops = [p for p in d["per_population"] if p != "World"]
    short = [p.replace("United States ", "US ") for p in pops]
    new = [d["per_population"][p] * 100 for p in pops]
    old = [d["legacy_per_population"][p] * 100 for p in pops]
    y = np.arange(len(pops))
    ax2.barh(y + 0.19, new, 0.36, color=ACCENT, label="designed panel (24 DR)")
    ax2.barh(y - 0.19, old, 0.36, color="#c2ccd6", label="legacy panel (15 DR)")
    for i, (a, b) in enumerate(zip(new, old)):
        ax2.text(a + 1, i + 0.19, f"{a:.0f}", va="center", fontsize=7, color=ACCENT)
        ax2.text(b + 1, i - 0.19, f"{b:.0f}", va="center", fontsize=7, color=MUTED)
    ax2.set_yticks(y, short, fontsize=8)
    ax2.set_xlim(0, 112)
    ax2.set_xlabel("DRB1 phenotypic coverage (%)")
    ax2.set_title("Coverage by population", pad=18)
    ax2.legend(frameon=False, fontsize=7.5, loc="lower center",
               bbox_to_anchor=(0.5, 1.0), ncol=2, handlelength=1.2)
    ax2.grid(axis="x", color=GRID, lw=0.6)
    ax2.set_axisbelow(True)
    fig.savefig(figures_path("fig1_panel_coverage.png"))
    plt.close(fig)


# --------------------------------------------------------------------------
def fig_binding_heatmap(target="AAVX_VHH"):
    cfg = load_config()
    seqs = read_fasta(data_path("sequences.fasta"))
    L = len(seqs[target])
    with open(results_path("m2_panel_alleles.txt")) as f:
        panel = [l.strip() for l in f if l.strip()]

    M = np.full((len(panel), L - 14), np.nan)
    aidx = {a: i for i, a in enumerate(panel)}
    for r in tsv("m3_binding_long.tsv"):
        if r["id"] != target or r["el_rank"] in ("", "None"):
            continue
        j = int(r["start"]) - 1
        if 0 <= j < M.shape[1]:
            M[aidx[r["allele"]], j] = float(r["el_rank"])

    # log-rank, capped at 20% so the colour range spends itself on binders
    V = np.clip(M, 0.01, 20.0)
    V = -np.log10(V)
    cmap = LinearSegmentedColormap.from_list(
        "rank", ["#f6f8fa", "#cfe0ee", "#7fb0d4", "#2f6f9f", "#1b3f5e"])

    fig, (ax, axd) = plt.subplots(2, 1, figsize=(11, 5.4), sharex=True,
                                  gridspec_kw={"height_ratios": [3.1, 1]})
    im = ax.imshow(V, aspect="auto", cmap=cmap, vmin=-1.3, vmax=2.0,
                   interpolation="nearest")
    ax.set_yticks(range(len(panel)),
                  [a.replace("HLA-DR", "DR") for a in panel], fontsize=6.2)
    ax.set_ylabel("HLA-DR molecule")
    ax.set_title(f"{target} — NetMHCIIpan EL %Rank landscape (15-mer scan)")
    cb = fig.colorbar(im, ax=ax, pad=0.012, fraction=0.022)
    cb.set_ticks([-1.3, 0, 0.7, 2.0])
    cb.set_ticklabels(["20%", "1%", "0.2%", "0.01%"], fontsize=6.5)
    cb.set_label("EL %Rank", fontsize=7)
    cb.outline.set_visible(False)

    # cluster bands
    for c in tsv("m5_clusters.tsv"):
        if c["id"] != target:
            continue
        s, e = int(c["start"]) - 1, int(c["end"]) - 1
        col = BAD if c["tolerance_class"] == "foreign" else \
              WARN if c["tolerance_class"] == "mixed" else GOOD
        ax.add_patch(plt.Rectangle((s - 0.5, -0.5), e - s + 1, len(panel),
                                   fill=False, ec=col, lw=1.5, zorder=5))
        ax.text((s + e) / 2, -1.4, f"{c['peak_core']}", color=col, fontsize=7,
                ha="center", fontweight="bold")

    # per-position promiscuity
    prom = np.nansum(M < 1.0, axis=0).astype(float)
    axd.fill_between(range(M.shape[1]), prom, color=ACCENT, alpha=0.75, lw=0)
    axd.set_ylabel("# DR with\nSB (rank<1%)", fontsize=7.5)
    axd.set_xlabel("15-mer start position in ligand (aa)")
    axd.grid(axis="y", color=GRID, lw=0.6)
    axd.set_axisbelow(True)
    fig.savefig(figures_path("fig2_binding_landscape.png"))
    plt.close(fig)


# --------------------------------------------------------------------------
def fig_ranking():
    rows = tsv("m6_calibrated_ranking.tsv")
    rows = sorted(rows, key=lambda r: float(r["pIRS"]))
    y = np.arange(len(rows))
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(11, 4.2),
                                  gridspec_kw={"width_ratios": [1.35, 1]})

    colors = [ROLE_COLOR.get(r["role"], MUTED) for r in rows]
    ax.barh(y, [float(r["pIRS"]) for r in rows], 0.62, color=colors)
    anchor = next(float(r["pIRS"]) for r in rows if r["id"] == "ProteinA_Z")
    ax.axvline(anchor, color=INK, ls="--", lw=1.1)
    ax.text(anchor, len(rows) - 0.2, " Protein A Z-domain\n benchmark",
            fontsize=7, color=INK, va="top")
    for i, r in enumerate(rows):
        ax.text(float(r["pIRS"]) + 0.06, i, f"{float(r['pIRS']):.2f}"
                + (f"   {float(r['fold_vs_ProteinA_Z']):.1f}×"
                   if r["risk_band"] != "n/a (control)" else ""),
                va="center", fontsize=7, color=INK)
    ax.set_yticks(y, [r["id"] for r in rows], fontsize=8)
    ax.set_xlabel("pIRS — population-weighted foreign DR epitope content / 100 aa")
    ax.set_title("Calibrated ranking (controls run in the same batch)")
    ax.set_xlim(0, max(float(r["pIRS"]) for r in rows) * 1.32)
    ax.grid(axis="x", color=GRID, lw=0.6)
    ax.set_axisbelow(True)

    # tolerance-filter effect
    raw = [float(r["pIRS_raw"]) for r in rows]
    filt = [float(r["pIRS"]) for r in rows]
    for i, r in enumerate(rows):
        ax2.plot([raw[i], filt[i]], [i, i], color=GRID, lw=2.2, zorder=1,
                 solid_capstyle="round")
    ax2.scatter(raw, y, s=22, color="#c2ccd6", zorder=3, label="no tolerance filter")
    ax2.scatter(filt, y, s=26, color=[ROLE_COLOR.get(r["role"], MUTED) for r in rows],
                zorder=4, label="human-self filtered")
    ax2.set_yticks(y, ["" for _ in rows])
    ax2.set_xlabel("pIRS")
    ax2.set_title("What the self/tolerance filter removes")
    ax2.legend(frameon=False, fontsize=7.5, loc="lower right")
    ax2.grid(axis="x", color=GRID, lw=0.6)
    ax2.set_axisbelow(True)
    fig.savefig(figures_path("fig3_calibrated_ranking.png"))
    plt.close(fig)


# --------------------------------------------------------------------------
def fig_tb(target="AAVX_VHH"):
    seqs = read_fasta(data_path("sequences.fasta"))
    L = len(seqs[target])
    with open(results_path("m7_bcell_per_residue.json")) as f:
        bc = json.load(f)[target]

    dens = np.zeros(L)
    for r in tsv("m5_epitopes.tsv"):
        if r["id"] != target:
            continue
        p = int(r["pos"]) - 1
        dens[max(p, 0):min(p + 9, L)] += float(r["pop_presenting"]) * float(r["tolerance_weight"])

    fig, (a1, a2) = plt.subplots(2, 1, figsize=(11, 3.9), sharex=True)
    a1.fill_between(range(1, L + 1), dens, color=BAD, alpha=0.8, lw=0)
    a1.set_ylabel("foreign DR\nepitope weight", fontsize=7.5)
    a1.set_title(f"{target} — T-helper epitope load vs B-cell epitope propensity")
    a2.fill_between(range(1, L + 1), bc, color=ACCENT, alpha=0.7, lw=0)
    a2.axhline(0.5, color=MUTED, ls="--", lw=1)
    a2.set_ylabel("BepiPred-2.0", fontsize=7.5)
    a2.set_xlabel("residue")
    for ax in (a1, a2):
        ax.grid(axis="y", color=GRID, lw=0.6)
        ax.set_axisbelow(True)

    for c in tsv("m7_tb_coincidence.tsv") if os.path.exists(results_path("m7_tb_coincidence.tsv")) else []:
        if c.get("id") != target:
            continue
        s, e = c["t_cluster"].split("-")
        for ax in (a1, a2):
            ax.axvspan(int(s), int(e), color=WARN, alpha=0.16, lw=0, zorder=0)
    fig.savefig(figures_path("fig4_tb_coincidence.png"))
    plt.close(fig)


# --------------------------------------------------------------------------
def fig_deimmunization():
    path = results_path("m9_deimmunization_scan.tsv")
    if not os.path.exists(path):
        return
    rows = tsv("m9_deimmunization_scan.tsv")
    wt = next(r for r in rows if r["variant"] == "WT")
    muts = [r for r in rows if r["variant"] != "WT"][:14][::-1]
    y = np.arange(len(muts))
    fig, ax = plt.subplots(figsize=(7.6, 4.2))
    cols = [GOOD if int(r["blosum62"]) >= 0 else WARN for r in muts]
    ax.barh(y, [float(r["pop_presenting"]) * 100 for r in muts], 0.62, color=cols)
    ax.axvline(float(wt["pop_presenting"]) * 100, color=BAD, ls="--", lw=1.2)
    ax.text(float(wt["pop_presenting"]) * 100, len(muts) - 0.3,
            f" wild type {float(wt['pop_presenting'])*100:.0f}%", color=BAD, fontsize=7.5,
            va="top")
    ax.set_yticks(y, [f"{r['variant']}  ({r['core']})" for r in muts], fontsize=7.5)
    ax.set_xlabel("% of weighted US/EU population predicted to present the epitope")
    ax.set_title("Anchor-position deimmunisation scan of the dominant epitope")
    ax.grid(axis="x", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=GOOD, label="BLOSUM62 ≥ 0 (conservative)"),
                       Patch(color=WARN, label="BLOSUM62 < 0")],
              frameon=False, fontsize=7.5, loc="lower right")
    fig.savefig(figures_path("fig5_deimmunization.png"))
    plt.close(fig)


if __name__ == "__main__":
    todo = sys.argv[1:] or ["panel", "heatmap", "ranking", "tb", "deimm"]
    fns = {"panel": fig_panel_coverage, "heatmap": fig_binding_heatmap,
           "ranking": fig_ranking, "tb": fig_tb, "deimm": fig_deimmunization}
    for t in todo:
        fns[t]()
        print(f"  {t} ok")
