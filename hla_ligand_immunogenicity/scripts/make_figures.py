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
    p = results_path(name)
    if not os.path.exists(p):
        return []
    with open(p) as f:
        if f.readline().startswith("(none)"):
            return []
        f.seek(0)
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
    ax2.barh(y + 0.19, new, 0.36, color=ACCENT, label=f"designed panel ({d['panel_size_total']} DR)")
    ax2.barh(y - 0.19, old, 0.36, color="#c2ccd6", label=f"legacy panel ({len(d['legacy_drb1_panel'])+4} DR)")
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
    C = np.zeros((len(panel), L - 14))          # consensus strong-binder calls
    aidx = {a: i for i, a in enumerate(panel)}
    for r in tsv("m3_binding_long.tsv"):
        if r["id"] != target or r["el_rank"] in ("", "None"):
            continue
        j = int(r["start"]) - 1
        if 0 <= j < M.shape[1]:
            M[aidx[r["allele"]], j] = float(r["el_rank"])
            C[aidx[r["allele"]], j] = r["call_consensus"] == "SB"

    # 15-mer frame that each epitope core was called from, so cluster boxes land
    # on the same x-axis as the heatmap (15-mer start, not core position)
    frame = {}
    for e in tsv("m5_epitopes.tsv"):
        if e["id"] != target:
            continue
        frame[e["core"]] = int(e["pos"]) - e["peptide"].find(e["core"])

    # log-rank, capped at 20% so the colour range spends itself on binders
    V = np.clip(M, 0.01, 20.0)
    V = -np.log10(V)
    cmap = LinearSegmentedColormap.from_list(
        "rank", ["#f6f8fa", "#cfe0ee", "#7fb0d4", "#2f6f9f", "#1b3f5e"])

    fig, (ax, axd) = plt.subplots(2, 1, figsize=(11, 5.6), sharex=True,
                                  gridspec_kw={"height_ratios": [3.1, 1]})
    im = ax.imshow(V, aspect="auto", cmap=cmap, vmin=-0.7, vmax=2.0,
                   interpolation="nearest")
    ax.set_yticks(range(len(panel)),
                  [a.replace("HLA-DR", "DR") for a in panel], fontsize=6.2)
    ax.set_ylabel("HLA-DR molecule")
    ax.set_title(f"{target} — NetMHCIIpan EL %Rank landscape (15-mer scan)", pad=22)
    cb = fig.colorbar(im, ax=ax, pad=0.012, fraction=0.022)
    cb.set_ticks([-0.7, 0, 0.7, 2.0])
    cb.set_ticklabels(["5%", "1%", "0.2%", "0.01%"], fontsize=6.5)
    cb.set_label("EL %Rank", fontsize=7)
    cb.outline.set_visible(False)

    # cluster bands; only clusters that survive the tolerance filter are named,
    # and labels alternate rows so adjacent clusters do not collide
    ax.set_ylim(len(panel) - 0.5, -3.0)
    labelled = 0
    eps = [e for e in tsv("m5_epitopes.tsv") if e["id"] == target]
    for c in tsv("m5_clusters.tsv"):
        if c["id"] != target:
            continue
        members = [e for e in eps
                   if int(c["start"]) <= int(e["pos"]) <= int(c["end"]) - 8]
        starts = [frame[e["core"]] for e in members] or [int(c["start"])]
        s, e = min(starts) - 1, max(starts) - 1
        col = BAD if c["tolerance_class"] == "foreign" else \
              WARN if c["tolerance_class"] == "mixed" else GOOD
        ax.add_patch(plt.Rectangle((s - 0.5, -0.5), e - s + 1, len(panel),
                                   fill=False, ec=col, lw=1.4, zorder=5))
        if c["tolerance_class"] != "all_tolerised" and int(c["union_sb_alleles"]) >= 2:
            ax.text((s + e) / 2, -2.4 + 1.1 * (labelled % 2), c["peak_core"],
                    color=col, fontsize=7.5, ha="center", va="center",
                    fontweight="bold", zorder=6)
            labelled += 1

    # per-position promiscuity, on the same consensus call the tables use
    prom = C.sum(axis=0)
    axd.fill_between(range(M.shape[1]), prom, color=ACCENT, alpha=0.75, lw=0,
                     step="mid")
    axd.set_ylabel("DR molecules per\n15-mer (consensus SB)", fontsize=7.5)
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
    ax.set_ylim(-0.65, len(rows) + 0.15)
    ax.annotate("Protein A Z-domain benchmark", (anchor, len(rows) - 0.28),
                textcoords="offset points", xytext=(5, 0), fontsize=7.5,
                color=INK, va="center", ha="left")
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
    from matplotlib.patches import Patch
    seen, handles = set(), []
    for r in rows[::-1]:
        if r["role"] not in seen:
            seen.add(r["role"])
            handles.append(Patch(color=ROLE_COLOR.get(r["role"], MUTED),
                                 label=r["role"].replace("_", " ")))
    ax.legend(handles=handles, frameon=False, fontsize=7, loc="lower right",
              handlelength=1.1, labelspacing=0.35)

    # tolerance-filter effect
    raw = [float(r["pIRS_raw"]) for r in rows]
    filt = [float(r["pIRS"]) for r in rows]
    for i, r in enumerate(rows):
        ax2.plot([raw[i], filt[i]], [i, i], color=GRID, lw=2.2, zorder=1,
                 solid_capstyle="round")
    ax2.scatter(raw, y, s=22, color="#c2ccd6", zorder=3,
                label="before self/tolerance filter")
    ax2.scatter(filt, y, s=26, color=[ROLE_COLOR.get(r["role"], MUTED) for r in rows],
                zorder=4)
    ax2.set_yticks(y, ["" for _ in rows])
    ax2.set_xlabel("pIRS")
    ax2.set_title("What the self/tolerance filter removes")
    from matplotlib.lines import Line2D
    ax2.legend(handles=[
        Line2D([], [], marker="o", ls="", ms=4.5, color="#c2ccd6",
               label="before self/tolerance filter"),
        Line2D([], [], marker="o", ls="", ms=5, color=INK,
               label="after filter (dot colour = role)")],
        frameon=False, fontsize=7, loc="lower right", labelspacing=0.35)
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
    x = np.arange(1, L + 1)
    b = np.array(bc)
    a2.plot(x, b, color=ACCENT, lw=1.3)
    a2.fill_between(x, 0.5, b, where=b >= 0.5, color=ACCENT, alpha=0.65, lw=0,
                    interpolate=True)
    a2.axhline(0.5, color=MUTED, ls="--", lw=1)
    a2.text(L, 0.505, "epitope threshold ", fontsize=6.5, color=MUTED,
            ha="right", va="bottom")
    a2.set_ylim(min(b.min(), 0.3) - 0.03, max(b.max(), 0.6) + 0.05)
    a2.set_ylabel("BepiPred-2.0", fontsize=7.5)
    a2.set_xlabel("residue")
    for ax in (a1, a2):
        ax.grid(axis="y", color=GRID, lw=0.6)
        ax.set_axisbelow(True)

    for c in tsv("m7_tb_coincidence.tsv"):
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
    muts = [r for r in rows if r["variant"] != "WT"]

    def key(r):
        return (float(r["pop_presenting"]), int(r["n_sb_alleles"]), -int(r["blosum62"]))

    # Best variant at each anchor pocket, in two classes: the best substitution
    # of any kind, and the best chemically conservative one. Keeping them apart
    # is the point - the knockouts and the designable changes are not the same
    # substitutions.
    pos = sorted({m["variant"].split("_")[0] for m in muts})
    picked, seen = [], set()
    for p in pos:
        at_p = [m for m in muts if m["variant"].startswith(p + "_")]
        for cls, subset in (("conservative", [m for m in at_p if int(m["blosum62"]) >= 0]),
                            ("any", at_p)):
            if not subset:
                continue
            b = min(subset, key=key)
            if b["variant"] in seen:
                continue
            seen.add(b["variant"])
            picked.append({**b, "cls": cls})
    picked.sort(key=lambda r: -float(r["pop_presenting"]))

    y = np.arange(len(picked))
    fig, ax = plt.subplots(figsize=(8.4, 0.34 * len(picked) + 1.9))
    cols = [GOOD if r["cls"] == "conservative" else WARN for r in picked]
    vals = [float(r["pop_presenting"]) * 100 for r in picked]
    ax.barh(y, vals, 0.62, color=cols)
    # a zero-length bar is invisible, so mark the complete knockouts explicitly
    zero = [i for i, v in enumerate(vals) if v == 0]
    ax.scatter([0] * len(zero), zero, marker="|", s=140,
               color=[cols[i] for i in zero], linewidths=2.2, zorder=3)
    wtp = float(wt["pop_presenting"]) * 100
    ax.axvline(wtp, color=BAD, ls="--", lw=1.2)
    ax.set_ylim(-0.7, len(picked) + 0.2)
    ax.annotate(f"wild type {wtp:.0f}%", (wtp, len(picked) - 0.25),
                textcoords="offset points", xytext=(6, 0), color=BAD, fontsize=7.5,
                va="center", ha="left")
    for i, r in enumerate(picked):
        bl = int(r["blosum62"])
        ax.text(vals[i] + 0.7, i,
                f"{vals[i]:.1f}%   {r['n_sb_alleles']} DR   "
                f"BLOSUM {bl:+d}" if bl else f"{vals[i]:.1f}%   {r['n_sb_alleles']} DR   BLOSUM 0",
                va="center", fontsize=7, color=INK)
    ax.set_yticks(y, [f"{r['variant']}   {r['core']}" for r in picked], fontsize=7.5)
    ax.set_xlim(0, max(wtp, 1) * 1.55)
    ax.set_xlabel("% of weighted US/EU population predicted to present the epitope")
    ax.set_title("Best substitution at each anchor pocket, by chemical conservatism", pad=20)
    ax.grid(axis="x", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=GOOD, label="BLOSUM62 ≥ 0 — conservative, designable"),
                       Patch(color=WARN, label="BLOSUM62 < 0 — chemically disruptive")],
              frameon=False, fontsize=7.5, loc="lower center",
              bbox_to_anchor=(0.5, 1.0), ncol=2, handlelength=1.2)
    fig.savefig(figures_path("fig5_deimmunization.png"))
    plt.close(fig)


if __name__ == "__main__":
    todo = sys.argv[1:] or ["panel", "heatmap", "ranking", "tb", "deimm"]
    fns = {"panel": fig_panel_coverage, "heatmap": fig_binding_heatmap,
           "ranking": fig_ranking, "tb": fig_tb, "deimm": fig_deimmunization}
    for t in todo:
        fns[t]()
        print(f"  {t} ok")
