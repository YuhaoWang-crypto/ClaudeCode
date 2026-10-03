"""Chinese-language versions of the campaign figures.

Same data, same validated palette and chrome as `figures.py`; only the text is
localised and type sizes are nudged for CJK glyphs, which need slightly more
room than Latin at the same point size.
"""

from __future__ import annotations

import json
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch

from . import cjk
from .figures import (BASELINE, CONTEXT, GRID, INK, MUTED, S1, S2, S3,
                      SURFACE)
from .figures import INK_2 as INK2

# Cascade stage names come from the English pipeline log; map them for display.
STAGE_ZH = {
    "valid, unique, un-gated": "合法、去重、通过性质门",
    "novel (max Tanimoto < 0.4)": "新颖性（最大 Tanimoto < 0.4）",
    "activity >= 0.5, QED >= 0.6": "活性 ≥ 0.5，QED ≥ 0.6",
    "PAINS-free": "无 PAINS 子结构",
    "ring sanity": "环系合理性",
    "SA_Score <= 4.5": "合成可及性 SA ≤ 4.5",
}


def _setup():
    cjk.register_matplotlib()
    plt.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
        "axes.edgecolor": BASELINE, "axes.labelcolor": INK2,
        "axes.titlecolor": INK, "text.color": INK,
        "xtick.color": MUTED, "ytick.color": MUTED,
        "xtick.labelcolor": INK2, "ytick.labelcolor": INK2,
        "grid.color": GRID, "grid.linewidth": 0.8, "axes.grid": True,
        "axes.spines.top": False, "axes.spines.right": False,
        "legend.frameon": False, "figure.dpi": 200, "savefig.dpi": 200,
        "savefig.facecolor": SURFACE, "axes.unicode_minus": False,
    })


def _clean(ax, title=None, xlabel=None, ylabel=None):
    if title:
        ax.set_title(title, fontsize=11.5, loc="left", pad=10, color=INK)
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=9.5)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=9.5)
    ax.set_axisbelow(True)
    ax.tick_params(length=0, labelsize=8.5)
    return ax


def fig_convergence(history, out):
    gens = [h["generation"] for h in history]
    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    series = [("种群最优", [h["best"] for h in history], S1),
              ("种群均值", [h["mean"] for h in history], S2),
              ("新颖分子最优（Tanimoto < 0.4）",
               [h["best_novel"] for h in history], S3)]
    for label, ys, c in series:
        ax.plot(gens, ys, color=c, linewidth=2.0, solid_capstyle="round",
                label=label)
    ends = sorted(((s[1][-1], s[2]) for s in series), key=lambda e: e[0])
    placed: list[float] = []
    for value, colour in ends:
        y = value if not placed else max(value, placed[-1] + 0.016)
        placed.append(y)
        ax.annotate(f"{value:.3f}", (gens[-1], y), textcoords="offset points",
                    xytext=(7, 0), fontsize=8.5, color=colour, va="center")
    _clean(ax, "遗传算法收敛曲线", "迭代代数", "目标函数得分")
    ax.set_xlim(1, max(gens) + 3.5)
    ax.set_ylim(0.55, 0.95)
    ax.legend(fontsize=8.5, loc="lower right", labelcolor=INK2)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)


def fig_activity_qed(all_scored, survivors, top10, out):
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    a = np.array([[r["activity"], r["qed"]] for r in all_scored
                  if not r.get("gated")])
    s = np.array([[r["activity"], r["qed"]] for r in survivors])
    t = np.array([[r["activity"], r["qed"]] for r in top10])
    ax.scatter(a[:, 0], a[:, 1], s=7, c=CONTEXT, alpha=0.45, linewidths=0,
               label=f"全部生成分子（n={len(a)}）", rasterized=True)
    ax.scatter(s[:, 0], s[:, 1], s=26, c=S1, alpha=0.85, linewidths=0.6,
               edgecolors=SURFACE, label=f"通过筛选级联（n={len(s)}）")
    ax.scatter(t[:, 0], t[:, 1], s=95, c=S2, linewidths=1.6,
               edgecolors=SURFACE, label=f"最终入选（n={len(t)}）", zorder=5)
    for r in top10:
        ax.annotate(r["design_id"].replace("GSK3B-DN-", ""),
                    (r["activity"], r["qed"]), fontsize=6.5, color=SURFACE,
                    ha="center", va="center", zorder=6)
    ax.axvline(0.5, color=MUTED, linewidth=0.9, linestyle=(0, (4, 3)), zorder=1)
    ax.axhline(0.6, color=MUTED, linewidth=0.9, linestyle=(0, (4, 3)), zorder=1)
    ax.annotate("活性阈值 0.50", (0.49, 0.03), fontsize=7.5, color=MUTED,
                rotation=90, ha="right", va="bottom")
    ax.annotate("QED 阈值 0.60", (0.015, 0.615), fontsize=7.5, color=MUTED,
                ha="left", va="bottom")
    _clean(ax, "代理活性与类药性分布", "GSK3B 代理模型活性评分", "QED 类药性")
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(0, 1.0)
    ax.legend(fontsize=8.5, loc="upper left", labelcolor=INK2)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)


def fig_cascade(log, out):
    labels = [STAGE_ZH.get(e["stage"], e["stage"]) for e in log]
    vals = [e["out"] for e in log]
    start = log[0]["in"]
    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    for i, v in enumerate(vals):
        ax.add_patch(FancyBboxPatch(
            (0, i - 0.34), max(v, start * 0.004), 0.68,
            boxstyle="round,pad=0,rounding_size=%f" % (start * 0.006),
            mutation_aspect=0.0012, linewidth=0, facecolor=S1))
        ax.annotate(f"{v:,}", (v, i), xytext=(8, 0), textcoords="offset points",
                    fontsize=9, color=INK2, va="center")
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=8.5)
    ax.invert_yaxis()
    ax.set_xlim(0, start * 1.14)
    ax.set_ylim(len(labels) - 0.5, -0.5)
    _clean(ax, f"筛选级联各阶段剩余分子数（起始 {start:,} 个）", "剩余分子数", None)
    ax.grid(axis="y", visible=False)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)


def fig_properties(gen_props, ref_props, out):
    panels = [("mw", "分子量（Da）", (100, 650)),
              ("logp", "脂水分配系数 cLogP", (-2, 8)),
              ("tpsa", "拓扑极性表面积 TPSA（Å$^2$）", (0, 200)),
              ("sa", "合成可及性 SA（1 易 – 10 难）", (1, 7))]
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.0))
    for ax, (key, label, rng) in zip(axes.ravel(), panels):
        g = np.array([p[key] for p in gen_props if p.get(key) is not None])
        r = np.array([p[key] for p in ref_props if p.get(key) is not None])
        bins = np.linspace(rng[0], rng[1], 36)
        ax.hist(g, bins=bins, color=S1, alpha=0.75, density=True,
                label="本次生成（通过级联）", linewidth=0)
        ax.hist(r, bins=bins, histtype="step", color=S2, linewidth=2.0,
                density=True, label="ChEMBL 已知活性分子")
        _clean(ax, None, label, "密度")
        ax.set_yticks([])
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, fontsize=8.5, loc="upper right",
               bbox_to_anchor=(0.995, 0.975), ncol=2, labelcolor=INK2)
    fig.suptitle("理化性质分布：设计分子 vs 已知活性分子", fontsize=11.5,
                 x=0.02, ha="left", y=0.995, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)


def fig_novelty(all_sims, surv_sims, threshold, out):
    fig, ax = plt.subplots(figsize=(7.2, 3.4))
    bins = np.linspace(0, 1, 51)
    ax.hist(all_sims, bins=bins, color=CONTEXT, linewidth=0,
            label=f"全部生成分子（n={len(all_sims)}）")
    ax.hist(surv_sims, bins=bins, color=S1, linewidth=0,
            label=f"通过筛选级联（n={len(surv_sims)}）")
    ax.axvline(threshold, color=S2, linewidth=2.0)
    ax.annotate(f"新颖性阈值 {threshold:.2f}",
                (threshold, ax.get_ylim()[1] * 0.92), fontsize=8.5, color=S2,
                ha="right", va="top", xytext=(-6, 0),
                textcoords="offset points")
    _clean(ax, "生成分子与最近已知活性分子的相似度分布",
           "与任一 ChEMBL GSK3B 活性分子的最大 Tanimoto 相似度（ECFP4）", "分子数")
    ax.legend(fontsize=8.5, loc="upper right", labelcolor=INK2)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)


def build_all(result_dir: str):
    """Regenerate figures 1-5 in Chinese; structure grids are language-neutral."""
    _setup()
    R = pathlib.Path(result_dir)
    figs = R / "figures_zh"
    figs.mkdir(parents=True, exist_ok=True)
    hist = json.load(open(R / "ga_history.json"))
    allsc = json.load(open(R / "ga_all_scored.json"))
    surv = json.load(open(R / "cascade_survivors.json"))
    top = json.load(open(R / "top10_designs.json"))
    log = json.load(open(R / "cascade_log.json"))
    fd = json.load(open(R / "fig_data.json"))

    fig_convergence(hist, figs / "fig1_convergence.png")
    fig_activity_qed(allsc, surv, top, figs / "fig2_activity_qed.png")
    fig_cascade(log, figs / "fig3_cascade.png")
    fig_properties(fd["gen_props"], fd["ref_props"], figs / "fig4_properties.png")
    fig_novelty(fd["all_sims"], [r["max_sim_to_known"] for r in surv],
                0.40, figs / "fig5_novelty.png")
    return figs
