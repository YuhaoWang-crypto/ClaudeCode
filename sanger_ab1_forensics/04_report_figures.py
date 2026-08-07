#!/usr/bin/env python3
"""为 PDF 报告生成两张图：波形对比图 与 dapA 位点示意图。"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from matplotlib.patches import FancyArrow, Rectangle
from Bio import SeqIO

CJK = "/usr/share/fonts/truetype/wqy/wqy-zenhei.ttc"
font_manager.fontManager.addfont(CJK)
plt.rcParams["font.family"] = "WenQuanYi Zen Hei"
plt.rcParams["axes.unicode_minus"] = False

INK, MUTED, BAD, GOOD = "#1f2933", "#7b8794", "#c0392b", "#1e7a5a"
COL = {"G": "#111111", "A": "#1e7a5a", "T": "#c0392b", "C": "#1f4e9c"}


def fig_traces(ab1, out):
    ab = SeqIO.read(ab1, "abi").annotations["abif_raw"]
    fwo = ab["FWO_1"].decode()
    ana = {fwo[i]: np.array(ab[f"DATA{9 + i}"], float) for i in range(4)}
    raw = {fwo[i]: np.array(ab[f"DATA{1 + i}"], float) for i in range(4)}
    seq, loc = ab["PBAS1"].decode(), np.array(ab["PLOC1"])

    fig, ax = plt.subplots(3, 1, figsize=(9.6, 8.2))

    a, s, e = ax[0], 1300, 1900
    for b in "GATC":
        a.plot(range(s, e), ana[b][s:e], COL[b], lw=1.3)
    for p, bb in zip(loc, seq):
        if s <= p < e:
            a.text(p, 15200, bb, ha="center", fontsize=6.5, color=COL[bb])
    a.set_ylim(-600, 16600)
    a.set_title("① 文件中的「分析道」DATA9-12（第 100–146 位碱基）\n"
                "峰高恒定、峰形完美对称、峰间距恒为 13 scan、基线硬性归零",
                fontsize=10, color=BAD, loc="left")

    a, s2 = ax[1], 6000
    for b in "GATC":
        a.plot(range(s2, s2 + (e - s)), raw[b][s2:s2 + (e - s)], COL[b], lw=1.3)
    a.set_title("② 同一文件的「原始道」DATA1-4（等宽窗口）\n"
                "峰高峰宽各不相同、基线有真实噪音 —— 这才是真实毛细管电泳信号",
                fontsize=10, color=GOOD, loc="left")

    a = ax[2]
    env_a = np.maximum.reduce([ana[b] for b in "GATC"])
    env_r = np.maximum.reduce([raw[b] for b in "GATC"])
    xa = np.linspace(0, 1, len(env_a))
    xr = np.linspace(0, 1, len(env_r))
    k = 120
    a.plot(xa, np.convolve(env_a, np.ones(k) / k, "same"), color=BAD, lw=1.6,
           label="分析道包络（735 碱基全程无衰减）")
    a.plot(xr, np.convolve(env_r, np.ones(k) / k, "same"), color=GOOD, lw=1.6,
           label="原始道包络（指数衰减，符合物理规律）")
    a.set_ylim(0, 13500)
    a.legend(fontsize=8.5, frameon=False, loc="upper right", ncol=2)
    a.set_xlabel("归一化迁移位置", fontsize=8.5)
    a.set_title("③ 全长信号包络对比", fontsize=10, color=INK, loc="left")

    for a in ax:
        a.tick_params(labelsize=7, colors=MUTED)
        for sp in a.spines.values():
            sp.set_color("#d8dee4")
        a.spines["top"].set_visible(False)
        a.spines["right"].set_visible(False)
    ax[0].set_xlabel("scan", fontsize=8.5)
    ax[1].set_xlabel("scan", fontsize=8.5)
    fig.tight_layout(h_pad=2.0)
    fig.savefig(out, dpi=200, facecolor="white")
    print("wrote", out)


def fig_junction(out):
    """野生型 vs 送检样品 的位点示意。

    上行按 MG1655 真实坐标绘制；下行把缺失的 878 bp 折叠掉——即 dapA 右侧
    (gcvR 一侧) 整体左移 878 bp，与 bamC 侧直接相接，这才是样品分子的真实构型。
    """
    GENES = [("purC", 2596905, 2597618), ("bamC", 2597831, 2598865),
             ("dapA", 2598882, 2599760), ("gcvR", 2599906, 2600478)]
    DEL_S, DEL_E, DEL_N = 2598882, 2599759, 878
    LO, HI = 2596600, 2600800
    fig, ax = plt.subplots(figsize=(9.6, 3.6))

    def gene(y, s, e, name, color):
        ax.add_patch(FancyArrow(e, y, s - e, 0, width=0.19, head_width=0.19,
                                head_length=(e - s) * 0.16, length_includes_head=True,
                                color=color, zorder=3))
        ax.text((s + e) / 2, y, name, ha="center", va="center", fontsize=8.5,
                color="white", zorder=4)

    # ── 上行：野生型 ────────────────────────────────────────────────
    ax.text(LO, 1.32, "野生型 E. coli  (K-12 MG1655, NC_000913.3)",
            fontsize=9.5, color=INK)
    ax.add_patch(Rectangle((LO, 0.965), HI - LO, 0.07, color="#c9d2da", zorder=1))
    for g, s, e in GENES:
        gene(1.0, s, e, g, BAD if g == "dapA" else "#41608f")
    ax.annotate("", xy=(DEL_S, 0.70), xytext=(DEL_E + 1, 0.70),
                arrowprops=dict(arrowstyle="<->", color=BAD, lw=1.3))
    ax.text((DEL_S + DEL_E) / 2, 0.76,
            f"缺失 {DEL_N} bp  ({DEL_S:,}–{DEL_E:,}) = dapA CDS 的 99.9%",
            ha="center", fontsize=9, color=BAD)

    # ── 下行：样品（缺失折叠，两侧直接相接）──────────────────────────
    shift = -DEL_N                      # gcvR 一侧左移
    ax.text(LO, 0.30, "送检样品 DAPA21  (本次 735 bp 读段所示构型)",
            fontsize=9.5, color=INK)
    ax.add_patch(Rectangle((LO, -0.035), (HI + shift) - LO, 0.07,
                           color="#c9d2da", zorder=1))
    gene(0.0, 2596905, 2597618, "purC", "#41608f")
    gene(0.0, 2597831, 2598865, "bamC", "#41608f")
    gene(0.0, 2599906 + shift, 2600478 + shift, "gcvR", "#41608f")

    J = DEL_S                            # 折叠后的连接点
    ax.plot([DEL_S, J], [0.62, 0.10], ls=":", lw=1, color=BAD)
    ax.plot([DEL_E + 1, J], [0.62, 0.10], ls=":", lw=1, color=BAD)
    ax.plot([J], [0.0], marker="v", ms=9, color=BAD, zorder=5)
    ax.text(J, -0.23, "无痕连接点：1 bp 微同源 “T”，无 FRT / loxP / 抗性盒",
            ha="center", fontsize=8.5, color=BAD)

    # 读段覆盖（折叠坐标下是一段连续 735 bp）
    r0, r1 = 2598517, 2600130 + shift
    ax.add_patch(Rectangle((r0, -0.44), r1 - r0, 0.055, color=MUTED, zorder=3))
    ax.text((r0 + r1) / 2, -0.56, "读段覆盖 735 bp（跨断点连续）",
            ha="center", fontsize=8.5, color=MUTED)

    ax.set_xlim(LO, HI)
    ax.set_ylim(-0.70, 1.52)
    ax.set_yticks([])
    ax.set_xlabel("MG1655 坐标 (bp) —— 下行按缺失折叠后绘制", fontsize=8.5, color=MUTED)
    ax.ticklabel_format(axis="x", style="plain", useOffset=False)
    ax.tick_params(axis="x", labelsize=7.5, colors=MUTED)
    for sp in ax.spines.values():
        sp.set_visible(False)
    fig.tight_layout()
    fig.savefig(out, dpi=200, facecolor="white")
    print("wrote", out)


if __name__ == "__main__":
    fig_traces(sys.argv[1], "results/fig1_trace_comparison.png")
    fig_junction("results/fig2_dapA_junction.png")
