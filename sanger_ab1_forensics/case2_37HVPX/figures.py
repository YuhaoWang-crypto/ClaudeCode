#!/usr/bin/env python3
"""案例 2 报告用图。

用法: python3 figures.py <37HVPX_1_D12.ab1>
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
from matplotlib import font_manager
from matplotlib.patches import Rectangle

font_manager.fontManager.addfont("/usr/share/fonts/truetype/wqy/wqy-zenhei.ttc")
plt.rcParams["font.family"] = "WenQuanYi Zen Hei"
plt.rcParams["axes.unicode_minus"] = False

INK, MUTED, BAD, GOOD, ACC = "#1f2933", "#7b8794", "#c0392b", "#1e7a5a", "#41608f"
COL = {"G": "#111111", "A": "#1e7a5a", "T": "#c0392b", "C": "#1f4e9c"}

# 500 bp 滑窗：(读段起点, 该窗平均QV, BLAST 一致性%, 是否落在自定义 CAR/VHH 区)
WINDOWS = [(1, 16.4, 100.0, False), (501, 10.2, 85.9, True), (1001, 9.8, 78.6, True),
           (1501, 9.1, 78.6, True), (2001, 9.0, 98.6, True), (2501, 9.0, 100.0, False),
           (3001, 9.8, 100.0, False), (3501, 10.5, 99.8, False), (4001, 11.2, 100.0, False)]


def fig_render(ab1, out):
    ab = SeqIO.read(ab1, "abi").annotations["abif_raw"]
    fwo = ab["FWO_1"].decode()
    ana = {fwo[i]: np.array(ab[f"DATA{9 + i}"], float) for i in range(4)}
    seq = ab["PBAS1"].decode().upper()
    loc = np.array(ab["PLOC1"])
    qv = np.array([ord(c) for c in ab["PCON1"].decode("latin-1")])
    amp = np.array([ana[b][p] for p, b in zip(loc, seq)])

    fig, ax = plt.subplots(1, 2, figsize=(9.0, 4.2),
                           gridspec_kw={"width_ratios": [1.35, 1]})

    a, s, e = ax[0], 800, 940
    for b in "GATC":
        a.step(range(s, e), ana[b][s:e], COL[b], lw=1.1, where="post")
    for p, bb in zip(loc, seq):
        if s <= p < e:
            a.text(p, -22, bb, ha="center", fontsize=6.5, color=COL[bb])
    a.set_xlabel("scan（每碱基固定 4 个采样点）", fontsize=10)
    a.set_ylabel("强度", fontsize=10)
    a.set_title("① 色谱图是方波渲染，不是电泳信号", fontsize=11.5, color=BAD, loc="left")

    a = ax[1]
    a.scatter(qv, amp, s=3, alpha=0.25, color=ACC, edgecolors="none")
    k = np.polyfit(qv, amp, 1)
    xs = np.array([qv.min(), qv.max()])
    a.plot(xs, np.polyval(k, xs), color=BAD, lw=1.3)
    r = float(np.corrcoef(amp, qv)[0, 1])
    a.text(0.04, 0.93, f"峰高 ≈ {k[0]:.1f} × QV + {k[1]:.0f}\nPearson r = {r:.4f}",
           transform=a.transAxes, fontsize=10.5, color=BAD, va="top")
    a.set_xlabel("该碱基的质量值 QV", fontsize=10)
    a.set_ylabel("峰高", fontsize=10)
    a.set_title("② 峰高完全由质量值决定", fontsize=11.5, color=BAD, loc="left")

    for a in ax:
        a.tick_params(labelsize=8.5, colors=MUTED)
        for sp in a.spines.values():
            sp.set_color("#d8dee4")
        a.spines["top"].set_visible(False)
        a.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(out, dpi=200, facecolor="white")
    print("wrote", out)


def fig_quality(out):
    fig, a = plt.subplots(figsize=(9.2, 4.5))
    x = np.arange(len(WINDOWS))
    implied = [100 * 10 ** (-q / 10) for _, q, _, _ in WINDOWS]
    observed = [100 - i for _, _, i, _ in WINDOWS]
    custom = [c for _, _, _, c in WINDOWS]
    w = 0.38

    a.bar(x - w / 2, implied, w, color=BAD, label="质量值隐含的错误率（QV → 10^(-QV/10)）")
    a.bar(x + w / 2, observed, w, color=ACC, label="实测与公共参考的不一致率")
    lo = min(i for i, c in enumerate(custom) if c) - 0.5
    hi = max(i for i, c in enumerate(custom) if c) + 0.5
    a.add_patch(Rectangle((lo, 0), hi - lo, 26, color="#f0d9a0", alpha=0.35, zorder=0))
    for i, (_, _, ident, c) in enumerate(WINDOWS):
        if not c:                              # 零高柱：画个短横线当锚点，免得标签像贴在红柱上
            a.plot([i + w / 2 - w / 2, i + w / 2 + w / 2], [0.06, 0.06], lw=2.2, color=ACC)
            a.text(i + w / 2, 0.9, f"一致 {ident:.0f}%", ha="center",
                   fontsize=9, color=GOOD)
    a.text(2.0, 24.2, "自定义 VHH / CAR 区：公共库无确切参考，\n低一致性不代表测序错误",
           ha="center", fontsize=10, color="#8a6d1f")

    a.set_xticks(x)
    a.set_xticklabels([f"{s}-{s + 499}" for s, _, _, _ in WINDOWS], fontsize=8.5)
    a.set_xlabel("读段上的 500 bp 窗口", fontsize=10)
    a.set_ylabel("错误率 / 不一致率 (%)", fontsize=10)
    a.set_ylim(0, 27)
    a.legend(fontsize=9.5, frameon=False, loc="upper right")
    a.set_title("质量值声称约 10% 错误率，但有参考的窗口实测 100% 一致",
                fontsize=12, color=INK, loc="left")
    a.tick_params(labelsize=8.5, colors=MUTED)
    for sp in a.spines.values():
        sp.set_color("#d8dee4")
    a.spines["top"].set_visible(False)
    a.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(out, dpi=200, facecolor="white")
    print("wrote", out)


def fig_map(out):
    L = 4580
    ELEM = [(1, 850, "HIV-1 env / RRE", ACC), (861, 899, "cPPT", "#6b7f9e"),
            (1176, 2350, "人 EF-1α 启动子（含内含子）", "#2f6f4f"),
            (2358, 2401, "MCS", MUTED), (2411, 3898, "CAR 阅读框 1488 nt", BAD),
            (3919, 4500, "WPRE", "#2f6f4f"), (4510, 4544, "3'ΔU3 LTR", ACC),
            (4545, 4580, "折返", "#c9773a")]
    DOM = [(1, 21, "CD8α\n信号肽", MUTED), (22, 142, "VHH1", BAD),
           (143, 147, "", "#c9c9c9"), (148, 260, "VHH2", BAD),
           (261, 272, "", "#c9c9c9"), (273, 317, "CD8α 铰链", "#6b7f9e"),
           (318, 341, "CD8α TM", "#41608f"), (342, 383, "4-1BB", "#2f6f4f"),
           (384, 495, "CD3ζ", "#1f4e9c")]

    fig, ax = plt.subplots(2, 1, figsize=(9.0, 4.8),
                           gridspec_kw={"height_ratios": [1, 1]})

    a = ax[0]
    a.add_patch(Rectangle((1, -0.055), L, 0.11, color="#dfe4ea", zorder=1))
    tier = 0
    for s, e, n, c in ELEM:
        a.add_patch(Rectangle((s, -0.20), e - s, 0.40, color=c, zorder=3))
        if e - s > 200:
            a.text((s + e) / 2, 0, n, ha="center", va="center", fontsize=9.5,
                   color="white", zorder=4)
        else:                                   # 窄元件：标签交替错开两层
            y = 0.28 + 0.16 * (tier % 2)
            a.plot([(s + e) / 2, (s + e) / 2], [0.20, y - 0.03], lw=0.6, color=c)
            a.text((s + e) / 2, y, n, ha="center", fontsize=9, color=c)
            tier += 1
    a.set_xlim(-60, L + 230)
    a.set_ylim(-0.55, 0.72)
    a.set_yticks([])
    a.set_xlabel("载体正义方向坐标 (bp)　— 全长 4580 bp", fontsize=10, color=MUTED)
    a.set_title("慢病毒转移载体表达盒", fontsize=12, color=INK, loc="left")

    a = ax[1]
    for s, e, n, c in DOM:
        a.add_patch(Rectangle((s, -0.20), e - s, 0.40, color=c, zorder=3))
        if n and e - s > 25:
            a.text((s + e) / 2, 0, n, ha="center", va="center", fontsize=8.8,
                   color="white", zorder=4)
        elif n:
            a.text((s + e) / 2, 0.30, n, ha="center", fontsize=8.5, color=c)
    a.set_xlim(-6, 501)
    a.set_ylim(-0.55, 0.62)
    a.set_yticks([])
    a.set_xlabel("CAR 蛋白氨基酸位置 — 全长 495 aa，无移码、无提前终止",
                 fontsize=10, color=MUTED)
    a.set_title("CAR 结构域组成（二代 BBζ，双价 VHH）", fontsize=12, color=INK, loc="left")

    for a in ax:
        a.tick_params(axis="x", labelsize=8.5, colors=MUTED)
        for sp in a.spines.values():
            sp.set_visible(False)
    fig.tight_layout(h_pad=2.4)
    fig.savefig(out, dpi=200, facecolor="white")
    print("wrote", out)


if __name__ == "__main__":
    fig_render(sys.argv[1], "fig1_render_evidence.png")
    fig_quality("fig2_quality_vs_accuracy.png")
    fig_map("fig3_vector_map.png")
