r"""
E12 - Schematic: every pattern the atlas systems can make, and what opens the gate.
===================================================================================

A drawing, not a measurement.  Twelve morphology classes, each with the
architecture that produces it, the condition that has to hold, and the real
protein groups that sit there.  Numbers in the captions are the computed ones
from E10/E11 and ``test_pathways``/``test_function``; the cartoons themselves
are illustrative.

Layout
------
row 1  what the valence gate decides:  1:1 complex / 1-D filament / 3-D network
row 2  what opens the switch:          membrane 2-D / template length / phosphorylation
row 3  many bodies and selectivity:    sorted droplets / core-shell / client vs excluded
row 4  the boundary:                   allostery (bolted on) / ageing (no) / motion (forbidden)
"""

import numpy as np
from common import ACC, BG, FG, GRID, save
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, Rectangle

plt.rcParams["font.sans-serif"] = ["WenQuanYi Zen Hei", "DejaVu Sans"]

RNG = np.random.default_rng(7)
P1, P2, P3 = ACC[0], ACC[1], ACC[2]      # green, magenta, blue
DIM = "#8b949e"


# --------------------------------------------------------------------------- #
#  drawing primitives
# --------------------------------------------------------------------------- #
def blob(a, x, y, r=0.06, c=P1, alpha=1.0, ec=None, lw=0.0, z=3):
    a.add_patch(Circle((x, y), r, facecolor=c, alpha=alpha, zorder=z,
                       edgecolor=ec or c, linewidth=lw))


def bond(a, p, q, c=DIM, lw=1.3, alpha=0.9, z=2):
    a.plot([p[0], q[0]], [p[1], q[1]], color=c, lw=lw, alpha=alpha, zorder=z,
           solid_capstyle="round")


def arrow(a, p, q, c=FG, lw=1.6, style="-|>", ms=9, z=6):
    a.add_patch(FancyArrowPatch(p, q, arrowstyle=style, mutation_scale=ms,
                                color=c, lw=lw, zorder=z,
                                shrinkA=0, shrinkB=0))


def mark(a, kind, x=0.945, y=0.925, s=0.030, tr=None):
    """Verdict badge drawn as a shape - no font glyph needed."""
    tr = tr or a.transAxes
    col = {"ok": ACC[0], "warn": ACC[3], "no": ACC[6]}[kind]
    if kind == "ok":                                   # filled disc
        a.add_patch(Circle((x, y), s, facecolor=col, transform=tr,
                           zorder=8, clip_on=False))
    elif kind == "warn":                               # triangle
        a.add_patch(plt.Polygon([[x, y + s], [x - s, y - s * 0.8],
                                 [x + s, y - s * 0.8]], facecolor=col,
                                transform=tr, zorder=8, clip_on=False))
    else:                                              # cross
        for dx in (+1, -1):
            a.plot([x - dx * s * 0.8, x + dx * s * 0.8], [y - s * 0.8, y + s * 0.8],
                   color=col, lw=2.4, transform=tr, zorder=8, clip_on=False,
                   solid_capstyle="round")
    return col


def frame(a, title, cond, systems, ok="ok"):
    """Common panel furniture: title, condition strip, systems strip."""
    a.set_xlim(0, 1); a.set_ylim(0, 1); a.set_aspect("equal")
    a.axis("off")
    a.text(0.5, 1.06, title, ha="center", va="bottom", fontsize=11.5,
           color=FG, transform=a.transAxes)
    col = mark(a, ok)
    a.text(0.02, 0.145, cond, ha="left", va="top", fontsize=8.2,
           color=col, transform=a.transAxes)
    a.text(0.02, 0.015, systems, ha="left", va="bottom", fontsize=8.2,
           color=DIM, transform=a.transAxes)
    a.add_patch(Rectangle((0.0, 0.19), 1.0, 0.79, transform=a.transAxes,
                          facecolor="none", edgecolor=GRID, lw=1.0, zorder=0))


def band(a):
    """Usable drawing area inside the frame (data coords)."""
    return 0.06, 0.94, 0.24, 0.93


# --------------------------------------------------------------------------- #
#  the twelve cartoons
# --------------------------------------------------------------------------- #
def p_complex(a):
    """1:1 complex - valence 1 x 1, nothing connects."""
    for x, y in [(0.22, 0.75), (0.60, 0.82), (0.36, 0.45), (0.75, 0.52),
                 (0.20, 0.34), (0.66, 0.33)]:
        blob(a, x, y, 0.055, P1)
        blob(a, x + 0.10, y, 0.055, P2)
        bond(a, (x, y), (x + 0.10, y), lw=2.2, c=ACC[3])


def p_filament(a):
    """1-D filament - one head, one tail per molecule."""
    for k, y0 in enumerate((0.75, 0.50)):
        n = 9
        t = np.linspace(0, 1, n)
        xs = 0.10 + 0.80 * t
        ys = y0 + 0.055 * np.sin(2 * np.pi * t * 1.4 + k)
        for i in range(n - 1):
            bond(a, (xs[i], ys[i]), (xs[i + 1], ys[i + 1]), lw=2.4, c=ACC[3])
        for i in range(n):
            blob(a, xs[i], ys[i], 0.045, P1 if k == 0 else P3)
    a.text(0.5, 0.31, "每个分子只有一头一尾 → 只能延长，不能分叉",
           ha="center", fontsize=7.6, color=DIM)


def p_network(a):
    """3-D percolating network."""
    pts = np.column_stack([RNG.uniform(0.12, 0.88, 22), RNG.uniform(0.32, 0.90, 22)])
    d = np.linalg.norm(pts[:, None] - pts[None, :], axis=-1)
    for i in range(len(pts)):
        for j in range(i + 1, len(pts)):
            if d[i, j] < 0.26:
                bond(a, pts[i], pts[j], lw=1.2, c=DIM, alpha=0.8)
    for k, (x, y) in enumerate(pts):
        blob(a, x, y, 0.042, P1 if k % 2 else P2)


def p_membrane(a):
    """2-D condensate on a membrane; dilute above."""
    a.add_patch(Rectangle((0.05, 0.33), 0.90, 0.055, facecolor=ACC[4],
                          alpha=0.45, zorder=1))
    a.text(0.5, 0.295, "膜", ha="center", fontsize=8, color=ACC[4])
    xs = np.linspace(0.20, 0.80, 9)
    for i, x in enumerate(xs):
        blob(a, x, 0.415, 0.040, P1)
        if i < len(xs) - 1:
            bond(a, (x, 0.415), (xs[i + 1], 0.415), lw=1.4)
        blob(a, x + 0.02, 0.50, 0.034, P2)
        bond(a, (x, 0.415), (x + 0.02, 0.50), lw=1.1)
        if i % 2 == 0:
            blob(a, x - 0.01, 0.575, 0.030, P3)
            bond(a, (x + 0.02, 0.50), (x - 0.01, 0.575), lw=1.1)
    for x, y in RNG.uniform([0.10, 0.72], [0.90, 0.90], (7, 2)):
        blob(a, x, y, 0.030, P2, alpha=0.5)
    a.text(0.5, 0.965, "胞浆里是散的", ha="center", fontsize=7.6, color=DIM)


def p_template(a):
    """Polymer length supplies the valence."""
    t = np.linspace(0, 1, 120)
    xs = 0.08 + 0.84 * t
    ys = 0.72 + 0.05 * np.sin(2 * np.pi * t * 2.2)
    a.plot(xs, ys, color=ACC[5], lw=2.4, zorder=2)
    idx = np.linspace(5, 114, 9).astype(int)
    prev = None
    for i in idx:
        blob(a, xs[i], ys[i] - 0.075, 0.040, P1)
        bond(a, (xs[i], ys[i]), (xs[i], ys[i] - 0.075), lw=1.1)
        if prev is not None:
            bond(a, (prev, ys[i] - 0.075), (xs[i], ys[i] - 0.075), lw=1.3, c=DIM)
        prev = xs[i]
    a.text(0.5, 0.80, "长模板：位点多 → 成网", ha="center", fontsize=7.6, color=ACC[0])
    a.plot([0.36, 0.64], [0.44, 0.44], color=ACC[5], lw=2.4)
    for x in (0.42, 0.58):
        blob(a, x, 0.375, 0.040, P1)
        bond(a, (x, 0.44), (x, 0.375), lw=1.1)
    a.text(0.5, 0.32, "短模板：价数不够 → 点不着", ha="center", fontsize=7.6, color=ACC[6])


def p_phospho(a):
    """Phosphorylation creates the valence."""
    a.add_patch(Rectangle((0.08, 0.74), 0.34, 0.045, facecolor=DIM, alpha=0.7))
    a.text(0.25, 0.685, "未磷酸化：无位点", ha="center", fontsize=7.6, color=ACC[6])
    for x, y in RNG.uniform([0.10, 0.83], [0.40, 0.90], (5, 2)):
        blob(a, x, y, 0.030, P2, alpha=0.55)
    a.add_patch(Rectangle((0.56, 0.74), 0.36, 0.045, facecolor=ACC[3], alpha=0.85))
    px = np.linspace(0.60, 0.88, 4)
    for x in px:
        blob(a, x, 0.7625, 0.022, FG, z=5)
        blob(a, x, 0.845, 0.038, P2)
        bond(a, (x, 0.7625), (x, 0.845), lw=1.2)
    for i in range(len(px) - 1):
        bond(a, (px[i], 0.845), (px[i + 1], 0.845), lw=1.3, c=DIM)
    a.text(0.74, 0.685, "磷酸化：价数 0 → n", ha="center", fontsize=7.6, color=ACC[0])
    arrow(a, (0.455, 0.78), (0.545, 0.78), c=FG, lw=1.4)
    a.text(0.5, 0.40, "细胞开关用的是价数，不是亲和力",
           ha="center", fontsize=8.4, color=FG)


def p_sorted(a):
    """Several immiscible bodies - the 'cells' pattern."""
    spec = [(0.26, 0.72, 0.135, P1), (0.68, 0.76, 0.115, P2),
            (0.30, 0.40, 0.105, P3), (0.70, 0.42, 0.125, ACC[3])]
    for cx, cy, r, col in spec:
        a.add_patch(Circle((cx, cy), r, facecolor=col, alpha=0.30, zorder=1))
        for _ in range(9):
            th, rr = RNG.uniform(0, 2 * np.pi), r * np.sqrt(RNG.uniform(0, 0.7))
            blob(a, cx + rr * np.cos(th), cy + rr * np.sin(th), 0.028, col)


def p_coreshell(a):
    """Core-shell / concentric layers."""
    a.add_patch(Circle((0.5, 0.60), 0.245, facecolor=P3, alpha=0.25, zorder=1))
    a.add_patch(Circle((0.5, 0.60), 0.125, facecolor=P2, alpha=0.42, zorder=2))
    for _ in range(14):
        th, rr = RNG.uniform(0, 2 * np.pi), RNG.uniform(0.145, 0.225)
        blob(a, 0.5 + rr * np.cos(th), 0.60 + rr * np.sin(th), 0.026, P3)
    for _ in range(7):
        th, rr = RNG.uniform(0, 2 * np.pi), RNG.uniform(0, 0.09)
        blob(a, 0.5 + rr * np.cos(th), 0.60 + rr * np.sin(th), 0.026, P2)
    a.text(0.5, 0.295, "谁包谁 = 界面张力，模型只能给排序",
           ha="center", fontsize=7.6, color=ACC[3])


def p_client(a):
    """Clients enriched, big inert molecules excluded by size."""
    a.add_patch(Circle((0.44, 0.62), 0.245, facecolor=P1, alpha=0.22, zorder=1))
    for _ in range(10):
        th, rr = RNG.uniform(0, 2 * np.pi), 0.235 * np.sqrt(RNG.uniform(0, 0.85))
        blob(a, 0.44 + rr * np.cos(th), 0.62 + rr * np.sin(th), 0.028, P1)
    for _ in range(6):
        th, rr = RNG.uniform(0, 2 * np.pi), 0.20 * np.sqrt(RNG.uniform(0, 0.8))
        blob(a, 0.44 + rr * np.cos(th), 0.62 + rr * np.sin(th), 0.020, ACC[3])
    for x, y in ((0.83, 0.83), (0.87, 0.58), (0.79, 0.36)):
        blob(a, x, y, 0.070, P2, alpha=0.85)
    a.text(0.5, 0.285, "小客体富集 2×10⁴　　大惰性分子 K = 0.0025",
           ha="center", fontsize=7.4, color=DIM)


def p_allostery(a):
    """A ligand that changes a PARTNER's valence."""
    a.add_patch(Circle((0.20, 0.78), 0.075, facecolor=DIM, alpha=0.8, zorder=3))
    a.text(0.20, 0.665, "自抑制", ha="center", fontsize=7.4, color=ACC[6])
    arrow(a, (0.30, 0.78), (0.40, 0.78), c=ACC[0], lw=1.4)
    a.text(0.35, 0.815, "CAPRIN1", ha="center", fontsize=7.2, color=ACC[0])
    cx = 0.56
    a.add_patch(Circle((cx, 0.78), 0.075, facecolor=P1, alpha=0.9, zorder=3))
    for k in range(4):
        th = np.pi * (0.15 + 0.55 * k / 3)
        blob(a, cx + 0.13 * np.cos(th), 0.78 + 0.13 * np.sin(th), 0.028, P2)
        bond(a, (cx, 0.78), (cx + 0.13 * np.cos(th), 0.78 + 0.13 * np.sin(th)), lw=1.1)
    a.text(cx, 0.635, "价数 0.8 → 4", ha="center", fontsize=7.4, color=ACC[0])
    arrow(a, (0.20, 0.68), (0.20, 0.52), c=ACC[6], lw=1.4)
    a.text(0.315, 0.575, "USP10 抢同一位点", fontsize=7.2, color=ACC[6])
    a.add_patch(Circle((0.20, 0.44), 0.075, facecolor=DIM, alpha=0.8, zorder=3))
    a.text(0.5, 0.325, "[USP10] > 4.5×[CAPRIN1] → 解聚",
           ha="center", fontsize=8.0, color=ACC[3])


def p_ageing(a):
    """Liquid -> gel -> fibre: kinetics, outside the model."""
    a.add_patch(Circle((0.18, 0.68), 0.105, facecolor=P3, alpha=0.35, zorder=1))
    a.text(0.18, 0.53, "液滴", ha="center", fontsize=7.6, color=DIM)
    pts = np.column_stack([RNG.uniform(0.42, 0.60, 10), RNG.uniform(0.60, 0.78, 10)])
    for i in range(len(pts)):
        for j in range(i + 1, len(pts)):
            if np.linalg.norm(pts[i] - pts[j]) < 0.11:
                bond(a, pts[i], pts[j], lw=1.1)
    for x, y in pts:
        blob(a, x, y, 0.025, ACC[3])
    a.text(0.51, 0.53, "凝胶", ha="center", fontsize=7.6, color=DIM)
    for k in range(4):
        y = 0.61 + k * 0.045
        a.plot([0.74, 0.94], [y, y], color=ACC[6], lw=2.0, solid_capstyle="round")
    a.text(0.84, 0.53, "纤维", ha="center", fontsize=7.6, color=DIM)
    arrow(a, (0.10, 0.42), (0.92, 0.42), c=FG, lw=1.3)
    a.text(0.5, 0.325, "时间 —— 平衡模型里没有这根轴",
           ha="center", fontsize=8.2, color=ACC[6])


def p_motion(a):
    """Directed motion: forbidden at equilibrium."""
    for k, (cx, cy) in enumerate(((0.26, 0.70), (0.50, 0.62), (0.74, 0.54))):
        a.add_patch(Circle((cx, cy), 0.095, facecolor=P1,
                           alpha=0.20 + 0.25 * k, zorder=1))
        for _ in range(5):
            th, rr = RNG.uniform(0, 2 * np.pi), 0.085 * np.sqrt(RNG.uniform(0, 0.8))
            blob(a, cx + rr * np.cos(th), cy + rr * np.sin(th), 0.025, P1)
    arrow(a, (0.20, 0.44), (0.86, 0.44), c=ACC[3], lw=1.8, ms=12)
    a.text(0.53, 0.375, "定向运动", ha="center", fontsize=8.4, color=ACC[3])
    a.text(0.5, 0.86, "必须有 ATP 循环，且 $k_{cat}[E] \\geq k_{off}$",
           ha="center", fontsize=8.4, color=ACC[3])
    a.text(0.5, 0.295, "任何 Kd 组合都做不到 —— 互易定理", ha="center",
           fontsize=8.0, color=ACC[6])


# --------------------------------------------------------------------------- #
PANELS = [
    (p_complex, "1:1 复合物", "价数 1×1：门关着", "PKA–AKAP（Kd 1 nM 也不成相）、Raf–MEK–ERK", "ok"),
    (p_filament, "一维纤维 / 链", "每类位点只有 1 个 → 分支因子 ≤ 1", "MAVS、ASC、RIPK–RHIM、DVL2–AXIN", "ok"),
    (p_network, "三维网络 / 液滴", "$(n_A\\!-\\!1)(n_B\\!-\\!1)>1$ 且位点浓度 ≥ Kd",
     "p62–polyUb、SPOP–DAXX、SynGAP–PSD-95、NPM1", "ok"),
    (p_membrane, "膜上二维凝聚", "二维化把局部位点浓度抬 ~133×",
     "LAT–Grb2–SOS1（需 18×）、Nephrin–NCK–N-WASP（需 21×）", "ok"),
    (p_template, "模板上的凝聚：长度 = 价数", "阈值长度随受体丰度反比移动",
     "cGAS–dsDNA（160–1592 bp）、G3BP1–mRNA、p62–Ub 链（3–10 个）", "ok"),
    (p_phospho, "磷酸化造价数", "$n = n_{max}\\cdot p$，$n_H \\approx 4$（无协同）",
     "LAT、Nephrin —— 修饰即开关", "ok"),
    (p_sorted, "多个互不相溶的体（cells）", "跨通路 $K_d$ > 组内 $K_d$ 的调和平均",
     "胞浆本底：5 个模块 → 5 个独立的体", "ok"),
    (p_coreshell, "核–壳 / 同心多层", "两个稠相共存；但包裹次序需要界面张力",
     "核仁 FBL 相 / NPM1 相", "warn"),
    (p_client, "客体富集 + 尺寸排阻", "尺寸本身就是排阻机制，不需特异识别",
     "LAT 留住 ZAP-70、排除 CD45（大）", "ok"),
    (p_allostery, "别构：配体改对方的价数", "平模型给出相反符号 → 需外挂别构层",
     "G3BP1–CAPRIN1–USP10", "warn"),
    (p_ageing, "液 → 凝胶 → 纤维（老化）", "动力学：成核、粗化、液–固转变",
     "FUS–RNA、ALS 突变加速转变", "no"),
    (p_motion, "会动的图案", "平衡态被定理禁止（$E=K+V$ 是 Lyapunov 函数）",
     "LAT 斑点向心运动 ← 肌动球蛋白逆向流", "no"),
]


def main():
    fig, ax = plt.subplots(4, 3, figsize=(15.5, 19.0))
    for a, (draw, title, cond, systems, ok) in zip(ax.ravel(), PANELS):
        draw(a)
        frame(a, title, cond, systems, ok)
    fig.suptitle("E12  真实蛋白体系能长出的 pattern：架构 → 形态 → 条件",
                 fontsize=17, y=0.997)
    for x, col, lab in ((0.315, ACC[0], "● 模型算得出"),
                        (0.500, ACC[3], "▲ 只能算一半"),
                        (0.690, ACC[6], "× 超出平衡框架")):
        fig.text(x, 0.9785, lab, ha="center", va="center", fontsize=10.5, color=col)
    fig.tight_layout(rect=[0, 0, 1, 0.972])
    save(fig, "e12_pattern_schematic.png")


if __name__ == "__main__":
    main()
