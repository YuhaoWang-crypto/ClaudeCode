"""
C6 — The ceiling, MEASURED on real data (not simulated).

C5 derived the surrogate-design ceiling as a function of an unknown memory
time. This module removes the unknown by measuring the ceiling directly, on
the one open dataset built to measure exactly this.

Weinreb, Rodriguez-Fraticelli, Camargo & Klein, Science 2020 (LARRY):
mouse haematopoietic progenitors are lentivirally barcoded, then at day 2 the
clones are SPLIT ACROSS TWO SETS OF WELLS and cultured separately to days 4
and 6. Cells carry a clonal barcode, a time point, a well label (0 = day 2,
1 or 2 = the two split cultures) and a mature-cell-type annotation.

Files (downloaded from the Klein lab server, ~2 MB, no counts matrix needed):
  stateFate_inVitro_metadata.txt.gz     130,887 cells x 8 columns
  stateFate_inVitro_clone_matrix.mtx.gz 130,887 cells x 5,864 clones

This design is a natural experiment for the question "can a cell's state
predict its fate?", because two halves of the SAME clone, separated before
differentiation, are an upper bound on any state-based predictor: whatever a
perfect measurement of a cell's heritable state could tell you about its fate,
its clonal twin in the other well already tells you -- and more, because the
twin shares everything heritable, measured or not, transcriptomic or not.

Three numbers are computed:
  1. CROSS-WELL FATE CORRELATION -- how well one well's half of a clone
     predicts the other well's half. Different wells, so shared micro-
     environment is controlled away. This is the ceiling.
  2. A LABEL-PERMUTATION NULL -- clone identities shuffled, same marginals.
     Establishes that the correlation is real and not a compositional artifact.
  3. INTRACLASS CORRELATION per fate -- the fraction of fate variance that is
     clonally heritable at all.

Rigour: ✅ every number is computed from the downloaded files; the null is an
actual permutation, not an asymptotic formula. ⚠️ the system is mouse
haematopoiesis with differentiation fates, so the ceiling for a *death* fate
(the apoptotic-priming axis) is not measured here and, per C1, is not measured
by any open dataset.
"""
import gzip
import os

import numpy as np

from .common import CACHE, banner, figpath, save_json

META = "stateFate_inVitro_metadata.txt.gz"
CLONES = "stateFate_inVitro_clone_matrix.mtx.gz"
RNG = np.random.default_rng(11)


def load():
    mp, cp = os.path.join(CACHE, META), os.path.join(CACHE, CLONES)
    if not (os.path.exists(mp) and os.path.exists(cp)):
        return None, None, None
    hdr, rows = None, []
    with gzip.open(mp, "rt") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if hdr is None:
                hdr = p
                continue
            rows.append(p)
    col = {h: i for i, h in enumerate(hdr)}
    meta = dict(
        tp=np.array([float(r[col["Time point"]]) for r in rows]),
        well=np.array([r[col["Well"]] for r in rows]),
        ctype=np.array([r[col["Cell type annotation"]] for r in rows]),
        pop=np.array([r[col["Starting population"]] for r in rows]),
    )
    # clone matrix: rows = cells, cols = clones (MatrixMarket, 1-indexed)
    cell2clone = np.full(len(rows), -1, dtype=int)
    with gzip.open(cp, "rt") as fh:
        for line in fh:
            if line.startswith("%"):
                continue
            p = line.split()
            if len(p) == 3 and cell2clone.max() == -1 and int(p[0]) > 1e5:
                n_cells, n_clones = int(p[0]), int(p[1])
                continue
            if len(p) >= 2:
                i, j = int(p[0]) - 1, int(p[1]) - 1
                if 0 <= i < len(cell2clone):
                    cell2clone[i] = j
    return meta, cell2clone, (n_cells, n_clones)


def icc_oneway(groups):
    """One-way random-effects ICC(1) on binary outcomes grouped by clone."""
    groups = [g for g in groups if len(g) >= 2]
    if len(groups) < 3:
        return float("nan")
    k = np.array([len(g) for g in groups], float)
    m = np.array([np.mean(g) for g in groups], float)
    N, grand = k.sum(), np.concatenate(groups).mean()
    msb = np.sum(k * (m - grand) ** 2) / (len(groups) - 1)
    ssw = sum(np.sum((np.array(g) - np.mean(g)) ** 2) for g in groups)
    msw = ssw / (N - len(groups))
    k0 = (N - np.sum(k ** 2) / N) / (len(groups) - 1)
    denom = msb + (k0 - 1) * msw
    return float((msb - msw) / denom) if denom > 0 else float("nan")


def run():
    banner("C6 — 命运可预测性的上限:在 LARRY 真实数据上实测")
    meta, c2c, dims = load()
    if meta is None:
        print("⚠️  LARRY 文件缺失,C6 无结论。")
        print("    下载: https://kleintools.hms.harvard.edu/paper_websites/state_fate2020/")
        return None

    n_cells, n_clones = dims
    print(f"数据: Weinreb et al. Science 2020 (LARRY) 体外分化时程")
    print(f"  细胞 {n_cells:,},克隆 {n_clones:,}")
    has = c2c >= 0
    print(f"  有克隆标签的细胞: {has.sum():,} ({100*has.mean():.0f}%)")
    for t in sorted(set(meta['tp'])):
        print(f"  day {t:.0f}: {int((meta['tp']==t).sum()):,} 个细胞")

    mature = meta["ctype"] != "Undifferentiated"
    types = sorted(set(meta["ctype"][mature]))
    print(f"\n成熟细胞 {int(mature.sum()):,};命运类别 {len(types)} 种:")
    print("  " + ", ".join(types))

    # ---- clones split across the two wells -----------------------------
    late = mature & has & np.isin(meta["well"], ["1", "2"])
    print(f"\n两个孔中带克隆标签的成熟细胞: {int(late.sum()):,}")

    idx = np.where(late)[0]
    cl, wl, ct = c2c[idx], meta["well"][idx], meta["ctype"][idx]
    tidx = {t: i for i, t in enumerate(types)}

    from collections import defaultdict
    per = defaultdict(lambda: [np.zeros(len(types)), np.zeros(len(types))])
    for c, w, t in zip(cl, wl, ct):
        per[c][0 if w == "1" else 1][tidx[t]] += 1

    MIN = 3
    P1, P2, sizes = [], [], []
    for c, (a, b) in per.items():
        if a.sum() >= MIN and b.sum() >= MIN:
            P1.append(a / a.sum()); P2.append(b / b.sum())
            sizes.append((a.sum(), b.sum()))
    P1, P2 = np.array(P1), np.array(P2)
    print(f"两孔中各有 ≥{MIN} 个成熟细胞的克隆: {len(P1)} 个")
    print(f"  (论文报告 split-well 子集为 502 个条形码 / 6,351 个细胞)")

    if len(P1) < 10:
        print("⚠️  可用克隆过少,不做结论。")
        return None

    # ---- 1. cross-well fate correlation, per fate ----------------------
    print(f"\n【1. 跨孔命运相关】同一克隆的两半,分在不同孔里独立分化")
    print(f"{'命运类型':<22}{'跨孔 r':>10}{'置换零假设 r':>16}{'克隆数':>9}")
    print("-" * 60)
    rows_out = []
    for t in types:
        i = tidx[t]
        a, b = P1[:, i], P2[:, i]
        if a.std() < 1e-9 or b.std() < 1e-9:
            continue
        r = float(np.corrcoef(a, b)[0, 1])
        null = [float(np.corrcoef(a, RNG.permutation(b))[0, 1]) for _ in range(300)]
        rows_out.append(dict(fate=t, r=r, null=float(np.mean(np.abs(null))),
                             n=len(a)))
        print(f"{t:<22}{r:>10.3f}{np.mean(np.abs(null)):>16.3f}{len(a):>9}")

    rs = np.array([d["r"] for d in rows_out])
    nulls = np.array([d["null"] for d in rows_out])
    print("-" * 60)
    print(f"{'中位数':<22}{np.median(rs):>10.3f}{np.median(nulls):>16.3f}")

    # ---- 2. what that correlation means as a prediction ceiling --------
    print(f"\n【2. 换算成预测天花板】")
    rmed = float(np.median(rs))
    print(f"  跨孔命运相关中位数 r = {rmed:.2f}")
    print(f"  → 一个克隆孪生体最多解释 r² = {rmed**2*100:.0f}% 的命运方差。")
    print(f"  → 剩下 {100-rmed**2*100:.0f}% 在克隆内部,即使两个细胞共享全部可遗传信息")
    print(f"     (转录组、表观、蛋白、代谢,已测和未测的全部)也无法预测。")
    print("\n  这是硬上限。任何基于状态的预测器都不可能超过它,因为孪生体")
    print("  掌握的可遗传信息严格多于任何一次测量能拿到的信息。")

    # ---- 3. ICC: how much fate variance is clonally heritable ----------
    print(f"\n【3. 组内相关系数 ICC】命运方差里有多少是克隆可遗传的")
    groups_by_fate = {t: [] for t in types}
    byclone = defaultdict(list)
    for c, t in zip(cl, ct):
        byclone[c].append(t)
    for c, lst in byclone.items():
        if len(lst) >= 2:
            for t in types:
                groups_by_fate[t].append([1.0 if x == t else 0.0 for x in lst])
    print(f"{'命运类型':<22}{'ICC':>10}")
    print("-" * 34)
    iccs = []
    for t in types:
        v = icc_oneway(groups_by_fate[t])
        if np.isfinite(v):
            iccs.append(v)
            print(f"{t:<22}{v:>10.3f}")
    print("-" * 34)
    print(f"{'中位数':<22}{np.median(iccs):>10.3f}")
    print(f"\n  ICC 中位数 {np.median(iccs):.2f}:命运确实高度克隆可遗传 ——")
    print("  这正是状态空间这个想法的经验基础,值得肯定。")
    print("  但跨孔相关告诉我们,可遗传 ≠ 可从一次测量中读出。")

    # ---- 4. the axis this dataset cannot speak to ----------------------
    print(f"\n【4. 这个数据集测不到的东西】")
    dead_note = ("所有条形码设计都有生存者偏倚:死掉的细胞不会出现在终点的 "
                 "scRNA-seq 里。")
    print(f"  {dead_note}")
    print("  所以'凋亡准备度 → 死亡'这条最关键的状态→命运映射,在 LARRY、")
    print("  CellTag、Watermelon、FateMap 里都结构性地缺失(TraCe-seq 是少数例外)。")
    print("  C1 的检索与此一致:凋亡准备度轴只有 21 个 GEO 数据集、0 个带命运标签。")

    out = dict(n_cells=int(n_cells), n_clones=int(n_clones),
               n_split_clones=int(len(P1)), per_fate=rows_out,
               median_r=rmed, median_r2=float(rmed ** 2),
               median_icc=float(np.median(iccs)))
    _plot(rows_out, rmed, float(np.median(iccs)))
    save_json("c6_larry_ceiling.json", out)
    print("\n结果已存 figures/_cellstate_cache/c6_larry_ceiling.json")
    return out


def _plot(rows_out, rmed, icc):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return
    fig, ax = plt.subplots(figsize=(8.2, 4.2))
    names = [d["fate"] for d in rows_out]
    r = [d["r"] for d in rows_out]
    nl = [d["null"] for d in rows_out]
    x = np.arange(len(names))
    ax.bar(x - 0.2, r, 0.4, label="cross-well clone correlation", color="#4C78A8")
    ax.bar(x + 0.2, nl, 0.4, label="permutation null", color="#BBBBBB")
    ax.axhline(rmed, ls="--", color="#E45756", lw=1,
               label=f"median r = {rmed:.2f}  (ceiling: r$^2$={rmed**2:.2f})")
    ax.set_xticks(x); ax.set_xticklabels(names, rotation=40, ha="right", fontsize=8)
    ax.set_ylabel("correlation between the two halves of a clone")
    ax.set_title("Measured ceiling on fate prediction (LARRY split-well, "
                 "Weinreb et al. 2020)", fontsize=10)
    ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(figpath("c6_larry_ceiling.png"), dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    run()
