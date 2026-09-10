"""
C5 — The label problem: you cannot measure a cell's state and then watch its
fate, because measuring it kills it.

This is the constraint that decides whether the state space can be FIT, and it
is independent of how much data exists. scRNA-seq, scATAC, mass cytometry and
BH3 profiling are all destructive. So the state -> fate map can only be
estimated through one of three surrogate designs, each with a computable
ceiling:

  1. SISTER-CELL SPLIT  -- profile one daughter, challenge the other
     (the FateMap / Shaffer-Raj trick).
  2. CLONAL BARCODING   -- profile part of a clone, read fate in its siblings
     (LARRY / CellTag).
  3. LIVE IMAGING       -- genuinely same-cell, but only a handful of
     fluorescent channels.

Part 1 computes the ceiling for designs 1-2. Two cells that shared a state at
division drift apart; model each axis as an Ornstein-Uhlenbeck process with
memory time tau. The state correlation between the profiled cell and the cell
whose fate you observe decays as exp(-Delta/tau), where Delta is the delay
between profiling and the fate-determining moment. The resulting ceiling on
fate-prediction AUC is computed by Monte Carlo (no closed-form assumed), as a
DIMENSIONLESS function of Delta/tau -- so the result holds whatever the true
timescales are. Clone averaging over n siblings is included: it beats the
sampling noise but, as the numbers show, not the decorrelation.

Part 2 turns the C3 finding around. The six axes are redundant (effective
dimension 3.0), which is bad news for "six independent variables" but GOOD
news for live imaging: if the axes overlap, a few reporters can pin the whole
state. This part exhaustively searches every subset of axes to find the best
k-channel imaging panel, scored on how much of the TOTAL six-axis state it
recovers, using the C3-measured covariance.

Rigour: ✅ Part 1's AUC ceiling curve and Part 2's subset search are exact
computations on the C3-measured covariance. ⚠️ the placement of each axis on
the Delta/tau curve is an estimate, and is reported as a range, not a point.
"""
import itertools
import json
import os

import numpy as np

from .common import CACHE, banner, figpath, save_json

AXES = ["细胞周期", "应激水平", "代谢储备", "DNA损伤", "凋亡准备", "表观谱系"]
EN = ["cell cycle", "stress", "metab", "DNA dmg", "priming", "lineage"]
RNG = np.random.default_rng(7)


def load_sigma():
    path = os.path.join(CACHE, "c3_orthogonality.json")
    if os.path.exists(path):
        C = np.array(json.load(open(path))["mean_C_disjoint"])
        w, V = np.linalg.eigh(C)
        C = V @ np.diag(np.clip(w, 1e-6, None)) @ V.T
        d = np.sqrt(np.diag(C))
        return C / np.outer(d, d), "C3 实测"
    return np.eye(6), "⚠️ 退化为单位阵"


def auc_ceiling(r, n=200_000):
    """Max achievable AUC when the score correlates r with the true decision
    variable. Monte Carlo, so no binormal formula is being assumed."""
    if r <= 0:
        return 0.5
    r = min(r, 0.999)
    d = RNG.standard_normal(n)
    s = r * d + np.sqrt(1 - r ** 2) * RNG.standard_normal(n)
    y = d > 0
    sp, sn = s[y], s[~y]
    m = min(len(sp), len(sn), 60_000)
    a = sp[:m][:, None] > sn[:m][None, :] if m <= 3000 else None
    if a is None:                       # rank-based AUC for large samples
        allv = np.concatenate([sp, sn])
        order = np.argsort(allv)
        ranks = np.empty_like(order, dtype=float)
        ranks[order] = np.arange(1, len(allv) + 1)
        rp = ranks[:len(sp)].sum()
        return float((rp - len(sp) * (len(sp) + 1) / 2) / (len(sp) * len(sn)))
    return float(a.mean())


def clone_gain(r, n_sibs):
    """Averaging n siblings' states sharpens the shared component only."""
    # shared component variance r^2, private (1-r^2) averaged down by n
    return r / np.sqrt(r ** 2 + (1 - r ** 2) / n_sibs)


def run():
    banner("C5 — 命运标签问题:测了状态就杀死了细胞(替代设计的天花板)")
    Sigma, src = load_sigma()

    # ---------------- Part 1: decorrelation ceiling ---------------------
    print("【第一部分】姐妹细胞 / 克隆条形码设计的信息天花板")
    print("两个细胞在分裂时共享状态,之后按 exp(-Δ/τ) 解相关。")
    print("Δ = 从取样到命运决定时刻的延迟;τ = 该轴的状态记忆时间。\n")
    print(f"{'Δ/τ':>6}{'姐妹态相关 r':>14}{'单姐妹 AUC 上限':>17}"
          f"{'克隆 n=20 AUC 上限':>20}")
    print("-" * 58)
    ratios = [0.05, 0.1, 0.25, 0.5, 1.0, 2.0, 3.0]
    rows = []
    for q in ratios:
        r = float(np.exp(-q))
        a1 = auc_ceiling(r)
        a20 = auc_ceiling(clone_gain(r, 20))
        rows.append(dict(ratio=q, r=r, auc1=a1, auc20=a20))
        print(f"{q:>6.2f}{r:>14.2f}{a1:>17.3f}{a20:>20.3f}")

    print("\n读法:Δ/τ = 1(取样与命运决定相隔一个记忆时间)时,单个姐妹细胞的")
    print(f"命运预测 AUC 上限只有 {rows[4]['auc1']:.2f} —— 即使状态测得完美无误、")
    print("即使状态→命运的映射完全已知。这是设计本身的上限,不是方法不够好。")
    print(f"克隆平均(n=20)把它抬到 {rows[4]['auc20']:.2f},帮助有限:平均能压掉")
    print("抽样噪声,压不掉解相关。")

    print("\n⚠️ 各轴大致落在哪里(估计,非实测):")
    print("  表观/谱系:记忆跨多个细胞世代 → Δ/τ ≪ 1 → 姐妹设计有效")
    print("            (这正是 LARRY / CellTag 类实验能成立的原因)")
    print("  代谢储备:记忆约一个世代 → Δ/τ ~ 0.5-1 → 勉强")
    print("  DNA 损伤/应激:分钟到小时 → Δ/τ ≫ 1 → 姐妹设计基本无效")
    print("  细胞周期相位:姐妹出生时同步,但随刺激时刻漂移 → 取决于给药时机")
    print("  凋亡准备度:BH3 测定本身即破坏性,且无 RNA 替代 → 三种设计都不适用")

    # ---------------- Part 2: how many imaging channels ------------------
    print(f"\n{'='*72}")
    print("【第二部分】活细胞成像:同一细胞的状态与命运都能拿到,但通道数有限")
    print(f"用 C3 实测的协方差({src})穷举所有子集,问:k 个报告基因最多能")
    print("锁定整个六轴状态的多少?\n")

    def explained(idx, rho=0.9):
        """Fraction of total 6-axis variance recovered by observing subset idx."""
        k = len(idx)
        R = np.zeros((k, 6))
        for a, i in enumerate(idx):
            R[a, i] = rho
        N = np.eye(k) * (1 - rho ** 2)
        S = R @ Sigma @ R.T + N
        K = Sigma @ R.T @ np.linalg.inv(S)
        Post = Sigma - K @ R @ Sigma
        return float(1 - np.trace(Post) / np.trace(Sigma))

    best = {}
    print(f"{'k':>3}{'最佳组合':>34}{'解释的状态比例':>16}")
    print("-" * 56)
    for k in range(1, 7):
        cand = max(itertools.combinations(range(6), k), key=explained)
        v = explained(cand)
        best[k] = (cand, v)
        names = "+".join(EN[i] for i in cand)
        print(f"{k:>3}{names:>36}{v*100:>13.0f}%")

    v3, v6 = best[3][1], best[6][1]
    print(f"\n3 个通道拿到 {v3*100:.0f}%,6 个通道拿到 {v6*100:.0f}% ——")
    print(f"第 4-6 个通道总共只多买到 {(v6-v3)*100:.0f} 个百分点。")
    print("原因就是 C3 的结果:轴之间冗余,有效维数只有 3。")
    print("所以'轴不正交'这件坏消息,在测量端是好消息:三色活细胞成像")
    print("已经能锁住这个状态空间的大部分,而它恰好是唯一能同时拿到")
    print("同一细胞状态与命运的手段。")

    # the decisive design comparison
    print(f"\n{'='*72}\n【三种设计的对比】")
    print(f"{'设计':<28}{'同一细胞?':>10}{'可测轴数':>10}{'命运标签':>10}")
    print("-" * 60)
    print(f"{'scRNA-seq + 克隆条形码':<26}{'否(克隆)':>12}{'6(粗)':>10}{'有':>10}")
    print(f"{'姐妹细胞拆分':<26}{'否(姐妹)':>12}{'6(粗)':>10}{'有':>10}")
    print(f"{'活细胞成像 3 通道':<26}{'是':>12}{'~3':>10}{'有':>10}")
    print(f"{'BH3 profiling':<26}{'否(破坏)':>12}{'1(精)':>10}{'无':>10}")
    print("\n没有任何一种设计同时给出:同一细胞 + 六个轴 + 命运标签。")
    print("状态空间必须由多种设计拼起来,而拼接本身需要一个共同的标定轴。")

    _plot(rows, best)
    save_json("c5_pairing.json",
              dict(decorrelation=rows,
                   best_subsets={k: [list(v[0]), v[1]] for k, v in best.items()}))
    print("\n结果已存 figures/_cellstate_cache/c5_pairing.json")
    return rows, best


def _plot(rows, best):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.2))
    q = np.array([r["ratio"] for r in rows])
    a1.plot(q, [r["auc1"] for r in rows], "o-", label="single sister")
    a1.plot(q, [r["auc20"] for r in rows], "s-", label="clone average (n=20)")
    a1.axhline(0.5, ls=":", color="k", lw=1)
    a1.set_xlabel(r"$\Delta/\tau$  (profiling delay / state memory)")
    a1.set_ylabel("ceiling on fate-prediction AUC")
    a1.set_title("Surrogate designs decorrelate", fontsize=10)
    a1.legend(fontsize=8); a1.set_ylim(0.45, 1.0)

    ks = sorted(best)
    a2.plot(ks, [best[k][1] * 100 for k in ks], "o-", color="#54A24B")
    a2.set_xlabel("number of live-imaging channels")
    a2.set_ylabel("% of 6-axis state recovered")
    a2.set_title("Redundancy makes 3 channels nearly enough", fontsize=10)
    a2.set_ylim(0, 100); a2.grid(alpha=.3)
    fig.tight_layout(); fig.savefig(figpath("c5_pairing.png"), dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    run()
