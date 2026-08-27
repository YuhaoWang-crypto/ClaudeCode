"""
E4 — "靶位点回避" (target avoidance) 的可检验性：一个基因长期【不获得】某个
     miRNA 的位点，能不能被统计上确认？需要多少物种？

这是方案里 Avoidance score 的严格版本。做法：把位点数 n 在系统树上跑一个
生灭 CTMC，比较两个模型的似然

    H0: 中性获得（速率 λ_g）
    H1: 获得受负选择（速率 λ_g · h(-4N_e s_a)）

数据 = "全部 N 个物种都没有该位点"。输出：达到 ΔlnL > 2（≈ Bayes factor 7.4）
所需的最少物种数，以及不同 UTR 长度下的判定边界。

关键洞察：短 UTR 的基因【本来就】大概率没有位点（Λ_liq 很小），所以"没有位点"
毫无信息量。回避只有在【中性期望位点数足够大】时才可检测 —— 这把 avoidance
scan 自动限制到长 3'UTR / 高 Λ_liq 的基因上。
"""
import numpy as np
from scipy.linalg import expm

from .e2_liquidity import lam_gain, lam_loss_per_site, liquidity, h, hs

NE = 1.0e4
N_MAX = 6


def generator(L, s_avoid=0.0, k=7, n_max=N_MAX, ne=NE):
    """位点数 n 的生灭 CTMC。s_avoid>0 表示"获得位点是有害的"。"""
    lg, ll = lam_gain(L, k), lam_loss_per_site(k)
    g = lg * hs(-4 * ne * s_avoid)          # 获得被压制
    l = ll * hs(+4 * ne * s_avoid)          # 丢失被加速
    Q = np.zeros((n_max + 1, n_max + 1))
    for n in range(n_max + 1):
        if n < n_max:
            Q[n, n + 1] = g
        if n > 0:
            Q[n, n - 1] = n * l
    np.fill_diagonal(Q, -Q.sum(axis=1))
    return Q


def stationary(Q):
    w, v = np.linalg.eig(Q.T)
    p = np.real(v[:, np.argmin(np.abs(w))])
    return p / p.sum()


def loglik_all_absent(L, s_avoid, n_species, branch_len):
    """星形树：根取稳态，各物种独立演化 branch_len，观测到全部为 0 的对数似然。"""
    Q = generator(L, s_avoid)
    P = expm(Q * branch_len)
    pi = stationary(Q)
    return np.log(np.sum(pi * P[:, 0] ** n_species))


def min_species(L, s_avoid, branch_len, thresh=2.0, cap=400):
    for n in range(2, cap + 1):
        d = loglik_all_absent(L, s_avoid, n, branch_len) - \
            loglik_all_absent(L, 0.0, n, branch_len)
        if d > thresh:
            return n, d
    return None, None


def run():
    print("=" * 78)
    print("E4 — 回避 (avoidance) 的统计可检验性")
    print("=" * 78)

    b = 0.25   # 单物种根到叶的中性分歧（哺乳动物量级）
    print(f"\n[1] 检出'回避'所需最少物种数（星形树, 根到叶分歧 b={b}, ΔlnL>2）")
    print(f"{'UTR 长度':>10} {'Λ_liq':>9} {'P0(中性)':>10} " +
          " ".join(f"{'s_a=%.0e' % s:>10}" for s in (1e-4, 1e-3, 1e-2)))
    for L in (1000, 2000, 5000, 10000, 20000, 50000):
        Q0 = generator(L, 0.0)
        p0 = stationary(Q0)[0]
        row = f"{L:>10} {liquidity(L):>9.3f} {p0:>10.3f} "
        for s_a in (1e-4, 1e-3, 1e-2):
            n, _ = min_species(L, s_a, b)
            row += f" {('%d' % n) if n else '>400':>10}"
        print(row)

    print("\n  → 读法：Λ_liq 越大（UTR 越长），'全都没有位点'越反常，需要的物种越少。")
    print("    L=2 kb 的基因即使真的在回避，也需要 ~20 个独立谱系才够。")
    print("  ⚠ 星形树是【功效上界】：真实系统树共享内部枝，独立信息量更少，")
    print("    实际所需物种数应按有效独立枝数（~树总长/b）折算后再放大。")

    print("\n[2] 固定 100 个物种时，能检出的最弱回避强度")
    print(f"{'UTR 长度':>10} {'ΔlnL @ s_a=1e-4':>18} {'@1e-3':>10} {'@1e-2':>10} {'判定':>12}")
    for L in (2000, 10000, 20000, 50000):
        ds = [loglik_all_absent(L, s, 100, b) - loglik_all_absent(L, 0.0, 100, b)
              for s in (1e-4, 1e-3, 1e-2)]
        verdict = "可检出" if max(ds) > 2 else "信息不足"
        print(f"{L:>10} {ds[0]:>18.2f} {ds[1]:>10.2f} {ds[2]:>10.2f} {verdict:>12}")

    print("\n[3] 与'边保守'的对称性")
    print("  同一套 CTMC 同时给出两种极端候选：")
    print("   · 保守边 (conserved edge)：观测位点数 >> 中性期望 → 正选择/纯化选择")
    print("   · 回避边 (avoided edge)  ：观测位点数 == 0 且 << 中性期望 → 负选择")
    print("  二者用【同一个 Λ_liq 零模型】校准，因此可以放进同一张 z-score 表，")
    print("  这是现有 TargetScan/conservation 流程做不到的（它们没有 gain 的零模型）。")
    return None


if __name__ == "__main__":
    run()
