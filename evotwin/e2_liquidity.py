"""
E2 — 位点流动性 (site liquidity)：把"靶位点 turnover"写成 3'UTR 上的
     标记点过程 (marked point process)，并分离出两件本来被混为一谈的事：

       (a) 边(edge)是否保留       —— 功能层
       (b) 位点坐标是否保留        —— 序列层（conservation scan 看的就是这个）

核心发现是一个无量纲数：

       Λ_liq  =  λ_gain / λ_loss(单位点)  =  L · 4^-K

即"中性稳态期望位点数"。它同时是补偿速度与丢失速度之比：
   Λ_liq >> 1  → site-liquid   ：位点可被就地快速替换，坐标不保守但边保守
   Λ_liq << 1  → site-crystalline：一旦丢失无法补偿，选择只能把坐标冻住

时间单位一律取【中性分歧度 d（每位点替换数）】，与代数/物种无关。
"""
import numpy as np
from scipy.linalg import expm

K = 7
NE = 1.0e4
D_HUMAN_MOUSE = 0.5          # 人鼠中性分歧 ~0.5 subs/site


def h(S):
    """相对中性的固定速率倍数（Kimura）：h(S)=S/(1-e^-S)，S=4N_e s。"""
    S = np.asarray(S, dtype=float)
    out = np.ones_like(S)
    nz = np.abs(S) > 1e-9
    out[nz] = S[nz] / (1 - np.exp(-S[nz]))
    return out if out.shape else float(out)


def hs(S):
    return float(h(np.array(S)))


def lam_gain(L, k=K):
    """每单位分歧度 d 的位点获得率：(3k·4^-k·L)·(1/3)。"""
    return L * 3 * k * (0.25 ** k) / 3.0


def lam_loss_per_site(k=K):
    """每单位 d 的单位点中性丢失率：k 个碱基任一替换即失效。"""
    return float(k)


def liquidity(L, k=K):
    return lam_gain(L, k) / lam_loss_per_site(k)


def build_generator_rates(lg, ll, alpha, n_star=1, n_max=4, ne=NE):
    """通用版：直接给定 λ_gain / λ_loss(单位点)（单位：每单位中性分歧度 d）。

    状态 = (n 个位点, 祖先位点是否仍在)。稳定化选择 w(n)=exp(-alpha (n-n*)^2)。
    """
    states = [(0, 0)] + [(n, o) for n in range(1, n_max + 1) for o in (0, 1)]
    idx = {s: i for i, s in enumerate(states)}
    Q = np.zeros((len(states), len(states)))

    def w(n):
        return np.exp(-alpha * (n - n_star) ** 2)

    for (n, o) in states:
        i = idx[(n, o)]
        if n < n_max:                                     # 获得一个新位点（坐标必为新）
            s = w(n + 1) / w(n) - 1.0
            Q[i, idx[(n + 1, o)]] += lg * hs(4 * ne * s)
        if n > 0:                                         # 丢失一个位点
            s = w(n - 1) / w(n) - 1.0
            rate = n * ll * hs(4 * ne * s)
            if o == 1:                                    # 丢的是祖先位点的概率 1/n
                Q[i, idx[(n - 1, 1)] if n - 1 > 0 else idx[(0, 0)]] += rate * (n - 1) / n
                Q[i, idx[(n - 1, 0)] if n - 1 > 0 else idx[(0, 0)]] += rate * (1 / n)
            else:
                Q[i, idx[(n - 1, 0)] if n - 1 > 0 else idx[(0, 0)]] += rate
    np.fill_diagonal(Q, -Q.sum(axis=1))
    return Q, states, idx


def build_generator(L, alpha, k=K, **kw):
    """便捷版：由 UTR 长度与 k-mer 长度推出中性速率（均匀 25% 组成假设）。"""
    return build_generator_rates(lam_gain(L, k), lam_loss_per_site(k), alpha, **kw)


def branch_stats_rates(lg, ll, alpha, d, **kw):
    """从 (n=1, 祖先位点在) 出发，走过分歧度 d 之后的边保留率 / 坐标保留率。"""
    Q, states, idx = build_generator_rates(lg, ll, alpha, **kw)
    p = np.zeros(len(states))
    p[idx[(1, 1)]] = 1.0
    pt = p @ expm(Q * d)
    edge = sum(pt[idx[(n, o)]] for (n, o) in states if n > 0)
    pos = sum(pt[idx[(n, o)]] for (n, o) in states if o == 1)
    return edge, pos


def branch_stats(L, alpha, d, k=K, **kw):
    return branch_stats_rates(lam_gain(L, k), lam_loss_per_site(k), alpha, d, **kw)


def run():
    print("=" * 78)
    print("E2 — 位点流动性 Λ_liq = L·4^-K：liquid（可补偿）vs crystalline（须冻结）")
    print("=" * 78)

    print("\n[1] 无量纲流动性数（中性）")
    print(f"{'UTR 长度':>10} {'Λ_liq(8mer,K=8)':>17} {'Λ_liq(7mer,K=7)':>17} "
          f"{'Λ_liq(6mer,K=6)':>17} {'判定(7mer)':>12}")
    for L in (500, 1000, 2000, 5000, 10000, 20000):
        l7 = liquidity(L, 7)
        verdict = "liquid" if l7 > 1 else ("marginal" if l7 > 0.3 else "crystalline")
        print(f"{L:>10} {liquidity(L,8):>17.3f} {l7:>17.3f} "
              f"{liquidity(L,6):>17.3f} {verdict:>12}")
    print("  → 7mer 边要进入 liquid 区需要 L > 4^7 = 16.4 kb 的 3'UTR；")
    print("    6mer 边只要 L > 4.1 kb。神经元基因的超长 3'UTR 正好落在这一侧。")

    print("\n[2] 走过人鼠分歧 d=0.5 后：边保留 vs 坐标保留（选择强度 α 扫描）")
    print(f"{'α':>8} {'4N_eα':>8} " + " ".join(f"{'L=%d' % L:>20}" for L in (2000, 20000)))
    print(f"{'':>8} {'':>8} " + " ".join(f"{'边保留/坐标保留':>20}" for _ in range(2)))
    for alpha in (0.0, 1e-5, 1e-4, 1e-3):
        row = f"{alpha:>8.0e} {4*NE*alpha:>8.1f} "
        for L in (2000, 20000):
            e, p = branch_stats(L, alpha, D_HUMAN_MOUSE)
            row += f" {e:>9.3f} / {p:<9.3f}"
        print(row)

    print("\n[3] 判读：conservation scan 的漏检率 = 1 - 坐标保留/边保留")
    print(f"{'α':>8} {'L':>7} {'边保留':>9} {'坐标保留':>10} {'漏检率':>9}")
    for alpha in (1e-5, 1e-4, 1e-3):
        for L in (2000, 20000):
            e, p = branch_stats(L, alpha, D_HUMAN_MOUSE)
            print(f"{alpha:>8.0e} {L:>7} {e:>9.3f} {p:>10.3f} {1 - p / e:>9.1%}")
    print("\n  → 这就是【site-liquid edge】：功能上高度保守、坐标上完全不保守。")
    print("    任何以'同源坐标 seed 保守'为筛选条件的方法，在长 UTR + 中等选择区")
    print("    会系统性丢掉这一整类真边 —— 而它们恰恰是剂量敏感基因的主力。")
    return None


if __name__ == "__main__":
    run()
