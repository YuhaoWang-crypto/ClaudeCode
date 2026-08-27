"""
E1 — 把 EvoTwin 里抽象的 U_eff 换成一个可以真正算出来的数：
     "一个 3'UTR 在一代里获得某个特定 miRNA family 的 7mer 靶位点的概率"。

关键点：U_eff 不是 μ × UTR 长度。真正的 functional mutational target size 是
"距离目标 7mer 只差一个碱基的位置数" —— 这是一个纯组合量，可以精确算出来。

L_eff(7mer, L) = L × (7mer 的单突变邻居数) × 4^-7
               = L × 21 × 4^-7
且每个这样的位置只有 1/3 的替换方向能真正生成该 7mer。

输出：U_gain、U_loss、Θ=2N_e·U、Λ_success、E[T_wait]，以及全基因组"调控创新速率"。
"""
import math

MU = 1.25e-8          # 人类生殖系每 bp 每代替换率 (Garcia 2026 ~1.30e-8, Kong ~1.2e-8)
NE = 1.0e4            # 人类长期有效群体大小量级
GEN_YEARS = 25.0
K = 7                 # 7mer-m8 seed


def n_one_off_positions(utr_len, k=K):
    """UTR 中"单个替换即可生成指定 k-mer"的期望位置数（随机 25% 组成近似）。

    k-mer 的单突变邻居有 3k 个；其中 k 个位置各有 3 个替换方向，
    但只有 1 个方向命中目标 → 有效速率 = (3k 个邻居序列) × 4^-k × (1/3 的方向命中率)。
    """
    n_neighbors = 3 * k
    p_neighbor = n_neighbors * (0.25 ** k)      # 该位置正好是"差一个碱基"的概率
    return utr_len * p_neighbor


def u_gain(utr_len, mu=MU, k=K):
    """每代获得 >=1 个新位点的概率（Poisson 强度）。"""
    return n_one_off_positions(utr_len, k) * (mu / 3.0)


def u_loss_per_site(mu=MU, k=K):
    """单个已有位点每代被破坏的概率：k 个碱基任一被替换即失效。"""
    return k * mu


def expected_neutral_sites(utr_len, k=K):
    """中性生灭稳态期望位点数 = λ_gain/λ_loss，应等于组成学期望 L×4^-k。"""
    return u_gain(utr_len, k=k) / u_loss_per_site(k=k)


def p_fix(s, ne=NE):
    """二倍体 Kimura 固定概率（新突变初始频率 1/(2N)）。"""
    if abs(s) < 1e-12:
        return 1.0 / (2 * ne)
    return (1 - math.exp(-2 * s)) / (1 - math.exp(-4 * ne * s))


def waiting_time(utr_len, s, ne=NE, mu=MU):
    ug = u_gain(utr_len, mu)
    theta = 2 * ne * ug
    lam = theta * p_fix(s, ne)
    return dict(U_gain=ug, Theta=theta, P_fix=p_fix(s, ne),
                Lambda=lam, T_wait_gen=1 / lam, T_wait_My=1 / lam * GEN_YEARS / 1e6)


def run():
    print("=" * 78)
    print("E1 — 靶位点获得的功能突变供给 U_eff（可计算版）")
    print("=" * 78)

    print("\n[1] functional mutational target size —— 纯组合量")
    print(f"{'UTR 长度':>10} {'单突变可达位置数':>18} {'中性稳态期望位点数':>20} {'组成学核对 L·4^-7':>18}")
    for L in (500, 1000, 2000, 5000, 10000):
        print(f"{L:>10} {n_one_off_positions(L):>18.3f} "
              f"{expected_neutral_sites(L):>20.4f} {L * 0.25 ** K:>18.4f}")
    print("  → 两列相等是内部一致性检验：生灭稳态 λ_g/λ_l 必须回到序列组成期望。")

    print("\n[2] U_gain / U_loss（每代，μ=1.25e-8）")
    for L in (1000, 2000, 5000):
        print(f"  L={L:>5}bp   U_gain={u_gain(L):.3e}/gen   "
              f"U_loss(单位点)={u_loss_per_site():.3e}/gen   "
              f"比值={u_gain(L)/u_loss_per_site():.4f}")

    print("\n[3] 等待一条【特定】新边建立的时间（L=2000bp）")
    print(f"{'s':>8} {'Θ=2N_eU':>12} {'P_fix':>10} {'Λ/gen':>12} {'T_wait(代)':>14} {'T_wait(百万年)':>14}")
    for s in (0.0, 1e-4, 1e-3, 1e-2):
        r = waiting_time(2000, s)
        print(f"{s:>8.0e} {r['Theta']:>12.3e} {r['P_fix']:>10.3e} "
              f"{r['Lambda']:>12.3e} {r['T_wait_gen']:>14.3e} {r['T_wait_My']:>14.1f}")

    print("\n[4] 反演：候选边里有多大比例真的被选择青睐？")
    # 一旦 4N_e·s >> 1，位点一旦固定就几乎不再被中性丢失（P_fix(-s)~0），
    # 所以在时间窗 T 内"有益候选边"的实现比例 ≈ 1 - exp(-Λ·T)。
    n_fam, n_genes, L, s = 300, 20000, 2000, 1e-3
    n_cand = n_fam * n_genes
    lam_edge = 2 * NE * u_gain(L) * p_fix(s)
    T_gen = 90e6 / GEN_YEARS                     # 人鼠分歧 ~90 My
    realized = 1 - math.exp(-lam_edge * T_gen)
    n_obs = 45000                                # Friedman 2009: ~45k 保守靶位点量级
    f_ben = n_obs / (realized * n_cand)
    print(f"  候选边总数        = {n_cand:.0f}   (300 family × 20k gene)")
    print(f"  单边 Λ            = {lam_edge:.3e}/gen → 单边平均等待 "
          f"{1/lam_edge*GEN_YEARS/1e6:.0f} 百万年")
    print(f"  90 My 内实现比例   = {realized:.2f}  (有益边中已被抓到的份额)")
    print(f"  观测保守位点数     ≈ {n_obs}")
    print(f"  → 推得【有益候选边比例 f_ben ≈ {f_ben:.3f}】，即约 {f_ben*100:.1f}% 的"
          f"(miRNA family × 基因) 组合是受选择的")
    print("\n  → 结论：单条【指定】边极度 mutation-limited（10^7-10^8 年量级），")
    print("    网络之所以看起来'布满'靶位点，是因为候选空间有 6×10^6 之大；")
    print("    而真正受选择的只有百分之一量级 —— 其余是随机 seed match（假阳性的来源）。")
    return dict(u_gain_2k=u_gain(2000), lam_edge=lam_edge, f_ben=f_ben)


if __name__ == "__main__":
    run()
