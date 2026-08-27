"""
E3 — 把"遗传图谱"真正接进模型：两位点（miRNA 位点 × 靶位点）的
     重组–上位性阈值，以及由此推出的【共依赖模块的最大遗传距离】。

原方案里 r 只出现在符号表。但对任何"单独有害/中性、组合才有利"(sign epistasis)
的调控共依赖，重组是主要的破坏力：稀有的 11 单倍型在 00 背景中每代以 ~r 的速率
被拆散。因此建立条件是

        s_c  >  r        （重组阈值）
        4N_e s_c > 1     （漂变阈值）

把 r 用重组图谱换算成物理距离（人类平均 ~1.2 cM/Mb），就得到一个可以直接
拿去检验的预言：共依赖模块的成员必须落在某个 kb 尺度的窗口内。

模型：单倍型两位点选择–重组确定性递推（标准配子模型）。
"""
import numpy as np

CM_PER_MB = 1.2          # 人类常染色体平均重组率
NE = 1.0e4


def iterate(x, w, r, n_gen):
    """x=[x00,x01,x10,x11]，w 同序；选择后重组。返回 x11 轨迹。"""
    traj = np.empty(n_gen + 1)
    traj[0] = x[3]
    for t in range(n_gen):
        xs = w * x
        xs /= xs.sum()
        D = xs[0] * xs[3] - xs[1] * xs[2]
        x = np.array([xs[0] - r * D, xs[1] + r * D, xs[2] + r * D, xs[3] - r * D])
        x = np.clip(x, 0.0, None)
        x /= x.sum()
        traj[t + 1] = x[3]
    return traj


def invasion_rate(s_single, s_combo, r, x0=1e-8, n_gen=50):
    """稀有 11 单倍型的每代对数增长率（数值）。解析上应为 s_combo - r。"""
    w = np.array([1.0, 1 + s_single, 1 + s_single, 1 + s_combo])
    x = np.array([1 - x0 - 2e-10, 1e-10, 1e-10, x0])
    traj = iterate(x, w, r, n_gen)
    return (np.log(traj[-1]) - np.log(traj[0])) / n_gen


def establishes(s_single, s_combo, r):
    return invasion_rate(s_single, s_combo, r) > 0


def critical_r(s_single, s_combo, lo=1e-8, hi=0.5, iters=60):
    """二分求 invasion_rate = 0 的 r。"""
    if invasion_rate(s_single, s_combo, lo) <= 0:
        return 0.0
    if invasion_rate(s_single, s_combo, hi) > 0:
        return hi
    for _ in range(iters):
        mid = np.sqrt(lo * hi)
        if invasion_rate(s_single, s_combo, mid) > 0:
            lo = mid
        else:
            hi = mid
    return np.sqrt(lo * hi)


def r_to_kb(r, cm_per_mb=CM_PER_MB):
    """r(每代重组概率) → 物理距离 kb。r ≈ cM/100，1 cM = 1/cm_per_mb Mb。"""
    return (r * 100.0 / cm_per_mb) * 1000.0


def run():
    print("=" * 78)
    print("E3 — 重组–上位性阈值：共依赖模块能相隔多远？")
    print("=" * 78)

    print("\n[1] 数值求出的临界重组率 r*（sign epistasis：单独有害、组合有利）")
    print(f"{'s_single':>10} {'s_combo':>10} {'r* (数值)':>12} {'s_combo(理论)':>14} "
          f"{'最大物理距离':>16}")
    for s_single in (-0.002, -0.0005, 0.0):
        for s_combo in (1e-4, 1e-3, 1e-2, 5e-2):
            rc = critical_r(s_single, s_combo)
            kb = r_to_kb(rc)
            dist = f"{kb:,.0f} kb" if kb < 1e6 else "无限制(同染色体外)"
            print(f"{s_single:>10.4f} {s_combo:>10.0e} {rc:>12.3e} {s_combo:>14.0e} "
                  f"{dist:>16}")
    print("  → r* 与 s_combo 数值一致：阈值就是 r < s_combo。")

    print("\n[2] 换算成'共依赖模块的最大跨度'（人类 1.2 cM/Mb）")
    print(f"{'上位性增益 s_ε':>14} {'r*':>10} {'最大跨度':>14} {'现实对照':>34}")
    ref = {1e-5: "miRNA 簇内间距 (<10 kb)",
           1e-4: "单个基因/内含子尺度",
           1e-3: "TAD 内 (~100 kb)",
           1e-2: "染色体臂局部 (~1 Mb)",
           1e-1: "整条染色体 / 无约束"}
    for se in (1e-5, 1e-4, 1e-3, 1e-2, 1e-1):
        kb = r_to_kb(min(se, 0.5))
        print(f"{se:>14.0e} {min(se,0.5):>10.0e} {kb:>11,.0f} kb {ref[se]:>34}")

    print("\n[3] 漂变阈值同时要满足 4N_e·s_ε > 1")
    print(f"{'s_ε':>10} {'4N_e s_ε':>12} {'漂变判定':>14} {'重组判定(r=0.5,跨染色体)':>26}")
    for se in (1e-5, 1e-4, 1e-3, 1e-2, 1e-1):
        S = 4 * NE * se
        print(f"{se:>10.0e} {S:>12.1f} {'可见' if S > 1 else '被漂变淹没':>14} "
              f"{'可维持' if se > 0.5 else '被重组拆散':>26}")

    print("\n[4] 由此得到的可检验预言")
    print("  (a) 严格 sign-epistatic 的调控共依赖，在有性群体里只能存在于")
    print(f"      r < s_ε 的窗口内。s_ε=1e-4 时窗口只有 {r_to_kb(1e-4):,.0f} kb —— ")
    print("      这正好落在 miRNA 簇 '<10 kb' 的经验定义上。")
    print("  (b) 跨染色体(r=0.5)的共依赖必须 s_ε>0.5（不现实），或每个成员单独有利。")
    print("      → 数据库里所有跨染色体的 miRNA–target'共进化'边，都不可能靠")
    print("        上位性维持，必须由各自独立的选择解释。这是一条硬性排除规则。")
    print("  (c) 因此：重组图谱冷点 / 倒位 / TAD 内部应当富集共依赖模块，")
    print("      而 miRNA 成簇本身就可以解读为'重组回避适应'。")
    return None


if __name__ == "__main__":
    run()
