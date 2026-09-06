"""
M3 — 多击中占据模型 (multi-hit occupancy)：从分子亲和力 Kd 到中和滴度 NT50

模型
----
病毒粒子表面有 n 个可被抗体结合的刺突/表位。抗体过量条件下，单个位点的
占据概率服从 Langmuir 等温式

    θ([Ab]) = [Ab] / ([Ab] + Kd)

若中和要求"至少 k 个位点被占据"，则残留感染性 = 被占据位点数 < k 的概率

    P_inf(θ) = Σ_{i=0}^{k-1} C(n,i) θ^i (1-θ)^(n-i)

令 P_inf = 1/2 解出 θ*，回代 Langmuir 得

    NT50 = Kd · θ* / (1 - θ*)

k = 1 时有闭式解：P_inf = (1-θ)^n = 1/2  =>  θ* = 1 - 2^(-1/n)
    => NT50 = Kd · (2^(1/n) - 1)

文献依据
--------
Klasse PJ, Sattentau QJ (2002) J Gen Virol 83:2091 — 占据式中和的综述
Klasse PJ (2014) Adv Biol 2014:157895 — 中和的化学计量学

⚠️ 参数 k (击中阈值) 无法从结构预测，必须由中和曲线拟合。
   n 可从冷冻电镜结构计数。Kd 必须实测 (SPR/BLI)。
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from scipy.optimize import brentq
from scipy.stats import binom


@dataclass(frozen=True)
class Virion:
    """一种病毒粒子的中和几何。"""

    name: str
    n_spikes: int  # 表面可结合位点数
    k_hits: int  # 中和所需的最少占据位点数
    source: str = ""  # n 的结构来源

    def __post_init__(self) -> None:
        if not 1 <= self.k_hits <= self.n_spikes:
            raise ValueError("需要 1 <= k <= n")


def occupancy(ab_conc: float, kd: float) -> float:
    """Langmuir 单位点占据率 θ。ab_conc 与 kd 同单位 (如 nM)。"""
    if ab_conc < 0 or kd <= 0:
        raise ValueError("浓度须非负、Kd 须为正")
    return ab_conc / (ab_conc + kd)


def residual_infectivity(theta: float, n: int, k: int) -> float:
    """P_inf(θ) = P(占据位点数 <= k-1)，即二项分布的 CDF。"""
    if not 0.0 <= theta <= 1.0:
        raise ValueError("θ 必须在 [0,1]")
    return float(binom.cdf(k - 1, n, theta))


def neutralization_curve(ab_concs, kd: float, virion: Virion) -> list[float]:
    """给定抗体浓度序列，返回残留感染性 (1 = 完全未中和, 0 = 完全中和)。"""
    return [
        residual_infectivity(occupancy(c, kd), virion.n_spikes, virion.k_hits)
        for c in ab_concs
    ]


def theta_star(n: int, k: int, frac: float = 0.5) -> float:
    """解 P_inf(θ) = frac 的 θ*。"""
    if k == 1:
        # 闭式：(1-θ)^n = frac
        return 1.0 - frac ** (1.0 / n)
    lo, hi = 1e-15, 1.0 - 1e-15
    f = lambda t: residual_infectivity(t, n, k) - frac
    if f(lo) * f(hi) > 0:
        raise ValueError("θ* 不在 (0,1) 内，检查 n/k")
    return brentq(f, lo, hi, xtol=1e-14, rtol=1e-15)


def nt50_from_kd(kd: float, virion: Virion) -> float:
    """核心接口：分子亲和力 Kd -> 中和滴度 NT50 (同单位)。"""
    ts = theta_star(virion.n_spikes, virion.k_hits, 0.5)
    return kd * ts / (1.0 - ts)


def nt50_closed_form_single_hit(kd: float, n: int) -> float:
    """k=1 的解析解，用于验证数值解。"""
    return kd * (2.0 ** (1.0 / n) - 1.0)


def kd_from_nt50(nt50: float, virion: Virion) -> float:
    """反向接口：由实测 NT50 反推分子 Kd (在给定 n,k 假设下)。"""
    ts = theta_star(virion.n_spikes, virion.k_hits, 0.5)
    return nt50 * (1.0 - ts) / ts


def amplification_factor(virion: Virion) -> float:
    """NT50 / Kd。<1 表示位点冗余把中和滴度推到低于 Kd 的浓度。"""
    return nt50_from_kd(1.0, virion)


# --- 实测 θ* 先验（本模块最重要的一段） ---------------------------------
#
# θ* 是"半中和时的单位点占据率"，等价于击中阈值 k 的连续化形式。
# 它**无法从结构预测**，只能由中和曲线拟合 —— 早期版本用 k=1（单击中）作默认，
# 这个默认没有任何数据依据，而且被实测**证伪**了。
#
# 外部工作（VC-CON）用分层贝叶斯模型
#     logit(θ*) ~ StudentT(ν, μ_virus + β_occupancy·assay, σ)
# 在 38 个**单价**抗体配对点上拟合（冠状病毒 25 纳米抗体 / HIV-1 5 Fab /
# 流感 8 混合），θ* = IC50/(IC50+Kd)，与本模块 NT50 = Kd·θ*/(1−θ*) 同一式。
#
#   病毒        n_spike   θ* 中位数   95% CI            隐含 k
#   冠状病毒        24     0.576    [0.345, 0.777]     ≈14 (58%)
#   HIV-1          14     0.586    [0.260, 0.850]     ≈9  (64%)
#   流感 A        375     0.620    [0.318, 0.876]     ≈233 (62%)
#   全局（新病毒）   —     0.596    [0.289, 0.839]     ≈0.6·n
#
# 结论：**守恒的是占据分数 θ*，不是绝对击中数 k**。k 随 n 等比例放大。
# 若改为假设 k 恒定，θ* 在三个病毒间会相差 16.3× [9.6, 29.6] —— 数据不支持。
#
# 对本模块的直接后果：k=1 会把 NT50/Kd 低估 12×(HIV) 到 336×(流感)。
# 真实的单价抗体 NT50 ≈ 1.4–1.6 × Kd，**略高于** Kd，而不是低几个数量级。
#
# ⚠️ 边界：这 38 个点全是**单价**结合物（纳米抗体 / Fab / 设计小蛋白）。
#    双价 IgG 的亲合力效应可能把 θ* 压低，本先验不覆盖那一档。
#    ⚠️ 95% CI 很宽（θ* 从 0.29 到 0.84），跨病毒离散 τ 中位数 0.522 ——
#    "守恒"是指三个病毒的可信区间高度重叠，不是指 θ* 被钉死在 0.6。

EMPIRICAL_THETA_STAR = {
    "coronavirus": {"median": 0.5764, "ci": (0.3451, 0.7773), "n_spike": 24},
    "HIV-1": {"median": 0.5863, "ci": (0.2604, 0.8495), "n_spike": 14},
    "influenza": {"median": 0.6200, "ci": (0.3183, 0.8763), "n_spike": 375},
    "global": {"median": 0.5961, "ci": (0.2886, 0.8387), "n_spike": None},
}
THETA_STAR_SOURCE = (
    "VC-CON 分层贝叶斯拟合，38 个单价抗体配对点；NUTS 4 链 × 6000 抽样，"
    "最大 R̂ 1.001。仅覆盖单价结合物。")


def k_from_theta_star(n: int, theta: float) -> int:
    """由 θ* 反解最接近的整数击中阈值 k（给定位点数 n）。"""
    if not 0.0 < theta < 1.0:
        raise ValueError("θ* 必须在 (0,1)")
    return min(range(1, n + 1), key=lambda k: abs(theta_star(n, k) - theta))


def nt50_from_kd_empirical(kd: float, virus: str = "global") -> dict:
    """用**实测 θ* 先验**直接由 Kd 得 NT50，并带 95% 区间。

    这条路径绕开 k —— 既然守恒的是 θ*，就不必先猜 k 再算回来。
    """
    if virus not in EMPIRICAL_THETA_STAR:
        raise ValueError(f"virus 须为 {sorted(EMPIRICAL_THETA_STAR)} 之一")
    e = EMPIRICAL_THETA_STAR[virus]
    f = lambda t: kd * t / (1.0 - t)
    return {"nt50": f(e["median"]),
            "ci": (f(e["ci"][0]), f(e["ci"][1])),
            "theta_star": e["median"], "source": THETA_STAR_SOURCE}


def fit_k_from_curve(ab_concs, residual, kd: float, n: int) -> int:
    """从实测中和曲线拟合击中阈值 k (n 与 Kd 已知)。

    k 是整数，直接穷举取残差平方和最小者 —— 这是 k 唯一诚实的求法。
    """
    best_k, best_sse = 1, math.inf
    for k in range(1, n + 1):
        sse = 0.0
        for c, r in zip(ab_concs, residual):
            pred = residual_infectivity(occupancy(c, kd), n, k)
            sse += (pred - r) ** 2
        if sse < best_sse:
            best_k, best_sse = k, sse
    return best_k


# --- 已知病毒的刺突计数与**实测**击中阈值 --------------------------------
# ⚠️ n 是"表面刺突数"，不等于"可被单克隆抗体同时结合的位点数"；
#    位阻与表位可及性会让有效 n 更小。这些值是上界。
#
# k_hits 现在取自 EMPIRICAL_THETA_STAR 反解的整数值，**不再是 k=1**。
# 早期版本默认 k=1（单击中），该默认无数据依据且与实测相差 12–336 倍。
KNOWN_VIRIONS = {
    "influenza_A": Virion("流感 A", n_spikes=375, k_hits=233,
                          source="HA 三聚体 375/粒子；k 由实测 θ*=0.620 反解 (VC-CON)"),
    "SARS_CoV_2": Virion("SARS-CoV-2", n_spikes=24, k_hits=14,
                         source="S 三聚体 24/粒子；k 由实测 θ*=0.576 反解 (VC-CON)"),
    "HIV_1": Virion("HIV-1", n_spikes=14, k_hits=9,
                    source="Env 三聚体 14/粒子；k 由实测 θ*=0.586 反解 (VC-CON)"),
    # 单击中假设保留一个显式条目，专门用来演示它错得有多离谱
    "influenza_A_singlehit": Virion("流感 A（k=1，已被证伪）", n_spikes=375, k_hits=1,
                                    source="⚠️ 仅作反例：实测 θ* 高 336 倍"),
}
