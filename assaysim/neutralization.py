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
    """NT50 / Kd。<1 表示多价效应把中和滴度推到远低于 Kd 的浓度。"""
    return nt50_from_kd(1.0, virion)


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


# --- 已知病毒的刺突计数 (结构文献) --------------------------------------
# ⚠️ n 是"表面刺突数"，不等于"可被单克隆抗体同时结合的位点数"；
#    位阻与表位可及性会让有效 n 更小。这些值是上界。
KNOWN_VIRIONS = {
    "influenza_A": Virion("流感 A", n_spikes=340, k_hits=1,
                          source="HA 三聚体 ~300-400/粒子 (Harris 2006 PNAS)"),
    "SARS_CoV_2": Virion("SARS-CoV-2", n_spikes=25, k_hits=1,
                         source="S 三聚体 ~24-40/粒子 (Ke 2020 Nature)"),
    "HIV_1": Virion("HIV-1", n_spikes=14, k_hits=1,
                    source="Env 三聚体 ~7-14/粒子 (Zhu 2006 Nature)"),
}
