"""
M4 — 靶细胞受限病毒动力学 ODE：单孔数字孪生

状态方程 (带潜隐期 eclipse 的四状态模型)

    dT/dt = -β T V
    dE/dt =  β T V - k E
    dI/dt =  k E - δ I
    dV/dt =  p_eff I - c V - β T V

其中 β T V 项同时从 V 中扣除 (吸附消耗病毒)。

基本再生数        R0 = β p T0 / (δ c)
细胞病变读数      CPE(t) = (E + I + D) / T0     D = 已裂解细胞累计
杀伤率            同 CPE

文献依据
--------
Baccam P et al. (2006) J Virol 80:7590 — 人流感 A 的靶细胞受限模型 (含 eclipse)
Perelson AS (2002) Nat Rev Immunol 2:28 — 该模型族的综述
Beauchemin CAA, Handel A (2011) BMC Public Health 11:S7 — 体外/体内参数综述

药物接口 (本模块的核心设计)
--------------------------
药物**不是**一个自由的抑制因子，而是通过明确的作用机制修改某一个速率常数：

    entry        进入抑制剂 / 中和抗体   ->  β  *= (1-ε)
    replication  聚合酶 / 逆转录酶抑制剂  ->  p  *= (1-ε)
    maturation   蛋白酶抑制剂            ->  产出病毒的感染性比例 *= (1-ε)
    release      NA 抑制剂               ->  释放速率 p *= (1-ε)，且已结合病毒不脱落

ε(D) = D^h / (D^h + EC50_mol^h)   —— EC50_mol 是**分子层**效力

关键结论 (validate.py 里有数值验证):
    **测出来的孔水平 EC50 ≠ 分子层 EC50_mol**。
    二者的比值随 MOI 与读板时间系统性变化 —— 这是模型的预言，不是实验噪声。
"""

from __future__ import annotations

import math
from dataclasses import dataclass, replace

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, curve_fit

from .neutralization import Virion, occupancy, residual_infectivity

MOA = ("entry", "replication", "maturation", "release", "neutralizing_ab")


@dataclass(frozen=True)
class CellVirusParams:
    """一个 病毒 x 细胞系 组合的动力学参数。

    ⚠️ 这些参数**不可跨组合迁移**，必须对每个组合拟合。
    时间单位统一为小时 (h)。
    """

    name: str
    beta: float  # 1/(virion*h)   感染速率
    k_eclipse: float  # 1/h            潜隐期 -> 产毒 转换
    delta: float  # 1/h            产毒细胞死亡
    p: float  # virion/(cell*h) 产毒速率
    c: float  # 1/h            游离病毒清除
    T0: float = 1e5  # cells/well     初始靶细胞
    source: str = ""

    @property
    def R0(self) -> float:
        """R0 = β p T0 / (δ (c + β T0))

        分母含 β·T0 是因为 V 方程里的吸附项 −βTV 也在消耗游离病毒：
        一个病毒粒子被细胞吸附后就不再能感染别的细胞。
        当 c >> β·T0 时退化为常见的 β p T0/(δc)。
        """
        return (self.beta * self.p * self.T0
                / (self.delta * (self.c + self.beta * self.T0)))

    def growth_rate(self) -> float:
        """指数期增长率 = 线性化系统的最大实特征值。

        在 (T0,0,0,0) 处对 (E,I,V) 线性化：

            J = [[-k,  0,        β T0      ],
                 [ k, -δ,        0         ],
                 [ 0,  p, -(c + β T0)      ]]

        特征方程: (λ+k)(λ+δ)(λ+c+β T0) = k·p·β·T0

        ⚠️ λ 是**渐近**增长率。接种后有一段瞬态（只有 V、还没有 E/I），
           要在瞬态之后取窗口才能测到 λ —— 见 validate.py。
        """
        k, d = self.k_eclipse, self.delta
        c = self.c + self.beta * self.T0  # 吸附损耗并入清除
        rhs = k * self.p * self.beta * self.T0
        f = lambda lam: (lam + k) * (lam + d) * (lam + c) - rhs
        if f(0.0) >= 0:  # R0 <= 1，无正根
            lo, hi = -min(k, d, c) + 1e-12, 0.0
            if f(lo) * f(hi) > 0:
                return 0.0
            return brentq(f, lo, hi, xtol=1e-14)
        hi = 1.0
        while f(hi) < 0:
            hi *= 2
            if hi > 1e6:
                raise RuntimeError("增长率求解发散")
        return brentq(f, 0.0, hi, xtol=1e-14)


@dataclass(frozen=True)
class Drug:
    """药物 = (作用机制, 分子层效力, Hill 系数)。"""

    name: str
    moa: str
    ec50_mol: float  # 分子层 EC50 (nM)，来自生化测定
    hill: float = 1.0

    def __post_init__(self) -> None:
        if self.moa not in MOA:
            raise ValueError(f"moa 必须是 {MOA} 之一")
        if self.ec50_mol <= 0:
            raise ValueError("EC50 须为正")

    def efficacy(self, dose_nM: float) -> float:
        """ε(D) ∈ [0,1)。"""
        if dose_nM <= 0:
            return 0.0
        x = (dose_nM / self.ec50_mol) ** self.hill
        return x / (1.0 + x)


@dataclass(frozen=True)
class Antibody:
    """中和抗体 = (Kd, 病毒几何)。通过 M3 占据模型进入 ODE。"""

    name: str
    kd_nM: float
    virion: Virion

    def efficacy(self, conc_nM: float) -> float:
        """被中和的病毒比例 = 1 - P_inf。"""
        th = occupancy(conc_nM, self.kd_nM)
        return 1.0 - residual_infectivity(th, self.virion.n_spikes, self.virion.k_hits)


def _apply(params: CellVirusParams, agent, dose: float) -> tuple[CellVirusParams, float]:
    """把药物作用施加到参数上。返回 (修改后的参数, 感染性折减因子)。"""
    if agent is None or dose <= 0:
        return params, 1.0
    eps = agent.efficacy(dose)
    if isinstance(agent, Antibody):
        return replace(params, beta=params.beta * (1.0 - eps)), 1.0
    if agent.moa == "entry":
        return replace(params, beta=params.beta * (1.0 - eps)), 1.0
    if agent.moa in ("replication", "release"):
        return replace(params, p=params.p * (1.0 - eps)), 1.0
    if agent.moa == "maturation":
        # 病毒仍被产出，但感染性下降：p 不变，感染性比例下降
        return params, (1.0 - eps)
    raise ValueError(agent.moa)


@dataclass
class Trajectory:
    t: np.ndarray
    T: np.ndarray
    E: np.ndarray
    I: np.ndarray
    V: np.ndarray
    D: np.ndarray  # 累计死亡细胞
    T0: float

    @property
    def cpe(self) -> np.ndarray:
        """细胞病变比例 = 1 - 存活未感染细胞比例。"""
        return 1.0 - self.T / self.T0

    def cpe_at(self, t_h: float) -> float:
        return float(np.interp(t_h, self.t, self.cpe))

    def peak_titer(self) -> float:
        return float(self.V.max())


def simulate(params: CellVirusParams, *, moi: float = 0.01, t_end_h: float = 96.0,
             agent=None, dose: float = 0.0, n_points: int = 400) -> Trajectory:
    """跑一个孔。moi = 初始 V / T0。"""
    par, infect_frac = _apply(params, agent, dose)
    V0 = moi * params.T0

    def rhs(t, y):
        T, E, I, V, D = y
        inf = par.beta * T * V
        return [
            -inf,
            inf - par.k_eclipse * E,
            par.k_eclipse * E - par.delta * I,
            par.p * I * infect_frac - par.c * V - inf,
            par.delta * I,
        ]

    sol = solve_ivp(rhs, (0.0, t_end_h), [params.T0, 0.0, 0.0, V0, 0.0],
                    t_eval=np.linspace(0.0, t_end_h, n_points),
                    method="LSODA", rtol=1e-9, atol=1e-9)
    if not sol.success:
        raise RuntimeError(f"积分失败: {sol.message}")
    return Trajectory(sol.t, *sol.y, T0=params.T0)


# --- 剂量-反应 -> 孔水平表观 EC50 ---------------------------------------

def _hill(d, bottom, top, ec50, h):
    return bottom + (top - bottom) / (1.0 + (ec50 / d) ** h)


def dose_response(params: CellVirusParams, agent, doses_nM, *, moi: float = 0.01,
                  readout_h: float = 72.0, t_end_h: float | None = None) -> np.ndarray:
    """返回每个剂量下 readout_h 时刻的 CPE (0..1)。"""
    t_end = t_end_h if t_end_h is not None else readout_h * 1.2
    return np.array([
        simulate(params, moi=moi, t_end_h=t_end, agent=agent, dose=float(d)).cpe_at(readout_h)
        for d in doses_nM
    ])


def apparent_ec50(params: CellVirusParams, agent, *, moi: float = 0.01,
                  readout_h: float = 72.0, n_doses: int = 24,
                  span: float = 4.0) -> dict:
    """孔水平表观 EC50：把模拟出的 CPE 抑制曲线当作真实实验数据来拟合。

    这就是湿实验真正测到的数, 与 agent 的分子层 EC50 不同。
    """
    centre = agent.ec50_mol if isinstance(agent, Drug) else agent.kd_nM
    doses = np.logspace(math.log10(centre) - span, math.log10(centre) + span, n_doses)

    cpe = dose_response(params, agent, doses, moi=moi, readout_h=readout_h)
    cpe0 = simulate(params, moi=moi, t_end_h=readout_h * 1.2).cpe_at(readout_h)
    if cpe0 <= 1e-9:
        raise RuntimeError("无药对照没有 CPE，无法定义保护率")

    protection = 1.0 - cpe / cpe0  # 0 = 无保护, 1 = 完全保护

    try:
        popt, _ = curve_fit(
            _hill, doses, protection,
            p0=[0.0, 1.0, centre, 1.0],
            bounds=([-0.2, 0.5, doses[0], 0.2], [0.2, 1.5, doses[-1], 8.0]),
            maxfev=40000,
        )
        ec50_app, hill_app = float(popt[2]), float(popt[3])
    except RuntimeError:
        ec50_app, hill_app = float("nan"), float("nan")

    return {
        "doses": doses,
        "protection": protection,
        "cpe": cpe,
        "cpe_untreated": cpe0,
        "ec50_apparent": ec50_app,
        "hill_apparent": hill_app,
        "ec50_molecular": centre,
        "shift_log10": math.log10(ec50_app / centre) if ec50_app == ec50_app else float("nan"),
    }


# --- 参考参数集 ----------------------------------------------------------
# ⚠️ 这些是文献量级的**起始猜测**，不是对任一具体 病毒x细胞 组合的拟合结果。
#    真实使用时必须用自家的生长曲线重新拟合 (见 calibrate.py)。
REFERENCE_PARAMS = {
    "influenza_MDCK": CellVirusParams(
        name="流感 A / MDCK",
        beta=3e-7, k_eclipse=1 / 6.0, delta=1 / 12.0, p=5.0, c=1 / 6.0, T0=1e5,
        source="量级取自 Baccam 2006 / Beauchemin 2011 体外流感参数区间",
    ),
    "generic_slow": CellVirusParams(
        name="慢速 CPE 通用",
        beta=1e-7, k_eclipse=1 / 10.0, delta=1 / 24.0, p=2.0, c=1 / 12.0, T0=1e5,
        source="通用慢生长病毒的示例参数",
    ),
}
