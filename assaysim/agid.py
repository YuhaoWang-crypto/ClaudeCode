"""
M5 — 琼脂免疫扩散 (AGID / Ouchterlony)：反应-扩散 PDE 与沉淀环

物理
----
抗原与抗体从两个孔向凝胶中扩散，在浓度比落入**等价带**的位置形成免疫复合物
晶格并沉淀成可见的沉淀线。

    ∂A/∂t = D_A ∂²A/∂x² − R
    ∂B/∂t = D_B ∂²B/∂x² − ν R
    ∂P/∂t = R                      (沉淀物不扩散)

沉淀速率取 Heidelberger–Kendall 的等价带行为：

    R = k_p · A · B · W(r),    r = A / (ν B)
    W(r) = exp( −(ln r / w)² )

W 在 r = 1 (等价) 处最大，抗原过量 (r≫1) 与抗体过量 (r≪1) 两侧都衰减 ——
这就是**后带 (postzone)** 与**前带 (prozone)** 效应，不是外加的假设，
而是窗函数的直接推论。

扩散系数
--------
Stokes–Einstein:  D = k_B T / (6 π η r_h)

⚠️ 凝胶中的 D 低于自由水中的值。琼脂糖对 IgG 大小的蛋白通常给出
   0.5–0.9 的阻滞因子 (取决于凝胶浓度)，本模块用 `gel_factor` 显式暴露，
   默认 1.0 (自由水)，不替使用者猜。

文献依据
--------
Heidelberger M, Kendall FE (1935) J Exp Med 62:697 — 定量沉淀曲线
Ouchterlony Ö (1958) Prog Allergy 5:1 — 双向免疫扩散
Crowle AJ (1973) Immunodiffusion, 2nd ed. — 扩散-沉淀几何
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import erfc

K_B = 1.380649e-23  # J/K
ETA_WATER_20C = 1.002e-3  # Pa·s
ETA_WATER_37C = 6.913e-4  # Pa·s


def stokes_einstein_D(r_h_nm: float, T_K: float = 293.15,
                      eta_Pa_s: float = ETA_WATER_20C,
                      gel_factor: float = 1.0) -> float:
    """Stokes–Einstein 扩散系数，返回 cm²/s。

    r_h_nm      流体力学半径 (nm)
    gel_factor  凝胶阻滞因子 (0<f<=1)；1.0 = 自由水
    """
    if r_h_nm <= 0 or T_K <= 0 or eta_Pa_s <= 0:
        raise ValueError("参数必须为正")
    if not 0 < gel_factor <= 1:
        raise ValueError("gel_factor 应在 (0,1]")
    D_m2 = K_B * T_K / (6.0 * math.pi * eta_Pa_s * r_h_nm * 1e-9)
    return D_m2 * 1e4 * gel_factor  # m²/s -> cm²/s


# 实测流体力学半径与扩散系数 (20 °C, 水)，用于验证 Stokes–Einstein 实现
# 数值为文献常用值，量级与彼此关系是被反复测量确认的
REFERENCE_PROTEINS = {
    "IgG": {"MW_kDa": 150.0, "r_h_nm": 5.3, "D_measured_cm2_s": 4.0e-7},
    "BSA": {"MW_kDa": 66.5, "r_h_nm": 3.5, "D_measured_cm2_s": 6.0e-7},
    "lysozyme": {"MW_kDa": 14.3, "r_h_nm": 1.9, "D_measured_cm2_s": 1.06e-6},
}


# --- 纯扩散 (用于验证求解器) --------------------------------------------

def diffuse_semi_infinite_analytic(x_cm, t_s: float, D: float, C0: float = 1.0):
    """半无限域、边界恒定浓度的解析解 C = C0·erfc(x / 2√(Dt))。"""
    x = np.asarray(x_cm, dtype=float)
    if t_s <= 0:
        return np.where(x <= 0, C0, 0.0)
    return C0 * erfc(x / (2.0 * math.sqrt(D * t_s)))


def diffuse_numeric(L_cm: float, n: int, t_s: float, D: float,
                    C0: float = 1.0) -> tuple[np.ndarray, np.ndarray]:
    """同一问题的数值解：x=0 处 Dirichlet C0，右端零通量，初值 0。"""
    dx = L_cm / (n - 1)
    x = np.linspace(0.0, L_cm, n)

    def rhs(t, C):
        lap = np.empty_like(C)
        lap[1:-1] = (C[2:] - 2 * C[1:-1] + C[:-2]) / dx**2
        lap[0] = 0.0  # Dirichlet，由下面强制
        lap[-1] = (2 * C[-2] - 2 * C[-1]) / dx**2  # 零通量
        d = D * lap
        d[0] = 0.0
        return d

    C = np.zeros(n)
    C[0] = C0
    sol = solve_ivp(rhs, (0.0, t_s), C, t_eval=[t_s], method="LSODA",
                    rtol=1e-10, atol=1e-12)
    if not sol.success:
        raise RuntimeError(sol.message)
    return x, sol.y[:, 0]


# --- 完整 AGID 模型 ------------------------------------------------------

@dataclass(frozen=True)
class AgidSetup:
    """一次双向免疫扩散实验的配置。"""

    L_cm: float = 1.0  # 两孔间距
    n_grid: int = 201
    D_ag: float = 6.0e-7  # 抗原扩散系数 cm²/s
    D_ab: float = 4.0e-7  # 抗体扩散系数 cm²/s (IgG 较大、较慢)
    A_well: float = 1.0  # 抗原孔浓度 (任意单位)
    B_well: float = 1.0  # 抗体孔浓度 (同单位)
    nu: float = 2.0  # 等价时每个抗原结合的抗体分子数
    k_p: float = 5.0  # 沉淀速率常数
    log_width: float = 0.9  # 等价带宽度 (ln r 的高斯宽度)

    @property
    def dx(self) -> float:
        return self.L_cm / (self.n_grid - 1)


@dataclass
class AgidResult:
    x: np.ndarray
    t: np.ndarray
    A: np.ndarray  # (nt, nx)
    B: np.ndarray
    P: np.ndarray
    setup: AgidSetup

    def precipitate_total(self) -> np.ndarray:
        """每个时刻的沉淀总量 (对 x 积分)。"""
        return np.trapezoid(self.P, self.x, axis=1)

    def band_position(self, idx: int = -1) -> float:
        """沉淀带位置 = 沉淀分布的质心 (cm，从抗原孔量起)。无沉淀返回 nan。"""
        p = self.P[idx]
        s = np.trapezoid(p, self.x)
        if s <= 1e-14:
            return float("nan")
        return float(np.trapezoid(p * self.x, self.x) / s)

    def band_sharpness(self, idx: int = -1) -> float:
        """沉淀带宽度 (标准差, cm)。越小线越锐利。"""
        p = self.P[idx]
        s = np.trapezoid(p, self.x)
        if s <= 1e-14:
            return float("nan")
        m = self.band_position(idx)
        var = np.trapezoid(p * (self.x - m) ** 2, self.x) / s
        return float(math.sqrt(max(var, 0.0)))


def equivalence_window(A, B, nu: float, log_width: float):
    """Heidelberger–Kendall 等价带窗函数 W ∈ (0,1]。

    ν = 晶格中每个抗原分子消耗的抗体分子数，因此"抗体需求量 / 抗体供给量"

        r = ν·A / B

    在 r = 1 处恰好化学计量配平 —— 窗函数的峰与化学计量等价点重合。
    r ≫ 1 抗原过量（后带），r ≪ 1 抗体过量（前带），两侧都不成晶格。
    """
    A = np.asarray(A, dtype=float)
    B = np.asarray(B, dtype=float)
    eps = 1e-300
    r = (nu * A + eps) / (B + eps)
    return np.exp(-((np.log(r) / log_width) ** 2))


def simulate_agid(s: AgidSetup, t_end_s: float = 24 * 3600.0,
                  n_out: int = 40) -> AgidResult:
    """跑一次双向免疫扩散。

    x=0 为抗原孔，x=L 为抗体孔，两端按恒浓度储库处理 (孔远大于凝胶通道)。
    """
    n, dx = s.n_grid, s.dx
    x = np.linspace(0.0, s.L_cm, n)

    def lap(C, left, right):
        out = np.empty_like(C)
        ext_l = 2 * left - C[0]  # Dirichlet 镜像
        ext_r = 2 * right - C[-1]
        out[0] = (C[1] - 2 * C[0] + ext_l) / dx**2
        out[1:-1] = (C[2:] - 2 * C[1:-1] + C[:-2]) / dx**2
        out[-1] = (ext_r - 2 * C[-1] + C[-2]) / dx**2
        return out

    def rhs(t, y):
        A, B, P = y[:n], y[n:2 * n], y[2 * n:]
        A = np.clip(A, 0.0, None)
        B = np.clip(B, 0.0, None)
        R = s.k_p * A * B * equivalence_window(A, B, s.nu, s.log_width)
        dA = s.D_ag * lap(A, s.A_well, 0.0) - R
        dB = s.D_ab * lap(B, 0.0, s.B_well) - s.nu * R
        return np.concatenate([dA, dB, R])

    y0 = np.concatenate([np.zeros(n), np.zeros(n), np.zeros(n)])
    y0[0] = s.A_well
    y0[2 * n - 1] = s.B_well

    t_eval = np.linspace(0.0, t_end_s, n_out)
    sol = solve_ivp(rhs, (0.0, t_end_s), y0, t_eval=t_eval, method="LSODA",
                    rtol=1e-8, atol=1e-11)
    if not sol.success:
        raise RuntimeError(f"AGID 积分失败: {sol.message}")

    Y = sol.y.T
    return AgidResult(x=x, t=sol.t, A=Y[:, :n], B=Y[:, n:2 * n],
                      P=Y[:, 2 * n:], setup=s)


def tube_precipitin(A_total: float, B_total: float, nu: float = 2.0,
                    k_p: float = 5.0, log_width: float = 0.9,
                    t_end: float = 2000.0) -> dict:
    """试管沉淀试验（封闭、混匀体系）—— 经典 Heidelberger–Kendall 构型。

    抗原与抗体**无论比例如何都会结合**；区别在于生成的复合物能否长成
    可沉淀的晶格。落在等价带内的那一份沉淀，其余留在溶液里：

        dA/dt = −k·A·B
        dB/dt = −ν·k·A·B
        dP/dt =     W ·k·A·B      (沉淀)
        dS/dt = (1−W)·k·A·B       (可溶复合物，不再回到游离态)

    可溶复合物这一支是关键：抗原大幅过量时抗体被迅速消耗进**不沉淀**的
    小复合物，沉淀量因此塌陷 —— 这才是后带效应的机制。若只留沉淀一个去路，
    给足时间总会缓慢累积，钟形曲线就出不来。

    与双向免疫扩散的区别在于这里**没有储库补料**。
    """
    def rhs(t, y):
        A, B = max(y[0], 0.0), max(y[1], 0.0)
        flux = k_p * A * B
        W = float(equivalence_window(np.array([A]), np.array([B]), nu, log_width)[0])
        return [-flux, -nu * flux, W * flux, (1.0 - W) * flux]

    sol = solve_ivp(rhs, (0.0, t_end), [A_total, B_total, 0.0, 0.0],
                    method="LSODA", rtol=1e-10, atol=1e-14)
    if not sol.success:
        raise RuntimeError(sol.message)
    return {"precipitate": float(sol.y[2, -1]), "soluble": float(sol.y[3, -1]),
            "A_left": float(sol.y[0, -1]), "B_left": float(sol.y[1, -1])}


def heidelberger_kendall_curve(ag_amounts, B_total: float = 1.0,
                               nu: float = 2.0, **kw):
    """固定抗体量、扫抗原量，返回经典的定量沉淀曲线（钟形）。

    化学计量等价点在 ν·A = B，即 A/B = 1/ν。
    """
    out = []
    for a in ag_amounts:
        r = tube_precipitin(float(a), B_total, nu, **kw)
        out.append({"A_total": float(a), "ratio": float(a) / B_total,
                    "equivalence_r": nu * float(a) / B_total, **r})
    return out


def titration_curve(base: AgidSetup, ag_concs, t_end_s: float = 24 * 3600.0):
    """棋盘滴定：固定抗体孔浓度，扫抗原孔浓度，返回每个点的沉淀总量。

    这是本模型**唯一真正省钱的用途** —— 预测出现清晰沉淀线的抗原浓度窗口，
    省掉真正的棋盘滴定。
    """
    import dataclasses

    out = []
    for a in ag_concs:
        s = dataclasses.replace(base, A_well=float(a))
        r = simulate_agid(s, t_end_s, n_out=6)
        out.append({
            "A_well": float(a),
            "ratio": float(a) / base.B_well,
            "precipitate": float(r.precipitate_total()[-1]),
            "band_x": r.band_position(),
            "band_sd": r.band_sharpness(),
        })
    return out
