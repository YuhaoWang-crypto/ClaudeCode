"""
M7 — 平板凝集：Smoluchowski 聚集方程与前带效应

物理
----
颗粒 (乳胶球 / 红细胞 / 细菌) 通过抗体交联成簇。簇尺寸分布服从
Smoluchowski 聚集方程

    dn_k/dt = ½ Σ_{i+j=k} K_ij n_i n_j  −  n_k Σ_i K_ik n_i

常核 K_ij = K 时有**闭式解**（用于验证数值解）：

    τ = K n₀ t / 2
    n_k(t) = n₀ · τ^(k−1) / (1+τ)^(k+1)
    N(t)   = Σ n_k = n₀ / (1+τ)
    Σ k·n_k = n₀                              (质量守恒)

前带效应从何而来
----------------
交联需要**一端已结合抗体、另一端仍是空表位**。设表位占据率 θ，
桥联概率正比于 θ(1−θ)，归一化后

    K_eff = K₀ · 4θ(1−θ),      θ = [Ab] / ([Ab] + Kd)

θ=0.5（即 [Ab] = Kd）时最大。抗体极大时 θ→1，所有表位都被占满、
再没有空位可供搭桥 —— 这就是**前带 (prozone) / 钩状效应**。
它是模型的推论，不是外加的假设。

文献依据
--------
von Smoluchowski M (1917) Z Phys Chem 92:129 — 聚集方程与常核解
Bell GI (1978) Science 200:618 — 细胞黏附中的受体-配体桥联
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp


# --- 闭式解 --------------------------------------------------------------

def smoluchowski_analytic(k_sizes, t_s: float, K: float, n0: float) -> np.ndarray:
    """常核聚集方程的精确解 n_k(t)。"""
    tau = K * n0 * t_s / 2.0
    k = np.asarray(k_sizes, dtype=float)
    return n0 * tau ** (k - 1.0) / (1.0 + tau) ** (k + 1.0)


def total_particles_analytic(t_s: float, K: float, n0: float) -> float:
    """N(t) = n₀/(1+τ) 的精确解。"""
    return n0 / (1.0 + K * n0 * t_s / 2.0)


# --- 桥联动力学 ----------------------------------------------------------

def epitope_occupancy(ab_conc: float, kd: float) -> float:
    """Langmuir 表位占据率 θ。"""
    if ab_conc < 0 or kd <= 0:
        raise ValueError("浓度须非负、Kd 须为正")
    return ab_conc / (ab_conc + kd)


def bridging_efficiency(theta: float) -> float:
    """桥联效率 4θ(1−θ) ∈ [0,1]，θ=0.5 时取 1。"""
    if not 0.0 <= theta <= 1.0:
        raise ValueError("θ 必须在 [0,1]")
    return 4.0 * theta * (1.0 - theta)


def effective_kernel(ab_conc: float, kd: float, K0: float) -> float:
    """给定抗体浓度下的有效聚集核。"""
    return K0 * bridging_efficiency(epitope_occupancy(ab_conc, kd))


# --- 数值解 --------------------------------------------------------------

@dataclass
class AggResult:
    t: np.ndarray
    n: np.ndarray  # (nt, kmax) 各簇尺寸的数浓度
    k: np.ndarray  # 簇尺寸 1..kmax
    n0: float

    @property
    def N_total(self) -> np.ndarray:
        return self.n.sum(axis=1)

    @property
    def mass(self) -> np.ndarray:
        return (self.n * self.k).sum(axis=1)

    @property
    def mass_loss(self) -> np.ndarray:
        """因截断而流失的质量分数。

        ⚠️ 这是**截断误差的直接度量**，不是物理。i+j > k_max 的产物被丢弃，
           所以质量会缓慢下降。只有在 mass_loss 足够小的时间窗内，
           数值解才可信 —— validate.py 会检查这一点。
        """
        return 1.0 - self.mass / self.n0

    def visible_fraction(self, min_size: int = 4) -> np.ndarray:
        """处于 >= min_size 的簇中的颗粒质量分数 —— 肉眼可见凝集的代理量。"""
        m = self.n[:, self.k >= min_size] * self.k[self.k >= min_size]
        return m.sum(axis=1) / self.n0

    def mean_cluster_size(self) -> np.ndarray:
        return self.mass / np.maximum(self.N_total, 1e-300)


def simulate_aggregation(K: float, n0: float, t_end_s: float,
                         k_max: int = 60, n_out: int = 60) -> AggResult:
    """常核 Smoluchowski 方程的数值解（尺寸截断在 k_max）。

        dn_k/dt = ½K Σ_{i+j=k} n_i n_j  −  K n_k Σ_i n_i

    i+j > k_max 的产物被丢弃，因此总质量会缓慢流失；用 `mass_loss`
    监控截断误差，只在其足够小的时间窗内使用结果。
    """
    if k_max < 4:
        raise ValueError("k_max 至少为 4")
    k = np.arange(1, k_max + 1)

    def rhs(t, n):
        n = np.clip(n, 0.0, None)
        # 生成项用卷积：conv[m] = Σ_{i+j=m+2} n_i n_j (1-indexed)
        conv = np.convolve(n, n)[: k_max - 1]
        d = -K * n * n.sum()
        d[1:] += 0.5 * K * conv
        return d

    n_init = np.zeros(k_max)
    n_init[0] = n0
    sol = solve_ivp(rhs, (0.0, t_end_s), n_init,
                    t_eval=np.linspace(0.0, t_end_s, n_out),
                    method="LSODA", rtol=1e-10, atol=1e-14)
    if not sol.success:
        raise RuntimeError(f"聚集方程积分失败: {sol.message}")
    return AggResult(t=sol.t, n=sol.y.T, k=k, n0=n0)


def prozone_curve(ab_concs, kd: float, K0: float, n0: float,
                  t_read_s: float, k_max: int = 40, min_size: int = 4):
    """抗体浓度滴定曲线：给出经典的钩状 (prozone) 形状。

    返回每个浓度下读数时刻的可见凝集分数与平均簇大小。
    """
    out = []
    for c in ab_concs:
        th = epitope_occupancy(float(c), kd)
        Keff = effective_kernel(float(c), kd, K0)
        if Keff <= 0:
            out.append({"ab": float(c), "theta": th, "K_eff": 0.0,
                        "visible": 0.0, "mean_size": 1.0})
            continue
        r = simulate_aggregation(Keff, n0, t_read_s, k_max=k_max, n_out=4)
        out.append({
            "ab": float(c), "theta": th, "K_eff": Keff,
            "visible": float(r.visible_fraction(min_size)[-1]),
            "mean_size": float(r.mean_cluster_size()[-1]),
        })
    return out


def optimal_ab_concentration(kd: float) -> float:
    """使桥联效率最大的抗体浓度：θ=0.5 ⇒ [Ab] = Kd。"""
    return kd
