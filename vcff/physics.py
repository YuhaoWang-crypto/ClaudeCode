"""力场的组合律。

这一层是「模型产品」与「数据库」的分界线。数据库告诉你 HNF1A 在 HepG2 的
Chronos 值是 −1.42; 力场告诉你: 你的 3 个化合物、各自的浓度和效价、72 小时、
在 HepG2 上叠加之后, 相对正常上皮的选择性窗口是多少倍。

四条组合律 (全部为 MODELED 假设, 未在本数据上拟合):

  L1 剂量→占据率      θ = C^h / (C^h + IC50^h)                     [Hill]
  L2 占据率→表型      响应 = θ × 完全敲低的单位响应向量             [线性缩放]
  L3 多化合物叠加      log V = −λ Σ θᵢ·压力ᵢ                        [Bliss 独立]
  L4 交叉项           log V ← γ(上位性) × log V                     [ComboMap 修正]

λ 是**约定锚定**: 取 ln2, 使「单个中位必需靶点 (Chronos = −1) 在 100% 占据、
72 h」恰好给出 0.5 的伪存活率。它没有拟合任何存活率数据。因此绝对百分比无意义,
只有**比值**、**剂量位移**、**相对加和预期的偏离**是可解释的 —— 与分子力场
只有相对能量可解释完全同理。
"""

from __future__ import annotations

import math
from dataclasses import dataclass

#: L3 的约定锚定常数。见模块 docstring。
LAMBDA_ANCHOR = math.log(2.0)

#: 暴露时长的参考点 (小时)。Chronos 本身是 ~21 天筛的稳态量,
#: 时间缩放是纯粹的建模假设, 幅度上不可信, 只用于让曲线随时长单调移动。
REFERENCE_EXPOSURE_H = 72.0

#: 上位性修正的阻尼系数。γ 被夹在 [0.5, 2.0], 避免单个小样本组合主导结果。
EPISTASIS_DAMPING = 0.5
EPISTASIS_CLIP = (0.5, 2.0)


def occupancy(conc_uM: float, ic50_uM: float, hill: float = 1.0) -> float:
    """L1 · Hill 占据率。客户参数进, 无量纲 [0,1] 出。"""
    if conc_uM <= 0:
        return 0.0
    if ic50_uM <= 0:
        return 1.0
    x = (conc_uM / ic50_uM) ** hill
    return x / (1.0 + x)


def kill_pressure(dep: float, theta: float) -> float:
    """L2 · 单个化合物在某个系上的杀伤压力。

    dep 是 Chronos gene effect (越负越必需)。只取负的部分:
    正的 gene effect 表示敲掉反而长得更好, 抑制剂不会因此杀细胞。
    """
    return max(0.0, -dep) * theta


def bliss_log_viability(
    pressures: list[float],
    exposure_h: float = REFERENCE_EXPOSURE_H,
    gamma: float = 1.0,
) -> float:
    """L3 + L4 · 多化合物叠加 -> log(伪存活率)。"""
    total = sum(pressures)
    scale = exposure_h / REFERENCE_EXPOSURE_H
    return -LAMBDA_ANCHOR * total * scale * gamma


def viability(log_v: float) -> float:
    return math.exp(log_v)


def epistasis_gamma(epi: float, epi_median: float) -> float:
    """L4 · 由 ComboMap 的上位性强度得到交叉项系数 γ。

    γ > 1 -> 超加和 (协同); γ < 1 -> 次加和 (缓冲/上位)。
    """
    if epi_median <= 0:
        return 1.0
    raw = 1.0 + EPISTASIS_DAMPING * (epi / epi_median - 1.0)
    lo, hi = EPISTASIS_CLIP
    return min(hi, max(lo, raw))


def combine_pathway_vectors(
    terms: list[tuple[float, list[str], list[float]]],
) -> tuple[list[str], list[float]]:
    """按占据率加权叠加多个 Hallmark 通路向量。

    terms 是 (权重, 通路名列表, 通路值列表) 的列表。不同 kernel 给的通路子集
    不同, 这里按名字取并集再相加。
    """
    acc: dict[str, float] = {}
    for w, names, vals in terms:
        for n, v in zip(names, vals):
            acc[n] = acc.get(n, 0.0) + w * v
    if not acc:
        return [], []
    items = sorted(acc.items(), key=lambda kv: kv[1])
    return [k for k, _ in items], [v for _, v in items]


@dataclass
class DoseScan:
    """在客户给的剂量组合上按公共倍数扫描。"""

    multipliers: list[float]
    viability: dict[str, list[float]]      # line -> V(m)

    def crossing(self, line: str, level: float) -> float | None:
        """求 V(line) 首次降到 level 以下时的剂量倍数 (对数插值)。"""
        ms, vs = self.multipliers, self.viability[line]
        for i in range(1, len(ms)):
            if vs[i] <= level <= vs[i - 1]:
                v0, v1 = vs[i - 1], vs[i]
                if v0 == v1:
                    return ms[i]
                f = (v0 - level) / (v0 - v1)
                return math.exp(
                    math.log(ms[i - 1]) + f * (math.log(ms[i]) - math.log(ms[i - 1]))
                )
        return None


def log_grid(lo: float, hi: float, n: int) -> list[float]:
    step = (math.log(hi) - math.log(lo)) / (n - 1)
    return [math.exp(math.log(lo) + i * step) for i in range(n)]
