"""诚实标注层.

力场类比的核心诚实性约束: **只有相对量有物理意义**。

分子力场给出的绝对能量依赖于任意的参考态, 只有能量差、相对排序、构象位移
是可解释的。本引擎完全一样:

  - 可解释: 细胞系之间的比值 (选择性窗口)、剂量位移 (IC50 相对多少倍)、
            组合相对加和预期的偏离、通路响应的相对排序。
  - 不可解释: 绝对存活率百分比、绝对 IC50、NOAEL、安全窗数值。

每一个输出数字都必须带一个 Evidence 等级。
"""

from __future__ import annotations

from dataclasses import dataclass, field, asdict
from enum import Enum
from typing import Any


class Evidence(str, Enum):
    """证据等级。"""

    #: 直接来自已发表实测数据 (DepMap Chronos / Replogle Perturb-seq /
    #: Norman 双扰动 / Tahoe-100M), 未经外推。
    ANCHORED = "ANCHORED"

    #: 实测量经过模型外推 (跨系共识、指纹→表型预测), 有报告过的验证指标。
    PREDICTED = "PREDICTED"

    #: 引擎自己引入的建模假设 (Hill 占据率、指数杀伤、λ 锚定、跨模态迁移)。
    #: 这些假设**没有**在本数据上拟合或验证过。
    MODELED = "MODELED"

    #: 数据覆盖不到, 引擎拒绝给数。
    UNSUPPORTED = "UNSUPPORTED"


#: 引擎明确拒绝输出的东西 —— 出现在每张结果卡上。
REFUSED_OUTPUTS = [
    "NOAEL / LOAEL 数值",
    "绝对 IC50 或 EC50 (效价是你的输入参数, 不是本引擎的预测)",
    "绝对细胞存活率百分比 (λ 为约定锚定, 非拟合)",
    "临床安全窗、给药剂量、毒理放行结论",
]


@dataclass
class Claim:
    """一个带证据等级的输出量。"""

    key: str
    label: str
    value: Any
    unit: str = ""
    evidence: Evidence = Evidence.MODELED
    basis: str = ""
    caveat: str = ""

    def to_dict(self) -> dict:
        d = asdict(self)
        d["evidence"] = self.evidence.value
        return d

    def __str__(self) -> str:
        mark = {
            Evidence.ANCHORED: "✅",
            Evidence.PREDICTED: "◐",
            Evidence.MODELED: "⚠️",
            Evidence.UNSUPPORTED: "⛔",
        }[self.evidence]
        v = self.value
        vs = f"{v:.4g}" if isinstance(v, float) else str(v)
        line = f"{mark} {self.label}: {vs}{(' ' + self.unit) if self.unit else ''}"
        if self.caveat:
            line += f"\n      ↳ {self.caveat}"
        return line


@dataclass
class Ledger:
    """一次求值中所有 Claim 的账本。"""

    claims: list[Claim] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    def add(self, claim: Claim) -> Claim:
        self.claims.append(claim)
        return claim

    def note(self, text: str) -> None:
        self.notes.append(text)

    def by_evidence(self, level: Evidence) -> list[Claim]:
        return [c for c in self.claims if c.evidence is level]

    def to_dict(self) -> dict:
        return {
            "claims": [c.to_dict() for c in self.claims],
            "notes": list(self.notes),
            "refused_outputs": list(REFUSED_OUTPUTS),
            "counts": {
                lvl.value: len(self.by_evidence(lvl)) for lvl in Evidence
            },
        }
