"""客户输入的实验规格 (AssaySpec) —— 引擎的输入适配接口。

这一层就是「力场的参数文件」: 模型不知道你的化合物有多强、你用多少浓度、
你测几小时。这些全部由客户提供; 模型只提供响应的形状。
"""

from __future__ import annotations

from dataclasses import dataclass, field, asdict
from typing import Any


#: 引擎支持的细胞上下文 (来自 TwinCell 的 4 个系)。
CONTEXTS = {
    "hepg2": "HepG2 · 肝细胞癌",
    "rpe1": "RPE1 · 正常视网膜色素上皮 (正常组织对照)",
    "k562": "K562 · 慢性髓系白血病",
    "jurkat": "Jurkat · T 细胞急性淋巴细胞白血病",
}

#: 必需性 (Chronos gene effect) 骨干层实际覆盖的两个系。
#: 毒性/选择性主轴永远在这一对上计算。
ESSENTIALITY_LINES = ("hepg2", "rpe1")


@dataclass
class Compound:
    """一个化合物条目。

    Parameters
    ----------
    name:
        化合物名 (会尝试映射到靶点), 或直接写靶基因符号。
    conc_uM:
        客户实际使用的浓度 (µM)。**客户参数, 非模型预测。**
    ic50_uM:
        该化合物对其靶点的效价 (µM)。**客户参数。** 引擎不预测效价;
        不填时用 1.0 µM 占位, 结果里会明确标注这是占位值。
    hill:
        Hill 系数, 默认 1.0 (无协同结合)。
    target:
        显式指定靶基因; 留空则从化合物名自动映射。
    modality:
        "inhibitor" (抑制, 等效于功能缺失) 或 "activator" (激活)。
        本引擎的必需性骨干只对功能缺失型有意义, activator 会被拒绝定量。
    """

    name: str
    conc_uM: float
    ic50_uM: float | None = None
    hill: float = 1.0
    target: str | None = None
    modality: str = "inhibitor"

    @property
    def ic50_is_placeholder(self) -> bool:
        return self.ic50_uM is None

    @property
    def effective_ic50(self) -> float:
        return 1.0 if self.ic50_uM is None else float(self.ic50_uM)

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class ReadoutRequest:
    """客户想要什么读出 —— 决定引擎装配哪些项。"""

    #: 剂量-响应曲线 + 选择性窗口 (需要必需性骨干层)
    viability_curve: bool = True
    #: 通路响应谱 (需要 PhenoMap 或 PerturbLens 覆盖)
    pathway_profile: bool = True
    #: 组合非加和性分析 (需要 >=2 个化合物)
    combination: bool = True
    #: 差异基因表 (需要 PhenoMap / PerturbLens 覆盖)
    top_genes: bool = True

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class AssaySpec:
    """一次完整的实验规格。

    这是引擎唯一的输入。同一批模型、不同的 AssaySpec ->
    完全不同的输出结构和数值, 这正是「模型产品」与「数据库」的区别。
    """

    #: 客户关心的细胞上下文
    context: str = "hepg2"
    #: 正常组织参照系 (用于算选择性窗口)
    normal_reference: str = "rpe1"
    compounds: list[Compound] = field(default_factory=list)
    #: 暴露时长 (小时)。**客户参数。**
    exposure_h: float = 72.0
    #: 想要的读出
    readout: ReadoutRequest = field(default_factory=ReadoutRequest)
    #: 客户标签 (只是元数据)
    customer: str = ""
    note: str = ""

    def validate(self) -> list[str]:
        errs: list[str] = []
        if self.context not in CONTEXTS:
            errs.append(
                f"未知细胞上下文 {self.context!r}; 可选: {', '.join(CONTEXTS)}"
            )
        if self.normal_reference not in CONTEXTS:
            errs.append(f"未知正常参照系 {self.normal_reference!r}")
        if not self.compounds:
            errs.append("至少需要一个化合物条目")
        for c in self.compounds:
            if c.conc_uM <= 0:
                errs.append(f"{c.name}: 浓度必须 > 0")
            if c.ic50_uM is not None and c.ic50_uM <= 0:
                errs.append(f"{c.name}: IC50 必须 > 0")
            if c.hill <= 0:
                errs.append(f"{c.name}: Hill 系数必须 > 0")
        if self.exposure_h <= 0:
            errs.append("暴露时长必须 > 0")
        return errs

    def to_dict(self) -> dict[str, Any]:
        return {
            "customer": self.customer,
            "note": self.note,
            "context": self.context,
            "normal_reference": self.normal_reference,
            "exposure_h": self.exposure_h,
            "compounds": [c.to_dict() for c in self.compounds],
            "readout": self.readout.to_dict(),
        }

    @classmethod
    def from_dict(cls, d: dict) -> "AssaySpec":
        return cls(
            context=d.get("context", "hepg2"),
            normal_reference=d.get("normal_reference", "rpe1"),
            compounds=[Compound(**c) for c in d.get("compounds", [])],
            exposure_h=float(d.get("exposure_h", 72.0)),
            readout=ReadoutRequest(**d.get("readout", {})),
            customer=d.get("customer", ""),
            note=d.get("note", ""),
        )
