"""VC-FF · Virtual Cell Force Field.

把 5 个虚拟细胞模型 (ToxSentinel / PerturbLens / PhenoMap / ComboMap / TwinCell)
从「查表式数据库」重构成「可装配的力场」:

    模型提供泛函形式 (响应方向、脆弱性地形、通路签名、耦合项)
    客户提供标量参数 (浓度、IC50、Hill 系数、暴露时间、细胞系、组合)
    引擎按组合律装配   ->  完全客制化的输出

核心 API::

    from vcff import ForceField, AssaySpec, Compound

    ff = ForceField()
    spec = AssaySpec(
        context="hepg2",
        compounds=[Compound(name="sorafenib", conc_uM=10.0, ic50_uM=3.0)],
    )
    card = ff.evaluate(spec)
    print(card.to_json())
"""

from .spec import AssaySpec, Compound, ReadoutRequest
from .engine import ForceField
from .honesty import Evidence, Claim

__all__ = [
    "AssaySpec",
    "Compound",
    "ReadoutRequest",
    "ForceField",
    "Evidence",
    "Claim",
]

__version__ = "0.1.0"
