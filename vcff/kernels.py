"""五个模型 kernel —— 力场的「项」。

每个 kernel 只声明三件事:

  provides   它能贡献力场的哪一项
  accepts    它的输入接口 (适配器要喂什么)
  coverage   它覆盖多大, 覆盖不到时怎么优雅降级

引擎不硬编码任何一个 kernel 的内部; 它只按接口装配。少一个 kernel,
输出就少一层, 而不是崩掉 —— 这是「可装配」的实际含义。
"""

from __future__ import annotations

import json
import re
import statistics
from dataclasses import dataclass
from functools import cached_property
from pathlib import Path
from typing import Any

DATA_DIR = Path(__file__).with_name("data")


def _norm(s: str) -> str:
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


def _load(name: str) -> dict:
    with open(DATA_DIR / f"{name}.json", encoding="utf-8") as fh:
        return json.load(fh)


#: 补充化合物→靶点映射。原始 Tahoe 映射表 (260 条) 漏掉了几个最常用的
#: 细胞毒药物, 这里人工补一小张表。每条都在结果里标注为「补充映射」,
#: 并对多靶点化合物给出警告 —— 单靶点近似对它们是最弱的一环。
SUPPLEMENTARY_TARGETS: dict[str, tuple[str, bool]] = {
    # normalized name -> (canonical gene, is_multi_target)
    "paclitaxel": ("TUBB", False),
    "docetaxel": ("TUBB", False),
    "vincristine": ("TUBB", False),
    "vinblastine": ("TUBB", False),
    "doxorubicin": ("TOP2A", False),
    "daunorubicin": ("TOP2A", False),
    "etoposide": ("TOP2A", False),
    "methotrexate": ("DHFR", False),
    "gemcitabine": ("RRM1", False),
    "hydroxyurea": ("RRM1", False),
    "sorafenib": ("BRAF", True),
    "sunitinibmalate": ("KDR", True),
    "dasatinib": ("ABL1", True),
    "imatinib": ("ABL1", True),
    "nilotinib": ("ABL1", True),
    "midostaurin": ("FLT3", True),
    "cabozantinib": ("MET", True),
    "lenvatinib": ("KDR", True),
}


# --------------------------------------------------------------------------
# K1 · 必需性骨干 (ToxSentinel)
# --------------------------------------------------------------------------


@dataclass
class EssentialityHit:
    gene: str
    dep: dict[str, float]          # {"hepg2": -1.23, "rpe1": -0.11}
    selectivity: float             # hepg2 - rpe1
    klass: str
    via_drug: str | None = None
    via_supplement: bool = False
    multi_target: bool = False


class EssentialityKernel:
    """VC-TOX · ToxSentinel。

    provides : 每个靶点在 HepG2 (肝) 与 RPE1 (正常上皮) 上的适应度效应
               (DepMap Chronos gene effect), 即力场的「必需性项」。
    accepts  : 基因符号 (17,787) 或化合物名 (260 个 Tahoe 化合物 -> 首个靶点)
    coverage : 覆盖最广的一层 —— 任意靶点都能进入毒性主轴。
               只覆盖 HepG2 / RPE1 两个系。
    """

    name = "ToxSentinel"
    code = "VC-TOX"
    provides = "essentiality"
    lines = ("hepg2", "rpe1")
    source = "DepMap Chronos (Tsherniak 2017 / Meyers 2017)"

    def __init__(self) -> None:
        self._d = _load("tox")

    @cached_property
    def n_genes(self) -> int:
        return len(self._d["genes"])

    @cached_property
    def n_drugs(self) -> int:
        return len(self._d["drug2target"])

    def resolve(self, query: str, explicit_target: str | None = None) -> EssentialityHit | None:
        """化合物名或基因符号 -> 必需性。返回 None 表示覆盖不到。"""
        d = self._d
        via = None
        supp = False
        multi = False
        key = _norm(explicit_target) if explicit_target else _norm(query)
        idx = d["norm2idx"].get(key)
        if idx is None and not explicit_target:
            q = _norm(query)
            tgt = d["drug2target"].get(q)
            if tgt is None and q in SUPPLEMENTARY_TARGETS:
                tgt, multi = SUPPLEMENTARY_TARGETS[q]
                supp = True
            if tgt is not None:
                idx = d["norm2idx"].get(_norm(tgt))
                via = tgt
        if idx is None:
            return None
        h = float(d["hepg2"][idx])
        r = float(d["rpe1"][idx])
        return EssentialityHit(
            gene=d["genes"][idx],
            dep={"hepg2": h, "rpe1": r},
            selectivity=h - r,
            klass=self.classify(h, r),
            via_drug=via,
            via_supplement=supp,
            multi_target=multi,
        )

    @staticmethod
    def classify(h: float, r: float) -> str:
        if h < -1 and r > -0.5:
            return "HepG2 选择性必需 (肝特异脆弱点)"
        if h < -1 and r < -1:
            return "两系共同必需 (广谱毒性风险)"
        if r < -1 and h > -0.5:
            return "RPE1 选择性必需 (正常组织脆弱点)"
        return "非必需 (两系依赖性均弱)"

    def suggest(self, prefix: str, limit: int = 12) -> list[str]:
        p = _norm(prefix)
        if not p:
            return []
        out = [g for g in self._d["genes"] if _norm(g).startswith(p)]
        return out[:limit]

    #: 已发表的诚实标注: 必需性与毒性通路激活只是弱耦合。
    reliability = (
        "必需性与毒性通路激活仅弱耦合 (Pearson 0.137) —— "
        "「肝选择性必需」不等于「特异性激活肝毒性通路」, 两者须独立评估。"
    )


# --------------------------------------------------------------------------
# K2 · 转录响应项 (PerturbLens)
# --------------------------------------------------------------------------


class PerturbationKernel:
    """VC-PRT · PerturbLens。

    provides : 单基因扰动后, 指定细胞系的全转录组响应方向 + Hallmark 通路签名,
               即力场的「转录应激项」。这是**单位响应向量**: 数据来自完全敲低,
               引擎按占据率线性缩放它。
    accepts  : (靶基因, 细胞系)
    coverage : 6 个基因 x 4 个系 = 24 个预计算组合。覆盖窄但深 (含实测比对)。
    """

    name = "PerturbLens"
    code = "VC-PRT"
    provides = "transcriptional_response"
    source = "Replogle 2022 Perturb-seq + 跨系共识模型 (Adduri 2025 STATE)"

    def __init__(self) -> None:
        self._d = _load("prt")

    @property
    def genes(self) -> list[str]:
        return list(self._d["genes"])

    @property
    def lines(self) -> list[str]:
        return list(self._d["keys"])

    def get(self, gene: str, line: str) -> dict[str, Any] | None:
        key = f"{gene}|{line}"
        rec = self._d["combos"].get(key)
        if rec is None:
            return None
        return {
            "gene": gene,
            "line": line,
            "live_r": rec["live_r"],
            "ceiling": rec["ceiling"],
            "frac_of_ceiling": rec["live_r"] / rec["ceiling"],
            "pw_names": rec["pw_names"],
            "pw_vals": rec["pw_vals"],
            "up": rec["up"],
            "dn": rec["dn"],
        }

    reliability = (
        "跨系共识 (其余 3 系平均) 加权模型; LCO 留一扰动 Pearson Δ "
        "K562 0.275 / HepG2 0.308 / Jurkat 0.279 / RPE1 0.318; "
        "DEG 方向一致性 0.865–0.894。K562 共识为 EXPLORATORY。"
    )


# --------------------------------------------------------------------------
# K3 · 化学表型项 (PhenoMap)
# --------------------------------------------------------------------------


class ChemPhenotypeKernel:
    """VC-PHE · PhenoMap。

    provides : 化合物 -> 全转录组表型 + Hallmark 通路签名 + MoA 最近邻,
               即力场的「化学表型项」。同样是单位响应向量 (Tahoe 筛选浓度)。
    accepts  : 化合物名
    coverage : 8 个预计算化合物。任意新 SMILES 需要跑 rdkit 指纹模型 (未内置)。
    """

    name = "PhenoMap"
    code = "VC-PHE"
    provides = "chemical_phenotype"
    source = "Tahoe-100M (Zhang 2025) + Morgan 指纹→表型模型"

    def __init__(self) -> None:
        self._d = _load("phe")

    @property
    def drugs(self) -> list[str]:
        return list(self._d["drugs"])

    def get(self, drug: str) -> dict[str, Any] | None:
        for k in self._d["drugs"]:
            if _norm(k) == _norm(drug):
                rec = dict(self._d["per"][k])
                rec["drug"] = k
                return rec
        return None

    reliability = (
        "基因级 LOO Pearson 中位 0.305 / 通路级 0.399; 80 个配对药物上 95% 优于 "
        "L1000 同药检索。已知 MoA 与预测 MoA 不一致时提示需实验确认。"
    )


# --------------------------------------------------------------------------
# K4 · 上位性耦合项 (ComboMap)
# --------------------------------------------------------------------------


@dataclass
class EpistasisHit:
    pair: str
    new_signal_x: float | None
    epistasis: float
    rank: int
    n_cells: int


class EpistasisKernel:
    """VC-CMB · ComboMap。

    provides : 双扰动的非加和性 (上位性), 即力场里唯一的**交叉项**。
               没有它, 组合就只能按 Bliss 独立叠加。
    accepts  : (靶基因A, 靶基因B), 顺序无关
    coverage : Norman 2019 K562 CRISPRa 的 126 个双扰动对。
               命中率低 —— 命中时给实测耦合, 未命中时引擎明确回落到加和并加宽误差带。
    """

    name = "ComboMap"
    code = "VC-CMB"
    provides = "epistasis"
    source = "Norman 2019 Science 双扰动 Perturb-seq (K562, CRISPRa)"

    def __init__(self) -> None:
        self._d = _load("cmb")
        self._index: dict[frozenset[str], int] = {}
        for i, p in enumerate(self._d["pairs"]):
            a, b = p.split("+")
            self._index[frozenset((_norm(a), _norm(b)))] = i

    @cached_property
    def epi_median(self) -> float:
        return statistics.median(self._d["epi"])

    @cached_property
    def epi_iqr(self) -> tuple[float, float]:
        vals = sorted(self._d["epi"])
        n = len(vals)
        return vals[n // 4], vals[(3 * n) // 4]

    @cached_property
    def frac_above_noise(self) -> float:
        """分母是全部 126 对; 缺测 (None) 的对不计入「超过噪声」。

        这样得到 115/126 = 91%, 与原报告的 91% 一致。若只在非缺测的对里算,
        会得到 100% —— 一个由缺失值造成的假象。
        """
        vals = self._d["non"]
        return sum(1 for v in vals if v is not None and v > 1.0) / len(vals)

    @cached_property
    def non_median(self) -> float:
        return statistics.median(v for v in self._d["non"] if v is not None)

    def get(self, gene_a: str, gene_b: str) -> EpistasisHit | None:
        i = self._index.get(frozenset((_norm(gene_a), _norm(gene_b))))
        if i is None:
            return None
        non = self._d["non"][i]
        return EpistasisHit(
            pair=self._d["pairs"][i],
            new_signal_x=None if non is None else float(non),
            epistasis=float(self._d["epi"][i]),
            rank=int(self._d["rank"][i]),
            n_cells=int(self._d["ncells"][i]),
        )

    @property
    def n_pairs(self) -> int:
        return len(self._d["pairs"])

    reliability = (
        "126 个双扰动中新信号中位为噪声的 2.72×, 91% 组合新信号 >1× 噪声; "
        "双扰动可解释方差约 52%。单组合仅数百细胞, 精确上位性数值需谨慎。"
    )


# --------------------------------------------------------------------------
# K5 · 上下文置信权重 (TwinCell)
# --------------------------------------------------------------------------


class ContextKernel:
    """VC-TWN · TwinCell。

    provides : 细胞系作为数字孪生的可信度 (STATE 性能 / oracle 上限),
               在力场里不改变均值, 只调制**误差带宽度**。
    accepts  : 细胞系
    coverage : 4 个系。
    """

    name = "TwinCell"
    code = "VC-TWN"
    provides = "context_confidence"
    source = "STATE vs 数据驱动 oracle 上限 (Adduri 2025)"

    def __init__(self) -> None:
        self._d = _load("twn")

    @property
    def lines(self) -> list[str]:
        return list(self._d["keys"])

    def get(self, line: str) -> dict[str, Any] | None:
        rec = self._d["per"].get(line)
        if rec is None:
            return None
        out = dict(rec)
        out["verdict"] = self.verdict(rec["frac"])
        out["note"] = re.sub(r"<[^>]+>", "", self._d["rel"][line])
        return out

    @staticmethod
    def verdict(frac: float) -> str:
        if frac >= 0.80:
            return "GO · 高置信"
        if frac >= 0.70:
            return "CAUTION · 需验证"
        return "GAP · 数据不足"

    reliability = "达 oracle 上限 ≥80% 判 GO, 70–80% 判 CAUTION。"


ALL_KERNELS = (
    EssentialityKernel,
    PerturbationKernel,
    ChemPhenotypeKernel,
    EpistasisKernel,
    ContextKernel,
)
