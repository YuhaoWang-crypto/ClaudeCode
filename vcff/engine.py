"""ForceField —— 装配器。

按 AssaySpec 决定装配哪些 kernel、以什么权重、能算到哪一层, 然后输出一张
带证据等级的 ResultCard。同一批 kernel + 不同 spec = 不同的输出结构。
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from typing import Any

from . import physics as ph
from .honesty import Claim, Evidence, Ledger
from .kernels import (
    ChemPhenotypeKernel,
    ContextKernel,
    EpistasisKernel,
    EssentialityKernel,
    PerturbationKernel,
)
from .spec import CONTEXTS, ESSENTIALITY_LINES, AssaySpec

#: 装配的层名。tier = 实际装配上的层数, 不是一个固定阶梯 ——
#: 客户可能只有通路层没有必需性层, 引擎要照实说。
LAYER_NAMES = {
    "essentiality": "必需性骨干 (VC-TOX)",
    "pathway": "通路/表型层 (VC-PHE/VC-PRT)",
    "epistasis": "实测上位性交叉项 (VC-CMB)",
}


@dataclass
class CompoundTrace:
    """单个化合物在装配过程中的完整轨迹 —— 让「装配」可见。"""

    name: str
    target: str | None
    resolved_via: str
    occupancy: float
    conc_uM: float
    ic50_uM: float
    ic50_is_placeholder: bool
    hill: float
    dep: dict[str, float] = field(default_factory=dict)
    dep_class: str = ""
    pressure: dict[str, float] = field(default_factory=dict)
    pathway_source: str = "none"
    kernels_used: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "name": self.name,
            "target": self.target,
            "resolved_via": self.resolved_via,
            "conc_uM": self.conc_uM,
            "ic50_uM": self.ic50_uM,
            "ic50_is_placeholder": self.ic50_is_placeholder,
            "hill": self.hill,
            "occupancy": round(self.occupancy, 4),
            "dep": {k: round(v, 4) for k, v in self.dep.items()},
            "dep_class": self.dep_class,
            "pressure": {k: round(v, 4) for k, v in self.pressure.items()},
            "pathway_source": self.pathway_source,
            "kernels_used": self.kernels_used,
            "warnings": self.warnings,
        }


@dataclass
class ResultCard:
    """结果卡 —— 可直接进客户流水线的 JSON。"""

    spec: dict
    tier: int
    tier_label: str
    compounds: list[CompoundTrace]
    viability_at_spec: dict[str, float]
    dose_curve: dict[str, Any]
    selectivity: dict[str, Any]
    pathway: dict[str, Any]
    combination: dict[str, Any]
    confidence: dict[str, Any]
    ledger: dict
    errors: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "engine": "VC-FF v0.1.0",
            "spec": self.spec,
            "coverage_tier": self.tier,
            "coverage_tier_label": self.tier_label,
            "compound_trace": [c.to_dict() for c in self.compounds],
            "pseudo_viability_at_spec_dose": self.viability_at_spec,
            "dose_response": self.dose_curve,
            "selectivity": self.selectivity,
            "pathway_response": self.pathway,
            "combination": self.combination,
            "confidence": self.confidence,
            "honesty_ledger": self.ledger,
            "errors": self.errors,
        }

    def to_json(self, indent: int = 2) -> str:
        return json.dumps(self.to_dict(), ensure_ascii=False, indent=indent)


class ForceField:
    """把 5 个 kernel 装配成一个可参数化的引擎。"""

    def __init__(self, kernels: dict | None = None) -> None:
        k = kernels or {}
        self.ess: EssentialityKernel = k.get("ess") or EssentialityKernel()
        self.prt: PerturbationKernel = k.get("prt") or PerturbationKernel()
        self.phe: ChemPhenotypeKernel = k.get("phe") or ChemPhenotypeKernel()
        self.cmb: EpistasisKernel = k.get("cmb") or EpistasisKernel()
        self.twn: ContextKernel = k.get("twn") or ContextKernel()

    # -- 接口自述 ---------------------------------------------------------

    def describe(self) -> dict:
        return {
            "engine": "VC-FF · Virtual Cell Force Field v0.1.0",
            "kernels": [
                {
                    "code": k.code,
                    "name": k.name,
                    "provides": k.provides,
                    "source": k.source,
                    "reliability": k.reliability,
                }
                for k in (self.ess, self.prt, self.phe, self.cmb, self.twn)
            ],
            "composition_laws": [
                "L1 剂量→占据率  θ = C^h/(C^h+IC50^h)          [Hill, MODELED]",
                "L2 占据率→表型  响应 = θ × 完全敲低单位向量    [线性缩放, MODELED]",
                "L3 多药叠加     log V = −λ Σ θᵢ·压力ᵢ          [Bliss 独立, MODELED]",
                "L4 交叉项       log V ← γ(上位性) × log V      [ComboMap 修正, MODELED]",
            ],
            "customer_parameters": [
                "浓度 conc_uM", "效价 ic50_uM", "Hill 系数", "暴露时长",
                "细胞上下文", "正常参照系", "组合设计",
            ],
            "coverage": {
                "essentiality_genes": self.ess.n_genes,
                "mapped_compounds": self.ess.n_drugs,
                "essentiality_lines": list(ESSENTIALITY_LINES),
                "perturbation_combos": len(self.prt.genes) * len(self.prt.lines),
                "chem_phenotype_drugs": len(self.phe.drugs),
                "epistasis_pairs": self.cmb.n_pairs,
                "contexts": list(CONTEXTS),
            },
        }

    # -- 主入口 -----------------------------------------------------------

    def evaluate(self, spec: AssaySpec) -> ResultCard:
        errs = spec.validate()
        led = Ledger()
        if errs:
            return ResultCard(
                spec=spec.to_dict(), tier=0, tier_label="T0 · 规格无效",
                compounds=[], viability_at_spec={}, dose_curve={},
                selectivity={}, pathway={}, combination={},
                confidence={}, ledger=led.to_dict(), errors=errs,
            )

        traces = [self._resolve(c, spec) for c in spec.compounds]
        tier, tier_label = self._tier(traces, spec)
        has_ess = any(t.dep for t in traces)

        via = self._viability_at_spec(traces, spec, led, has_ess)
        curve = (self._dose_curve(traces, spec, led)
                 if spec.readout.viability_curve and has_ess else {})
        selec = self._selectivity(curve, traces, spec, led) if curve else {}
        pw = self._pathway(traces, spec, led) if spec.readout.pathway_profile else {}
        combo = (self._combination(traces, spec, led)
                 if spec.readout.combination else {})
        conf = self._confidence(traces, spec, tier, tier_label, led)

        self._global_notes(led, traces, spec, has_ess)

        return ResultCard(
            spec=spec.to_dict(), tier=tier, tier_label=tier_label,
            compounds=traces, viability_at_spec=via, dose_curve=curve,
            selectivity=selec, pathway=pw, combination=combo,
            confidence=conf, ledger=led.to_dict(), errors=[],
        )

    # -- 装配步骤 ---------------------------------------------------------

    def _resolve(self, c, spec: AssaySpec) -> CompoundTrace:
        theta = ph.occupancy(c.conc_uM, c.effective_ic50, c.hill)
        tr = CompoundTrace(
            name=c.name, target=c.target, resolved_via="",
            occupancy=theta, conc_uM=c.conc_uM,
            ic50_uM=c.effective_ic50,
            ic50_is_placeholder=c.ic50_is_placeholder,
            hill=c.hill,
        )
        if c.ic50_is_placeholder:
            tr.warnings.append(
                "IC50 未提供, 使用 1.0 µM 占位。效价是你的输入参数, 引擎不预测效价; "
                "占据率与由它推出的一切数值都随该占位值线性漂移。"
            )
        if c.modality != "inhibitor":
            tr.warnings.append(
                f"modality={c.modality}: 必需性骨干只对功能缺失型 (抑制/敲低) 有意义, "
                "激活型化合物的杀伤压力不予定量。"
            )

        hit = self.ess.resolve(c.name, c.target)
        if hit is None:
            tr.resolved_via = "未命中"
            tr.warnings.append(
                f"「{c.name}」既不在 17,787 基因表中, 也无法映射到已知化合物靶点; "
                "请直接指定 target=基因符号。该化合物不进入毒性主轴。"
            )
            return tr

        tr.target = hit.gene
        if hit.via_drug and hit.via_supplement:
            tr.resolved_via = f"补充映射 (人工整理) → {hit.via_drug}"
        elif hit.via_drug:
            tr.resolved_via = f"Tahoe 化合物→靶点映射 → {hit.via_drug}"
        elif c.target:
            tr.resolved_via = "显式指定靶点"
        else:
            tr.resolved_via = "基因符号直接命中"
        if hit.multi_target:
            tr.warnings.append(
                f"{c.name} 是多靶点化合物, 这里只按主靶点 {hit.gene} 建模。"
                "单靶点近似是本次装配最弱的一环, 实际毒性可能由未建模的次级靶点主导。"
            )
        tr.dep = hit.dep
        tr.dep_class = hit.klass
        tr.kernels_used.append("VC-TOX")
        if c.modality == "inhibitor":
            tr.pressure = {
                ln: ph.kill_pressure(hit.dep[ln], theta) for ln in ESSENTIALITY_LINES
            }
        else:
            tr.pressure = {ln: 0.0 for ln in ESSENTIALITY_LINES}

        if self.phe.get(c.name) is not None:
            tr.pathway_source = "VC-PHE"
            tr.kernels_used.append("VC-PHE")
        elif self.prt.get(hit.gene, spec.context) is not None:
            tr.pathway_source = "VC-PRT"
            tr.kernels_used.append("VC-PRT")
        else:
            tr.warnings.append(
                f"通路层未覆盖: 「{c.name}」不在 PhenoMap 的 8 个化合物中, "
                f"靶点 {hit.gene} 也不在 PerturbLens 的 6 个基因中。本化合物只贡献必需性项。"
            )
        return tr

    def _tier(self, traces: list[CompoundTrace], spec: AssaySpec) -> tuple[int, str]:
        layers: list[str] = []
        if any(t.dep for t in traces):
            layers.append(LAYER_NAMES["essentiality"])
        if any(t.pathway_source != "none" for t in traces):
            layers.append(LAYER_NAMES["pathway"])
        if self._gamma(traces)[1] is not None:
            layers.append(LAYER_NAMES["epistasis"])
        if not layers:
            return 0, "T0 · 无任何层可装配 (引擎拒绝定量)"
        label = f"T{len(layers)} · " + " + ".join(layers)
        if LAYER_NAMES["essentiality"] not in layers:
            label += "  [无必需性骨干 → 不出剂量-响应与选择性窗口]"
        return len(layers), label

    def _gamma(self, traces: list[CompoundTrace]) -> tuple[float, Any]:
        tg = [t.target for t in traces if t.target]
        for i in range(len(tg)):
            for j in range(i + 1, len(tg)):
                hit = self.cmb.get(tg[i], tg[j])
                if hit is not None:
                    return ph.epistasis_gamma(hit.epistasis, self.cmb.epi_median), hit
        return 1.0, None

    def _viability_at_spec(
        self, traces, spec, led: Ledger, has_ess: bool
    ) -> dict[str, float]:
        if not has_ess:
            led.add(Claim(
                key="pseudo_viability", label="伪存活率 / 剂量-响应", value="不可用",
                evidence=Evidence.UNSUPPORTED,
                basis="没有任何化合物解析到 Chronos 依赖性数据",
                caveat="注意: 引擎在这里返回「不可用」而不是 1.0。"
                       "「算不出杀伤」与「无毒」是两回事。",
            ))
            return {}
        gamma, hit = self._gamma(traces)
        out = {}
        for ln in ESSENTIALITY_LINES:
            press = [t.pressure.get(ln, 0.0) for t in traces]
            lv = ph.bliss_log_viability(press, spec.exposure_h, gamma)
            out[ln] = round(ph.viability(lv), 4)
        led.add(Claim(
            key="pseudo_viability", label="规格剂量下的伪存活率",
            value=out, evidence=Evidence.MODELED,
            basis="Chronos gene effect (✅实测) 经 Hill 占据率 + 指数杀伤 (⚠️假设) 组合",
            caveat="λ 为约定锚定 (ln2), 未拟合任何存活率数据。绝对百分比无意义; "
                   "只有系间比值与剂量位移可解释。",
        ))
        return out

    def _dose_curve(self, traces, spec, led: Ledger) -> dict:
        gamma, _ = self._gamma(traces)
        ms = ph.log_grid(1e-3, 1e3, 121)
        curves: dict[str, list[float]] = {}
        for ln in ESSENTIALITY_LINES:
            vals = []
            for m in ms:
                press = []
                for t, c in zip(traces, spec.compounds):
                    if not t.dep:
                        continue
                    th = ph.occupancy(c.conc_uM * m, c.effective_ic50, c.hill)
                    press.append(ph.kill_pressure(t.dep[ln], th)
                                 if c.modality == "inhibitor" else 0.0)
                vals.append(ph.viability(
                    ph.bliss_log_viability(press, spec.exposure_h, gamma)))
            curves[ln] = [round(v, 5) for v in vals]
        return {
            "x_label": "剂量倍数 (× 你规格里的浓度组合)",
            "multipliers": [round(m, 6) for m in ms],
            "viability": curves,
            "gamma_applied": round(gamma, 4),
        }

    def _selectivity(self, curve: dict, traces, spec, led: Ledger) -> dict:
        scan = ph.DoseScan(curve["multipliers"], curve["viability"])
        gamma = curve["gamma_applied"]

        # 主指标: 等效剂量比。在**同一个抑制水平**上比两条曲线的剂量。
        # 两系完全相同时恰好等于 1.0 —— 这是「无选择性」的正确零点。
        # 取两条曲线都能达到的最高抑制水平, 因为强抑制常常压根到不了。
        # 比值随抑制水平变化是真实药理 (占据率趋饱和时正常系需要不成比例的剂量),
        # 所以固定几个水平各报一次, 不同规格之间才可比。
        by_level: dict[str, float | None] = {}
        iso_ratio = iso_level = m_t = m_n = None
        for lvl in (0.5, 0.3, 0.2, 0.1):
            a = scan.crossing("hepg2", 1.0 - lvl)
            b = scan.crossing("rpe1", 1.0 - lvl)
            by_level[f"{lvl:.0%}"] = round(b / a, 3) if (a and b) else None
            if a and b and iso_ratio is None:
                iso_ratio, iso_level, m_t, m_n = b / a, lvl, a, b

        # 次指标 (更保守, 治疗指数式): 正常系刚开始受损 ÷ 肝癌系达半数抑制。
        m50_t = scan.crossing("hepg2", 0.5)
        m90_n = scan.crossing("rpe1", 0.9)
        window = (m90_n / m50_t) if (m50_t and m90_n) else None

        # 饱和下限: 占据率 → 1 时的伪存活率。这是「靠 on-target 必需性
        # 最多能杀到哪」的解析上限, 也是曲线不跨阈值时唯一诚实的说法。
        floor, sat_press = {}, {}
        for ln in ESSENTIALITY_LINES:
            p = [max(0.0, -t.dep[ln]) for t in traces if t.dep]
            sat_press[ln] = round(sum(p), 4)
            floor[ln] = round(
                ph.viability(ph.bliss_log_viability(p, spec.exposure_h, gamma)), 4)

        ratio = (sat_press["hepg2"] / sat_press["rpe1"]
                 if sat_press["rpe1"] > 0 else None)
        diff = round(sat_press["hepg2"] - sat_press["rpe1"], 4)

        if iso_ratio:
            if iso_ratio > 1.2:
                interp = (f"等效剂量比 {iso_ratio:.2f}× —— 正常系要 {iso_ratio:.2f} 倍"
                          f"剂量才达到同样的 {iso_level:.0%} 抑制, 存在选择性")
            elif iso_ratio < 0.83:
                interp = (f"等效剂量比 {iso_ratio:.2f}× < 1 —— 正常系比肝癌系更敏感, "
                          f"选择性方向是反的")
            else:
                interp = (f"等效剂量比 {iso_ratio:.2f}× ≈ 1 —— 两系敏感性基本相同, "
                          f"无选择性 (广谱杀伤)")
            reason = None
        elif floor["hepg2"] > 0.5:
            interp = "等效剂量比无法定义 —— 见 saturation_floor"
            reason = (
                f"即使占据率饱和 (剂量 →∞), 靶点必需性最多把 HepG2 伪存活率压到 "
                f"{floor['hepg2']:.3f}, 到不了 50%。含义: **这些靶点的 on-target "
                f"必需性不足以解释 50% 的细胞杀伤**。若你在实验里确实看到 ≥50% 杀伤, "
                f"那部分毒性来自必需性以外的机制 —— 脱靶、活性代谢物、非依赖性应激 —— "
                f"必需性轴对它没有预测力, 应另设对照。"
            )
        else:
            interp = "等效剂量比无法定义 —— 有一条曲线在扫描范围内未达 10% 抑制"
            reason = ("在 10⁻³–10³ × 规格剂量内两条曲线没有共同可比的抑制水平; "
                      "选择性的判据落在扫描范围之外。")

        led.add(Claim(
            key="iso_effect_ratio",
            label=(f"等效剂量比 (正常系 ÷ 肝癌系, 同为 {iso_level:.0%} 抑制)"
                   if iso_level else "等效剂量比"),
            value=(round(iso_ratio, 3) if iso_ratio else "无法定义"),
            unit="×" if iso_ratio else "", evidence=Evidence.MODELED,
            basis="HepG2 vs RPE1 Chronos 差值 (✅实测) 经组合律传播",
            caveat="对称指标: 两系敏感性相同时恰为 1.0。它是**比值**, 是曲线类量里"
                   "最可解释的; 但它不是治疗指数, 不能替代体内 PK/PD 与毒理。",
        ))
        led.add(Claim(
            key="saturation_pressure_ratio",
            label="饱和杀伤压力比 (HepG2 ÷ RPE1)",
            value=(round(ratio, 3) if ratio else "RPE1 无压力, 比值发散"),
            unit="×" if ratio else "", evidence=Evidence.ANCHORED,
            basis="两个 Chronos gene effect 的直接比值, 不经过任何剂量假设",
            caveat="这是本卡上唯一**完全不依赖建模假设**的选择性量 —— "
                   "它与浓度、IC50、Hill、暴露时长都无关。",
        ))
        return {
            "target_line": "hepg2", "normal_line": "rpe1",
            "iso_effect_ratio_x": (round(iso_ratio, 3) if iso_ratio else None),
            "iso_effect_level": iso_level,
            "iso_effect_ratio_by_level": by_level,
            "iso_dose_target": (round(m_t, 4) if m_t else None),
            "iso_dose_normal": (round(m_n, 4) if m_n else None),
            "dose_multiplier_target_IC50": (round(m50_t, 4) if m50_t else None),
            "dose_multiplier_normal_10pct": (round(m90_n, 4) if m90_n else None),
            "conservative_window_x": (round(window, 3) if window else None),
            "saturation_floor": floor,
            "saturation_pressure": sat_press,
            "saturation_pressure_ratio": (round(ratio, 3) if ratio else None),
            "saturation_pressure_diff": diff,
            "interpretation": interp,
            "why_undefined": reason,
        }

    def _pathway(self, traces, spec, led: Ledger) -> dict:
        terms, srcs = [], []
        for t, c in zip(traces, spec.compounds):
            rec = self.phe.get(c.name)
            if rec is not None:
                terms.append((t.occupancy, rec["pw_names"], rec["pw_vals"]))
                srcs.append(f"{c.name} → VC-PHE (LOO r={rec['loo_r']:.3f})")
                continue
            if t.target:
                rec2 = self.prt.get(t.target, spec.context)
                if rec2 is not None:
                    terms.append((t.occupancy, rec2["pw_names"], rec2["pw_vals"]))
                    srcs.append(
                        f"{c.name} → VC-PRT {t.target}@{spec.context} "
                        f"(live r={rec2['live_r']:.3f}/上限 {rec2['ceiling']:.3f})"
                    )
        if not terms:
            led.add(Claim(
                key="pathway", label="通路响应谱", value="不可用",
                evidence=Evidence.UNSUPPORTED,
                basis="无任何化合物落入 PhenoMap(8 化合物) 或 PerturbLens(6 基因×4 系) 覆盖",
                caveat="引擎在这里拒绝外推, 而不是给一个看起来合理的数。",
            ))
            return {"available": False, "sources": [], "names": [], "values": []}

        names, vals = ph.combine_pathway_vectors(terms)
        led.add(Claim(
            key="pathway", label="通路响应谱 (占据率加权叠加)",
            value=f"{len(names)} 条 Hallmark", evidence=Evidence.PREDICTED,
            basis="; ".join(srcs),
            caveat="单位响应向量来自完全敲低/筛选浓度, 按占据率线性缩放是 ⚠️ 假设; "
                   "多化合物按通路名相加没有考虑通路间串扰。",
        ))
        return {
            "available": True, "sources": srcs,
            "names": names, "values": [round(v, 5) for v in vals],
            "top_down": names[:5], "top_up": names[-5:][::-1],
        }

    def _combination(self, traces, spec, led: Ledger) -> dict:
        n = len([t for t in traces if t.dep])
        if n < 2:
            return {"applicable": False,
                    "reason": "少于 2 个可定量化合物, 无组合项"}
        gamma, hit = self._gamma(traces)
        press = {ln: [t.pressure.get(ln, 0.0) for t in traces]
                 for ln in ESSENTIALITY_LINES}
        add = {ln: round(ph.viability(ph.bliss_log_viability(p, spec.exposure_h, 1.0)), 4)
               for ln, p in press.items()}
        cor = {ln: round(ph.viability(ph.bliss_log_viability(p, spec.exposure_h, gamma)), 4)
               for ln, p in press.items()}

        if hit is not None:
            led.add(Claim(
                key="epistasis", label=f"实测上位性交叉项 ({hit.pair})",
                value=hit.epistasis, evidence=Evidence.ANCHORED,
                basis=f"Norman 2019 双扰动 Perturb-seq, 排名 {hit.rank}/{self.cmb.n_pairs}, "
                      f"{hit.n_cells} 细胞",
                caveat="原始测量在 K562 CRISPRa 遗传扰动上; 迁移到你的细胞系与化学抑制 "
                       "是 ⚠️ 外推, 且单组合仅数百细胞。",
            ))
            return {
                "applicable": True, "evidence": "measured",
                "pair": hit.pair, "epistasis": hit.epistasis,
                "new_signal_x": hit.new_signal_x, "rank": hit.rank,
                "n_cells": hit.n_cells, "gamma": round(gamma, 4),
                "bliss_additive": add, "epistasis_corrected": cor,
                "deviation_x": {
                    ln: round(cor[ln] / add[ln], 3) if add[ln] > 0 else None
                    for ln in add
                },
            }

        lo, hi = self.cmb.epi_iqr
        g_lo = ph.epistasis_gamma(lo, self.cmb.epi_median)
        g_hi = ph.epistasis_gamma(hi, self.cmb.epi_median)
        band = {
            ln: sorted([
                round(ph.viability(ph.bliss_log_viability(p, spec.exposure_h, g_lo)), 4),
                round(ph.viability(ph.bliss_log_viability(p, spec.exposure_h, g_hi)), 4),
            ]) for ln, p in press.items()
        }
        led.add(Claim(
            key="epistasis", label="上位性交叉项", value="无实测证据, 回落到 Bliss 加和",
            evidence=Evidence.MODELED,
            basis=f"该靶点对不在 ComboMap 的 {self.cmb.n_pairs} 对里",
            caveat=f"Norman 2019 中 {self.cmb.frac_above_noise*100:.0f}% 的组合新信号 "
                   f">1× 噪声 (中位 {self.cmb.non_median:.2f}×), 即加和预期系统性偏保守; "
                   "误差带按该经验分布的四分位给出, 而非声称加和成立。",
        ))
        return {
            "applicable": True, "evidence": "no_measurement_fallback_additive",
            "pair": None, "gamma": 1.0,
            "bliss_additive": add, "epistasis_corrected": add,
            "empirical_band": band,
            "empirical_prior": {
                "n_pairs": self.cmb.n_pairs,
                "frac_new_signal_above_noise": round(self.cmb.frac_above_noise, 3),
                "median_new_signal_x": round(self.cmb.non_median, 3),
            },
        }

    def _confidence(self, traces, spec, tier: int, tier_label: str, led: Ledger) -> dict:
        ctx = self.twn.get(spec.context) or {}
        nrm = self.twn.get(spec.normal_reference) or {}
        frac = float(ctx.get("frac", 0.0))
        cov = sum(1 for t in traces if t.dep) / max(1, len(traces))
        pw_cov = sum(1 for t in traces if t.pathway_source != "none") / max(1, len(traces))
        placeholder = any(t.ic50_is_placeholder for t in traces)

        score = frac * (0.5 + 0.3 * cov + 0.2 * pw_cov)
        if placeholder:
            score *= 0.6
        led.add(Claim(
            key="confidence", label="综合置信度",
            value=round(score, 3), evidence=Evidence.MODELED,
            basis=f"TwinCell {spec.context} 达 oracle 上限 {frac:.1%} (✅) × 覆盖度",
            caveat="这是一个用于排序的启发式分数, 不是校准过的概率。",
        ))
        return {
            "context": spec.context, "context_label": CONTEXTS.get(spec.context, ""),
            "twin_readiness": ctx.get("verdict"), "frac_of_oracle": frac,
            "context_note": ctx.get("note", ""),
            "normal_reference_readiness": nrm.get("verdict"),
            "essentiality_coverage": round(cov, 3),
            "pathway_coverage": round(pw_cov, 3),
            "ic50_placeholder_used": placeholder,
            "composite_score": round(score, 3),
            "tier": tier, "tier_label": tier_label,
        }

    def _global_notes(self, led: Ledger, traces, spec, has_ess: bool) -> None:
        if has_ess and all(
            all(v >= 0 for v in t.dep.values()) for t in traces if t.dep
        ):
            led.note(
                "所有已解析靶点在 HepG2 与 RPE1 上的 Chronos 均为非负 —— "
                "敲除它们不损伤 (甚至促进) 生长。必需性轴在此规格下不产生任何杀伤信号, "
                "伪存活率恒为 1 并**不代表无毒**, 只代表这条轴对该组合无预测力。"
            )
        if spec.context not in ESSENTIALITY_LINES:
            led.note(
                f"毒性/选择性主轴固定在 HepG2 ↔ RPE1 上 —— Chronos 依赖性数据只覆盖这两个系。"
                f"你选的上下文 {spec.context} 只参与通路层与置信度加权, 不改变主轴数值。"
            )
        led.note(
            "力场类比的边界: 与分子力场只有相对能量可解释一样, 本引擎只有"
            "系间比值、剂量位移、相对加和预期的偏离是可解释的。"
        )
        led.note(self.ess.reliability)
        if any(t.pathway_source == "VC-PRT" for t in traces) and spec.context == "k562":
            led.note("K562 的跨系共识被原报告标为 EXPLORATORY, 通路层外推需额外验证。")
