"""不变量测试。

测的不是「数值对不对」(没有金标准可对), 而是**引擎的行为约束**:
对称性、单调性、拒答边界、以及与已发表统计量的一致性。
这些是把力场当产品用时真正要保证的东西。

    python -m vcff.tests
"""

from __future__ import annotations

import math
import sys

from .engine import ForceField
from .kernels import EpistasisKernel
from .physics import occupancy
from .spec import AssaySpec, Compound

FAILS: list[str] = []


def check(name: str, cond: bool, detail: str = "") -> None:
    if cond:
        print(f"  ok   {name}")
    else:
        print(f"  FAIL {name} {detail}")
        FAILS.append(name)


def approx(a: float, b: float, tol: float = 1e-6) -> bool:
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def main() -> int:
    ff = ForceField()
    print("VC-FF 不变量测试\n")

    # --- L1 Hill 占据率 ---------------------------------------------------
    print("[L1 占据率]")
    check("C=IC50 时 θ=0.5", approx(occupancy(1.0, 1.0), 0.5))
    check("θ 有界 [0,1]", all(0.0 <= occupancy(c, 1.0, 2.0) <= 1.0
                             for c in (0, 1e-9, 1e-3, 1, 1e3, 1e9)))
    check("θ 对浓度单调不减",
          all(occupancy(c, 1.0) <= occupancy(c * 1.01, 1.0)
              for c in (0.01, 0.1, 1, 10, 100)))
    check("θ 只依赖 C/IC50 的比",
          approx(occupancy(2.0, 4.0), occupancy(0.5, 1.0)))

    # --- 对称性: 两系相同 -> 等效剂量比 = 1 -------------------------------
    print("\n[选择性零点]")
    # PSMB5 在两系上 Chronos 几乎相同 (−1.847 / −1.873)
    card = ff.evaluate(AssaySpec(
        context="hepg2", compounds=[
            Compound(name="Bortezomib", conc_uM=0.5, ic50_uM=0.1)]))
    r = card.selectivity["iso_effect_ratio_x"]
    check("两系依赖性几乎相同时, 等效剂量比 ≈ 1", 0.9 < r < 1.1, f"(得到 {r})")

    # --- 无量纲不变量: 压力比不受剂量假设影响 -----------------------------
    print("\n[剂量假设不变性]")
    base = dict(context="hepg2", normal_reference="rpe1")
    a = ff.evaluate(AssaySpec(**base, exposure_h=72.0, compounds=[
        Compound(name="Pemetrexed", conc_uM=0.5, ic50_uM=0.1)]))
    b = ff.evaluate(AssaySpec(**base, exposure_h=144.0, compounds=[
        Compound(name="Pemetrexed", conc_uM=50.0, ic50_uM=3.0, hill=2.0)]))
    check("饱和压力比与浓度/IC50/Hill/时长全部无关",
          approx(a.selectivity["saturation_pressure_ratio"],
                 b.selectivity["saturation_pressure_ratio"]),
          f"({a.selectivity['saturation_pressure_ratio']} vs "
          f"{b.selectivity['saturation_pressure_ratio']})")
    check("同药不同剂量确实给出不同的伪存活率",
          a.viability_at_spec["hepg2"] != b.viability_at_spec["hepg2"])

    # --- IC50 缩放: 曲线整体平移 -----------------------------------------
    print("\n[剂量-响应平移]")
    c1 = ff.evaluate(AssaySpec(**base, compounds=[
        Compound(name="Bortezomib", conc_uM=1.0, ic50_uM=0.1)]))
    c2 = ff.evaluate(AssaySpec(**base, compounds=[
        Compound(name="Bortezomib", conc_uM=1.0, ic50_uM=1.0)]))
    m1 = c1.selectivity["dose_multiplier_target_IC50"]
    m2 = c2.selectivity["dose_multiplier_target_IC50"]
    check("IC50 放大 10 倍, 达半效所需剂量倍数也放大 10 倍",
          approx(m2 / m1, 10.0, tol=0.02), f"({m2/m1:.3f})")

    # --- 拒答边界 ---------------------------------------------------------
    print("\n[拒答边界]")
    none_card = ff.evaluate(AssaySpec(**base, compounds=[
        Compound(name="完全查不到的化合物", conc_uM=1.0, ic50_uM=1.0)]))
    check("无覆盖时不输出伪存活率 (而不是给 1.0)",
          none_card.viability_at_spec == {})
    check("无覆盖时不输出剂量曲线", none_card.dose_curve == {})
    check("无覆盖时装配层数为 0", none_card.tier == 0)
    lv = [c for c in none_card.ledger["claims"] if c["key"] == "pseudo_viability"]
    check("无覆盖时给出 UNSUPPORTED 证据等级",
          bool(lv) and lv[0]["evidence"] == "UNSUPPORTED")

    nopw = ff.evaluate(AssaySpec(**base, compounds=[
        Compound(name="cmpd", target="PLK4", conc_uM=1.0, ic50_uM=0.2)]))
    check("通路层无覆盖时明确标 available=False",
          nopw.pathway["available"] is False)

    # --- 正向依赖性: 不产生杀伤, 但要说清楚 -------------------------------
    print("\n[非必需靶点]")
    pos = ff.evaluate(AssaySpec(**base, compounds=[
        Compound(name="cmpd", target="CDKN1A", conc_uM=2.0, ic50_uM=0.5)]))
    check("敲除促生长的靶点 -> 伪存活率为 1",
          approx(pos.viability_at_spec["hepg2"], 1.0))
    check("并且账本里明确说「不代表无毒」",
          any("不代表无毒" in n for n in pos.ledger["notes"]))

    # --- 组合项 -----------------------------------------------------------
    print("\n[组合与上位性]")
    ek = EpistasisKernel()
    check("上位性查询与靶点顺序无关",
          ek.get("PLK4", "STIL") is not None and ek.get("STIL", "PLK4") is not None)
    check("未测过的对返回 None", ek.get("PLK4", "TYMS") is None)
    check("91% 的组合新信号超过噪声底 (与原报告一致)",
          approx(round(ek.frac_above_noise, 2), 0.91, tol=0.02),
          f"({ek.frac_above_noise:.3f})")

    hit = ff.evaluate(AssaySpec(**base, compounds=[
        Compound(name="a", target="PLK4", conc_uM=1.0, ic50_uM=0.2),
        Compound(name="b", target="STIL", conc_uM=1.0, ic50_uM=0.3)]))
    check("命中实测对时证据等级为 measured",
          hit.combination["evidence"] == "measured")
    check("命中时上位性 Claim 标为 ANCHORED",
          any(c["key"] == "epistasis" and c["evidence"] == "ANCHORED"
              for c in hit.ledger["claims"]))

    miss = ff.evaluate(AssaySpec(**base, compounds=[
        Compound(name="Doxorubicin", conc_uM=0.5, ic50_uM=0.05),
        Compound(name="Pemetrexed", conc_uM=0.5, ic50_uM=0.1)]))
    check("未命中时回落到加和并给经验区间",
          miss.combination["evidence"] == "no_measurement_fallback_additive"
          and "empirical_band" in miss.combination)
    check("未命中时上位性 Claim 标为 MODELED",
          any(c["key"] == "epistasis" and c["evidence"] == "MODELED"
              for c in miss.ledger["claims"]))
    check("加和预期不高于任一单药 (Bliss 单调性)",
          miss.combination["bliss_additive"]["hepg2"] <= min(
              ff.evaluate(AssaySpec(**base, compounds=[
                  Compound(name="Doxorubicin", conc_uM=0.5, ic50_uM=0.05)]
              )).viability_at_spec["hepg2"],
              ff.evaluate(AssaySpec(**base, compounds=[
                  Compound(name="Pemetrexed", conc_uM=0.5, ic50_uM=0.1)]
              )).viability_at_spec["hepg2"]) + 1e-9)

    # --- 上下文 -----------------------------------------------------------
    print("\n[上下文]")
    jur = ff.evaluate(AssaySpec(context="jurkat", normal_reference="rpe1",
                                compounds=[Compound(name="a", target="PLK4",
                                                    conc_uM=1.0, ic50_uM=0.2)]))
    hep = ff.evaluate(AssaySpec(context="hepg2", normal_reference="rpe1",
                                compounds=[Compound(name="a", target="PLK4",
                                                    conc_uM=1.0, ic50_uM=0.2)]))
    check("换上下文不改变毒性主轴数值 (主轴只覆盖 HepG2/RPE1)",
          jur.viability_at_spec == hep.viability_at_spec)
    check("换到主轴不覆盖的系时账本要说明",
          any("只参与通路层" in n for n in jur.ledger["notes"]))
    check("置信度随上下文变化",
          jur.confidence["frac_of_oracle"] != hep.confidence["frac_of_oracle"])

    # --- IC50 占位 --------------------------------------------------------
    print("\n[占位参数]")
    ph_card = ff.evaluate(AssaySpec(**base, compounds=[
        Compound(name="Bortezomib", conc_uM=1.0)]))
    check("缺 IC50 时标记占位", ph_card.confidence["ic50_placeholder_used"] is True)
    check("缺 IC50 时置信度被扣分",
          ph_card.confidence["composite_score"] < c1.confidence["composite_score"])

    # --- 规格校验 ---------------------------------------------------------
    print("\n[规格校验]")
    bad = ff.evaluate(AssaySpec(context="不存在的系", compounds=[]))
    check("非法规格返回错误而不是数字", len(bad.errors) >= 2)

    # --- 序列化 -----------------------------------------------------------
    print("\n[序列化]")
    import json
    js = hit.to_json()
    check("结果卡可 JSON 往返", isinstance(json.loads(js), dict))
    check("结果卡含诚实账本", "honesty_ledger" in json.loads(js))
    check("结果卡含拒绝输出清单",
          len(json.loads(js)["honesty_ledger"]["refused_outputs"]) == 4)

    print(f"\n{'全部通过' if not FAILS else str(len(FAILS)) + ' 项失败: ' + ', '.join(FAILS)}")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
