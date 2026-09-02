"""命令行接口 —— 力场的第三种适配器 (库 / JSON / CLI)。

    python -m vcff.cli describe
    python -m vcff.cli scenarios
    python -m vcff.cli run A_hepatotox_single
    python -m vcff.cli run A_hepatotox_single --json
    echo '{"context":"hepg2","compounds":[...]}' | python -m vcff.cli stdin
"""

from __future__ import annotations

import argparse
import json
import sys

from .engine import ForceField
from .honesty import Evidence
from .scenarios import SCENARIOS
from .spec import AssaySpec

MARK = {"ANCHORED": "✅", "PREDICTED": "◐", "MODELED": "⚠️", "UNSUPPORTED": "⛔"}


def _fmt(v, n=4):
    return f"{v:.{n}g}" if isinstance(v, float) else str(v)


def render(card) -> str:
    d = card.to_dict()
    L: list[str] = []
    add = L.append
    sp = d["spec"]
    add("=" * 74)
    add(f"VC-FF 结果卡 · {sp['customer'] or '(未命名客户)'}")
    if sp["note"]:
        add(f"  {sp['note']}")
    add("=" * 74)

    if d["errors"]:
        add("规格错误:")
        for e in d["errors"]:
            add(f"  ⛔ {e}")
        return "\n".join(L)

    add(f"\n【输入规格】上下文 {sp['context']} · 正常参照 {sp['normal_reference']} "
        f"· 暴露 {sp['exposure_h']:g} h")
    for c in sp["compounds"]:
        ic = "占位 1.0" if c["ic50_uM"] is None else f"{c['ic50_uM']:g}"
        add(f"  · {c['name']:<20} {c['conc_uM']:>8g} µM   IC50 {ic:>8} µM   "
            f"h={c['hill']:g}")

    add(f"\n【装配层级】{d['coverage_tier_label']}")
    for t in d["compound_trace"]:
        add(f"  · {t['name']} → 靶点 {t['target'] or '—'} ({t['resolved_via']})")
        add(f"      占据率 θ={t['occupancy']:.3f}   kernel: "
            f"{'+'.join(t['kernels_used']) or '无'}   通路来源: {t['pathway_source']}")
        if t["dep"]:
            add(f"      Chronos  HepG2 {t['dep']['hepg2']:+.3f} / "
                f"RPE1 {t['dep']['rpe1']:+.3f}   →  {t['dep_class']}")
            add(f"      杀伤压力 HepG2 {t['pressure']['hepg2']:.3f} / "
                f"RPE1 {t['pressure']['rpe1']:.3f}")
        for w in t["warnings"]:
            add(f"      ⚠️  {w}")

    v = d["pseudo_viability_at_spec_dose"]
    if v:
        add("\n【规格剂量下的伪存活率】(⚠️ 绝对值无意义, 看比值)")
        for k, x in v.items():
            add(f"  {k:<8} {x:.3f}")

    s = d["selectivity"]
    if s:
        add("\n【选择性】")
        if s.get("iso_effect_level"):
            lv = s["iso_effect_level"]
            add(f"  同为 {lv:.0%} 抑制:  肝癌系 {_fmt(s['iso_dose_target'])} × 剂量  vs  "
                f"正常系 {_fmt(s['iso_dose_normal'])} × 剂量")
        add(f"  {s['interpretation']}")
        bl = s.get("iso_effect_ratio_by_level") or {}
        if bl:
            add("  各抑制水平上的等效剂量比: " +
                "  ".join(f"{k}→{_fmt(v)}×" if v else f"{k}→不可达" for k, v in bl.items()))
        add(f"  保守指标 (正常系 10% 抑制 ÷ 肝癌系 50% 抑制) = "
            f"{_fmt(s['conservative_window_x'])} ×")
        if s.get("why_undefined"):
            add(f"  ↳ {s['why_undefined']}")
        add(f"  饱和下限 (剂量→∞): HepG2 {s['saturation_floor']['hepg2']:.3f} / "
            f"RPE1 {s['saturation_floor']['rpe1']:.3f}")
        add(f"  ✅ 饱和杀伤压力比 HepG2÷RPE1 = {_fmt(s['saturation_pressure_ratio'])} × "
            f"(不依赖任何剂量假设)")

    c = d["combination"]
    if c.get("applicable"):
        add("\n【组合项】")
        if c["evidence"] == "measured":
            add(f"  ✅ 命中 ComboMap 实测对 {c['pair']} "
                f"(排名 {c['rank']}, {c['n_cells']} 细胞)")
            add(f"     上位性 {c['epistasis']:.3f} → 交叉项 γ={c['gamma']:.3f}")
            for ln in c["bliss_additive"]:
                add(f"     {ln:<8} 加和预期 {c['bliss_additive'][ln]:.3f}  →  "
                    f"上位性修正 {c['epistasis_corrected'][ln]:.3f}  "
                    f"(偏离 {c['deviation_x'][ln]}×)")
        else:
            p = c["empirical_prior"]
            add(f"  ⚠️ 该靶点对不在 ComboMap 的 {p['n_pairs']} 对里 → 回落到 Bliss 加和")
            add(f"     经验先验: {p['frac_new_signal_above_noise']:.0%} 的实测组合"
                f"新信号 >1× 噪声 (中位 {p['median_new_signal_x']:.2f}×)")
            for ln, b in c["empirical_band"].items():
                add(f"     {ln:<8} 加和 {c['bliss_additive'][ln]:.3f}   "
                    f"经验区间 [{b[0]:.3f}, {b[1]:.3f}]")
    elif c:
        add(f"\n【组合项】不适用 — {c.get('reason')}")

    p = d["pathway_response"]
    if p:
        add("\n【通路响应谱】")
        if not p["available"]:
            add("  ⛔ 不可用 — 无化合物落入通路层覆盖, 引擎拒绝外推")
        else:
            for srcs in p["sources"]:
                add(f"  来源: {srcs}")
            add(f"  最受抑制: {', '.join(p['top_down'])}")
            add(f"  最被动员: {', '.join(p['top_up'])}")

    cf = d["confidence"]
    add("\n【置信度】")
    add(f"  上下文 {cf['context_label']}  数字孪生就绪度 {cf['twin_readiness']} "
        f"(达 oracle 上限 {cf['frac_of_oracle']:.1%})")
    add(f"  必需性覆盖 {cf['essentiality_coverage']:.0%} · "
        f"通路覆盖 {cf['pathway_coverage']:.0%} · "
        f"IC50 占位 {'是' if cf['ic50_placeholder_used'] else '否'}")
    add(f"  综合分数 {cf['composite_score']:.3f} (⚠️ 排序用启发式, 非校准概率)")

    led = d["honesty_ledger"]
    add("\n【诚实账本】")
    for cl in led["claims"]:
        add(f"  {MARK[cl['evidence']]} {cl['label']}: {_fmt(cl['value'])}{cl['unit']}")
        if cl["basis"]:
            add(f"       依据: {cl['basis']}")
        if cl["caveat"]:
            add(f"       边界: {cl['caveat']}")
    for n in led["notes"]:
        add(f"  · {n}")
    add("\n  引擎明确拒绝输出:")
    for r in led["refused_outputs"]:
        add(f"    ⛔ {r}")
    add("=" * 74)
    return "\n".join(L)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="vcff", description="VC-FF 虚拟细胞力场引擎")
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("describe", help="打印引擎接口自述")
    sub.add_parser("scenarios", help="列出内置客户场景")
    r = sub.add_parser("run", help="运行一个内置场景")
    r.add_argument("name")
    r.add_argument("--json", action="store_true")
    s = sub.add_parser("stdin", help="从 stdin 读 AssaySpec JSON")
    s.add_argument("--json", action="store_true")
    a = ap.parse_args(argv)

    ff = ForceField()
    if a.cmd == "describe":
        print(json.dumps(ff.describe(), ensure_ascii=False, indent=2))
        return 0
    if a.cmd == "scenarios":
        for k, v in SCENARIOS.items():
            print(f"{k:<32} {v.customer} — {v.note}")
        return 0
    if a.cmd == "run":
        if a.name not in SCENARIOS:
            print(f"未知场景 {a.name}; 可选: {', '.join(SCENARIOS)}", file=sys.stderr)
            return 2
        card = ff.evaluate(SCENARIOS[a.name])
    else:
        card = ff.evaluate(AssaySpec.from_dict(json.load(sys.stdin)))

    print(card.to_json() if a.json else render(card))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
