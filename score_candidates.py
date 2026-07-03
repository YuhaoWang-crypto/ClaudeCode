#!/usr/bin/env python3
"""
score_candidates.py — 读取 candidates.yaml + proto_variants_il23r_tyk2.tsv，
计算每个候选的优先级分数、规划 PROTO/ESM-2 变异打分调用（骨架）、
并产出可审计的追踪表 (candidates_tracking.csv)。

设计原则
- 纯 stdlib + PyYAML；无网络即可 dry-run。
- PROTO / ESM-2 / Boltz 的真实接口是"骨架"：默认 dry-run 只打印计划；
  接入真实服务时实现 EngineClient 子类的 score() 即可。
- 优先级 = f(证据级别, 新颖性, 可及性, 矛盾证据惩罚, 可计算性)。

用法
  python3 score_candidates.py                 # 打分 + 打印排名 + 写 tracking csv
  python3 score_candidates.py --plan-variants # 额外打印每个变异的引擎调用计划
  python3 score_candidates.py --run-engines   # 尝试真实调用(未实现则 NotImplementedError)
"""
from __future__ import annotations
import argparse, csv, sys, os
from dataclasses import dataclass, field

try:
    import yaml
except ImportError:
    sys.exit("需要 PyYAML：pip install pyyaml")

HERE = os.path.dirname(os.path.abspath(__file__))
CANDIDATES = os.path.join(HERE, "candidates.yaml")
VARIANTS = os.path.join(HERE, "proto_variants_il23r_tyk2.tsv")
TRACKING_OUT = os.path.join(HERE, "candidates_tracking.csv")

# ---------- 优先级评分权重（可调） ----------
NOVELTY_BONUS = {           # 新颖性：差异化价值
    "novel_high_value": 2.0,
    "semi_novel_as_combination": 1.0,
    "known_anchor_for_scRNA_approach": 0.6,
    "known_mechanism_make_blood_assay": 0.5,
    "known_extend_to_multiclass_serum_panel": 0.5,
    "emerging_2025": 1.2,
    "known_but_validation_conflicted": 0.2,
}
LANDING_BONUS = {"high": 1.5, "medium": 0.8, "low_medium": 0.4, "low": 0.2}


@dataclass
class Candidate:
    cid: str
    symbol: list
    mechanism_class: str
    grade: int
    novelty: str
    landing: str
    n_supporting: int
    n_contradicting: int
    has_gap: bool
    engines: list
    paired_drug: list
    validation_cohort: list
    priority: float = 0.0
    status: str = "candidate"

    def score(self) -> float:
        # 证据级别是主轴；但 grade=1 且 high-novelty 仍可高优先（差异化机会）
        s = self.grade * 2.0
        s += NOVELTY_BONUS.get(self.novelty, 0.0)
        s += LANDING_BONUS.get(self.landing, 0.0)
        s += 0.4 * self.n_supporting
        s -= 0.9 * self.n_contradicting          # 矛盾证据惩罚（对抗性验证）
        s += 0.3 * len([e for e in self.engines if e])  # 可计算性
        # gap 候选：若是 novel_high_value，给探索性加成（值得压算力）
        if self.has_gap and self.novelty == "novel_high_value":
            s += 1.0
        self.priority = round(s, 2)
        return self.priority


def load_candidates(path: str) -> list[Candidate]:
    with open(path) as f:
        doc = yaml.safe_load(f)
    out = []
    for c in doc.get("candidates", []):
        ev = c.get("evidence", {}) or {}
        supporting = ev.get("supporting") or []
        contradicting = ev.get("contradicting") or []
        has_gap = "gap_note" in ev
        cn = c.get("compute_nextstep", {}) or {}
        engines = [k for k, v in cn.items() if v and v != "na"]
        assay = c.get("assayability", {}) or {}
        out.append(Candidate(
            cid=c["candidate_id"],
            symbol=c.get("symbol") if isinstance(c.get("symbol"), list) else [c.get("symbol")],
            mechanism_class=c.get("mechanism_class", "?"),
            grade=int(c.get("evidence_grade", 1)),
            novelty=c.get("novelty", ""),
            landing=str(assay.get("landing", "medium")),
            n_supporting=len(supporting),
            n_contradicting=len(contradicting),
            has_gap=has_gap,
            engines=engines,
            paired_drug=c.get("paired_drug") or [],
            validation_cohort=c.get("validation_cohort") or [],
        ))
    return out


def load_variants(path: str) -> list[dict]:
    if not os.path.exists(path):
        return []
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


# ---------- 引擎接口骨架（接真实服务时实现 score） ----------
class EngineClient:
    name = "engine"

    def score(self, variant: dict) -> dict:
        raise NotImplementedError(
            f"{self.name}.score() 未实现 —— 接入真实服务后在此调用 API 并返回打分。"
        )

    def plan(self, variant: dict) -> str:
        raise NotImplementedError


class Esm2VepClient(EngineClient):
    """编码错义变异有害性（伪对数似然 / 零样本）。"""
    name = "esm2_vep"

    def plan(self, v: dict) -> str:
        return (f"[ESM-2 VEP] {v['gene']} {v.get('hgvs_protein','?')} "
                f"→ 计算 wt vs mut 序列的 masked-LM 伪对数似然差 (Δpll)")


class ProtoClient(EngineClient):
    """PROTO：核酸/调控层。Evo2 + AlphaGenome (表达)、Enformer/Borzoi (染色质)、SpliceAI/Pangolin (剪接)。"""
    name = "proto"

    def plan(self, v: dict) -> str:
        tasks = v.get("proto_tasks", "")
        return (f"[PROTO] {v['gene']} {v['rsid']} ({v.get('region','?')}) "
                f"→ 运行: {tasks}")


class BoltzClient(EngineClient):
    name = "boltz"

    def plan(self, v: dict) -> str:
        return f"[Boltz] {v['gene']} {v.get('hgvs_protein','?')} → wt/mut 结构+口袋比较"


def route_engine(v: dict) -> EngineClient:
    """按变异类型选主引擎。"""
    consequence = (v.get("consequence") or "").lower()
    if "missense" in consequence:
        return Esm2VepClient()
    return ProtoClient()


def plan_variants(variants: list[dict], run: bool = False) -> None:
    print("\n=== 变异 → 引擎调用计划 ===")
    for v in variants:
        eng = route_engine(v)
        print("  " + eng.plan(v)
              + (f"   [rsid_status={v.get('rsid_status')}]" if v.get("rsid_status") else ""))
        # 编码变异额外走 SpliceAI（PROTO）和可选 Boltz
        if "missense" in (v.get("consequence") or "").lower():
            print("     ↳ 附加 " + ProtoClient().plan(v))
        if run:
            try:
                eng.score(v)
            except NotImplementedError as e:
                print(f"     ⚠ {e}")


def write_tracking(cands: list[Candidate], path: str) -> None:
    cols = ["priority", "candidate_id", "mechanism_class", "symbol", "grade",
            "novelty", "landing", "supporting", "contradicting", "engines",
            "paired_drug", "validation_cohort", "status"]
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(cols)
        for c in sorted(cands, key=lambda x: x.priority, reverse=True):
            w.writerow([c.priority, c.cid, c.mechanism_class, "|".join(map(str, c.symbol)),
                        c.grade, c.novelty, c.landing, c.n_supporting, c.n_contradicting,
                        ";".join(c.engines), "|".join(c.paired_drug),
                        "|".join(c.validation_cohort), c.status])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--plan-variants", action="store_true", help="打印变异→引擎计划")
    ap.add_argument("--run-engines", action="store_true", help="尝试真实调用引擎(未实现则报错)")
    args = ap.parse_args()

    cands = load_candidates(CANDIDATES)
    for c in cands:
        c.score()

    print("=== 候选优先级排名（priority 降序）===")
    print(f"{'PRIOR':>6}  {'ID':<10} {'CLASS':<12} {'GR':>2} {'SYMBOL'}")
    for c in sorted(cands, key=lambda x: x.priority, reverse=True):
        flag = "  ⚠矛盾" if c.n_contradicting else ""
        gap = "  ★缺口机会" if (c.has_gap and c.novelty == "novel_high_value") else ""
        print(f"{c.priority:>6}  {c.cid:<10} {c.mechanism_class:<12} {c.grade:>2} "
              f"{'/'.join(map(str,c.symbol))}{flag}{gap}")

    write_tracking(cands, TRACKING_OUT)
    print(f"\n✓ 追踪表已写入 {os.path.relpath(TRACKING_OUT, HERE)}")

    variants = load_variants(VARIANTS)
    print(f"✓ 载入 {len(variants)} 个待打分变异（{VARIANTS.split('/')[-1]}）")
    if args.plan_variants or args.run_engines:
        plan_variants(variants, run=args.run_engines)
    else:
        print("  （加 --plan-variants 查看每个变异的引擎调用计划）")


if __name__ == "__main__":
    main()
