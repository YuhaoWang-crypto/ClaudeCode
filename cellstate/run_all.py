"""
Run the whole cell-state-space feasibility pipeline (C1 -> C7).

The question it answers: can a six-variable cell-state space -- cell cycle,
stress, metabolic reserve, DNA damage, apoptotic priming, epigenetic/lineage
state -- be BUILT out of existing open data?

Each module returns a measured number rather than an opinion, and later
modules consume the earlier ones' outputs, so the final budget in C7 is an
arithmetic consequence of C1-C6 rather than a separate claim.

    python3 -m cellstate.run_all          # everything
    python3 -m cellstate.c3_orthogonality # or any single module

Network is used for live index queries and data downloads; everything is
cached under figures/_cellstate_cache/ so re-runs are offline and identical.
Modules degrade explicitly (printing a warning and drawing no conclusion)
rather than substituting invented numbers when a download fails.
"""
import traceback

from . import (c1_inventory, c2_axis_anchors, c3_orthogonality,
               c4_observability, c5_pairing, c6_larry_ceiling, c7_budget)

STEPS = [
    ("C1  六轴开源数据盘点(GEO / BioModels / Europe PMC 实时检索)", c1_inventory),
    ("C2  每个轴锚定到真实策展模型(官方 SBML + libRoadRunner)", c2_axis_anchors),
    ("C3  正交性:真实人类表达数据 + 基因集重叠对照", c3_orthogonality),
    ("C4  可观测性:六个隐藏状态能否从组学恢复(闭式解)", c4_observability),
    ("C5  配对问题:替代设计的天花板 + 成像通道数", c5_pairing),
    ("C6  实测天花板:LARRY split-well 跨孔命运相关", c6_larry_ceiling),
    ("C7  综合预算:三个约束相乘 + 可执行建议", c7_budget),
]


def main():
    done, failed = [], []
    for title, mod in STEPS:
        print(f"\n\n{'#'*78}\n#  {title}\n{'#'*78}")
        try:
            mod.run()
            done.append(title)
        except Exception:
            failed.append(title)
            print(f"\n⚠️  该模块失败(不影响其它模块):")
            traceback.print_exc()

    print(f"\n\n{'='*78}\n完成 {len(done)}/{len(STEPS)} 个模块")
    for t in done:
        print(f"  ✅ {t}")
    for t in failed:
        print(f"  ⚠️  {t}")
    print("\n图见 figures/c*.png,数值结果见 figures/_cellstate_cache/*.json,")
    print("完整解读见 REPORT_CELLSTATE.md。")


if __name__ == "__main__":
    main()
