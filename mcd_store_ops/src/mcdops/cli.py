"""命令行入口。

    python -m mcdops.cli brief                 班前会简报（终端）
    python -m mcdops.cli report                生成 output/daily_report.md
    python -m mcdops.cli dashboard             生成 output/dashboard.html
    python -m mcdops.cli matrix <指标> <维度>   打印一张归因矩阵
    python -m mcdops.cli cross <指标> <维A> <维B>  打印一张交叉矩阵
    python -m mcdops.cli reports                打印每日报表清单
    python -m mcdops.cli dump                  导出模拟数据为 CSV
    python -m mcdops.cli validate              校验指标字典与数据字典是否自洽
    python -m mcdops.cli all                   跑全套并写 output/

公共参数：--date YYYY-MM-DD --days N --seed N --scenario clean|realistic
"""

from __future__ import annotations

import argparse
import csv
import sys
from datetime import date, datetime
from pathlib import Path

from . import dashboard, matrix as mx, report, schema
from .metrics import CALCULATORS, MetricEngine, load_metrics, load_profile
from .report import fmt
from .simulate import simulate_period

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "output"


def _engine(args):
    profile = load_profile()
    end = datetime.strptime(args.date, "%Y-%m-%d").date() if args.date else date(2026, 8, 18)
    data = simulate_period(profile, end, days=args.days, seed=args.seed, scenario=args.scenario)
    return MetricEngine(data, profile, load_metrics()), end


# ── 子命令 ───────────────────────────────────────────────────────────────
def cmd_brief(args):
    eng, d = _engine(args)
    print(report.build_console_brief(eng, d))


def cmd_report(args):
    eng, d = _engine(args)
    OUT.mkdir(exist_ok=True)
    path = OUT / "daily_report.md"
    path.write_text(report.build_daily_report(eng, d, args.days), encoding="utf-8")
    print(f"→ {path}")


def cmd_dashboard(args):
    eng, d = _engine(args)
    OUT.mkdir(exist_ok=True)
    path = OUT / "dashboard.html"
    path.write_text(dashboard.build_dashboard(eng, d, args.days), encoding="utf-8")
    print(f"→ {path}")


def cmd_matrix(args):
    eng, d = _engine(args)
    a = mx.attribution(eng, d, args.metric, args.dim)
    print(f"\n{a['metric_cn']} × {a['dim_cn']}   ({d})")
    print(f"全店 {fmt(a['store_value'], a['unit'])}  目标 {fmt(a['target'], a['unit'])}"
          f"  贡献度口径：{a['contribution_basis']}\n")
    print(f"{a['dim_cn']:<16s}{'值':>12s}{'状态':>6s}{a['gap_label']:>10s}{'拖累度':>10s}{'量占比':>8s}")
    print("─" * 64)
    for r in a["rows"]:
        if r["value"] is None:
            continue
        print(f"{str(r['dim_value']):<16s}{fmt(r['value'], a['unit']):>12s}"
              f"{report.STATUS_ICON[r['status']]:>5s}"
              f"{report.fmt_gap(r['gap']):>10s}"
              f"{r['contribution']:>10.2f}{r['weight']*100:>7.0f}%")
    print()


def cmd_cross(args):
    eng, d = _engine(args)
    c = mx.cross(eng, d, args.metric, args.dim_a, args.dim_b)
    print(f"\n{c['metric_cn']}   {c['dim_a_cn']} × {c['dim_b_cn']}   ({d})\n")
    if c["note"]:
        print(f"  ⚠ {c['note']}\n")
        return 0
    print(f"{'':<14s}" + "".join(f"{str(b):>14s}" for b in c["b_values"]))
    print("─" * (14 + 14 * len(c["b_values"])))
    for row in c["cells"]:
        print(f"{str(row[0]['a']):<14s}" + "".join(
            f"{fmt(cell['value'], c['unit']) + report.STATUS_ICON[cell['status']]:>14s}"
            for cell in row))
    if c["worst"]:
        w = c["worst"][0]
        print(f"\n最拖后腿：{w['a']} × {w['b']} = {fmt(w['value'], c['unit'])}"
              f"（拖累度 {w['contribution']:.2f}）\n")


def cmd_reports(_args):
    defs = load_metrics()
    print("\n【每日报表清单】按一天的时间轴\n")
    total = 0
    for r in report.DAILY_REPORTS:
        total += r["minutes"]
        print(f"{r['code']}  {r['name']}")
        print(f"      时间：{r['when']}（{r['minutes']} 分钟） · 谁看：{mx.OWNER_LABEL[r['owner']]}")
        print(f"      问题：{r['question']}")
        print(f"      指标：{'、'.join(defs[m]['name_cn'] for m in r['metrics'])}")
        print(f"      来源：{'、'.join(r['tables'])}")
        print(f"      动作：{r['action']}\n")
    print(f"合计 {total} 分钟/天\n")
    print("【每周报表】\n")
    for r in report.WEEKLY_REPORTS:
        print(f"{r['code']}  {r['name']}（{mx.OWNER_LABEL[r['owner']]}）— {r['question']}")
        print(f"      指标：{'、'.join(defs[m]['name_cn'] for m in r['metrics'])}\n")


def cmd_dump(args):
    eng, _d = _engine(args)
    out = OUT / "data"
    out.mkdir(parents=True, exist_ok=True)
    for table, rows in sorted(eng.data.items()):
        if not rows:
            continue
        cols = list(rows[0])
        with (out / f"{table}.csv").open("w", newline="", encoding="utf-8-sig") as fh:
            w = csv.DictWriter(fh, fieldnames=cols)
            w.writeheader()
            w.writerows(rows)
        print(f"→ {out / (table + '.csv')}  ({len(rows)} 行)")


def cmd_validate(args):
    """指标字典 ↔ 数据字典 ↔ 计算实现，三者必须对得上。"""
    defs = load_metrics()
    known = schema.table_names()
    # 模拟器造的中间聚合表：数据字典里以 agg_metric_* 描述其角色，这里单独放行
    generated = {"agg_labor_interval", "agg_labor_day", "agg_food_cost_day"}
    errs, warns = [], []

    for mid, m in defs.items():
        if mid not in CALCULATORS:
            errs.append(f"指标 {mid} 在 metrics.json 里有定义，但 metrics.py 没有实现公式")
        for t in m.get("source_tables", []):
            if t not in known and t not in generated:
                errs.append(f"指标 {mid} 引用了数据字典里不存在的表 {t}")
        if m["direction"] == "higher_better" and not (m["amber"] <= m["green"]):
            errs.append(f"指标 {mid} 阈值方向错：higher_better 应满足 amber <= green")
        if m["direction"] in ("lower_better", "neutral") and not (m["green"] <= m["amber"]):
            errs.append(f"指标 {mid} 阈值方向错：lower_better 应满足 green <= amber")
        if m.get("aggregation") not in ("sum", "ratio"):
            errs.append(f"指标 {mid} 缺少 aggregation（sum / ratio）")

    for mid in CALCULATORS:
        if mid not in defs:
            warns.append(f"metrics.py 实现了 {mid}，但 metrics.json 里没有它 —— 它不会出现在任何看板上")

    # 每张报表引用的指标必须存在
    for r in report.DAILY_REPORTS + report.WEEKLY_REPORTS:
        for m in r["metrics"]:
            if m not in defs:
                errs.append(f"报表 {r['code']} 引用了不存在的指标 {m}")

    # 每个指标至少要被一张报表看到，否则就是"没人看的指标"
    watched = {m for r in report.DAILY_REPORTS + report.WEEKLY_REPORTS for m in r["metrics"]}
    for mid in defs:
        if mid not in watched:
            warns.append(f"指标 {defs[mid]['name_cn']} ({mid}) 不在任何一张报表里 —— 没人会看到它")

    # 数据能不能真的把每个指标算出来
    eng, d = _engine(args)
    for mid in defs:
        try:
            v = eng.compute(d, [mid]).get(mid)
        except Exception as exc:                      # noqa: BLE001
            errs.append(f"指标 {mid} 计算抛异常：{exc!r}")
            continue
        if v is None:
            warns.append(f"指标 {defs[mid]['name_cn']} ({mid}) 在 {d} 算不出值（数据不足）")

    print(f"\n指标 {len(defs)} 项 · 数据表 {len(known)} 张 · 报表 "
          f"{len(report.DAILY_REPORTS)} 日 + {len(report.WEEKLY_REPORTS)} 周\n")
    for w in warns:
        print(f"  ⚠ {w}")
    for e in errs:
        print(f"  ✗ {e}")
    if not errs and not warns:
        print("  ✓ 全部自洽")
    elif not errs:
        print(f"\n  ✓ 无错误（{len(warns)} 条提醒）")
    print()
    return 1 if errs else 0


def cmd_docs(_args):
    """把数据字典和报表清单从代码里导出成文档，保证文档不会和实现走散。"""
    docs = ROOT / "docs"
    docs.mkdir(exist_ok=True)
    defs = load_metrics()

    # ── 01 数据字典
    L = ["# 01 · 一家门店每天在用的表（数据字典）", "",
         "> 本文件由 `python -m mcdops.cli docs` 从 `src/mcdops/schema.py` 生成，不要手改。",
         "",
         f"共 {len(schema.TABLES)} 张表，分 {len(schema.by_domain())} 个业务域。",
         "门店的现实是：这些表分散在 6~8 个互不相通的系统里，",
         "把它们对齐到同一个 `store_id` + `biz_date`，是做任何报表的前提工作，",
         "也是最容易出错的地方（尤其是跨零点的夜宵订单该算哪一天）。", ""]
    for domain, tables in schema.by_domain().items():
        if not tables:
            continue
        L += [f"## {domain}", ""]
        for t in tables:
            L += [f"### `{t['name']}` — {t['cn']}", "",
                  f"- **层级** {t['layer']} · **来源系统** {t['owner_system']} · "
                  f"**刷新** {t['refresh']}",
                  f"- **粒度** {t['grain']} · **主键** {', '.join(t['pk'])} · "
                  f"**日增行数** {t['rows_per_day']}",
                  f"- **为什么需要它**：{t['why']}", "",
                  "| 字段 | 类型 | 说明 |", "|---|---|---|"]
            for name, typ, note in t["fields"]:
                L.append(f"| `{name}` | {typ} | {note} |")
            L.append("")
    (docs / "01_source_tables.md").write_text("\n".join(L), encoding="utf-8")
    print(f"→ {docs / '01_source_tables.md'}")

    # ── 02 报表清单
    L = ["# 02 · 店长每日报表清单 (Daily Report List)", "",
         "> 本文件由 `python -m mcdops.cli docs` 从 `src/mcdops/report.py` 生成，不要手改。", "",
         "一张报表如果回答不了「看完要做什么」，它就不该出现在这个清单里。",
         "所以每一行的最后一栏是动作，不是描述。", ""]
    total = sum(r["minutes"] for r in report.DAILY_REPORTS)
    L += [f"全天 {len(report.DAILY_REPORTS)} 张日报表，合计 **{total} 分钟**。",
          "超过这个数，报表就开始和营运抢时间了。", "",
          "## 一天的时间轴", "",
          "| 编号 | 报表 | 时间 | 分钟 | 谁看 | 回答什么问题 |",
          "|---|---|---|--:|:--:|---|"]
    for r in report.DAILY_REPORTS:
        L.append(f"| {r['code']} | {r['name']} | {r['when']} | {r['minutes']} | "
                 f"{mx.OWNER_LABEL[r['owner']]} | {r['question']} |")
    L += ["", "## 每张报表的明细", ""]
    for r in report.DAILY_REPORTS:
        L += [f"### {r['code']} · {r['name']}", "",
              f"- **什么时候看**：{r['when']}，约 {r['minutes']} 分钟",
              f"- **谁看**：{mx.OWNER_LABEL[r['owner']]}",
              f"- **回答什么问题**：{r['question']}",
              f"- **看哪些指标**：" + "、".join(
                  f"{defs[m]['name_cn']}（目标 {defs[m]['target']}{defs[m]['unit']}）"
                  for m in r["metrics"]),
              f"- **数据来自**：" + "、".join(f"`{t}`" for t in r["tables"]),
              f"- **看完做什么**：{r['action']}", ""]
    L += ["## 每周报表", "",
          "| 编号 | 报表 | 谁看 | 回答什么问题 | 指标 |", "|---|---|:--:|---|---|"]
    for r in report.WEEKLY_REPORTS:
        L.append(f"| {r['code']} | {r['name']} | {mx.OWNER_LABEL[r['owner']]} | {r['question']} | "
                 + "、".join(defs[m]["name_cn"] for m in r["metrics"]) + " |")
    L.append("")
    (docs / "02_daily_report_list.md").write_text("\n".join(L), encoding="utf-8")
    print(f"→ {docs / '02_daily_report_list.md'}")


def cmd_all(args):
    cmd_report(args)
    cmd_dashboard(args)
    eng, d = _engine(args)
    print()
    print(report.build_console_brief(eng, d))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="mcdops", description="门店日运营报表 / 指标矩阵 / 看板")
    ap.add_argument("--date", default="2026-08-18", help="营业日 YYYY-MM-DD")
    ap.add_argument("--days", type=int, default=14, help="模拟多少天（趋势用）")
    ap.add_argument("--seed", type=int, default=20260819)
    ap.add_argument("--scenario", default="realistic", choices=["realistic", "clean"])

    # 同样的参数也允许写在子命令后面（`brief --date ...`），SUPPRESS 保证不覆盖顶层已给的值
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--date", default=argparse.SUPPRESS)
    common.add_argument("--days", type=int, default=argparse.SUPPRESS)
    common.add_argument("--seed", type=int, default=argparse.SUPPRESS)
    common.add_argument("--scenario", default=argparse.SUPPRESS, choices=["realistic", "clean"])

    sub = ap.add_subparsers(dest="cmd", required=True)

    for name, fn, helptext in (
        ("brief", cmd_brief, "班前会简报（终端）"),
        ("report", cmd_report, "生成 Markdown 日报"),
        ("dashboard", cmd_dashboard, "生成 HTML 看板"),
        ("reports", cmd_reports, "打印每日/每周报表清单"),
        ("dump", cmd_dump, "导出模拟数据为 CSV"),
        ("validate", cmd_validate, "校验指标/数据/报表三者自洽"),
        ("docs", cmd_docs, "从代码生成 docs/01 与 docs/02"),
        ("all", cmd_all, "跑全套"),
    ):
        sp = sub.add_parser(name, help=helptext, parents=[common])
        sp.set_defaults(func=fn)

    sp = sub.add_parser("matrix", help="打印归因矩阵：指标 × 维度", parents=[common])
    sp.add_argument("metric")
    sp.add_argument("dim")
    sp.set_defaults(func=cmd_matrix)

    sp = sub.add_parser("cross", help="打印交叉矩阵：指标 × 维A × 维B", parents=[common])
    sp.add_argument("metric")
    sp.add_argument("dim_a")
    sp.add_argument("dim_b")
    sp.set_defaults(func=cmd_cross)

    args = ap.parse_args(argv)
    return args.func(args) or 0


if __name__ == "__main__":
    sys.exit(main())
