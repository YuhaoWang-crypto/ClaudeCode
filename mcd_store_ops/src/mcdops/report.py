"""店长每日报表包：报表清单 (report list) + 生成日报 Markdown。

DAILY_REPORTS 是这个 package 的骨架——它就是"一个店长每天要看的报表清单"，
每张表都写清了：什么时候看、看几分钟、回答什么问题、看哪些指标、数据从哪来、
看完要做什么动作。没有最后一栏的报表，都是浪费时间的报表。
"""

from __future__ import annotations

from datetime import date, timedelta

from . import matrix as mx
from .metrics import MetricEngine

STATUS_ICON = {"green": "🟢", "amber": "🟡", "red": "🔴", "na": "⚪"}
SPARK = "▁▂▃▄▅▆▇█"


# ══ 每日报表清单 ══════════════════════════════════════════════════════════
DAILY_REPORTS = [
    {
        "code": "R01", "name": "昨日日结一页纸 (DSR)", "when": "07:00 开店前", "minutes": 5,
        "owner": "RGM",
        "question": "昨天到底做成什么样？有没有需要今天补的窟窿？",
        "metrics": ["net_sales", "sales_vs_forecast", "guest_count", "avg_check",
                    "food_cost_pct", "labor_cost_pct", "cash_over_short"],
        "tables": ["fct_daily_sales_summary", "agg_food_cost_day", "agg_labor_day",
                   "fct_cash_declaration"],
        "action": "只有一个动作：把偏差最大的那一项，写进今天班前会的第一句话。",
    },
    {
        "code": "R02", "name": "今日销售预测 × 排班对照表", "when": "07:10", "minutes": 8,
        "owner": "RGM",
        "question": "今天的人，够不够撑住今天的生意？缺口在哪个 15 分钟？",
        "metrics": ["sales_vs_forecast", "labor_schedule_variance", "splh", "attendance_rate"],
        "tables": ["dim_sales_forecast", "fct_schedule", "fct_timecard"],
        "action": "高峰前 30 分钟人力缺口 > 1 人 → 现在就打电话调班，不要等到高峰再救。",
    },
    {
        "code": "R03", "name": "开店食安与设备检查表", "when": "07:20 开店前", "minutes": 10,
        "owner": "SM",
        "question": "今天能不能安全地开门？",
        "metrics": ["food_safety_score", "critical_temp_fail", "equipment_downtime"],
        "tables": ["fct_food_safety_check", "fct_temperature_log", "fct_equipment_ticket",
                   "dim_equipment"],
        "action": "关键项任意一项 fail = 停止该产品出品并立即整改。这张表没有'先记着'这个选项。",
    },
    {
        "code": "R04", "name": "库存盘点与订货建议", "when": "07:40", "minutes": 8,
        "owner": "RGM",
        "question": "今天会不会断货？会不会又订多了明天扔掉？",
        "metrics": ["food_cost_pct", "waste_pct"],
        "tables": ["fct_inventory_count", "dim_sales_forecast", "fct_goods_receipt", "dim_item"],
        "action": "关键品缺口 → 调整备货计划；连续三天盘点差异同向 → 查配方/打料/收货。",
    },
    {
        "code": "R05", "name": "实时销售 vs 预测（日内滚动）", "when": "每小时 + 高峰前", "minutes": 2,
        "owner": "SM",
        "question": "今天走在预算线上吗？差在哪个渠道？",
        "metrics": ["net_sales", "sales_vs_forecast", "guest_count", "avg_check",
                    "items_per_check", "delivery_share"],
        "tables": ["fct_pos_transaction", "dim_sales_forecast"],
        "action": "落后 > 5% 且原因在客单价 → 立刻做附加销售动员；原因在客流 → 改动的是外场和外送曝光。",
    },
    {
        "code": "R06", "name": "服务速度日内看板", "when": "高峰进行中（11:30 / 18:00）", "minutes": 3,
        "owner": "SM",
        "question": "现在慢在哪一段？点单、厨房、还是交付？",
        "metrics": ["dt_total_time", "dt_window_time", "kds_prep_time", "front_otd",
                    "speed_hit_rate", "delivery_handover_time"],
        "tables": ["fct_dt_timer", "fct_kds_order", "fct_delivery_order"],
        "action": "窗口时长高 → 后厨备餐没跟上，加人到理货；点单段高 → 加点单员或开第二车道。",
    },
    {
        "code": "R07", "name": "人力投放与 SPLH 追踪", "when": "高峰后（14:00 / 20:30）", "minutes": 5,
        "owner": "RGM",
        "question": "这一波高峰，人是花值了还是花冤了？",
        "metrics": ["splh", "labor_cost_pct", "labor_schedule_variance", "attendance_rate", "late_rate"],
        "tables": ["fct_timecard", "fct_schedule", "fct_pos_transaction"],
        "action": "SPLH 低 + 速度也差 = 排班结构错（不是人不够，是人站错岗）→ 调岗位不是加人头。",
    },
    {
        "code": "R08", "name": "废弃与保温批次报表", "when": "14:00 / 21:00 交接班", "minutes": 5,
        "owner": "SM",
        "question": "今天扔了多少？为什么扔？",
        "metrics": ["waste_pct", "holding_time_violation", "food_cost_pct"],
        "tables": ["fct_waste_log", "fct_holding_batch", "fct_pos_line_item"],
        "action": "成品废弃高 = 备货节奏问题（改批量大小）；原料废弃高 = 先进先出与订货量问题。",
    },
    {
        "code": "R09", "name": "外送平台运营报表", "when": "15:00", "minutes": 5,
        "owner": "RGM",
        "question": "外送这四分之一的生意，是在赚钱还是在赔口碑？",
        "metrics": ["delivery_handover_time", "delivery_rating", "delivery_share", "order_accuracy"],
        "tables": ["fct_delivery_order", "fct_order_issue", "fct_guest_survey"],
        "action": "出餐时长超标 → 高峰给外送单独留一个打包位；差评集中在漏餐 → 上双人核对。",
    },
    {
        "code": "R10", "name": "顾客反馈与差错清单", "when": "16:00", "minutes": 5,
        "owner": "SM",
        "question": "顾客今天在抱怨什么？是不是同一件事在重复发生？",
        "metrics": ["csat", "complaints_per_1k", "order_accuracy", "cleanliness_score"],
        "tables": ["fct_guest_survey", "fct_order_issue", "fct_travel_path_check"],
        "action": "同一个 comment_tag 一周出现 3 次以上 → 它不是个案，是流程缺陷，进周复盘。",
    },
    {
        "code": "R11", "name": "闭店日结：销售 / 成本 / 现金", "when": "23:30 闭店", "minutes": 12,
        "owner": "RGM",
        "question": "今天的账对不对？成本控住了吗？",
        "metrics": ["net_sales", "food_cost_pct", "waste_pct", "labor_cost_pct", "splh",
                    "cash_over_short"],
        "tables": ["fct_pos_transaction", "fct_payment", "agg_food_cost_day", "agg_labor_day",
                   "fct_cash_declaration", "fct_daily_sales_summary"],
        "action": "现金短溢超阈值 → 当班复核录像与交接记录，当天闭环，不隔夜。",
    },
    {
        "code": "R12", "name": "次日准备与交接单", "when": "闭店前", "minutes": 5,
        "owner": "SM",
        "question": "明天早班一进门，需要知道哪三件事？",
        "metrics": ["equipment_downtime", "station_certification", "attendance_rate"],
        "tables": ["fct_equipment_ticket", "fct_schedule", "fct_training_record"],
        "action": "交接单只写三件事：坏了什么、缺了谁、明天几点最忙。写多了没人看。",
    },
]

WEEKLY_REPORTS = [
    {"code": "W01", "name": "周 P&L 与可控成本复盘", "owner": "RGM",
     "metrics": ["food_cost_pct", "labor_cost_pct", "waste_pct", "net_sales"],
     "question": "这周的钱花在哪，下周能省回来多少？"},
    {"code": "W02", "name": "周人员报表：排班效率 / 认证 / 流失", "owner": "RGM",
     "metrics": ["crew_turnover_30d", "station_certification", "attendance_rate", "splh"],
     "question": "队伍在变强还是在流血？"},
    {"code": "W03", "name": "周食安与设备健康", "owner": "SM",
     "metrics": ["food_safety_score", "critical_temp_fail", "equipment_downtime", "cleanliness_score"],
     "question": "有没有反复出现的同一个隐患？"},
    {"code": "W04", "name": "周顾客体验与活动复盘", "owner": "RGM",
     "metrics": ["csat", "delivery_rating", "complaints_per_1k", "avg_check", "upsell_rate"],
     "question": "顾客的评价，和我们内部的数说的是同一件事吗？"},
]


# ══ 格式化 ═══════════════════════════════════════════════════════════════
def fmt(v: float | None, unit: str) -> str:
    if v is None:
        return "—"
    if unit == "CNY":
        return f"{v:,.1f}" if abs(v) < 200 else f"{v:,.0f}"
    if unit == "%":
        return f"{v:.1f}%"
    if unit == "秒":
        return f"{v:.0f}s"
    if unit == "笔":
        return f"{v:,.0f}"
    if unit == "分":
        return f"{v:.1f}"
    if unit in ("次", "分钟"):
        return f"{v:.0f}"
    if unit == "CNY/工时":
        return f"{v:,.0f}"
    return f"{v:.2f}"


def fmt_gap(v: float | None) -> str:
    """差值：小数量级（百分点）多给一位，否则一位就够。"""
    if v is None:
        return "—"
    return f"{v:+.2f}" if abs(v) < 10 else f"{v:+.1f}"


def sparkline(values: list[float | None]) -> str:
    vals = [v for v in values if v is not None]
    if len(vals) < 2:
        return ""
    lo, hi = min(vals), max(vals)
    if hi == lo:
        return SPARK[3] * len(values)
    out = []
    for v in values:
        if v is None:
            out.append(" ")
        else:
            out.append(SPARK[min(7, int((v - lo) / (hi - lo) * 7.99))])
    return "".join(out)


def _heat(status: str) -> str:
    return {"green": "·", "amber": "▲", "red": "█", "na": " "}[status]


# ══ 日报生成 ═════════════════════════════════════════════════════════════
def build_daily_report(engine: MetricEngine, biz_date: date, trend_days: int = 14) -> str:
    p = engine.profile
    card = engine.scorecard(biz_date)
    by_id = {r["id"]: r for r in card}
    fam = mx.family_summary(engine, biz_date)
    pri = mx.priority(engine, biz_date)

    L: list[str] = []
    A = L.append

    A(f"# {p['store_name']} · 每日运营日报")
    A("")
    A(f"**营业日** {biz_date}（{'一二三四五六日'[biz_date.weekday()]}） · "
      f"**门店** {p['store_id']} · **店型** {p['format']}")
    A("")

    # ── 一句话结论
    reds = [r for r in card if r["status"] == "red"]
    ambers = [r for r in card if r["status"] == "amber"]
    ns = by_id["net_sales"]
    verdict = (f"净销售 {fmt(ns['value'], 'CNY')}（目标 {fmt(ns['target'], 'CNY')}，"
               f"达成 {fmt(by_id['sales_vs_forecast']['value'], '%')}）；"
               f"{len(reds)} 项红灯、{len(ambers)} 项黄灯。")
    A(f"> **今天一句话**：{verdict}")
    if pri["top3"]:
        A(f"> **今天先做这件事**：{pri['top3'][0]['name_cn']} —— {pri['top3'][0]['quadrant']}")
    A("")

    # ── 六大家族概览
    A("## 一、看板顶栏：六个家族的红黄绿")
    A("")
    A("| 家族 | 状态 | 绿/黄/红 | 最该看的一项 | 实际 | 目标 |")
    A("|---|:--:|:--:|---|--:|--:|")
    for f in fam:
        c = f["counts"]
        h = f["headline"]
        A(f"| {f['family_cn']} | {STATUS_ICON[f['status']]} | "
          f"{c['green']}/{c['amber']}/{c['red']} | {h['name_cn']} | "
          f"{fmt(h['value'], h['unit'])} | {fmt(h['target'], h['unit'])} |")
    A("")

    # ── 今日三件事（优先级矩阵 + 自动归因）
    A("## 二、今天的三件事（优先级矩阵 M4 + 归因矩阵 M1）")
    A("")
    if not pri["top3"]:
        A("今天全绿。把时间用在训练和顾客身上。")
    for i, it in enumerate(pri["top3"], 1):
        A(f"### {i}. {STATUS_ICON[it['status']]} {it['name_cn']}："
          f"{fmt(it['value'], it['unit'])}（目标 {fmt(it['target'], it['unit'])}，"
          f"差 {it['gap_pct']:+.1f}%）")
        A("")
        A(f"- **象限**：{it['quadrant']} · {it['quadrant_note']}")
        A(f"- **责任人**：{mx.OWNER_LABEL[it['owner']]} · **行动卡**：{it['playbook']} · "
          f"**行动分**：{it['score']}（超线 {it['severity']:.2f} 个带宽 × 影响 {it['impact']} "
          f"× 可控 {it['controllability']}"
          f"{' × 红灯加权 2.0' if it['status'] == 'red' else ''}）")
        cells = mx.worst_cells(engine, biz_date, it["id"], top=3)
        if cells:
            A(f"- **问题最集中在**：")
            for c in cells:
                A(f"  - {c['dim_cn']} = **{c['dim_value']}** → "
                  f"{fmt(c['value'], it['unit'])} {STATUS_ICON[c['status']]}"
                  f"（拖累度 {c['contribution']:.2f}，量占比 {c['weight']*100:.0f}%）")
        wins = mx.problem_windows(engine, biz_date, it["id"]) if engine.defs[it["id"]].get("grain") == "15min" else []
        if wins:
            w = wins[0]
            A(f"- **问题时段**：{w['from']} → {w['to']}（连续 {w['slots']} 个 15 分钟"
              f"{'，含红灯' if w['has_red'] else ''}，最差 {fmt(w['worst'], it['unit'])}）")
        A("")

    # ── 完整记分卡
    A("## 三、完整记分卡（32 项指标）")
    A("")
    A(f"| 家族 | 指标 | 实际 | 目标 | 达标 | 状态 | 责任人 | 近{trend_days}天 |")
    A("|---|---|--:|--:|--:|:--:|:--:|---|")
    for r in card:
        trend = engine.trend(r["id"], biz_date, trend_days)
        A(f"| {mx.FAMILY_LABEL[r['family']]} | {r['name_cn']} | "
          f"{fmt(r['value'], r['unit'])} | {fmt(r['target'], r['unit'])} | "
          f"{('%+.1f%%' % r['gap_pct']) if r['gap_pct'] is not None else '—'} | "
          f"{STATUS_ICON[r['status']]} | {mx.OWNER_LABEL[r['owner']]} | "
          f"`{sparkline([v for _d, v in trend])}` |")
    A("")

    # ── 时段矩阵
    A("## 四、时段矩阵 M2：一天里哪几个 15 分钟崩了")
    A("")
    A("`·` 绿灯　`▲` 黄灯　`█` 红灯")
    A("")
    grid_ids = ["sales_vs_forecast", "speed_hit_rate", "dt_total_time", "kds_prep_time", "splh"]
    g = mx.interval_grid(engine, biz_date, grid_ids)
    ivs = g["intervals"]
    hours = sorted({iv[:2] for iv in ivs})
    A("| 指标 | " + " | ".join(f"{h}" for h in hours) + " |")
    A("|---|" + "|".join([":--:"] * len(hours)) + "|")
    for mid in grid_ids:
        cells = {c["interval15"]: c for c in g["grid"][mid]}
        row = []
        for h in hours:
            row.append("".join(_heat(cells[iv]["status"]) for iv in ivs if iv[:2] == h))
        A(f"| {engine.defs[mid]['name_cn']} | " + " | ".join(row) + " |")
    A("")
    for mid in grid_ids:
        for w in mx.problem_windows(engine, biz_date, mid)[:1]:
            A(f"- **{engine.defs[mid]['name_cn']}** 连续告警：{w['from']} → {w['to']}"
              f"（{w['slots']} 个时段，均值 {fmt(w['avg'], engine.defs[mid]['unit'])}）")
    A("")

    # ── 交叉矩阵
    A("## 五、交叉矩阵 M3：时段 × 渠道，钉到具体格子")
    A("")
    for mid, da, db in (("speed_hit_rate", "daypart", "channel"),
                        ("net_sales", "daypart", "channel")):
        c = mx.cross(engine, biz_date, mid, da, db)
        A(f"**{c['metric_cn']}**（{c['dim_a_cn']} × {c['dim_b_cn']}）")
        A("")
        head = [ch for ch in c["b_values"]]
        A("| " + c["dim_a_cn"] + " | " + " | ".join(
            next(x["name_cn"] for x in p["channels"] if x["id"] == h) if db == "channel" else str(h)
            for h in head) + " |")
        A("|---|" + "|".join(["--:"] * len(head)) + "|")
        dp_cn = {d["id"]: d["name_cn"] for d in p["dayparts"]}
        for row in c["cells"]:
            label = dp_cn.get(row[0]["a"], row[0]["a"]) if da == "daypart" else row[0]["a"]
            A(f"| {label} | " + " | ".join(
                f"{fmt(cell['value'], c['unit'])}{STATUS_ICON[cell['status']]}" for cell in row) + " |")
        A("")
        if c["worst"]:
            w = c["worst"][0]
            a_cn = dp_cn.get(w["a"], w["a"])
            b_cn = next((x["name_cn"] for x in p["channels"] if x["id"] == w["b"]), w["b"])
            A(f"→ 最拖后腿的格子：**{a_cn} × {b_cn}** = {fmt(w['value'], c['unit'])}"
              f"（拖累度 {w['contribution']:.1f}）")
            A("")

    # ── 归因矩阵明细
    A("## 六、归因矩阵 M1：红黄灯指标的下钻明细")
    A("")
    for it in pri["items"][:5]:
        m = engine.defs[it["id"]]
        for dim in m.get("dims", [])[:2]:
            a = mx.attribution(engine, biz_date, it["id"], dim)
            rows = [r for r in a["rows"] if r["value"] is not None]
            if not rows:
                continue
            A(f"**{a['metric_cn']} × {a['dim_cn']}**"
              f"（全店 {fmt(a['store_value'], a['unit'])}，贡献度口径：{a['contribution_basis']}）")
            A("")
            A(f"| {a['dim_cn']} | 值 | 状态 | {a['gap_label']} | 拖累度 | 量占比 |")
            A("|---|--:|:--:|--:|--:|--:|")
            for r in rows:
                A(f"| {r['dim_value']} | {fmt(r['value'], a['unit'])} | "
                  f"{STATUS_ICON[r['status']]} | "
                  f"{fmt_gap(r['gap'])} | "
                  f"{r['contribution']:.2f} | {r['weight']*100:.0f}% |")
            A("")

    # ── 报表清单
    A("## 七、今天要看的报表清单（按时间轴）")
    A("")
    A("| 编号 | 报表 | 什么时候 | 几分钟 | 谁看 | 回答什么问题 | 看完做什么 |")
    A("|---|---|---|--:|:--:|---|---|")
    for r in DAILY_REPORTS:
        A(f"| {r['code']} | {r['name']} | {r['when']} | {r['minutes']} | "
          f"{mx.OWNER_LABEL[r['owner']]} | {r['question']} | {r['action']} |")
    A("")
    total = sum(r["minutes"] for r in DAILY_REPORTS)
    A(f"> 全部加起来 **{total} 分钟**。超过这个数，报表就开始和营运抢时间了。")
    A("")

    A("---")
    A("")
    A(f"*由 mcd_store_ops 生成 · 数据为模拟（scenario 见 config）· "
      f"口径以 `config/metrics.json` 为准*")
    return "\n".join(L)


def build_console_brief(engine: MetricEngine, biz_date: date) -> str:
    """班前会用的 20 行摘要，直接打在终端上。"""
    card = engine.scorecard(biz_date)
    by_id = {r["id"]: r for r in card}
    pri = mx.priority(engine, biz_date)
    fam = mx.family_summary(engine, biz_date)
    p = engine.profile

    L = [
        "═" * 68,
        f" {p['store_name']}  班前会简报   营业日 {biz_date}",
        "═" * 68,
        f" 净销售 {fmt(by_id['net_sales']['value'], 'CNY')}"
        f"   达成 {fmt(by_id['sales_vs_forecast']['value'], '%')}"
        f"   交易 {fmt(by_id['guest_count']['value'], '笔')}"
        f"   客单 {fmt(by_id['avg_check']['value'], 'CNY')}",
        f" 食品成本 {fmt(by_id['food_cost_pct']['value'], '%')}"
        f"   人工 {fmt(by_id['labor_cost_pct']['value'], '%')}"
        f"   SPLH {fmt(by_id['splh']['value'], 'CNY/工时')}"
        f"   废弃 {fmt(by_id['waste_pct']['value'], '%')}",
        "─" * 68,
        " " + "   ".join(f"{STATUS_ICON[f['status']]}{f['family_cn']}" for f in fam),
        "─" * 68,
        " 今天的三件事：",
    ]
    for i, it in enumerate(pri["top3"], 1):
        L.append(f"  {i}. {STATUS_ICON[it['status']]} {it['name_cn']} "
                 f"{fmt(it['value'], it['unit'])} (目标 {fmt(it['target'], it['unit'])}, "
                 f"{it['gap_pct']:+.1f}%)  →  {it['quadrant']} · {mx.OWNER_LABEL[it['owner']]} · {it['playbook']}")
        cells = mx.worst_cells(engine, biz_date, it["id"], top=2)
        for c in cells:
            L.append(f"       └ {c['dim_cn']} {c['dim_value']}: {fmt(c['value'], it['unit'])} "
                     f"(拖累 {c['contribution']:.2f})")
    if not pri["top3"]:
        L.append("  今天全绿。")
    L.append("═" * 68)
    return "\n".join(L)
