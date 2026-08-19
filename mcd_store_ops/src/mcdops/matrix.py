"""矩阵引擎：把一列一列的报表数字，变成能直接指向动作的二维表。

报表回答"多少"，矩阵回答"哪里"和"先做哪个"。本模块提供五种矩阵，
它们是店长真正会用的五个动作，不是五种画法：

  M1 归因矩阵  attribution()   一个指标 × 一个维度  → 问题出在哪个渠道/时段/品类
  M2 时段矩阵  interval_grid() 多个指标 × 15分钟    → 一天里哪几段崩了（热力图）
  M3 交叉矩阵  cross()         维度A × 维度B        → 钉到"晚市的得来速"这种具体格子
  M4 优先级矩阵 priority()     影响度 × 可控度      → 今天这 3 件事先做哪件
  M5 责任矩阵  ownership()     指标 × 责任人 × 频率 → 谁在什么节奏上盯什么（静态配置）

贡献度 (contribution) 是 M1 的核心，也是最容易做错的地方：
不能直接按指标值排序（外送废弃率再高，只占 3% 的量也不值得先动），
要按"这一格拖了全店多少"排序。三种算法见 _contribution()。
"""

from __future__ import annotations

from datetime import date

from .metrics import (MetricEngine, Slice, gap_pct, gap_to_target, severity, status)

# 已知维度的展示顺序（其余按字母序）
DIM_ORDER_KEYS = {
    "daypart": lambda p: [d["id"] for d in p["dayparts"]],
    "channel": lambda p: [c["id"] for c in p["channels"]],
    "category": lambda p: [c["id"] for c in p["categories"]],
    "station": lambda p: p["stations"],
    "platform": lambda p: p["delivery_platforms"],
}

DIM_LABEL = {
    "daypart": "时段", "channel": "渠道", "category": "品类", "station": "岗位",
    "platform": "外卖平台", "waste_type": "废弃类型", "zone": "区域", "shift": "班次",
    "interval15": "15分钟", "equipment_id": "设备", "cashier_id": "收银员",
    "issue_type": "差错类型", "tenure_band": "司龄段",
}

FAMILY_LABEL = {
    "SALES": "生意", "SPEED": "速度", "COST": "成本",
    "QSC": "品质与食安", "PEOPLE": "人员", "GUEST": "顾客",
}

STATUS_RANK = {"red": 0, "amber": 1, "green": 2, "na": 3}


def dim_values(engine: MetricEngine, biz_date, dim_name: str) -> list:
    """某维度在当天出现过的所有取值，按业务顺序排。"""
    if dim_name in DIM_ORDER_KEYS:
        return DIM_ORDER_KEYS[dim_name](engine.profile)
    seen = set()
    bd = str(biz_date)
    for rows in engine.data.values():
        if not rows or dim_name not in rows[0]:
            continue
        for r in rows:
            if r.get("biz_date", bd) == bd and r.get(dim_name) is not None:
                seen.add(r[dim_name])
    return sorted(seen)


# 订单头没有的维度（外卖平台、岗位、设备、区域…），退化到各自的事实表按行数取权重
_WEIGHT_FALLBACK = [
    "fct_delivery_order", "fct_kds_order", "fct_dt_timer", "fct_waste_log",
    "fct_holding_batch", "fct_timecard", "fct_schedule", "fct_guest_survey",
    "fct_order_issue", "fct_temperature_log", "fct_travel_path_check",
    "fct_food_safety_check", "fct_cash_declaration", "dim_employee",
]


def _slice_weight(engine: MetricEngine, biz_date, dim_name, value, extra=None) -> float:
    """该切片在全店里占多大的量。率型偏差乘上它，才是对全店的真实影响。"""
    s = Slice(engine.data, biz_date, engine.profile, dim_name, value, extra)
    part = sum(t["net_amount"] for t in s.rows("fct_pos_transaction"))
    whole = sum(t["net_amount"] for t in s.day_rows("fct_pos_transaction"))
    if part and whole:
        return part / whole
    if dim_name == "category":
        lines = s.rows("fct_pos_line_item")
        allv = sum(l["line_net"] for l in s.day_rows("fct_pos_line_item"))
        if allv:
            return sum(l["line_net"] for l in lines) / allv
    for table in _WEIGHT_FALLBACK:
        day = s.day_rows(table)
        if day and dim_name in day[0]:
            return len(s.rows(table)) / len(day)
    return 0.0


def _forecast_pair(engine: MetricEngine, biz_date, metric_id, dim_name, value,
                   extra=None) -> tuple[float, float] | None:
    """量型指标在该切片上的 (预测, 实际)。"""
    field = {"net_sales": "forecast_sales", "guest_count": "forecast_gc"}.get(metric_id)
    if not field:
        return None
    s = Slice(engine.data, biz_date, engine.profile, dim_name, value, extra)
    fc = sum(f[field] for f in s.rows("dim_sales_forecast"))
    if not fc:
        return None
    actual = (sum(t["net_amount"] for t in s.rows("fct_pos_transaction"))
              if metric_id == "net_sales" else len(s.rows("fct_pos_transaction")))
    return fc, actual


def _forecast_gap(engine: MetricEngine, biz_date, metric_id, dim_name, value, extra=None) -> float | None:
    """量型指标的缺口 = 预测 − 实际（正数 = 少做了这么多）。"""
    pair = _forecast_pair(engine, biz_date, metric_id, dim_name, value, extra)
    return None if pair is None else pair[0] - pair[1]


# 量型指标在切片上的达成率红黄线（与 sales_vs_forecast 保持一致）
_ACHIEVEMENT_BANDS = (98.0, 93.0)


def scaled_metric(metric: dict, weight: float) -> dict:
    """按该切片占的量，把全店阈值缩成"这一格应得的那份预算"。

    全店废弃率目标 0.8%，午市占 34% 的销售 → 午市的废弃预算是 0.27pp。
    午市实际 0.54pp，是它应得份额的两倍——该判红，而不是因为 0.54 < 0.8 判绿。
    各格的预算加起来正好等于全店目标，所以这样判色仍然是自洽的。
    """
    return {**metric, "target": metric["target"] * weight,
            "green": metric["green"] * weight, "amber": metric["amber"] * weight}


def cell_status(engine, biz_date, metric, val, dim_name, value, extra=None) -> str:
    """矩阵格子的红黄绿。切片上不能直接套全店阈值，两类指标各有各的换算。

    - 量型（净销售、交易笔数）：晚市 4,761 元对着全店目标 62,000 永远是红的，
      整张矩阵一片红等于没有信息 → 改判"对该切片预测的达成率"。
    - 全店口径的率型（废弃率、千单投诉数）：切片值天生小于全店目标，一片绿
      同样没有信息 → 阈值按量占比缩放，见 scaled_metric()。
    - 其余率型（速度、人工率、客单价）：目标本来就与切片大小无关，照常判。
    """
    if val is None:
        return "na"
    if metric.get("aggregation") == "sum":
        pair = _forecast_pair(engine, biz_date, metric["id"], dim_name, value, extra)
        if pair:
            fc, actual = pair
            pct = actual / fc * 100
            green, amber = _ACHIEVEMENT_BANDS
            return "green" if pct >= green else ("amber" if pct >= amber else "red")
    if metric.get("denominator") == "store_day" and dim_name is not None:
        w = _slice_weight(engine, biz_date, dim_name, value, extra)
        if w > 0:
            return status(scaled_metric(metric, w), val)
    return status(metric, val)


def _contribution(engine, biz_date, metric, val, dim_name, value, extra=None) -> tuple[float, str]:
    """这一格拖了全店多少。返回 (贡献值, 口径说明)。正数 = 拖后腿。"""
    mid = metric["id"]
    if metric.get("aggregation") == "sum":
        g = _forecast_gap(engine, biz_date, mid, dim_name, value, extra)
        if g is not None:
            return g, "对预测的缺口"
        return (val or 0.0), "绝对量"
    if metric.get("denominator") == "store_day":
        # 已经是"占全店"的口径，本身就可加
        return (val or 0.0), "占全店百分点"
    g = gap_to_target(metric, val)
    if g is None:
        return 0.0, "无数据"
    w = _slice_weight(engine, biz_date, dim_name, value, extra)
    return -g * w, "偏差 × 销售权重"


# ── M1 归因矩阵 ───────────────────────────────────────────────────────────
def attribution(engine: MetricEngine, biz_date, metric_id: str, dim_name: str) -> dict:
    m = engine.defs[metric_id]
    is_volume = m.get("aggregation") == "sum"
    rows = []
    basis = ""
    for v in dim_values(engine, biz_date, dim_name):
        val = engine.compute(biz_date, [metric_id], dim_name, v).get(metric_id)
        contrib, basis = _contribution(engine, biz_date, m, val, dim_name, v)
        w = _slice_weight(engine, biz_date, dim_name, v)
        pair = _forecast_pair(engine, biz_date, metric_id, dim_name, v) if is_volume else None
        if pair:
            # 量型指标在切片上的"差"是对该切片预测的差，不是对全店目标的差
            gap, gp = pair[1] - pair[0], (pair[1] - pair[0]) / pair[0] * 100
        elif m.get("denominator") == "store_day" and w > 0:
            # 全店口径的率型：和"这一格应得的那份预算"比，不和全店目标比
            sm = scaled_metric(m, w)
            gap, gp = gap_to_target(sm, val), gap_pct(sm, val)
        else:
            gap, gp = gap_to_target(m, val), gap_pct(m, val)
        rows.append({
            "dim_value": v, "value": val,
            "status": cell_status(engine, biz_date, m, val, dim_name, v),
            "gap": gap, "gap_pct": gp,
            "contribution": contrib, "weight": w,
        })
    rows.sort(key=lambda r: -(r["contribution"] or 0))
    return {
        "kind": "attribution", "metric": metric_id, "metric_cn": m["name_cn"],
        "unit": m["unit"], "dim": dim_name, "dim_cn": DIM_LABEL.get(dim_name, dim_name),
        "contribution_basis": basis if rows else "",
        "gap_label": "对预测" if basis == "对预测的缺口" else "对目标",
        "store_value": engine.compute(biz_date, [metric_id]).get(metric_id),
        "target": m["target"], "rows": rows,
    }


def worst_cells(engine: MetricEngine, biz_date, metric_id: str, top: int = 3) -> list[dict]:
    """跨该指标所有可下钻维度，找出最拖后腿的几格。"""
    m = engine.defs[metric_id]
    cells = []
    for dim in m.get("dims", []):
        a = attribution(engine, biz_date, metric_id, dim)
        for r in a["rows"]:
            if r["value"] is None:
                continue
            cells.append({**r, "dim": dim, "dim_cn": a["dim_cn"]})
    cells.sort(key=lambda r: -(r["contribution"] or 0))
    return cells[:top]


# ── M2 时段矩阵（热力图） ─────────────────────────────────────────────────
def interval_grid(engine: MetricEngine, biz_date, metric_ids: list[str] | None = None) -> dict:
    ids = metric_ids or [mid for mid, m in engine.defs.items() if m.get("grain") == "15min"]
    intervals = dim_values(engine, biz_date, "interval15")
    intervals = sorted(intervals)
    grid = {}
    for mid in ids:
        m = engine.defs[mid]
        row = []
        for iv in intervals:
            val = engine.compute(biz_date, [mid], "interval15", iv).get(mid)
            row.append({"interval15": iv, "value": val,
                        "status": cell_status(engine, biz_date, m, val, "interval15", iv)})
        grid[mid] = row
    return {"kind": "interval_grid", "intervals": intervals, "metrics": ids, "grid": grid}


def problem_windows(engine: MetricEngine, biz_date, metric_id: str, min_run: int = 2) -> list[dict]:
    """连续 min_run 个 15 分钟都不是绿灯的时间段——单点飘红是噪声，连片才是问题。"""
    m = engine.defs[metric_id]
    cells = interval_grid(engine, biz_date, [metric_id])["grid"][metric_id]
    runs, cur = [], []
    for c in cells:
        if c["status"] in ("amber", "red"):
            cur.append(c)
        else:
            if len(cur) >= min_run:
                runs.append(cur)
            cur = []
    if len(cur) >= min_run:
        runs.append(cur)
    out = []
    for r in runs:
        vals = [c["value"] for c in r if c["value"] is not None]
        out.append({
            "from": r[0]["interval15"], "to": r[-1]["interval15"],
            "slots": len(r), "worst": max(vals) if m["direction"] == "lower_better" else min(vals),
            "avg": sum(vals) / len(vals) if vals else None,
            "has_red": any(c["status"] == "red" for c in r),
        })
    out.sort(key=lambda r: -r["slots"])
    return out


# ── M3 交叉矩阵 ──────────────────────────────────────────────────────────
def cross(engine: MetricEngine, biz_date, metric_id: str, dim_a: str, dim_b: str) -> dict:
    """维度A × 维度B 的二维表。定位到具体格子用它。"""
    m = engine.defs[metric_id]
    unsupported = [d for d in (dim_a, dim_b) if d not in m.get("dims", [])]
    note = ("" if not unsupported else
            f"指标「{m['name_cn']}」不支持按 "
            f"{'、'.join(DIM_LABEL.get(d, d) for d in unsupported)} 下钻"
            f"（可用维度：{'、'.join(DIM_LABEL.get(d, d) for d in m.get('dims', [])) or '无'}）")
    a_vals = dim_values(engine, biz_date, dim_a)
    b_vals = dim_values(engine, biz_date, dim_b)
    cells = []
    for av in a_vals:
        row = []
        for bv in b_vals:
            val = engine.compute(biz_date, [metric_id], dim_b, bv, {dim_a: av}).get(metric_id)
            contrib, _ = _contribution(engine, biz_date, m, val, dim_b, bv, {dim_a: av})
            row.append({"a": av, "b": bv, "value": val,
                        "status": cell_status(engine, biz_date, m, val, dim_b, bv, {dim_a: av}),
                        "contribution": contrib})
        cells.append(row)
    flat = [c for row in cells for c in row if c["value"] is not None]
    flat.sort(key=lambda c: -(c["contribution"] or 0))
    return {
        "kind": "cross", "metric": metric_id, "metric_cn": m["name_cn"], "unit": m["unit"],
        "dim_a": dim_a, "dim_a_cn": DIM_LABEL.get(dim_a, dim_a), "a_values": a_vals,
        "dim_b": dim_b, "dim_b_cn": DIM_LABEL.get(dim_b, dim_b), "b_values": b_vals,
        "cells": cells, "worst": flat[:5], "note": note,
    }


# ── M4 优先级矩阵 ────────────────────────────────────────────────────────
QUADRANTS = {
    ("high", "high"): ("Q1 今天就做", "影响大、你控得住 —— 这就是今天的班前会内容"),
    ("high", "low"): ("Q2 上报求援", "影响大但不在你手上 —— 找营运经理/供应商/总部，别自己耗"),
    ("low", "high"): ("Q3 顺手做掉", "影响有限但一动就好 —— 交给值班经理清单化处理"),
    ("low", "low"): ("Q4 记录观察", "别在今天的会上讲 —— 记进周复盘，看趋势再说"),
}


def priority(engine: MetricEngine, biz_date, impact_cut: float = 0.7,
             control_cut: float = 0.7, include_green: bool = False) -> dict:
    """把当日非绿灯指标放进 影响度 × 可控度 四象限，并给出行动分。"""
    card = engine.scorecard(biz_date)
    items = []
    for r in card:
        if not include_green and r["status"] in ("green", "na"):
            continue
        m = engine.defs[r["id"]]
        # 严重度用"越过绿线几个黄灯带宽"，不用相对偏差——见 metrics.severity 的说明
        sev = r["severity"]
        # 红灯的紧迫度是黄灯的两倍：同样超出一个带宽，越过红线的那个先处理
        urgency = 2.0 if r["status"] == "red" else 1.0
        score = sev * m["impact"] * m["controllability"] * urgency
        q = QUADRANTS[("high" if m["impact"] >= impact_cut else "low",
                       "high" if m["controllability"] >= control_cut else "low")]
        items.append({
            **r, "impact": m["impact"], "controllability": m["controllability"],
            "score": round(score, 2), "quadrant": q[0], "quadrant_note": q[1],
            "family_cn": FAMILY_LABEL[m["family"]],
        })
    items.sort(key=lambda x: -x["score"])
    buckets: dict[str, list] = {q[0]: [] for q in QUADRANTS.values()}
    for it in items:
        buckets[it["quadrant"]].append(it)
    return {"kind": "priority", "items": items, "quadrants": buckets,
            "top3": items[:3], "cuts": {"impact": impact_cut, "control": control_cut}}


# ── M5 责任矩阵 ──────────────────────────────────────────────────────────
OWNER_LABEL = {"RGM": "店长", "SM": "值班经理", "TL": "组长", "TRN": "训练员", "AREA": "营运经理"}
CADENCE_LABEL = {"daily": "每日", "weekly": "每周", "monthly": "每月"}


def ownership(engine: MetricEngine) -> dict:
    """静态矩阵：谁 × 什么节奏 × 盯哪些指标。用来做岗位职责和会议议程。"""
    grid: dict[str, dict[str, list[str]]] = {}
    for mid, m in engine.defs.items():
        grid.setdefault(m["owner"], {}).setdefault(m["cadence"], []).append(m["name_cn"])
    return {
        "kind": "ownership", "owners": list(OWNER_LABEL), "cadences": list(CADENCE_LABEL),
        "owner_label": OWNER_LABEL, "cadence_label": CADENCE_LABEL, "grid": grid,
    }


def family_summary(engine: MetricEngine, biz_date) -> list[dict]:
    """六大家族的健康度概览，看板最上面那一排。"""
    card = engine.scorecard(biz_date)
    out = []
    for fam, label in FAMILY_LABEL.items():
        rows = [r for r in card if r["family"] == fam]
        if not rows:
            continue
        counts = {"green": 0, "amber": 0, "red": 0, "na": 0}
        for r in rows:
            counts[r["status"]] += 1
        worst = min(rows, key=lambda r: (STATUS_RANK[r["status"]], -(abs(r["gap_pct"] or 0))))
        out.append({
            "family": fam, "family_cn": label, "counts": counts, "total": len(rows),
            "status": "red" if counts["red"] else ("amber" if counts["amber"] else "green"),
            "headline": worst,
        })
    return out
