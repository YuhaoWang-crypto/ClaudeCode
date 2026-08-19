"""指标引擎：把行级数据算成指标，并按任意维度切片。

两个设计决定值得说清楚（`docs/03_metric_matrix.md` 展开讲）：

1. 一个指标只写一次公式，靠 `Slice` 传进来的切片决定它算的是全店、某个时段、
   某个渠道还是某个岗位。这样"矩阵"不是另一套逻辑，只是同一个公式跑 N 遍。

2. 率型指标的分母分两种口径，写在 metrics.json 的 `denominator` 字段里：
   - slice     （默认）分母取切片自己的量 → 回答"这一块自己好不好"
                例：晚市人工成本率 = 晚市人工 / 晚市销售
   - store_day  分母取全店当日的量 → 回答"这一块拖了全店多少"，好处是可加、可排序
                例：午市废弃率 0.6pp + 晚市 0.4pp = 全店 1.0pp
   归因矩阵要的是第二种；单块体检要的是第一种。混用会得出互相矛盾的结论。
"""

from __future__ import annotations

import json
from datetime import date, timedelta
from pathlib import Path

CONFIG_DIR = Path(__file__).resolve().parents[2] / "config"


# ── 配置装载 ──────────────────────────────────────────────────────────────
def load_metrics(path: Path | None = None) -> dict[str, dict]:
    p = path or CONFIG_DIR / "metrics.json"
    raw = json.loads(p.read_text(encoding="utf-8"))
    return {m["id"]: m for m in raw["metrics"]}


def load_profile(path: Path | None = None) -> dict:
    p = path or CONFIG_DIR / "store_profile.json"
    return json.loads(p.read_text(encoding="utf-8"))


# ── 切片 ─────────────────────────────────────────────────────────────────
class Slice:
    """某一天的数据，可选再按一个维度取值收窄。"""

    def __init__(self, data: dict[str, list[dict]], biz_date, profile: dict,
                 dim_name: str | None = None, dim_value=None,
                 extra_filters: dict | None = None):
        self.data = data
        self.biz_date = str(biz_date)
        self.profile = profile
        self.dim_name = dim_name
        self.dim_value = dim_value
        # 交叉矩阵用：在主维度之外再钉住一个维度（如 daypart=Dinner × channel=DriveThru）
        self.extra_filters = extra_filters or {}
        self._day: dict[str, list[dict]] = {}
        self._cut: dict[str, list[dict]] = {}

    def day_rows(self, table: str) -> list[dict]:
        """只按营业日过滤（维度表不含 biz_date，则原样返回）。"""
        if table not in self._day:
            rows = self.data.get(table, [])
            if rows and "biz_date" in rows[0]:
                rows = [r for r in rows if r["biz_date"] == self.biz_date]
            self._day[table] = rows
        return self._day[table]

    def rows(self, table: str) -> list[dict]:
        """按营业日 + 当前维度过滤。表里没有这个维度字段 → 视为不适用，返回空。"""
        if table not in self._cut:
            rows = self.day_rows(table)
            conds = dict(self.extra_filters)
            if self.dim_name is not None:
                conds[self.dim_name] = self.dim_value
            for k, v in conds.items():
                if rows and k not in rows[0]:
                    rows = []
                    break
                rows = [r for r in rows if r.get(k) == v]
            self._cut[table] = rows
        return self._cut[table]

    def sub(self, dim_name: str | None, dim_value=None, extra_filters: dict | None = None) -> "Slice":
        return Slice(self.data, self.biz_date, self.profile, dim_name, dim_value, extra_filters)


# ── 小工具 ───────────────────────────────────────────────────────────────
def _avg(xs) -> float | None:
    xs = list(xs)
    return sum(xs) / len(xs) if xs else None


def _div(a, b) -> float | None:
    return None if not b else a / b


def _cost_ratio(s: Slice) -> dict[str, float]:
    return {c["id"]: c["cost_ratio"] for c in s.profile["categories"]}


def _net_sales(s: Slice) -> float:
    """切片净销售。品类维度下订单头没有品类，改从订单行取。"""
    if s.dim_name == "category":
        return sum(l["line_net"] for l in s.rows("fct_pos_line_item"))
    return sum(t["net_amount"] for t in s.rows("fct_pos_transaction"))


def _store_day_net_sales(s: Slice) -> float:
    return sum(t["net_amount"] for t in s.day_rows("fct_pos_transaction"))


def _denominator_sales(s: Slice, metric: dict) -> float:
    if metric.get("denominator") == "store_day":
        return _store_day_net_sales(s)
    return _net_sales(s)


def _gc(s: Slice) -> int:
    txns = s.rows("fct_pos_transaction")
    if txns:
        return len(txns)
    if s.dim_name == "category":
        return len({l["order_id"] for l in s.rows("fct_pos_line_item")})
    return 0


def _labor_hours(s: Slice) -> float | None:
    """切片实际工时。不同维度只能从不同的表拿到。"""
    if s.dim_name is None:
        rows = s.day_rows("agg_labor_day")
        return rows[0]["actual_hours"] if rows else None
    if s.dim_name == "station":
        cards = s.rows("fct_timecard")
        return sum(c["actual_hours"] for c in cards) if cards else None
    rows = s.rows("agg_labor_interval")
    return sum(r["actual_crew"] for r in rows) * 0.25 if rows else None


def _wage(s: Slice) -> float:
    return s.profile["labor"]["avg_wage_per_hour"]


# ── 指标公式 ─────────────────────────────────────────────────────────────
def m_net_sales(s, m):
    return _net_sales(s) or None


def m_sales_vs_forecast(s, m):
    fc = sum(f["forecast_sales"] for f in s.rows("dim_sales_forecast"))
    return _div(_net_sales(s) * 100, fc)


def m_guest_count(s, m):
    return _gc(s) or None


def m_avg_check(s, m):
    return _div(_net_sales(s), _gc(s))


def m_items_per_check(s, m):
    lines = s.rows("fct_pos_line_item")
    return _div(sum(l["qty"] for l in lines), _gc(s))


def m_upsell_rate(s, m):
    manual = [t for t in s.rows("fct_pos_transaction") if t["cashier_id"] != "KIOSK"]
    return _div(sum(1 for t in manual if t["has_upsell"]) * 100, len(manual))


def m_delivery_share(s, m):
    txns = s.rows("fct_pos_transaction")
    total = sum(t["net_amount"] for t in txns)
    deliv = sum(t["net_amount"] for t in txns if t["channel"] == "Delivery")
    return _div(deliv * 100, total)


def m_dt_total_time(s, m):
    return _avg(d["total_seconds"] for d in s.rows("fct_dt_timer"))


def m_dt_window_time(s, m):
    return _avg(d["window_seconds"] for d in s.rows("fct_dt_timer"))


def m_kds_prep_time(s, m):
    return _avg(k["prep_seconds"] for k in s.rows("fct_kds_order"))


def m_front_otd(s, m):
    rows = [k for k in s.rows("fct_kds_order")
            if k["channel"] in ("FrontCounter", "Kiosk", "MobileApp")]
    return _avg(k["service_seconds"] for k in rows)


def m_delivery_handover_time(s, m):
    return _avg(d["handover_seconds"] for d in s.rows("fct_delivery_order"))


def m_speed_hit_rate(s, m):
    rows = s.rows("fct_kds_order")
    return _div(sum(1 for k in rows if k["on_time"]) * 100, len(rows))


def m_order_accuracy(s, m):
    gc = _gc(s)
    n_issue = len(s.rows("fct_order_issue"))
    return None if not gc else (1 - n_issue / gc) * 100


def m_food_cost_pct(s, m):
    denom = _denominator_sales(s, m)
    if s.dim_name == "category":
        ratio = _cost_ratio(s)
        theo = sum(l["line_net"] * ratio[l["category"]] for l in s.rows("fct_pos_line_item"))
        waste = sum(w["cost"] for w in s.rows("fct_waste_log"))
        return _div((theo + waste) * 100, denom)
    rows = s.day_rows("agg_food_cost_day")
    if not rows or s.dim_name is not None:
        return None
    return _div(rows[0]["actual_cost"] * 100, denom)


def m_waste_pct(s, m):
    waste = sum(w["cost"] for w in s.rows("fct_waste_log"))
    return _div(waste * 100, _denominator_sales(s, m))


def m_labor_cost_pct(s, m):
    hours = _labor_hours(s)
    if hours is None:
        return None
    return _div(hours * _wage(s) * 100, _net_sales(s))


def m_splh(s, m):
    hours = _labor_hours(s)
    return None if not hours else _div(_net_sales(s), hours)


def m_labor_schedule_variance(s, m):
    if s.dim_name in (None, "daypart", "interval15"):
        rows = s.rows("agg_labor_interval") if s.dim_name else s.day_rows("agg_labor_interval")
        if not rows:
            return None
        planned = sum(r["planned_crew"] for r in rows)
        actual = sum(r["actual_crew"] for r in rows)
        return _div(abs(actual - planned) * 100, planned)
    sched = sum(x["planned_hours"] for x in s.rows("fct_schedule"))
    act = sum(x["actual_hours"] for x in s.rows("fct_timecard"))
    return _div(abs(act - sched) * 100, sched)


def m_cash_over_short(s, m):
    rows = s.rows("fct_cash_declaration")
    return sum(abs(r["variance"]) for r in rows) if rows else None


def m_food_safety_score(s, m):
    rows = s.rows("fct_food_safety_check")
    if not rows:
        return None
    by_shift: dict[str, float] = {}
    for r in rows:
        pen = 0.0
        if r["result"] == "fail":
            pen = 8.0 if r["is_critical"] else 2.0
        by_shift[r["shift"]] = by_shift.get(r["shift"], 0.0) + pen
    return _avg(max(0.0, 100.0 - p) for p in by_shift.values())


def m_critical_temp_fail(s, m):
    rows = s.rows("fct_temperature_log")
    return float(sum(1 for r in rows if not r["in_range"])) if rows else None


def m_holding_time_violation(s, m):
    rows = s.rows("fct_holding_batch")
    made = sum(r["qty_made"] for r in rows)
    return _div(sum(r["qty_discard"] for r in rows) * 100, made)


def m_cleanliness_score(s, m):
    return _avg(r["score"] for r in s.rows("fct_travel_path_check"))


def m_equipment_downtime(s, m):
    rows = s.rows("fct_equipment_ticket") if s.dim_name else s.day_rows("fct_equipment_ticket")
    return float(sum(r["downtime_minutes"] for r in rows))


def m_attendance_rate(s, m):
    sched = s.rows("fct_schedule")
    cards = s.rows("fct_timecard")
    return _div(len(cards) * 100, len(sched))


def m_late_rate(s, m):
    cards = s.rows("fct_timecard")
    return _div(sum(1 for c in cards if c["late_minutes"] > 5) * 100, len(cards))


def m_station_certification(s, m):
    sched = s.rows("fct_schedule")
    return _div(sum(1 for x in sched if x["certified"]) * 100, len(sched))


def m_crew_turnover_30d(s, m):
    emps = s.rows("dim_employee")
    if not emps:
        return None
    today = date.fromisoformat(s.biz_date)
    cutoff = today - timedelta(days=30)
    left = [e for e in emps if e["term_date"] and cutoff <= date.fromisoformat(e["term_date"]) <= today]
    opening = [e for e in emps
               if date.fromisoformat(e["hire_date"]) <= cutoff
               and (not e["term_date"] or date.fromisoformat(e["term_date"]) >= cutoff)]
    return _div(len(left) * 100, len(opening))


def m_csat(s, m):
    rows = [r for r in s.rows("fct_guest_survey") if r["source"] == "receipt"]
    avg = _avg(r["score"] for r in rows)
    return None if avg is None else avg * 20


def m_delivery_rating(s, m):
    return _avg(d["rating"] for d in s.rows("fct_delivery_order"))


def m_complaints_per_1k(s, m):
    issues = s.rows("fct_order_issue")
    gc_total = len(s.day_rows("fct_pos_transaction")) if m.get("denominator") == "store_day" else _gc(s)
    return _div(len(issues) * 1000, gc_total)


CALCULATORS = {
    "net_sales": m_net_sales, "sales_vs_forecast": m_sales_vs_forecast,
    "guest_count": m_guest_count, "avg_check": m_avg_check,
    "items_per_check": m_items_per_check, "upsell_rate": m_upsell_rate,
    "delivery_share": m_delivery_share,
    "dt_total_time": m_dt_total_time, "dt_window_time": m_dt_window_time,
    "kds_prep_time": m_kds_prep_time, "front_otd": m_front_otd,
    "delivery_handover_time": m_delivery_handover_time,
    "speed_hit_rate": m_speed_hit_rate, "order_accuracy": m_order_accuracy,
    "food_cost_pct": m_food_cost_pct, "waste_pct": m_waste_pct,
    "labor_cost_pct": m_labor_cost_pct, "splh": m_splh,
    "labor_schedule_variance": m_labor_schedule_variance, "cash_over_short": m_cash_over_short,
    "food_safety_score": m_food_safety_score, "critical_temp_fail": m_critical_temp_fail,
    "holding_time_violation": m_holding_time_violation, "cleanliness_score": m_cleanliness_score,
    "equipment_downtime": m_equipment_downtime,
    "attendance_rate": m_attendance_rate, "late_rate": m_late_rate,
    "station_certification": m_station_certification, "crew_turnover_30d": m_crew_turnover_30d,
    "csat": m_csat, "delivery_rating": m_delivery_rating,
    "complaints_per_1k": m_complaints_per_1k,
}


# ── 红黄绿判定 ────────────────────────────────────────────────────────────
def status(metric: dict, value: float | None) -> str:
    if value is None:
        return "na"
    if metric["direction"] == "higher_better":
        if value >= metric["green"]:
            return "green"
        return "amber" if value >= metric["amber"] else "red"
    # lower_better 和 neutral 都是"越过上界越糟"
    if value <= metric["green"]:
        return "green"
    return "amber" if value <= metric["amber"] else "red"


def gap_to_target(metric: dict, value: float | None) -> float | None:
    """带符号的达标差：正数=好于目标。"""
    if value is None:
        return None
    d = value - metric["target"]
    return d if metric["direction"] == "higher_better" else -d


def severity(metric: dict, value: float | None) -> float:
    """越过绿线多少个"黄灯带宽"。0 = 还在绿区，1.0 = 正好踩在红线上。

    不能用相对偏差 (gap%) 来排优先级：迟到率目标 2%、实际 10%，相对偏差 -400%，
    会把一个小问题顶到当天第一。带宽标准化让不同量纲、不同基数的指标可比。
    """
    if value is None:
        return 0.0
    g, a = metric["green"], metric["amber"]
    band = abs(a - g) or 1.0
    d = (g - value) if metric["direction"] == "higher_better" else (value - g)
    return max(0.0, d / band)


def gap_pct(metric: dict, value: float | None) -> float | None:
    g = gap_to_target(metric, value)
    t = metric["target"]
    if g is None or t == 0:
        return None
    return g / abs(t) * 100


# ── 引擎 ─────────────────────────────────────────────────────────────────
class MetricEngine:
    def __init__(self, data: dict[str, list[dict]], profile: dict, defs: dict[str, dict] | None = None):
        self.data = data
        self.profile = profile
        self.defs = defs or load_metrics()

    def compute(self, biz_date, metric_ids=None, dim_name=None, dim_value=None,
                extra_filters: dict | None = None) -> dict[str, float | None]:
        s = Slice(self.data, biz_date, self.profile, dim_name, dim_value, extra_filters)
        ids = metric_ids or list(self.defs)
        out = {}
        for mid in ids:
            m = self.defs[mid]
            if dim_name and dim_name not in m.get("dims", []) and dim_name != "interval15":
                continue
            if dim_name == "interval15" and m.get("grain") != "15min":
                continue
            fn = CALCULATORS.get(mid)
            out[mid] = fn(s, m) if fn else None
        return out

    def scorecard(self, biz_date, metric_ids=None) -> list[dict]:
        """全店当日记分卡，按家族排序。"""
        vals = self.compute(biz_date, metric_ids)
        rows = []
        for mid, v in vals.items():
            m = self.defs[mid]
            rows.append({
                "id": mid, "name_cn": m["name_cn"], "family": m["family"], "unit": m["unit"],
                "value": v, "target": m["target"], "status": status(m, v),
                "gap": gap_to_target(m, v), "gap_pct": gap_pct(m, v),
                "severity": severity(m, v),
                "owner": m["owner"], "playbook": m.get("playbook"),
            })
        order = {"SALES": 0, "SPEED": 1, "COST": 2, "QSC": 3, "PEOPLE": 4, "GUEST": 5}
        rows.sort(key=lambda r: (order.get(r["family"], 9), r["id"]))
        return rows

    def trend(self, metric_id: str, end_date: date, days: int = 14) -> list[tuple[str, float | None]]:
        out = []
        for i in range(days - 1, -1, -1):
            d = end_date - timedelta(days=i)
            out.append((str(d), self.compute(d, [metric_id]).get(metric_id)))
        return out
