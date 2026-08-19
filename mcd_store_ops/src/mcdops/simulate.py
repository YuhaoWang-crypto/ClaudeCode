"""门店一天的运营数据模拟器。

生成 schema.py 里那些 fct_* 表的行级数据（不是拍脑袋的汇总值），
让下游的指标引擎、矩阵、看板都跑在真实形状的数据上：
双峰客流、渠道结构随时段漂移、人力紧张时出餐变慢、慢了顾客就打低分。

模拟器故意在"今天"埋了几个可被矩阵发现的问题（见 SCENARIOS），
否则报表全绿，就演示不出"怎么用矩阵找问题"。

纯标准库，无第三方依赖。
"""

from __future__ import annotations

import math
import random
from datetime import date, datetime, timedelta

# ── 日内形状参数 ──────────────────────────────────────────────────────────
# 每个时段的客流峰值位置与胖瘦（小时制）
DAYPART_SHAPE = {
    "Breakfast": (7.75, 1.10),
    "Lunch": (12.10, 0.85),
    "Afternoon": (15.30, 1.20),
    "Dinner": (18.40, 1.00),
    "LateNight": (21.80, 1.00),
}

# 渠道在各时段的相对偏好（乘数，之后按时段归一化）
CHANNEL_TILT = {
    "Breakfast": {"FrontCounter": 1.20, "Kiosk": 0.80, "DriveThru": 1.60, "MobileApp": 1.10, "Delivery": 0.50},
    "Lunch":     {"FrontCounter": 1.00, "Kiosk": 1.15, "DriveThru": 0.90, "MobileApp": 1.10, "Delivery": 0.95},
    "Afternoon": {"FrontCounter": 0.90, "Kiosk": 1.00, "DriveThru": 1.00, "MobileApp": 1.00, "Delivery": 1.20},
    "Dinner":    {"FrontCounter": 0.95, "Kiosk": 1.10, "DriveThru": 1.00, "MobileApp": 0.95, "Delivery": 1.05},
    "LateNight": {"FrontCounter": 0.70, "Kiosk": 0.90, "DriveThru": 1.20, "MobileApp": 0.80, "Delivery": 1.70},
}

# 渠道客单价系数：点餐机会推加购、外送有起送门槛、小程序自取多为单人餐
CHANNEL_CHECK_FACTOR = {
    "FrontCounter": 1.00, "Kiosk": 1.08, "DriveThru": 1.05, "MobileApp": 0.95, "Delivery": 1.25,
}

# 星期几因子（周一=0）
DOW_FACTOR = [0.92, 0.90, 0.94, 0.99, 1.14, 1.22, 1.12]

# 一名员工在 15 分钟里"舒服"能处理的订单数（含后厨/理货/清洁/打包所有岗位）。
# 超过这个负载，出餐时间就开始非线性上升。
LOAD_REF = 1.85

SCENARIOS = {
    "clean": {},
    # 默认剧本：晚市人手不足 + 午市备货过量 + 一次冷链异常 + 外送打包环节拖后腿
    "realistic": {
        "understaff_windows": [(17.0, 19.75, 2.5)],   # (起, 止, 少几个人)
        "overprep_dayparts": {"Lunch": 1.9},           # 废弃倍数
        "cold_chain_incident": True,
        "delivery_pack_drag": 1.42,                    # 外送出餐时长倍数
    },
}


def _norm(rng: random.Random, mu: float, sigma: float, lo: float | None = None) -> float:
    v = rng.gauss(mu, sigma)
    return v if lo is None else max(lo, v)


def _pick(rng: random.Random, weights: dict[str, float]) -> str:
    total = sum(weights.values())
    r = rng.random() * total
    acc = 0.0
    for k, w in weights.items():
        acc += w
        if r <= acc:
            return k
    return list(weights)[-1]


class StoreDaySimulator:
    def __init__(self, profile: dict, seed: int = 20260819, scenario: str = "realistic"):
        self.p = profile
        self.seed = seed
        self.scn = SCENARIOS[scenario]
        self.scenario_name = scenario
        self.dayparts = {d["id"]: d for d in profile["dayparts"]}
        self.channels = {c["id"]: c for c in profile["channels"]}
        self.categories = profile["categories"]

    # ── 时间骨架 ─────────────────────────────────────────────────────────
    def intervals(self) -> list[tuple[str, float, str]]:
        """[(interval15, 小时浮点, daypart), ...]"""
        out = []
        h = self.p["open_hour"] * 4
        end = self.p["close_hour"] * 4
        while h < end:
            hour = h / 4.0
            key = f"{int(hour):02d}:{int(round((hour % 1) * 60)):02d}"
            dp = self._daypart(hour)
            if dp:
                out.append((key, hour, dp))
            h += 1
        return out

    def _daypart(self, hour: float) -> str | None:
        for d in self.p["dayparts"]:
            if d["start"] <= hour < d["end"]:
                return d["id"]
        return None

    def _channel_mix(self, daypart: str) -> dict[str, float]:
        """该时段的渠道结构（已归一化到 1）。"""
        raw = {c: self.channels[c]["share"] * CHANNEL_TILT[daypart][c] for c in self.channels}
        s = sum(raw.values())
        return {k: v / s for k, v in raw.items()}

    def interval_weights(self) -> dict[str, float]:
        """每个 15 分钟时段占全天交易的比重。"""
        raw: dict[str, float] = {}
        for key, hour, dp in self.intervals():
            peak, sigma = DAYPART_SHAPE[dp]
            raw[key] = math.exp(-(((hour + 0.125) - peak) ** 2) / (2 * sigma ** 2))
        # 按 daypart 的销售份额重新配平
        by_dp: dict[str, float] = {}
        for key, _h, dp in self.intervals():
            by_dp[dp] = by_dp.get(dp, 0.0) + raw[key]
        out = {}
        for key, _h, dp in self.intervals():
            out[key] = raw[key] / by_dp[dp] * self.dayparts[dp]["sales_share"]
        return out

    # ── 人力 ────────────────────────────────────────────────────────────
    def _understaff(self, hour: float) -> float:
        for lo, hi, n in self.scn.get("understaff_windows", []):
            if lo <= hour < hi:
                return n
        return 0.0

    def crew_curve(self, biz_date: date, gc_by_interval: dict[str, int]) -> dict[str, dict]:
        """每个时段的计划/实际在岗人数。计划跟着预测走，实际会掉人。"""
        rng = random.Random(self.seed + biz_date.toordinal() * 7 + 3)
        out = {}
        for key, hour, dp in self.intervals():
            gc = gc_by_interval.get(key, 0)
            planned = max(5.0, round(4.2 + gc / 2.4, 1))
            actual = planned - self._understaff(hour) - max(0.0, _norm(rng, 0.12, 0.40, 0.0))
            out[key] = {
                "interval15": key, "hour": hour, "daypart": dp,
                "planned_crew": round(planned, 2),
                "actual_crew": round(max(3.0, actual), 2),
            }
        return out

    # ── 主流程 ──────────────────────────────────────────────────────────
    def simulate_day(self, biz_date: date) -> dict[str, list[dict]]:
        rng = random.Random(self.seed + biz_date.toordinal())
        sid = self.p["store_id"]
        base_gc = self.p["baseline"]["guest_count"]
        base_check = self.p["baseline"]["net_sales"] / base_gc

        dow = DOW_FACTOR[biz_date.weekday()]
        day_noise = _norm(rng, 1.0, 0.045, 0.75)
        day_gc = int(base_gc * dow * day_noise)

        weights = self.interval_weights()
        wsum = sum(weights.values())
        gc_by_interval = {k: max(0, int(round(day_gc * w / wsum))) for k, w in weights.items()}

        # 预测（T-1 生成）：围绕真实值有偏差，午市系统性偏高 → 备货过量
        forecast: list[dict] = []
        frng = random.Random(self.seed + biz_date.toordinal() * 13)
        for key, hour, dp in self.intervals():
            bias = 1.06 if (dp == "Lunch" and self.scn.get("overprep_dayparts")) else 1.0
            f_gc = gc_by_interval[key] * bias * _norm(frng, 1.0, 0.09, 0.4)
            mix = self._channel_mix(dp)
            for ch, share in mix.items():
                gc_ch = f_gc * share
                forecast.append({
                    "store_id": sid, "biz_date": str(biz_date), "interval15": key, "daypart": dp,
                    "channel": ch,
                    "forecast_gc": round(gc_ch, 2),
                    # 0.94 = 预测模型已经学过折扣与退款，口径对齐"净销售"
                    "forecast_sales": round(gc_ch * base_check * CHANNEL_CHECK_FACTOR[ch] * 0.94
                                            * _norm(frng, 1.0, 0.04, 0.6), 2),
                    "forecast_version": "v_t-1",
                })

        crew = self.crew_curve(biz_date, gc_by_interval)

        txns, lines, kds, dt, deliveries, issues = [], [], [], [], [], []
        dt_seq = 0
        order_no = 0

        for key, hour, dp in self.intervals():
            n = gc_by_interval[key]
            if n == 0:
                continue
            load = n / max(1.0, crew[key]["actual_crew"])
            # 压力系数：人均负载超过 LOAD_REF 后，出餐时间非线性上升
            stress = 1.0 + 1.60 * max(0.0, load / LOAD_REF - 1.0)

            mix = self._channel_mix(dp)

            for _ in range(n):
                order_no += 1
                oid = f"{biz_date:%Y%m%d}-{order_no:05d}"
                ch = _pick(rng, mix)
                minute = rng.random() * 15
                ts = datetime.combine(biz_date, datetime.min.time()) + timedelta(hours=hour, minutes=minute)

                # ── 金额与商品行
                n_items = max(1, int(round(_norm(rng, 3.2 * CHANNEL_CHECK_FACTOR[ch] ** 0.6, 1.1, 1))))
                gross = 0.0
                is_upsell = False
                for ln in range(n_items):
                    cat = _pick(rng, {c["id"]: c["share"] for c in self.categories})
                    unit = _norm(rng, base_check * 0.305 * CHANNEL_CHECK_FACTOR[ch] ** 0.4,
                                 base_check * 0.10, 4.0)
                    upsell_p = {"Kiosk": 0.13, "FrontCounter": 0.155, "DriveThru": 0.135,
                                "MobileApp": 0.08, "Delivery": 0.06}[ch]
                    # 剧本外的常态：高峰时段人工点单顾不上推荐
                    if stress > 1.15 and ch in ("FrontCounter", "DriveThru"):
                        upsell_p *= 0.55
                    up = ln > 0 and rng.random() < upsell_p
                    is_upsell = is_upsell or up
                    gross += unit
                    lines.append({
                        "store_id": sid, "biz_date": str(biz_date), "order_id": oid, "line_no": ln + 1,
                        "sku_id": f"SKU-{cat[:3].upper()}-{rng.randint(1, 40):03d}",
                        "category": cat, "qty": 1, "unit_price": round(unit, 2),
                        "line_net": round(unit, 2), "is_upsell": up,
                        # 冗余维度：数仓里常见的做法，省掉看板每次回头 join 订单头
                        "interval15": key, "daypart": dp, "channel": ch,
                    })

                discount = round(gross * (0.12 if rng.random() < 0.34 else 0.0) * _norm(rng, 1.0, 0.25, 0.0), 2)
                refund = round(gross * 0.5, 2) if rng.random() < 0.004 else 0.0
                net = round(max(0.0, gross - discount - refund), 2)

                txns.append({
                    "store_id": sid, "biz_date": str(biz_date), "order_id": oid,
                    "order_ts": ts.isoformat(timespec="seconds"), "interval15": key, "daypart": dp,
                    "channel": ch, "gross_amount": round(gross, 2), "discount_amount": discount,
                    "refund_amount": refund, "net_amount": net,
                    "pay_type": _pick(rng, {"WeChat": 0.42, "Alipay": 0.27, "Card": 0.08,
                                            "Cash": 0.06, "Platform": 0.17}),
                    "cashier_id": f"E{rng.randint(101, 138)}" if ch in ("FrontCounter", "DriveThru") else "KIOSK",
                    "status": "refunded" if refund else "completed",
                    "has_upsell": is_upsell,
                })

                # ── 后厨制作时长（KDS 进单 → 敲单）
                prep = _norm(rng, 62 * stress + 5.0 * n_items, 15, 25)
                # 端到端服务时长：各渠道口径不同，这里统一算出来供达标率使用
                service = prep + _norm(rng, 38 * stress, 12, 8)   # 前台/点餐机/小程序：制作 + 叫号交付

                # ── 得来速分段计时
                if ch == "DriveThru":
                    dt_seq += 1
                    order_sec = _norm(rng, 42 * (1 + 0.25 * (stress - 1)), 9, 15)
                    queue_sec = _norm(rng, 46 * stress ** 1.6, 16, 5)
                    window_sec = _norm(rng, 52 * stress, 13, 18)
                    service = order_sec + queue_sec + window_sec
                    dt.append({
                        "store_id": sid, "biz_date": str(biz_date), "car_sequence": dt_seq,
                        "interval15": key, "daypart": dp, "channel": ch, "lane": "A",
                        "order_seconds": round(order_sec, 1),
                        "queue_seconds": round(queue_sec, 1),
                        "window_seconds": round(window_sec, 1),
                        "total_seconds": round(service, 1),
                        "pull_forward_flag": window_sec > 110 and rng.random() < 0.3,
                    })

                # ── 外送
                if ch == "Delivery":
                    drag = self.scn.get("delivery_pack_drag", 1.0)
                    handover = _norm(rng, 262 * stress * drag, 62, 90)
                    service = handover
                    plat = _pick(rng, {"Meituan": 0.5, "Eleme": 0.33, "McDeliveryApp": 0.17})
                    late = handover > 420
                    deliveries.append({
                        "store_id": sid, "biz_date": str(biz_date), "platform": plat,
                        "platform_order_id": f"{plat[:2].upper()}{order_no:06d}", "pos_order_id": oid,
                        "interval15": key, "daypart": dp, "channel": ch,
                        "handover_seconds": round(handover, 1),
                        "commission_amount": round(net * 0.20, 2),
                        "rating": self._delivery_rating(rng, late),
                        "issue_flag": "超时" if late and rng.random() < 0.45 else "无",
                    })

                kds.append({
                    "store_id": sid, "biz_date": str(biz_date), "order_id": oid, "interval15": key,
                    "daypart": dp, "channel": ch, "station": "Assembly",
                    "prep_seconds": round(prep, 1),
                    "service_seconds": round(service, 1),
                    "speed_target_sec": self.channels[ch]["speed_target_sec"],
                    "on_time": service <= self.channels[ch]["speed_target_sec"],
                    "item_count": n_items,
                    "recall_flag": rng.random() < 0.006,
                })

                # ── 差错事件：压力越大越容易出错
                if rng.random() < 0.0011 * (1 + 1.8 * (stress - 1)):
                    issues.append({
                        "store_id": sid, "biz_date": str(biz_date), "issue_id": f"IS-{biz_date:%m%d}-{len(issues)+1:03d}",
                        "order_id": oid, "interval15": key, "daypart": dp, "channel": ch,
                        "source": _pick(rng, {"店内退换": 0.45, "平台售后": 0.35, "顾客热线": 0.2}),
                        "issue_type": _pick(rng, {"漏餐": 0.32, "错餐": 0.26, "超时": 0.22,
                                                  "温度": 0.12, "异物": 0.02, "态度": 0.06}),
                        "amount": round(_norm(rng, 22, 12, 5), 2),
                    })

        day = {
            "fct_pos_transaction": txns,
            "fct_pos_line_item": lines,
            "fct_kds_order": kds,
            "fct_dt_timer": dt,
            "fct_delivery_order": deliveries,
            "fct_order_issue": issues,
            "dim_sales_forecast": forecast,
        }
        day.update(self._people(biz_date, crew))
        day.update(self._inventory(biz_date, txns, lines, rng))
        day.update(self._quality(biz_date, rng))
        day.update(self._guest(biz_date, txns, kds, dt, deliveries, issues, rng))
        day.update(self._cash(biz_date, txns, rng))
        return day

    # ── 人员 ────────────────────────────────────────────────────────────
    def _people(self, biz_date: date, crew: dict) -> dict:
        rng = random.Random(self.seed + biz_date.toordinal() * 11)
        sid = self.p["store_id"]
        planned_hours = sum(c["planned_crew"] for c in crew.values()) * 0.25
        actual_hours = sum(c["actual_crew"] for c in crew.values()) * 0.25
        mgr = self.p["labor"]["manager_hours_per_day"]

        # 把连续的人力曲线摊回到"班次"上，便于按人复盘
        shifts, cards = [], []
        n_shifts = max(12, int(round(actual_hours / 6.0)))
        for i in range(n_shifts):
            emp = f"E{100 + i:03d}"
            start = self.p["open_hour"] + rng.choice([0, 2, 4, 5, 6, 8, 9, 11, 12, 13, 14])
            length = rng.choice([4.0, 5.0, 6.0, 7.5, 8.0])
            station = rng.choice(self.p["stations"])
            late = 0
            if rng.random() < 0.09:
                late = rng.choice([6, 8, 12, 18, 25])
            no_show = rng.random() < 0.035
            # 该员工在所站岗位是否已认证（真实来源是 fct_training_record join fct_schedule）
            cert_base = {"Kitchen_Fry": 0.78, "DriveThru_Window": 0.80, "Delivery_Pack": 0.82}.get(station, 0.94)
            shifts.append({
                "store_id": sid, "biz_date": str(biz_date), "sched_id": f"S-{biz_date:%m%d}-{i:02d}",
                "emp_id": emp, "station": station, "start_hour": start, "planned_hours": length,
                "daypart": self._daypart(min(start, self.p["close_hour"] - 0.1)) or "Lunch",
                "certified": rng.random() < cert_base,
            })
            if not no_show:
                cards.append({
                    "store_id": sid, "biz_date": str(biz_date), "emp_id": emp, "station": station,
                    "start_hour": start, "actual_hours": round(length + _norm(rng, 0.05, 0.28), 2),
                    "late_minutes": late,
                    "daypart": self._daypart(min(start, self.p["close_hour"] - 0.1)) or "Lunch",
                })

        return {
            "fct_schedule": shifts,
            "fct_timecard": cards,
            "agg_labor_interval": [
                {"store_id": sid, "biz_date": str(biz_date), **c} for c in crew.values()
            ],
            "agg_labor_day": [{
                "store_id": sid, "biz_date": str(biz_date),
                "planned_hours": round(planned_hours + mgr, 2),
                "actual_hours": round(actual_hours + mgr, 2),
                "avg_wage": self.p["labor"]["avg_wage_per_hour"],
            }],
        }

    # ── 库存 / 废弃 / 成本 ───────────────────────────────────────────────
    def _inventory(self, biz_date: date, txns: list, lines: list, rng: random.Random) -> dict:
        sid = self.p["store_id"]
        cost_ratio = {c["id"]: c["cost_ratio"] for c in self.categories}
        theo = sum(l["line_net"] * cost_ratio[l["category"]] for l in lines)

        waste, batches = [], []
        overprep = self.scn.get("overprep_dayparts", {})
        by_dp_sales: dict[str, float] = {}
        for t in txns:
            by_dp_sales[t["daypart"]] = by_dp_sales.get(t["daypart"], 0.0) + t["net_amount"]

        for dp, sales in by_dp_sales.items():
            mult = overprep.get(dp, 1.0)
            for cat in ("Burger", "Chicken", "Fries", "Beverage", "McCafe", "Dessert"):
                base_rate = {"Burger": 0.00225, "Chicken": 0.00204, "Fries": 0.00156,
                             "Beverage": 0.00054, "McCafe": 0.00090, "Dessert": 0.00180}[cat]
                for wtype, share in (("complete", 0.5), ("raw", 0.3), ("holding", 0.2)):
                    amt = sales * base_rate * share * mult * _norm(rng, 1.0, 0.3, 0.1)
                    if amt < 0.5:
                        continue
                    waste.append({
                        "store_id": sid, "biz_date": str(biz_date), "daypart": dp,
                        "waste_id": f"W-{biz_date:%m%d}-{len(waste)+1:03d}",
                        "waste_type": wtype, "category": cat, "cost": round(amt, 2),
                        "reason": {"complete": "做多了", "raw": "过期/掉落", "holding": "保温超时"}[wtype],
                        "logged_by": f"E{rng.randint(101, 138)}",
                    })
            # 保温批次
            for cat in ("Burger", "Chicken", "Fries"):
                made = max(1, int(sales * 0.012 * _norm(rng, 1.0, 0.15, 0.4)))
                discard = int(made * 0.022 * mult * _norm(rng, 1.0, 0.4, 0.1))
                batches.append({
                    "store_id": sid, "biz_date": str(biz_date), "daypart": dp, "category": cat,
                    "batch_id": f"B-{biz_date:%m%d}-{len(batches)+1:03d}",
                    "qty_made": made, "qty_discard": discard,
                    "triggered_by": "预测建议" if rng.random() < 0.7 else "人工判断",
                })

        waste_cost = sum(w["cost"] for w in waste)
        # 盘点差异：偷跑的成本（少记废弃、多打料、收货短少）
        variance = max(0.0, _norm(rng, theo * 0.011, theo * 0.006, 0.0))

        return {
            "fct_waste_log": waste,
            "fct_holding_batch": batches,
            "agg_food_cost_day": [{
                "store_id": sid, "biz_date": str(biz_date),
                "theoretical_cost": round(theo, 2),
                "waste_cost": round(waste_cost, 2),
                "inventory_variance_cost": round(variance, 2),
                "actual_cost": round(theo + waste_cost + variance, 2),
            }],
        }

    # ── 食安 / 品质 / 设备 ───────────────────────────────────────────────
    def _quality(self, biz_date: date, rng: random.Random) -> dict:
        sid = self.p["store_id"]
        checks, temps, walks, tickets = [], [], [], []

        for shift in ("早班", "中班", "晚班"):
            for i in range(18):
                critical = i < 6
                fail = rng.random() < (0.012 if critical else 0.045)
                checks.append({
                    "store_id": sid, "biz_date": str(biz_date), "shift": shift,
                    "check_id": f"FS-{biz_date:%m%d}-{shift}", "item_code": f"FS-{i+1:02d}",
                    "is_critical": critical, "result": "fail" if fail else "pass",
                    "corrective_action": "立即整改并复检" if fail else "",
                })

        incident = self.scn.get("cold_chain_incident") and biz_date.toordinal() % 6 == 0
        for eq in self.p["equipment"]:
            rng_c = eq.get("temp_range_c")
            if not rng_c:
                continue
            lo, hi = rng_c
            mid = (lo + hi) / 2
            for slot in range(int((self.p["close_hour"] - self.p["open_hour"]) * 2)):
                t = _norm(rng, mid, (hi - lo) * 0.16)
                if incident and eq["id"] == "WALKIN-C" and 14 <= slot <= 17:
                    t = hi + _norm(rng, 2.4, 0.8, 0.3)
                temps.append({
                    "store_id": sid, "biz_date": str(biz_date), "equipment_id": eq["id"],
                    "hour": round(self.p["open_hour"] + slot * 0.5, 2),
                    "temp_c": round(t, 2), "in_range": lo <= t <= hi, "source": "iot",
                })

        for shift in ("早班", "中班", "晚班"):
            for zone in ("大堂", "洗手间", "厨房", "外围", "得来速车道"):
                base = {"大堂": 93, "洗手间": 88, "厨房": 95, "外围": 91, "得来速车道": 92}[zone]
                walks.append({
                    "store_id": sid, "biz_date": str(biz_date), "shift": shift, "zone": zone,
                    "score": max(50, min(100, round(_norm(rng, base, 6)))),
                })

        if rng.random() < 0.22:
            eq = rng.choice([e for e in self.p["equipment"] if e.get("critical")])
            down = int(_norm(rng, 55, 35, 10))
            tickets.append({
                "store_id": sid, "biz_date": str(biz_date), "ticket_id": f"TK-{biz_date:%m%d}-1",
                "equipment_id": eq["id"], "downtime_minutes": down,
                "severity": "降级" if down < 90 else "停售", "vendor": "区域服务商",
            })

        return {
            "fct_food_safety_check": checks,
            "fct_temperature_log": temps,
            "fct_travel_path_check": walks,
            "fct_equipment_ticket": tickets,
        }

    # ── 顾客反馈：分数由当天真实体验反推 ──────────────────────────────────
    def _delivery_rating(self, rng: random.Random, late: bool) -> int:
        if late:
            return _pick_int(rng, {5: 0.30, 4: 0.28, 3: 0.20, 2: 0.12, 1: 0.10})
        return _pick_int(rng, {5: 0.74, 4: 0.18, 3: 0.05, 2: 0.02, 1: 0.01})

    def _guest(self, biz_date, txns, kds, dt, deliveries, issues, rng) -> dict:
        sid = self.p["store_id"]
        prep_by_dp: dict[str, list[float]] = {}
        for k in kds:
            prep_by_dp.setdefault(k["daypart"], []).append(k["prep_seconds"])
        issue_orders = {i["order_id"] for i in issues}

        surveys = []
        for t in txns:
            if rng.random() > 0.022:
                continue
            prep = prep_by_dp.get(t["daypart"], [90])
            avg_prep = sum(prep) / len(prep)
            slow_penalty = max(0.0, (avg_prep - 100) / 45)
            score = 4.55 - slow_penalty * 0.9 - (1.9 if t["order_id"] in issue_orders else 0)
            score = max(1, min(5, int(round(_norm(rng, score, 0.55)))))
            surveys.append({
                "store_id": sid, "biz_date": str(biz_date), "survey_id": f"SV-{biz_date:%m%d}-{len(surveys)+1:03d}",
                "source": "receipt", "channel": t["channel"], "daypart": t["daypart"],
                "score": score,
                "comment_tag": "速度" if slow_penalty > 0.4 else rng.choice(["口味", "干净", "态度", "速度"]),
            })
        for d in deliveries:
            if rng.random() < 0.16:
                surveys.append({
                    "store_id": sid, "biz_date": str(biz_date),
                    "survey_id": f"SV-{biz_date:%m%d}-D{len(surveys)+1:03d}",
                    "source": d["platform"].lower(), "channel": "Delivery", "daypart": d["daypart"],
                    "score": d["rating"], "comment_tag": "超时" if d["issue_flag"] != "无" else "口味",
                })
        return {"fct_guest_survey": surveys}

    # ── 现金 ────────────────────────────────────────────────────────────
    def _cash(self, biz_date: date, txns: list, rng: random.Random) -> dict:
        sid = self.p["store_id"]
        cash_sales = sum(t["net_amount"] for t in txns if t["pay_type"] == "Cash")
        rows = []
        for i, shift in enumerate(("早班", "中班", "晚班")):
            sys_amt = cash_sales / 3
            var = _norm(rng, 0.0, 9.0) if rng.random() > 0.12 else _norm(rng, -28.0, 18.0)
            rows.append({
                "store_id": sid, "biz_date": str(biz_date), "shift": shift,
                "declaration_id": f"CD-{biz_date:%m%d}-{i+1}",
                "system_amount": round(sys_amt, 2),
                "declared_amount": round(sys_amt + var, 2),
                "variance": round(var, 2),
                "cashier_id": f"E{rng.randint(101, 138)}",
            })
        return {"fct_cash_declaration": rows}


def _pick_int(rng: random.Random, weights: dict[int, float]) -> int:
    total = sum(weights.values())
    r = rng.random() * total
    acc = 0.0
    for k, w in weights.items():
        acc += w
        if r <= acc:
            return k
    return list(weights)[-1]


def _roster(profile: dict, end_date: date, seed: int) -> list[dict]:
    """员工花名册（维度表，不随日期增长）。用于 30 日流失率。"""
    rng = random.Random(seed + 999)
    rows = []
    for i in range(72):
        emp = f"E{100 + i:03d}"
        tenure_days = int(abs(rng.gauss(0, 1)) * 420) + 15
        hire = end_date - timedelta(days=tenure_days)
        term = None
        # 短工龄的人离职率高得多——这是餐饮用工的真实形状
        p_term = 0.34 if tenure_days < 90 else (0.12 if tenure_days < 365 else 0.05)
        if rng.random() < p_term:
            term = end_date - timedelta(days=rng.randint(1, 75))
            if term < hire:
                term = None
        rows.append({
            "store_id": profile["store_id"], "emp_id": emp,
            "role": "店长" if i == 0 else ("值班经理" if i < 4 else ("组长" if i < 10 else "员工")),
            "hire_date": str(hire), "term_date": str(term) if term else None,
            "tenure_days": tenure_days,
            "tenure_band": "<90天" if tenure_days < 90 else ("90天-1年" if tenure_days < 365 else ">1年"),
            "wage_per_hour": round(profile["labor"]["avg_wage_per_hour"] * rng.uniform(0.85, 1.35), 2),
            "employment_type": "全职" if i < 22 else "兼职",
        })
    return rows


def simulate_period(profile: dict, end_date: date, days: int = 14,
                    seed: int = 20260819, scenario: str = "realistic") -> dict[str, list[dict]]:
    """生成 [end_date - days + 1, end_date] 的全部数据，表名 → 行列表。"""
    sim = StoreDaySimulator(profile, seed=seed, scenario=scenario)
    merged: dict[str, list[dict]] = {"dim_employee": _roster(profile, end_date, seed)}
    for i in range(days - 1, -1, -1):
        d = end_date - timedelta(days=i)
        for table, rows in sim.simulate_day(d).items():
            merged.setdefault(table, []).extend(rows)
    return merged
