"""生成一张自包含的 HTML 门店看板（无外部依赖，双击即可打开）。

看板的分层原则见 docs/04_dashboard_spec.md，一句话版本：

  第 0 屏（3 秒）  今天是红是绿 —— 六个家族的灯 + 四个数
  第 1 屏（30 秒） 今天做什么   —— 优先级矩阵的三件事，带下钻结论
  第 2 屏（3 分钟） 为什么       —— 时段矩阵 + 交叉矩阵
  第 3 屏（复盘）  全量          —— 记分卡 + 14 天趋势

一个屏里塞不下的指标，不是"再加一行"，是"往下一层放"。
"""

from __future__ import annotations

import html
from datetime import date

from . import matrix as mx
from .metrics import MetricEngine
from .report import fmt

COLORS = {"green": "var(--ok)", "amber": "var(--warn)", "red": "var(--bad)", "na": "var(--muted)"}

CSS = """
:root{
  --bg:#f6f7f9; --panel:#fff; --ink:#16191d; --muted:#8b93a1; --line:#e3e6ea;
  --ok:#1f9d55; --warn:#d98315; --bad:#d4372c;
  --ok-bg:#e8f5ee; --warn-bg:#fdf1e0; --bad-bg:#fbe9e7;
  --accent:#b8232f;
}
@media (prefers-color-scheme: dark){
  :root:not([data-theme="light"]){
    --bg:#0f1114; --panel:#171a1f; --ink:#e8eaed; --muted:#8b93a1; --line:#272c33;
    --ok:#4ade80; --warn:#fbbf24; --bad:#f87171;
    --ok-bg:#132a1d; --warn-bg:#2e2410; --bad-bg:#33191a;
    --accent:#ef4b58;
  }
}
:root[data-theme="dark"]{
  --bg:#0f1114; --panel:#171a1f; --ink:#e8eaed; --muted:#8b93a1; --line:#272c33;
  --ok:#4ade80; --warn:#fbbf24; --bad:#f87171;
  --ok-bg:#132a1d; --warn-bg:#2e2410; --bad-bg:#33191a; --accent:#ef4b58;
}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);
  font:15px/1.55 -apple-system,BlinkMacSystemFont,"PingFang SC","Microsoft YaHei",sans-serif;}
.wrap{max-width:1180px;margin:0 auto;padding:24px 18px 64px}
header{display:flex;flex-wrap:wrap;align-items:baseline;gap:12px;margin-bottom:6px}
h1{font-size:21px;margin:0;letter-spacing:-.2px}
h2{font-size:15px;margin:32px 0 12px;color:var(--muted);font-weight:600;
   text-transform:uppercase;letter-spacing:.08em}
.sub{color:var(--muted);font-size:13px}
.grid{display:grid;gap:12px}
.kpis{grid-template-columns:repeat(auto-fit,minmax(150px,1fr))}
.fams{grid-template-columns:repeat(auto-fit,minmax(150px,1fr))}
.card{background:var(--panel);border:1px solid var(--line);border-radius:10px;padding:14px 16px}
.card .lab{font-size:12px;color:var(--muted);margin-bottom:6px}
.card .val{font-size:26px;font-weight:650;letter-spacing:-.5px;line-height:1.1}
.card .foot{font-size:12px;color:var(--muted);margin-top:6px}
.fam{border-left:4px solid var(--line)}
.fam.green{border-left-color:var(--ok)} .fam.amber{border-left-color:var(--warn)}
.fam.red{border-left-color:var(--bad)}
.pill{display:inline-block;padding:1px 8px;border-radius:99px;font-size:11px;font-weight:600}
.pill.green{background:var(--ok-bg);color:var(--ok)}
.pill.amber{background:var(--warn-bg);color:var(--warn)}
.pill.red{background:var(--bad-bg);color:var(--bad)}
.pill.na{background:var(--line);color:var(--muted)}
.todo{background:var(--panel);border:1px solid var(--line);border-left:4px solid var(--bad);
      border-radius:10px;padding:14px 16px;margin-bottom:10px}
.todo.amber{border-left-color:var(--warn)}
.todo h3{margin:0 0 4px;font-size:16px}
.todo .meta{font-size:12px;color:var(--muted);margin-bottom:8px}
.todo ul{margin:6px 0 0;padding-left:18px;font-size:13.5px}
.tblwrap{overflow-x:auto;background:var(--panel);border:1px solid var(--line);border-radius:10px}
table{border-collapse:collapse;width:100%;font-size:13px}
th,td{padding:7px 11px;text-align:right;white-space:nowrap;border-bottom:1px solid var(--line)}
th:first-child,td:first-child{text-align:left}
th{font-weight:600;color:var(--muted);font-size:11.5px;text-transform:uppercase;
   letter-spacing:.05em;position:sticky;top:0;background:var(--panel)}
tr:last-child td{border-bottom:none}
td.g{color:var(--ok)} td.a{color:var(--warn)} td.r{color:var(--bad)}
.heat td{padding:3px 5px;text-align:center;font-size:11px;border:1px solid var(--bg)}
.heat td.h-green{background:var(--ok-bg)} .heat td.h-amber{background:var(--warn-bg)}
.heat td.h-red{background:var(--bad-bg);color:var(--bad);font-weight:700}
.heat td.h-na{background:transparent;color:var(--muted)}
.note{color:var(--muted);font-size:12.5px;margin:8px 2px 0}
footer{margin-top:44px;color:var(--muted);font-size:12px;border-top:1px solid var(--line);padding-top:14px}
"""


def _e(x) -> str:
    return html.escape(str(x))


def _spark_svg(values, direction="higher_better", w=110, h=22) -> str:
    vals = [v for v in values if v is not None]
    if len(vals) < 2:
        return ""
    lo, hi = min(vals), max(vals)
    rng = (hi - lo) or 1
    pts = []
    n = len(values)
    for i, v in enumerate(values):
        if v is None:
            continue
        x = i / (n - 1) * (w - 2) + 1
        y = h - 2 - (v - lo) / rng * (h - 4)
        pts.append(f"{x:.1f},{y:.1f}")
    return (f'<svg width="{w}" height="{h}" viewBox="0 0 {w} {h}" aria-hidden="true">'
            f'<polyline points="{" ".join(pts)}" fill="none" stroke="currentColor" '
            f'stroke-width="1.4" stroke-opacity=".65" stroke-linejoin="round"/></svg>')


def build_dashboard(engine: MetricEngine, biz_date: date, trend_days: int = 14) -> str:
    p = engine.profile
    card = engine.scorecard(biz_date)
    by_id = {r["id"]: r for r in card}
    fams = mx.family_summary(engine, biz_date)
    pri = mx.priority(engine, biz_date)
    ch_cn = {c["id"]: c["name_cn"] for c in p["channels"]}
    dp_cn = {d["id"]: d["name_cn"] for d in p["dayparts"]}

    H: list[str] = []
    A = H.append
    A("<title>门店日运营看板</title>")
    A(f"<style>{CSS}</style>")
    A('<div class="wrap">')
    A(f'<header><h1>{_e(p["store_name"])} · 日运营看板</h1>'
      f'<span class="sub">营业日 {biz_date}（周{"一二三四五六日"[biz_date.weekday()]}） · '
      f'{_e(p["store_id"])} · {_e(p["format"])}</span></header>')

    # ── 第 0 屏：四个数 + 六盏灯
    A('<h2>第 0 屏 · 3 秒钟：今天是红是绿</h2>')
    A('<div class="grid kpis">')
    for mid, lab in (("net_sales", "净销售"), ("sales_vs_forecast", "销售达成"),
                     ("guest_count", "交易笔数"), ("avg_check", "客单价"),
                     ("food_cost_pct", "食品成本率"), ("labor_cost_pct", "人工成本率"),
                     ("splh", "人时销售 SPLH"), ("speed_hit_rate", "速度达标率")):
        r = by_id[mid]
        trend = [v for _d, v in engine.trend(mid, biz_date, trend_days)]
        A(f'<div class="card"><div class="lab">{lab}</div>'
          f'<div class="val" style="color:{COLORS[r["status"]]}">{fmt(r["value"], r["unit"])}</div>'
          f'<div class="foot">目标 {fmt(r["target"], r["unit"])} · '
          f'{("%+.1f%%" % r["gap_pct"]) if r["gap_pct"] is not None else "—"}'
          f'<span style="float:right;color:var(--muted)">{_spark_svg(trend, w=64, h=16)}</span></div></div>')
    A('</div>')

    A('<div class="grid fams" style="margin-top:12px">')
    for f in fams:
        c = f["counts"]
        h = f["headline"]
        A(f'<div class="card fam {f["status"]}"><div class="lab">{f["family_cn"]}'
          f'<span style="float:right">'
          f'<span class="pill green">{c["green"]}</span> '
          f'<span class="pill amber">{c["amber"]}</span> '
          f'<span class="pill red">{c["red"]}</span></span></div>'
          f'<div style="font-size:13.5px;margin-top:4px">最该看：{_e(h["name_cn"])} '
          f'<b style="color:{COLORS[h["status"]]}">{fmt(h["value"], h["unit"])}</b></div>'
          f'<div class="foot">目标 {fmt(h["target"], h["unit"])} · {mx.OWNER_LABEL[h["owner"]]}</div></div>')
    A('</div>')

    # ── 第 1 屏：今天的三件事
    A('<h2>第 1 屏 · 30 秒：今天做什么（优先级矩阵 M4）</h2>')
    if not pri["top3"]:
        A('<div class="card">今天全绿。把时间用在训练和顾客身上。</div>')
    for i, it in enumerate(pri["top3"], 1):
        cells = mx.worst_cells(engine, biz_date, it["id"], top=3)
        A(f'<div class="todo {it["status"]}">')
        A(f'<h3>{i}. {_e(it["name_cn"])} '
          f'<span style="color:{COLORS[it["status"]]}">{fmt(it["value"], it["unit"])}</span> '
          f'<span class="sub">/ 目标 {fmt(it["target"], it["unit"])}</span></h3>')
        A(f'<div class="meta">{_e(it["quadrant"])} · {mx.OWNER_LABEL[it["owner"]]} · '
          f'行动卡 {_e(it["playbook"])} · 行动分 {it["score"]}'
          f'（超线 {it["severity"]:.2f} 个带宽 × 影响 {it["impact"]} × 可控 {it["controllability"]}'
          f'{" × 红灯 2×" if it["status"] == "red" else ""}）</div>')
        A(f'<div style="font-size:13px;color:var(--muted)">{_e(it["quadrant_note"])}</div>')
        if cells:
            A('<ul>')
            for c in cells:
                A(f'<li>{_e(c["dim_cn"])} <b>{_e(dp_cn.get(c["dim_value"], ch_cn.get(c["dim_value"], c["dim_value"])))}</b>'
                  f' → <span style="color:{COLORS[c["status"]]}">{fmt(c["value"], it["unit"])}</span>'
                  f'　拖累度 {c["contribution"]:.2f}，量占比 {c["weight"]*100:.0f}%</li>')
            A('</ul>')
        if engine.defs[it["id"]].get("grain") == "15min":
            for w in mx.problem_windows(engine, biz_date, it["id"])[:1]:
                A(f'<div class="note">连续告警窗口：<b>{w["from"]} → {w["to"]}</b>'
                  f'（{w["slots"]} 个 15 分钟，均值 {fmt(w["avg"], it["unit"])}）</div>')
        A('</div>')

    # ── 第 2 屏：时段矩阵
    A('<h2>第 2 屏 · 3 分钟：为什么（时段矩阵 M2）</h2>')
    grid_ids = ["sales_vs_forecast", "speed_hit_rate", "dt_total_time", "kds_prep_time", "splh"]
    g = mx.interval_grid(engine, biz_date, grid_ids)
    ivs = g["intervals"]
    A('<div class="tblwrap"><table class="heat"><thead><tr><th>指标</th>')
    for iv in ivs:
        A(f'<th style="text-align:center;font-size:9.5px;padding:3px 2px">'
          f'{iv[3:] == "00" and iv[:2] or ""}</th>')
    A('</tr></thead><tbody>')
    for mid in grid_ids:
        m = engine.defs[mid]
        A(f'<tr><td style="text-align:left;white-space:nowrap">{_e(m["name_cn"])}</td>')
        for c in g["grid"][mid]:
            title = f'{c["interval15"]} {fmt(c["value"], m["unit"])}'
            A(f'<td class="h-{c["status"]}" title="{_e(title)}">'
              f'{"" if c["status"] == "green" else ("!" if c["status"] == "red" else "·")}</td>')
        A('</tr>')
    A('</tbody></table></div>')
    A('<p class="note">空白=达标　·=黄灯　<b>!</b>=红灯　（鼠标悬停看具体数值）'
      '单点飘红是噪声，连成一片才是问题。</p>')

    # ── 交叉矩阵
    A('<h2>第 2 屏 · 交叉矩阵 M3：时段 × 渠道</h2>')
    for mid in ("speed_hit_rate", "net_sales"):
        c = mx.cross(engine, biz_date, mid, "daypart", "channel")
        A(f'<p class="note" style="margin-top:14px"><b>{_e(c["metric_cn"])}</b></p>')
        A('<div class="tblwrap"><table><thead><tr><th>时段</th>')
        for b in c["b_values"]:
            A(f'<th>{_e(ch_cn.get(b, b))}</th>')
        A('</tr></thead><tbody>')
        for row in c["cells"]:
            A(f'<tr><td>{_e(dp_cn.get(row[0]["a"], row[0]["a"]))}</td>')
            for cell in row:
                cls = {"green": "g", "amber": "a", "red": "r", "na": ""}[cell["status"]]
                A(f'<td class="{cls}">{fmt(cell["value"], c["unit"])}</td>')
            A('</tr>')
        A('</tbody></table></div>')
        if c["worst"]:
            w = c["worst"][0]
            A(f'<p class="note">最拖后腿的格子：<b>{_e(dp_cn.get(w["a"], w["a"]))} × '
              f'{_e(ch_cn.get(w["b"], w["b"]))}</b> = {fmt(w["value"], c["unit"])}'
              f'（拖累度 {w["contribution"]:.1f}）</p>')

    # ── 第 3 屏：完整记分卡
    A('<h2>第 3 屏 · 复盘：完整记分卡</h2>')
    A('<div class="tblwrap"><table><thead><tr>'
      '<th>家族</th><th>指标</th><th>实际</th><th>目标</th><th>达标</th>'
      '<th>超线带宽</th><th>状态</th><th>责任人</th><th>行动卡</th>'
      f'<th>近 {trend_days} 天</th></tr></thead><tbody>')
    label = {"green": "达标", "amber": "关注", "red": "告警", "na": "无数据"}
    for r in card:
        trend = [v for _d, v in engine.trend(r["id"], biz_date, trend_days)]
        cls = {"green": "g", "amber": "a", "red": "r", "na": ""}[r["status"]]
        A(f'<tr><td>{mx.FAMILY_LABEL[r["family"]]}</td><td>{_e(r["name_cn"])}</td>'
          f'<td class="{cls}"><b>{fmt(r["value"], r["unit"])}</b></td>'
          f'<td>{fmt(r["target"], r["unit"])}</td>'
          f'<td>{("%+.1f%%" % r["gap_pct"]) if r["gap_pct"] is not None else "—"}</td>'
          f'<td>{r["severity"]:.2f}</td>'
          f'<td><span class="pill {r["status"]}">{label[r["status"]]}</span></td>'
          f'<td>{mx.OWNER_LABEL[r["owner"]]}</td><td>{_e(r["playbook"] or "")}</td>'
          f'<td style="color:{COLORS[r["status"]]}">{_spark_svg(trend)}</td></tr>')
    A('</tbody></table></div>')

    # ── 责任矩阵
    own = mx.ownership(engine)
    A('<h2>责任矩阵 M5：谁 × 什么节奏 × 盯什么</h2>')
    A('<div class="tblwrap"><table><thead><tr><th>责任人</th>'
      + "".join(f'<th style="text-align:left">{mx.CADENCE_LABEL[c]}</th>' for c in own["cadences"])
      + '</tr></thead><tbody>')
    for o in own["owners"]:
        if o not in own["grid"]:
            continue
        A(f'<tr><td>{mx.OWNER_LABEL[o]}</td>')
        for c in own["cadences"]:
            items = own["grid"][o].get(c, [])
            A(f'<td style="text-align:left;white-space:normal;font-size:12.5px">'
              f'{_e("、".join(items)) if items else "—"}</td>')
        A('</tr>')
    A('</tbody></table></div>')

    A(f'<footer>由 <code>mcd_store_ops</code> 生成 · 指标口径见 <code>config/metrics.json</code> · '
      f'数据为模拟数据，用于演示报表与矩阵的结构</footer>')
    A('</div>')
    return "\n".join(H)
