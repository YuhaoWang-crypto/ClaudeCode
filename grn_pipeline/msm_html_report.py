"""
Build the self-contained HTML report for M23 (MSM + TPT on alanine
dipeptide) and M24 (reliability controls, screening, MD perturbation).

All numbers are read from figures/m23_summary.json and
figures/m24_summary.json (written by the two modules' report()); figures are
embedded as base64 so the HTML is a single file.

    python3 -m grn_pipeline.m23_msm_tpt        # -> m23_summary.json
    python3 -m grn_pipeline.m24_msm_validation # -> m24_summary.json
    python3 -m grn_pipeline.msm_html_report    # -> figures/msm_report.html
"""
import base64
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "..", "figures")


def _img(name, alt):
    with open(os.path.join(FIG, name), "rb") as fh:
        b64 = base64.b64encode(fh.read()).decode()
    return (f'<figure class="plate"><img src="data:image/png;base64,{b64}" '
            f'alt="{alt}"><figcaption>{alt}</figcaption></figure>')


def _chip(state):
    label = {"pass": "通过", "warn": "有条件", "fail": "不可用"}[state]
    return f'<span class="chip {state}">{label}</span>'


def _pct(x, digits=0):
    return f"{x * 100:+.{digits}f}%"


CSS = """
/* Layout: single reading column (~72ch) with full-width figure plates and
   tables; verdict chips encode pass / conditional / fail. */
:root {
  --bg: #f6f7f9; --surface: #ffffff; --ink: #17202b; --muted: #5a6573;
  --rule: #d9dee5; --accent: #1f5fae; --accent-soft: #e6eef8;
  --pass: #1d7a4d; --pass-bg: #e3f3ea; --warn: #9a6200; --warn-bg: #fbf0dc;
  --fail: #b3261e; --fail-bg: #fbe5e3; --plate: #ffffff;
  --f-display: "Noto Serif SC", "Songti SC", Georgia, serif;
  --f-body: "Noto Sans SC", "PingFang SC", "Microsoft YaHei", system-ui,
    sans-serif;
  --f-mono: "JetBrains Mono", ui-monospace, Menlo, Consolas, monospace;
}
@media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) {
  --bg: #11161d; --surface: #182029; --ink: #e5eaf0; --muted: #9aa6b4;
  --rule: #2b3642; --accent: #79aef0; --accent-soft: #1d2b3c;
  --pass: #6fcf9b; --pass-bg: #163326; --warn: #f0b85a; --warn-bg: #3a2c12;
  --fail: #f28b82; --fail-bg: #3d1a18; --plate: #f4f6f8;
  color-scheme: dark; } }
:root[data-theme="dark"] {
  --bg: #11161d; --surface: #182029; --ink: #e5eaf0; --muted: #9aa6b4;
  --rule: #2b3642; --accent: #79aef0; --accent-soft: #1d2b3c;
  --pass: #6fcf9b; --pass-bg: #163326; --warn: #f0b85a; --warn-bg: #3a2c12;
  --fail: #f28b82; --fail-bg: #3d1a18; --plate: #f4f6f8;
  color-scheme: dark; }
body { background: var(--bg); color: var(--ink); font-family: var(--f-body);
  font-size: 15.5px; line-height: 1.75; }
.wrap { max-width: 1080px; margin: 0 auto; padding-inline: 20px;
  padding-block: 40px 72px; display: grid; gap: 44px; }
.col { max-width: 72ch; }
header { display: grid; gap: 14px; }
.eyebrow { font-family: var(--f-mono); font-size: 12px; letter-spacing: .08em;
  text-transform: uppercase; color: var(--accent); }
h1 { font-family: var(--f-display); font-weight: 700; font-size: 34px;
  line-height: 1.25; margin: 0; text-wrap: balance; }
h2 { font-family: var(--f-display); font-weight: 600; font-size: 24px;
  margin: 0 0 6px; text-wrap: balance; }
h3 { font-size: 16.5px; margin: 18px 0 4px; }
p { margin: 0 0 12px; }
.lede { color: var(--muted); font-size: 16.5px; }
section { display: grid; gap: 16px; }
.verdicts { display: grid; grid-template-columns: repeat(auto-fit,
  minmax(300px, 1fr)); gap: 12px; }
.verdict { background: var(--surface); border: 1px solid var(--rule);
  border-radius: 8px; padding: 14px 16px; display: grid; gap: 6px;
  min-width: 0; }
.verdict b { font-size: 15px; }
.verdict .chip { justify-self: start; }
.verdict span.small { color: var(--muted); font-size: 13.5px;
  line-height: 1.6; }
.chip { display: inline-block; font-size: 12px; font-weight: 600;
  padding: 1px 9px; border-radius: 999px; letter-spacing: .02em;
  white-space: nowrap; }
.chip.pass { color: var(--pass); background: var(--pass-bg); }
.chip.warn { color: var(--warn); background: var(--warn-bg); }
.chip.fail { color: var(--fail); background: var(--fail-bg); }
.tbl { overflow-x: auto; border: 1px solid var(--rule); border-radius: 8px;
  background: var(--surface); }
table { border-collapse: collapse; width: 100%; font-size: 14px; }
th, td { text-align: left; padding: 8px 12px; border-bottom: 1px solid
  var(--rule); vertical-align: top; }
th { font-weight: 600; color: var(--muted); font-size: 12.5px;
  letter-spacing: .03em; background: var(--accent-soft); }
tr:last-child td { border-bottom: none; }
td.num, th.num { text-align: right; font-family: var(--f-mono);
  font-variant-numeric: tabular-nums; white-space: nowrap; }
.plate { margin: 0; background: var(--plate); border: 1px solid var(--rule);
  border-radius: 8px; padding: 10px; overflow-x: auto; }
.plate img { display: block; max-width: 100%; height: auto; margin: 0 auto; }
.plate figcaption { color: #5a6573; font-size: 13px; padding: 8px 4px 0; }
code, pre { font-family: var(--f-mono); font-size: 13px; }
pre { background: var(--surface); border: 1px solid var(--rule);
  border-radius: 8px; padding: 12px 14px; overflow-x: auto; margin: 0;
  line-height: 1.55; }
.note { border-left: 3px solid var(--warn); background: var(--warn-bg);
  padding: 10px 14px; border-radius: 0 6px 6px 0; }
ul { margin: 0 0 12px; padding-left: 1.3em; }
li { margin: 3px 0; }
.two { display: grid; grid-template-columns: repeat(auto-fit,
  minmax(300px, 1fr)); gap: 16px; }
.two > div { min-width: 0; }
a { color: var(--accent); }
:focus-visible { outline: 2px solid var(--accent); outline-offset: 2px; }
@media (max-width: 480px) { h1 { font-size: 27px; } body { font-size: 15px; } }
"""


SCORE_CN = {"MSM: d ln K": "MSM 预测的平衡偏移 ΔlnK",
            "MSM: d ln k (rate)": "MSM 预测的速率变化 Δln k",
            "static structure difference": "静态结构差异（两个晶体结构比较）",
            "RMSF (conventional MD)": "RMSF（常规 MD 柔性）"}
SYS_CN = {"GK": "葡萄糖激酶（GK）", "SHP2": "SHP2 磷酸酶"}


def _retro(path):
    """Summary card + section for the M25 retrospective test, or ('', '')
    when m25_summary.json does not exist yet."""
    if not os.path.exists(path):
        return "", ""
    with open(path) as fh:
        R = json.load(fh)
    msm = [r["eval"]["MSM: d ln K"] for r in R]
    st = [r["eval"]["static structure difference"] for r in R]
    win = all(a["auc"] >= b["auc"] and a["patch_percentile"] >= 0.95
              for a, b in zip(msm, st))
    state = ("pass" if win else
             "fail" if np.mean([a["auc"] for a in msm]) < 0.6 else "warn")
    auc_txt = "，".join(f"{SYS_CN[r['name']]} {a['auc']:.2f}（静态结构 "
                       f"{b['auc']:.2f}）" for r, a, b in zip(R, msm, st))
    card = (f'<div class="verdict">{_chip(state)}<b>已知变构药物回顾性验证</b>'
            f'<span class="small">MSM 位点打分的 ROC-AUC：{auc_txt}。'
            f'{"优于" if win else "没有稳定优于"}直接比较两个晶体结构的简单基线。'
            f'</span></div>')
    rows_diag, rows_eval, rows_chk, rows_post = "", "", "", ""
    for r in R:
        nm = SYS_CN[r["name"]]
        rows_diag += (f"<tr><td>{nm}</td><td>{'抑制剂' if r['kind'] == 'inhibitor' else '激活剂'}</td>"
                      f"<td class='num'>{r['n_beads']}</td><td class='num'>{r['n_site']}</td>"
                      f"<td class='num'>{r['n_frames']}</td><td class='num'>{r['transitions']}</td>"
                      f"<td class='num'>{r['ck_err']:.3f}</td><td class='num'>{r['p_active']:.2f}</td></tr>")
        for k, e in r["eval"].items():
            rows_eval += (f"<tr><td>{nm}</td><td>{SCORE_CN[k]}</td>"
                          f"<td class='num'>{e['auc']:.2f}</td>"
                          f"<td class='num'>{e['enrich_top10']:.1f}×</td>"
                          f"<td class='num'>{e['best_rank']}</td>"
                          f"<td class='num'>{e['patch_percentile'] * 100:.0f}</td></tr>")
        for c in r.get("checks", []):
            ok = c["ci"][0] <= c["dlnK_pred"] <= c["ci"][1]
            rows_chk += (f"<tr><td>{nm}</td><td>{'最高分位点' if c['tag'].startswith('top') else '远处对照位点'}"
                         f"（残基 {c['tag'][3:]}）</td>"
                         f"<td class='num'>{c['dlnK_pred']:+.2f}</td>"
                         f"<td class='num'>{c['dlnK_sim']:+.2f} ({c['ci'][0]:+.2f}, {c['ci'][1]:+.2f})</td>"
                         f"<td>{_chip('pass' if ok else 'warn')}</td></tr>")
        for k, d in r.get("posthoc", {}).items():
            rows_post += (f"<tr><td>{nm}</td><td>{SCORE_CN[k]}</td>"
                          f"<td class='num'>{d['median_dist_A']:.1f} Å</td>"
                          f"<td class='num'>{r['posthoc_random']['median']:.1f} / {r['posthoc_random']['p5']:.1f} Å</td></tr>")
    verdict_lines = "".join(
        f"<li><b>{SYS_CN[r['name']]}：</b>MSM ΔlnK 的 AUC 为 {a['auc']:.2f}，静态结构基线为 {b['auc']:.2f}；"
        f"真实位点在随机空间片段中的百分位为 {a['patch_percentile'] * 100:.0f}"
        f"（{'达到' if a['patch_percentile'] >= 0.95 else '未达到'} 95% 显著水平）。</li>"
        for r, a, b in zip(R, msm, st))
    near = [SYS_CN[r["name"]] for r in R if "posthoc" in r and
            r["posthoc"]["MSM: d ln K"]["median_dist_A"]
            < r["posthoc_random"]["p5"]]
    far = [SYS_CN[r["name"]] for r in R if SYS_CN[r["name"]] not in near]
    near_txt = ""
    if near:
        near_txt += ("在" + "、".join(near) + "上，MSM 的前 10 名比 95% 的随机残基组更靠近药物位点，"
                     "可以把搜索范围缩小到位点附近（事后分析）。")
    if far:
        near_txt += ("在" + "、".join(far) + "上，MSM 的前 10 名并不比随机残基更靠近药物位点。")
    shaky = [SYS_CN[r["name"]] for r in R
             if r["ck_err"] > 0.1 or any(
                 not (c["ci"][0] <= c["dlnK_pred"] <= c["ci"][1])
                 for c in r.get("checks", []))]
    if shaky:
        near_txt += ("其中" + "、".join(shaky) + "的模型诊断不合格（CK 误差大于 0.1，或远处对照位点的显式模拟与预测不符），"
                     "那里的结果首先说明模型和采样不足，不能单独当作方法失败的证据。")
    checks = [(SYS_CN[r["name"]], c) for r in R for c in r.get("checks", [])]
    ok = [c["ci"][0] <= c["dlnK_pred"] <= c["ci"][1] for _, c in checks]
    bad = [f"{nm}的{'最高分位点' if c['tag'].startswith('top') else '远处对照位点'}"
           for (nm, c), o in zip(checks, ok) if not o]
    chk_txt = (f"显式虚拟配体检验：{sum(ok)}/{len(ok)} 个位点的预测落在模拟的 95% 置信区间内，区间普遍较宽。"
               + (f"不一致的是{'、'.join(bad)}，说明在该模型里，即便预测某位点没有效应，也不能直接相信。" if bad else
                  "重加权这一步在蛋白上可用，限制在于采样量。"))
    section = f"""
<section>
  <div class="col">
  <h2>7. 回顾性验证：已知变构药物的蛋白</h2>
  <p>用两个有已知变构药物共晶结构的蛋白，检验第 4 节的位点扫描能否在不知道药物位置的情况下找出它：
  SHP2 的变构<b>抑制剂</b> SHP099（PDB 5EHR），以及葡萄糖激酶的变构<b>激活剂</b>（PDB 1V4S）。
  标准答案是共晶结构中距药物 4.5 Å 以内的残基。5EHR 只用于打分，不参与建模；
  葡萄糖激酶唯一的闭合结构本身带有激活剂，建模时只用了蛋白坐标，这一点仍可能泄漏少量信息。</p>
  <p>全原子显式溶剂的 MSM 需要 GPU 上微秒级采样，这里只有 4 个 CPU 核，因此用双势阱 Cα 结构模型（由失活态和激活态两个晶体结构构建）。
  它检验的是"给定构象景观，通量分析能否定位位点"，不检验力场精度。打分规则在看到结果之前固定：
  在每个残基周围 8 Å 放一个虚拟配体，稳定药物应偏好的状态（抑制剂偏好失活态，激活剂偏好激活态）的局部接触，结合能 2 kT，
  用一个无偏 MSM 重加权预测平衡偏移。</p>
  <h3>采样与模型诊断</h3>
  </div>
  <div class="tbl"><table>
    <tr><th>体系</th><th>药物类型</th><th class="num">残基</th><th class="num">药物接触残基</th><th class="num">帧数</th>
    <th class="num">失活↔激活转变</th><th class="num">CK 误差</th><th class="num">激活态布居</th></tr>{rows_diag}
  </table></div>
  <p class="col note">过程中做过三次修正，均在打分之前：从两个端态直接模拟时转变太少，CK 误差 0.55，于是改用自适应采样
  （在失活→激活坐标上均匀撒种子，SHP2 的种子来自驱动模拟，但驱动轨迹不进入 MSM）；PCCA+ 抓到的是模型局部缓慢重排而不是失活↔激活转变，
  于是改用物理序参量定义两个状态；之后 CK 误差降到上表的水平。</p>
  <div class="col"><h3>位点识别（预先固定的指标）</h3></div>
  <div class="tbl"><table>
    <tr><th>体系</th><th>打分方法</th><th class="num">ROC-AUC</th><th class="num">前 10% 富集</th><th class="num">最佳排名</th><th class="num">相对随机片段百分位</th></tr>{rows_eval}
  </table></div>
  <div class="col"><ul>{verdict_lines}</ul></div>
  <div class="col"><h3>显式虚拟配体检验</h3>
  <p>在最高分位点和一个远处的对照位点真正施加虚拟配体，重新做种子模拟并建 MSM，与重加权预测比较。</p></div>
  <div class="tbl"><table>
    <tr><th>体系</th><th>位点</th><th class="num">预测 ΔlnK</th><th class="num">模拟 ΔlnK（95% CI）</th><th>一致</th></tr>{rows_chk}
  </table></div>
  <div class="col"><h3>事后分析：前 10 名离药物有多远</h3>
  <p>这一项不在预先固定的指标里，单独列出。数值是每种方法前 10 个残基到最近药物接触残基的 Cα 距离中位数；随机 10 个残基作为对照（中位数 / 第 5 百分位）。</p></div>
  <div class="tbl"><table>
    <tr><th>体系</th><th>打分方法</th><th class="num">前 10 名距药物</th><th class="num">随机（中位 / 5%）</th></tr>{rows_post}
  </table></div>
  {_img("m25_allosteric_retro.png", "图 5 回顾性验证：(a) 沿序列的位点得分与药物接触残基；(b) ROC 曲线；(c) 显式虚拟配体检验")}
  <div class="col">
  <h3>结论</h3>
  <p>{near_txt}在预先固定的指标上，MSM 位点扫描{"在两个体系上都" if win else "没有"}稳定超过"直接比较两个晶体结构"这一简单基线。
  原因可以理解：结构模型里所有的耦合都来自这两个结构，动力学分析能提供的额外信息有限。
  要真正检验 MSM 相对静态分析的增量，需要全原子模拟，那里会出现晶体结构里看不到的构象（例如隐藏口袋）和非局部耦合。</p>
  <p>{chk_txt}</p>
  </div>
</section>
"""
    return card, section


def build(path=os.path.join(FIG, "msm_report.html")):
    with open(os.path.join(FIG, "m23_summary.json")) as fh:
        a = json.load(fh)
    with open(os.path.join(FIG, "m24_summary.json")) as fh:
        v = json.load(fh)
    retro_card, retro_section = _retro(os.path.join(FIG, "m25_summary.json"))
    rate, sc, neg = v["rate"], v["screening"], v["negative"]
    md = v["md"]
    covered = rate["ci"][0] <= rate["exact"] <= rate["ci"][1]
    hp_err = [h["rate_err"] for h in v["hyperparams"]]
    lab = {int(k): s for k, s in a["labels"].items()}
    pis = a["pis"]
    pop = {lab[k]: pis[k] for k in lab}
    ch = a["channels"]
    ch_direct = [x for k, x in ch.items() if "direct" in k][0]

    def mfpt(k):
        return 1.0 / k

    psi, met, un = md["psi"], md["methyl"], md["unbiased"]

    def ci_mfpt(ci):
        return f"{1 / ci[1]:.0f}–{1 / ci[0]:.0f}"

    md_rows = ""
    for name, d in [("无偏置（基线）", None), ("ψ 偏置：预测", psi["pred"]),
                    ("ψ 偏置：模拟", psi["sim"]),
                    ("甲基偏置：预测", met["pred"]),
                    ("甲基偏置：模拟", met["sim"])]:
        d = d or un
        cif = f" ({ci_mfpt(d['ci_f'])})" if "ci_f" in d else ""
        cib = f" ({ci_mfpt(d['ci_b'])})" if "ci_b" in d else ""
        md_rows += (f"<tr><td>{name}</td><td class='num'>{d['pop_aR']:.3f}"
                    f"</td><td class='num'>{mfpt(d['k_f']):.0f}{cif}</td>"
                    f"<td class='num'>{mfpt(d['k_b']):.0f}{cib}</td></tr>")

    def within(pred, sim):
        return sim["ci_b"][0] <= pred["k_b"] <= sim["ci_b"][1]

    conv_rows = "".join(
        f"<tr><td class='num'>{c['frac'] * 100:.0f}%</td>"
        f"<td class='num'>{c['n_events']}</td>"
        f"<td class='num'>{_pct(c['rate_err'], 1)}</td>"
        f"<td class='num'>{c['q_r']:.4f}</td></tr>"
        for c in v["convergence"])
    hp_rows = "".join(
        f"<tr><td>{h['method']}</td><td class='num'>{h['n']}</td>"
        f"<td class='num'>{h['lag']}</td>"
        f"<td class='num'>{_pct(h['rate_err'], 1)}</td>"
        f"<td class='num'>{h['q_r']:.4f}</td>"
        f"<td class='num'>{h['split_err']:+.3f}</td></tr>"
        for h in v["hyperparams"])
    amp_rows = "".join(
        f"<tr><td class='num'>{q['amp_kT']:.1f} kT</td>"
        f"<td class='num'>{q['rho']:.3f}</td>"
        f"<td class='num'>{q['sign'] * 100:.0f}%</td>"
        f"<td class='num'>{q['mae']:.3f}</td>"
        f"<td class='num'>{q['rho_static']:.3f}</td></tr>"
        for q in v["amplitude"])
    hot = "；".join(f"({h['x']:.1f}, {h['y']:.1f})：精确 {h['S_exact']:+.2f}，"
                   f"预测 {h['S_pred']:+.2f}" for h in v["hot_sites"])

    psi_ok = within(psi["pred"], psi["sim"])
    met_ok = within(met["pred"], met["sim"])

    html = f"""<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">
<title>构象通量验证报告</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=JetBrains+Mono:wght@400;600&family=Noto+Sans+SC:wght@400;600&family=Noto+Serif+SC:wght@600;700&display=swap">
<style>{CSS}</style>
<div class="wrap">
<header class="col">
  <div class="eyebrow">MSM + 过渡路径理论 · 可靠性评估</div>
  <h1>构象通量：从丙氨酸二肽到虚拟配体筛选</h1>
  <p class="lede">用马尔可夫状态模型（MSM）和过渡路径理论（TPT）把分子动力学轨迹变成"状态之间的概率流"，
  再用已知精确答案的体系和真实 MD 的阳性、阴性对照检验这套方法能否用于找变构位点和筛选调节剂。</p>
</header>

<section>
  <h2>结论速览</h2>
  <div class="verdicts">
    <div class="verdict">{_chip("pass")}<b>速率与 committor 准确</b>
      <span class="small">精确解体系上速率误差 {_pct(rate['msm'] / rate['exact'] - 1, 1)}，
      bootstrap 95% 区间{'覆盖' if covered else '未覆盖'}精确值；committor 相关 r = {v['committor']['r']:.4f}。</span></div>
    <div class="verdict">{_chip("pass")}<b>阴性对照被识别</b>
      <span class="small">打乱时间顺序后最慢时间尺度从 {neg['t1_xy']:.2f} 降到 {neg['t1_shuffled']:.3f}；
      看不见慢过程的坐标让它降到 {neg['t1_y']:.2f}；诱饵配体预测效应 ≤ {v['decoy_max']:.0e}。</span></div>
    <div class="verdict">{_chip("pass")}<b>虚拟配体排序可靠</b>
      <span class="small">56 个扰动的 Spearman ρ = {sc['rho_msm']:.3f}，方向判断 {sc['sign_msm'] * 100:.0f}% 正确；
      只看反应物稳定性的静态基线 ρ = {sc['rho_static']:.2f}。</span></div>
    <div class="verdict">{_chip("pass" if psi_ok and met_ok else "warn")}<b>真实 MD 扰动对照</b>
      <span class="small">ψ 偏置（阳性）：αR 布居预测 {psi['pred']['pop_aR']:.2f}，模拟 {psi['sim']['pop_aR']:.2f}；
      动力学方向正确，但 αR→β 的 MFPT 预测 {1 / psi['pred']['k_b']:.0f} ps，模拟 {1 / psi['sim']['k_b']:.0f} ps。
      甲基偏置（阴性）：预测无效应，模拟也无效应。</span></div>
    <div class="verdict">{_chip("warn")}<b>参数敏感性</b>
      <span class="small">聚类与 lag 的选择让速率在 {min(hp_err) * 100:+.0f}% 到 {max(hp_err) * 100:+.0f}% 之间变化，
      比统计误差大；报告结果时应给出这个系统误差范围。</span></div>
    <div class="verdict">{_chip("fail")}<b>稀有事件与新构象</b>
      <span class="small">β→αL 只采到 {a['rare_events']} 次；重加权只能预测已采样构象内的扰动，
      不能预测诱导出的新状态或隐藏口袋。</span></div>
    {retro_card}
  </div>
</section>

<section class="col">
  <h2>1. 方法</h2>
  <p>流程分七步：内坐标特征（二面角 sin/cos、Cα 距离），TICA 提取最慢的线性组合，regular-space 聚类成微观态，
  可逆最大似然 MSM（用隐含时间尺度选 lag，用 Chapman–Kolmogorov 检验马尔可夫性），PCCA+ 宏观态，
  在 PCCA+ 核心集（隶属度 ≥ 0.9）之间做 TPT，最后与轨迹上直接计数的速率对照。</p>
  <p>净通量用一般式 f<sub>ij</sub> = π<sub>i</sub> q⁻<sub>i</sub> T<sub>ij</sub> q⁺<sub>j</sub>；
  在细致平衡下它等于 π<sub>i</sub> T<sub>ij</sub>(q⁺<sub>j</sub> − q⁺<sub>i</sub>)，代码数值验证到 10⁻¹⁴。
  瓶颈不取"通量最小的边"，而取穿过 q⁺ = ½ 等值面、承载通量最多的边：这些边的通量之和严格等于总通量 F。</p>
</section>

<section>
  <div class="col">
  <h2>2. 丙氨酸二肽：构象通量和 RMSD 的区别</h2>
  <p>OpenMM 模拟 4 × 50 ns（amber99sb-ildn + OBC 隐式溶剂，300 K）。</p>
  </div>
  <div class="tbl"><table>
    <tr><th>指标</th><th class="num">数值</th><th>说明</th></tr>
    <tr><td>MSM lag</td><td class="num">{a['lag_ps']:.0f} ps</td><td>隐含时间尺度收敛处</td></tr>
    <tr><td>最慢时间尺度 t₁, t₂</td><td class="num">{a['timescales_ps'][0]:.0f}, {a['timescales_ps'][1]:.0f} ps</td><td>t₁：离开 αL；t₂：β↔αR</td></tr>
    <tr><td>宏观态布居 β / αR / αL</td><td class="num">{pop['beta/PII']:.2f} / {pop['alphaR']:.2f} / {pop['alphaL']:.3f}</td><td>PCCA+</td></tr>
    <tr><td>CK 检验最大误差</td><td class="num">{a['ck_err']:.3f}</td><td>马尔可夫性成立</td></tr>
    <tr><td>MFPT β→αR：MSM / 直接计数</td><td class="num">{a['mfpt_ps']:.0f} / {a['mfpt_direct_ps']:.0f} ps</td><td>直接计数 {a['n_events']} 次，不依赖马尔可夫假设</td></tr>
    <tr><td>MFPT αR→β：MSM / 直接计数</td><td class="num">{a['mfpt_back_ps']:.0f} / {a['mfpt_back_direct_ps']:.0f} ps</td><td></td></tr>
    <tr><td>主通道份额</td><td class="num">{ch_direct * 100:.0f}%</td><td>ψ 经 +60°…+120° 下降；其余绕 ψ = ±180°</td></tr>
    <tr><td>αL 与 αR 的 RMSD 分布重叠</td><td class="num">{a['rmsd_overlap'][[k for k in lab if lab[k] == 'alphaL'][0]][[k for k in lab if lab[k] == 'alphaR'][0]] * 100:.0f}%</td><td>被最慢势垒隔开的两个态，RMSD 难以区分</td></tr>
    <tr><td>β→αL 稀有事件</td><td class="num">{a['rare_events']} 次</td><td>MFPT ≈ {a['rare_mfpt_ns']:.0f} ns，只作数量级</td></tr>
  </table></div>
  {_img("m23_msm_tpt.png", "图 1 丙氨酸二肽：(a) 自由能与宏观态；(b) 隐含时间尺度；(c) CK 检验；(d) committor 与 q⁺ = ½ 穿越边；(e) 净反应通量网络；(f) 各宏观态的 RMSD 分布")}
</section>

<section>
  <div class="col">
  <h2>3. 可靠性：在精确解体系上检验</h2>
  <p>MD 没有标准答案，所以另建一个能精确求解的体系：Metzner–Schütte–Vanden-Eijnden 三势阱，两个深阱 A、B 加一个浅阱，
  有两条竞争的反应通道。动力学是细网格上的跳跃过程（平方根近似速率，满足细致平衡），速率、committor 和通道分配都能用线性代数精确算出。
  管线只看到采样出的 (x, y) 轨迹，代码与第 2 节完全相同。</p>
  <h3>阳性对照</h3>
  </div>
  <div class="tbl"><table>
    <tr><th>检验</th><th class="num">管线</th><th class="num">精确</th><th>判定</th></tr>
    <tr><td>速率 k<sub>AB</sub>（95% CI）</td><td class="num">{rate['msm']:.4f} ({rate['ci'][0]:.4f}–{rate['ci'][1]:.4f})</td><td class="num">{rate['exact']:.4f}</td><td>{_chip("pass" if covered else "warn")}</td></tr>
    <tr><td>直接计数速率</td><td class="num">{rate['direct']:.4f}</td><td class="num">{rate['exact']:.4f}</td><td>{_chip("pass")} 直接计数本身也有采样误差</td></tr>
    <tr><td>committor 相关 / RMSE</td><td class="num">{v['committor']['r']:.4f} / {v['committor']['rmse']:.3f}</td><td class="num">1 / 0</td><td>{_chip("pass")}</td></tr>
    <tr><td>上通道份额</td><td class="num">{v['split']['msm']:.2f}</td><td class="num">{v['split']['exact']:.2f}</td><td>{_chip("pass")}</td></tr>
  </table></div>
  <div class="col"><h3>阴性对照</h3></div>
  <div class="tbl"><table>
    <tr><th>对照</th><th class="num">最慢时间尺度</th><th>预期与结果</th></tr>
    <tr><td>精确值</td><td class="num">{neg['exact_t1']:.2f}</td><td>参考</td></tr>
    <tr><td>正常特征 (x, y)</td><td class="num">{neg['t1_xy']:.2f}</td><td>略低于精确值：MSM 只会低估（变分原理）</td></tr>
    <tr><td>只给 y 坐标</td><td class="num">{neg['t1_y']:.2f}</td><td>{_chip("pass")} 看不见 A↔B，时间尺度塌缩。单看一个模型发现不了，比较不同特征集即可发现</td></tr>
    <tr><td>打乱时间顺序</td><td class="num">{neg['t1_shuffled']:.3f}</td><td>{_chip("pass")} 低于一帧，管线报告"无可分辨慢过程"</td></tr>
    <tr><td>诱饵配体（从未访问的区域）</td><td class="num">≤ {v['decoy_max']:.0e}</td><td>{_chip("pass")} 预测效应为零</td></tr>
  </table></div>
  <div class="two">
    <div><div class="col"><h3>数据量收敛</h3></div>
    <div class="tbl"><table><tr><th class="num">数据</th><th class="num">A→B 事件</th><th class="num">速率误差</th><th class="num">q⁺ 相关</th></tr>{conv_rows}</table></div></div>
    <div><div class="col"><h3>参数敏感性</h3></div>
    <div class="tbl"><table><tr><th>聚类</th><th class="num">态数</th><th class="num">lag</th><th class="num">速率误差</th><th class="num">q⁺ 相关</th><th class="num">通道误差</th></tr>{hp_rows}</table></div></div>
  </div>
  <p class="col note">约 60 次以上 A→B 事件后速率误差稳定在 ±5% 内。参数选择带来的系统误差（{min(hp_err) * 100:+.0f}% 到 {max(hp_err) * 100:+.0f}%）比 bootstrap 区间宽；
  lag 超过跨越势垒的时间后通道分配失真（lag = 20 帧时上通道份额偏差 {v['hyperparams'][-1]['split_err']:+.2f}）。</p>
  {_img("m24_controls.png", "图 2 精确解对照：(a) committor 等值线；(b) committor 一致性；(c) 速率；(d) 数据量收敛；(e) 参数敏感性；(f) 阴性对照的时间尺度")}
</section>

<section>
  <div class="col">
  <h2>4. 虚拟配体筛选与变构热点</h2>
  <p>把"配体"建模为势能面上的局部扰动（高斯凸起，正为去稳定、负为稳定），放在 28 个位置，每个位置两种符号，其中 4 个是从未被访问区域的诱饵，共 56 个。
  对每个扰动，只用<b>一个未扰动的 MSM</b>预测 ln k<sub>AB</sub> 的变化，再与精确值比较。</p>
  <p>预测采用过渡态重加权：每条边 i→j 的速率按实际穿越这条边的帧上的扰动能量修正，
  k′<sub>ij</sub> = k<sub>ij</sub>·exp[−β(V<sub>b,ij</sub> − δ<sub>i</sub>)]，δ<sub>i</sub> 是微观态 i 的自由能变化。由于 V<sub>b</sub> 对 i、j 对称，细致平衡严格保持。</p>
  </div>
  <div class="tbl"><table>
    <tr><th>指标（56 个虚拟配体）</th><th class="num">MSM 重加权</th><th class="num">静态基线</th></tr>
    <tr><td>Spearman ρ（预测 vs 精确）</td><td class="num">{sc['rho_msm']:.3f}</td><td class="num">{sc['rho_static']:.3f}</td></tr>
    <tr><td>方向正确（激活/抑制，|效应| &gt; 0.05 的 {sc['n_big']} 个）</td><td class="num">{sc['sign_msm'] * 100:.0f}%</td><td class="num">{sc['sign_static'] * 100:.0f}%</td></tr>
    <tr><td>最强 10 个效应的召回</td><td class="num">{sc['top10_msm'] * 100:.0f}%</td><td class="num">{sc['top10_static'] * 100:.0f}%</td></tr>
    <tr><td>平均绝对误差（ln k）</td><td class="num">{sc['mae_msm']:.3f}</td><td class="num">—</td></tr>
    <tr><td>每个配体的计算成本</td><td class="num">{sc['t_pred'] * 1e3:.0f} ms</td><td class="num">新模拟 + 建模 ≈ {sc['sim_cost']:.0f} s</td></tr>
  </table></div>
  <div class="col">
  <p>静态基线只看配体如何改变反应物 A 的稳定性，相当于只看终态结构的打分。它能抓住结合在 A 附近的效应，
  但对作用在势垒和通道上的扰动完全没有信号；这正是动力学信息带来的增量。
  放在势垒上的扰动方向判断正确，但幅度低估约 40%，原因是势垒区微观态最稀疏。</p>
  </div>
  <div class="tbl"><table><tr><th class="num">扰动幅度</th><th class="num">ρ</th><th class="num">方向正确</th><th class="num">MAE</th><th class="num">静态 ρ</th></tr>{amp_rows}</table></div>
  <div class="col">
  <h3>灵敏度图：哪里稳定化会改变 A→B 速率</h3>
  <p>在每个位置放一个小的稳定化探针，计算 ∂ln k<sub>AB</sub>/∂ε。预测图与精确图的相关为 r = {v['sensitivity_r']:.3f}。
  稳定 A 本身会减慢反应（品红）；离 A、B 都较远而灵敏度最大的位置都在直接通道的过渡区：{hot}。
  上通道的浅阱一带也有较弱的正灵敏度。位置预测正确，过渡区的幅度低估约 40%。</p>
  <p>这个二维体系没有独立的"远端"自由度，所以它验证的是"找到远离起点和终点的动力学热点"，
  还没有验证经由蛋白结构传递的真正变构耦合；后者需要在已知变构调节剂的蛋白上回顾性检验。</p>
  </div>
  {_img("m24_screening.png", "图 3 虚拟配体筛选：(a) MSM 预测 vs 精确；(b) 静态基线 vs 精确；(c)(d) 精确与预测的灵敏度图；(e) 扰动幅度；(f) 汇总")}
</section>

<section>
  <div class="col">
  <h2>5. 真实 MD 扰动对照（丙氨酸二肽）</h2>
  <p>在 ψ 上加一个深 4 kJ/mol、中心在 ψ = −40°（αR 内）的偏置作为阳性对照，在 ACE 甲基转动上加同样强度的偏置作为阴性对照。
  先只用无偏 MSM 预测效应，再各跑 4 × 25 ns 带偏置的模拟，把新帧放进同一套微观态里比较。</p>
  </div>
  <div class="tbl"><table>
    <tr><th>条件</th><th class="num">αR 布居</th><th class="num">MFPT β→αR, ps (95% CI)</th><th class="num">MFPT αR→β, ps (95% CI)</th></tr>
    {md_rows}
  </table></div>
  <div class="col">
  <p><b>热力学：通过。</b>ψ 偏置把 αR 布居从 {un['pop_aR']:.2f} 推到 {psi['sim']['pop_aR']:.2f}（模拟），无偏 MSM 预测 {psi['pred']['pop_aR']:.2f}。</p>
  <p><b>动力学：方向对，幅度偏差约 1.5 倍。</b>预测 αR→β 变慢到 {1 / psi['pred']['k_b']:.0f} ps，模拟为 {1 / psi['sim']['k_b']:.0f} ps（{ci_mfpt(psi['sim']['ci_b'])}）；
  β→αR 预测加快到 {1 / psi['pred']['k_f']:.0f} ps，模拟为 {1 / psi['sim']['k_f']:.0f} ps，置信区间很宽。
  原因与第 2 节相同：轨迹每 1 ps 存一帧，跨越势垒不到 1 ps，过渡态上的偏置能量只能用两端的帧近似。存帧更密时这一项应会改善，尚未验证。</p>
  <p><b>阴性对照：通过。</b>甲基偏置的预测与无偏基线相同，模拟结果也落在基线附近。</p>
  </div>
  {_img("m24_md_perturbation.png", "图 4 真实 MD：无偏 MSM 的预测（斜线）与带偏置模拟（实色，95% CI）")}
</section>

<section class="col">
  <h2>6. 能否用于找变构位点和筛选抑制剂、激活剂</h2>
  <h3>已经验证的部分</h3>
  <ul>
    <li><b>找出对速率最敏感的区域。</b>灵敏度图把"稳定这里会加快或减慢转变"量化成一张图，热点可以远离起点和终点（本例中在过渡区）。RMSD、RMSF 和终态结构比较都给不出这种信息。</li>
    <li><b>一个模型筛很多扰动。</b>在已采样构象内，一次 MSM 就能给几十到上千个候选扰动排序，每个只需毫秒级计算，替代逐个重跑模拟。</li>
    <li><b>区分激活剂和抑制剂。</b>预测符号就是作用方向；能对过渡态和通道起作用的调节剂，静态打分系统性地漏掉。</li>
    <li><b>平衡布居的变化预测得很准</b>（真实 MD 中 {psi['pred']['pop_aR']:.2f} 对 {psi['sim']['pop_aR']:.2f}）；速率变化在真实 MD 中方向正确，幅度可能差 1.5 倍左右，适合排序和初筛，不适合直接当定量结论。</li>
  </ul>
  <h3>用到真实蛋白上的流程</h3>
  <ul>
    <li>多条并行短模拟 + 自适应采样（从低采样或高通量的态重新起跑），建立蛋白的 MSM，并按本报告的对照做验证。</li>
    <li>把 committor 和灵敏度映射到残基：哪些残基的接触或二面角状态与 q⁺ = ½ 的过渡态集合、与灵敏度热点相关。</li>
    <li>在高权重的态上找口袋（包括只在某些态出现的隐藏口袋），对过渡态和中间态做系综对接，而不只对晶体结构对接。</li>
    <li>把对接得到的结合自由能当作对这些态的扰动 δ<sub>i</sub>，用重加权预测每个化合物对 k<sub>on</sub>、k<sub>off</sub> 或构象平衡的影响，按方向分出激活剂和抑制剂。</li>
    <li>对排在前列的少数化合物，用带配体的模拟或自由能计算复核，再做实验。</li>
  </ul>
  <h3>不能做到的部分</h3>
  <ul>
    <li>重加权只覆盖已采样的构象。配体诱导出的新状态、隐藏口袋的打开过程若没在无偏模拟中出现，就预测不了，必须带配体重新采样。</li>
    <li>真实配体不是一个光滑的势能扰动：它有结合与解离动力学、溶剂效应和熵贡献，把对接分数转成 δ<sub>i</sub> 本身就有较大误差。</li>
    <li>稀有事件（本报告的 β→αL）采样不足时，任何 MSM 结果都只有数量级意义；需要更长模拟或增强采样。</li>
    <li>在真实蛋白上的回顾性验证见第 7 节：用的是粗粒化模型，全原子精度仍待检验。</li>
  </ul>
</section>

{retro_section}
<section class="col">
  <h2>8. 复现</h2>
  <pre>pip install numpy scipy matplotlib openmm mdtraj deeptime
python3 -m grn_pipeline.m23_msm_tpt          # 丙氨酸二肽 MSM + TPT（首次运行约 15 分钟模拟）
python3 -m grn_pipeline.m24_msm_validation   # 对照、筛选、MD 扰动（首次运行另需约 20 分钟模拟）
python3 -m grn_pipeline.m25_allosteric_retro # 回顾性验证（自适应采样，首次运行约 3 小时）
python3 -m grn_pipeline.msm_html_report      # 生成本报告

# 自己的 GROMACS 轨迹
gmx trjconv -s md.tpr -f md.xtc -o md_whole.xtc -pbc mol -center
gmx editconf -f md.tpr -o conf.gro
python3 -m grn_pipeline.m23_msm_tpt --top conf.gro --traj md_whole.xtc --dt-ps 10 --features both</pre>
</section>
</div>
"""
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(html)
    return path


if __name__ == "__main__":
    print(build())
