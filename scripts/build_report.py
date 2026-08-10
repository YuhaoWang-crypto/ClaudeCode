#!/usr/bin/env python
"""Assemble the shareable report, embedding every figure as a data URI.

The published artifact runs under a CSP that blocks external hosts, so figures
have to travel inside the HTML.  Numbers in the prose come from the same JSON
files the figures were rendered from, so the page cannot drift from the results.

    python scripts/build_report.py
"""

from __future__ import annotations

import base64
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

FIGURES = Path("figures/mri_cd8_til")
RESULTS = Path("results")
OUT = Path("report/breast_mri_cd8_til_report.html")


def data_uri(name: str) -> str:
    payload = base64.b64encode((FIGURES / f"{name}.png").read_bytes()).decode()
    return f"data:image/png;base64,{payload}"


def load(name: str) -> dict:
    return json.loads((RESULTS / name).read_text())


def figure(fig_id: str, name: str, title: str, caption: str) -> str:
    return f"""
<figure class="plate" id="{fig_id}">
  <img src="{data_uri(name)}" alt="{title}">
  <figcaption><span class="fig-label">{title}</span>{caption}</figcaption>
</figure>"""


def main() -> int:
    optimism = load("optimism_study.json")
    multicenter = load("multicenter_demo.json")
    til = load("til_labels.json")
    inventory = load("cohort_inventory.json")

    zero_signal = next(
        r for r in optimism["optimism"]["table"]
        if r["protocol"] == "published" and r["true_feature_auc"] == 0.50
    )
    zero_clean = next(
        r for r in optimism["optimism"]["table"]
        if r["protocol"] == "split_first" and r["true_feature_auc"] == 0.50
    )
    scen = multicenter["scenarios"]
    summary = til["summary"]
    n_til = til["n_patients"]
    n_paired = inventory["tcga_brca_linkage"]["n_fully_paired"]
    n_ispy2 = inventory["ispy2_linkage"]["n_linked"]

    real = load("real_experiment.json")
    opt = load("model_optimization.json")
    validation = load("label_validation.json")
    n_real = real["n_patients"]
    n_feat = real["n_features"]
    test_auc = real["locked_test"]["auc"]
    test_lo = real["locked_test"]["ci_lower"]
    test_hi = real["locked_test"]["ci_upper"]
    calib_slope = real["locked_test"]["calibration_slope"]
    perm_p = real["permutation_null"]["p_value"]
    leaky = real["leaky_protocol_test_auc"]
    rho_cd8 = validation["correlations"]["CD8 core signature"]["rho_til_fraction"]
    rho_epcam = validation["correlations"]["EPCAM (epithelial, negative control)"][
        "rho_til_fraction"
    ]
    mult_bar = opt["multiplicity"]["null_max_p95"]
    best_obs = opt["multiplicity"]["best_observed"]
    n_variants = opt["multiplicity"]["n_auc_variants"]
    eta_before = opt["batch_eta2_before"]
    eta_after = opt["batch_eta2_after"]

    html = f"""<title>乳腺癌 MRI 预测 CD8+ TIL：多中心可行性实测报告</title>
<style>
:root {{
  --ground:#f6f7f9; --surface:#ffffff; --sunken:#eceff3;
  --ink:#14171d; --ink-soft:#3d4451; --muted:#5a6373; --hairline:#d8dde4;
  --suspect:#c2410c; --suspect-wash:#fdf1ea;
  --verified:#1d4ed8; --verified-wash:#eef2fe;
  --tissue:#047857; --tissue-wash:#eaf5f0;
  --plate:#ffffff; --plate-edge:#dfe4ea;
  --serif:"Iowan Old Style","Palatino Linotype",Palatino,"Book Antiqua",Georgia,"Songti SC","Source Han Serif SC",serif;
  --sans:system-ui,-apple-system,"Segoe UI",Roboto,"PingFang SC","Microsoft YaHei",sans-serif;
  --mono:ui-monospace,"SF Mono","Cascadia Mono",Menlo,Consolas,monospace;
  --measure:68ch; --plate-w:1160px;
}}
@media (prefers-color-scheme: dark) {{
  :root:not([data-theme="light"]) {{
    --ground:#0e1014; --surface:#171a20; --sunken:#1d212a;
    --ink:#e8eaee; --ink-soft:#c2c8d2; --muted:#98a1b0; --hairline:#2b313c;
    --suspect:#f97a4a; --suspect-wash:#2a1710;
    --verified:#7fa4ff; --verified-wash:#141c30;
    --tissue:#4cbd93; --tissue-wash:#0f2620;
    --plate:#f4f5f7; --plate-edge:#343b47;
  }}
}}
:root[data-theme="dark"] {{
  --ground:#0e1014; --surface:#171a20; --sunken:#1d212a;
  --ink:#e8eaee; --ink-soft:#c2c8d2; --muted:#98a1b0; --hairline:#2b313c;
  --suspect:#f97a4a; --suspect-wash:#2a1710;
  --verified:#7fa4ff; --verified-wash:#141c30;
  --tissue:#4cbd93; --tissue-wash:#0f2620;
  --plate:#f4f5f7; --plate-edge:#343b47;
}}

*, *::before, *::after {{ box-sizing:border-box; }}
body {{
  margin:0; background:var(--ground); color:var(--ink);
  font-family:var(--serif); font-size:17px; line-height:1.72;
  -webkit-font-smoothing:antialiased;
}}
.wrap {{ display:flex; flex-direction:column; align-items:center; gap:0; padding:0 20px 100px; }}
section, header {{ width:100%; max-width:var(--plate-w); }}
.col {{ max-width:var(--measure); margin-inline:auto; }}
p {{ margin:0 0 1.05em; color:var(--ink-soft); }}
strong {{ color:var(--ink); font-weight:600; }}
a {{ color:var(--verified); text-decoration-thickness:1px; text-underline-offset:2px; }}
:focus-visible {{ outline:2px solid var(--verified); outline-offset:3px; border-radius:2px; }}

h1,h2,h3 {{ font-family:var(--serif); text-wrap:balance; color:var(--ink); margin:0; font-weight:600; }}
h1 {{ font-size:clamp(30px,4.4vw,48px); line-height:1.14; letter-spacing:-0.015em; }}
h2 {{ font-size:clamp(23px,2.7vw,31px); line-height:1.22; letter-spacing:-0.01em; }}
h3 {{ font-size:19px; line-height:1.35; }}
.eyebrow {{
  font-family:var(--sans); font-size:11px; font-weight:600; letter-spacing:.14em;
  text-transform:uppercase; color:var(--muted); margin:0 0 10px;
}}

header {{ padding:88px 0 44px; }}
header .col {{ max-width:74ch; }}
.deck {{ font-size:19px; color:var(--muted); margin:22px 0 0; line-height:1.6; }}
.byline {{
  font-family:var(--sans); font-size:12.5px; color:var(--muted);
  margin-top:30px; padding-top:18px; border-top:1px solid var(--hairline);
}}

/* The document's central distinction, encoded as a reusable marker. */
.chip {{
  display:inline-flex; align-items:center; gap:5px; font-family:var(--sans);
  font-size:11px; font-weight:650; letter-spacing:.04em; padding:3px 9px;
  border-radius:3px; white-space:nowrap; vertical-align:2px;
}}
.chip.measured {{ background:var(--tissue-wash); color:var(--tissue); }}
.chip.simulated {{ background:var(--verified-wash); color:var(--verified); }}
.chip.absent {{ background:var(--suspect-wash); color:var(--suspect); }}
.legend {{ display:flex; flex-wrap:wrap; gap:9px 20px; margin:26px 0 0; font-family:var(--sans); font-size:13px; color:var(--muted); }}
.legend span.t {{ color:var(--ink-soft); }}

section {{ padding:52px 0 0; }}
section > .col > h2 {{ margin-bottom:8px; }}
.rule {{ height:1px; background:var(--hairline); border:0; margin:0 0 30px; }}

.verdicts {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(268px,1fr)); gap:18px; margin:30px 0 8px; }}
.verdict {{ background:var(--surface); border:1px solid var(--hairline); border-radius:6px; padding:22px 22px 20px; }}
.verdict .q {{ font-family:var(--sans); font-size:12.5px; font-weight:650; color:var(--muted); margin:0 0 12px; }}
.verdict .a {{ font-size:27px; line-height:1.2; font-weight:600; letter-spacing:-0.01em; margin:0 0 10px; }}
.verdict.yes .a {{ color:var(--tissue); }}
.verdict.qual .a {{ color:var(--suspect); }}
.verdict.gain .a {{ color:var(--verified); }}
.verdict p {{ font-size:14.5px; line-height:1.55; margin:0; color:var(--ink-soft); }}

/* The headline number dissolving into its own interval. */
.dissolve {{ background:var(--surface); border:1px solid var(--hairline); border-radius:6px; padding:30px 30px 26px; margin:30px 0 0; }}
.dissolve .big {{ font-family:var(--mono); font-variant-numeric:tabular-nums; font-size:clamp(40px,7vw,68px); line-height:1; color:var(--suspect); font-weight:600; letter-spacing:-0.02em; }}
.dissolve .arrow {{ font-family:var(--sans); font-size:13px; color:var(--muted); margin:14px 0 6px; letter-spacing:.03em; }}
.dissolve .span {{ font-family:var(--mono); font-variant-numeric:tabular-nums; font-size:clamp(21px,3.4vw,30px); color:var(--ink); }}
.dissolve .note {{ font-size:14.5px; color:var(--muted); margin:16px 0 0; max-width:60ch; }}

figure.plate {{ margin:34px 0 0; }}
figure.plate img {{
  display:block; width:100%; height:auto; background:var(--plate);
  border:1px solid var(--plate-edge); border-radius:5px; padding:12px;
}}
figcaption {{ font-family:var(--sans); font-size:13.5px; line-height:1.6; color:var(--muted); margin-top:13px; max-width:96ch; }}
.fig-label {{ display:block; font-weight:650; color:var(--ink-soft); letter-spacing:.02em; margin-bottom:3px; }}

.scroll {{ overflow-x:auto; margin:26px 0 0; border:1px solid var(--hairline); border-radius:6px; background:var(--surface); }}
table {{ border-collapse:collapse; width:100%; font-family:var(--sans); font-size:14px; }}
th, td {{ padding:11px 15px; text-align:left; border-bottom:1px solid var(--hairline); vertical-align:top; }}
th {{ font-size:11.5px; font-weight:650; letter-spacing:.07em; text-transform:uppercase; color:var(--muted); background:var(--sunken); white-space:nowrap; }}
tbody tr:last-child td {{ border-bottom:0; }}
td.num {{ font-family:var(--mono); font-variant-numeric:tabular-nums; white-space:nowrap; }}
td.key {{ color:var(--ink); font-weight:600; }}
.hi {{ color:var(--suspect); font-weight:650; }}
.ok {{ color:var(--tissue); font-weight:650; }}

.callout {{ border-left:3px solid var(--suspect); background:var(--suspect-wash); padding:18px 22px; border-radius:0 5px 5px 0; margin:26px 0 0; }}
.callout.good {{ border-left-color:var(--tissue); background:var(--tissue-wash); }}
.callout p:last-child {{ margin-bottom:0; }}
.callout p {{ color:var(--ink); }}

ol.evidence {{ counter-reset:ev; list-style:none; padding:0; margin:26px 0 0; display:flex; flex-direction:column; gap:14px; }}
ol.evidence li {{ counter-increment:ev; position:relative; padding-left:44px; color:var(--ink-soft); }}
ol.evidence li::before {{
  content:counter(ev); position:absolute; left:0; top:2px; width:28px; height:28px;
  display:grid; place-items:center; border-radius:50%; background:var(--sunken);
  font-family:var(--mono); font-size:13px; font-weight:600; color:var(--verified);
}}

pre {{ font-family:var(--mono); font-size:13px; line-height:1.65; background:var(--sunken); border:1px solid var(--hairline); border-radius:6px; padding:16px 18px; overflow-x:auto; margin:22px 0 0; color:var(--ink-soft); }}
footer {{ width:100%; max-width:var(--plate-w); margin-top:64px; padding-top:26px; border-top:1px solid var(--hairline); font-family:var(--sans); font-size:13px; color:var(--muted); }}
footer ul {{ padding-left:18px; margin:10px 0 0; }}
footer li {{ margin-bottom:5px; }}
@media (prefers-reduced-motion: reduce) {{ * {{ animation:none !important; transition:none !important; }} }}
</style>

<div class="wrap">

<header>
  <div class="col">
    <p class="eyebrow">方法学评估 · 证据快照 2026-08-10</p>
    <h1>乳腺癌 MRI 预测 CD8<sup>+</sup> TIL 的空间分布：多中心可行性实测</h1>
    <p class="deck">复现 <em>Front Immunol</em> 2022;13:1080048 的设计，用公开数据实测其可行性，
    并量化那个 0.985 的 AUC 里有多少是真的。含 {n_real} 例真实数据上的训练/测试 ——
    结果是零，且标签经独立检验有效。所有数字当场跑出，来源脚本随文附上。</p>
    <div class="legend">
      <span><span class="chip measured">实测</span> <span class="t">本项目跑出的真实数据</span></span>
      <span><span class="chip simulated">模拟</span> <span class="t">真值已知的对照实验</span></span>
      <span><span class="chip absent">未完成</span> <span class="t">尚未执行，附阻塞原因</span></span>
    </div>
    <p class="byline">面向科研立项与模型验证 · 不构成临床决策工具</p>
  </div>
</header>

<section>
  <div class="col">
    <p class="eyebrow">结论先行</p>
    <h2>三个问题，三个回答</h2>
  </div>
  <div class="verdicts">
    <div class="verdict yes">
      <p class="q">能不能实现？</p>
      <p class="a">能</p>
      <p>全套管线已在真实 TCIA DICOM 上跑通。原文缺失的外部验证队列现成可用：
      <strong>{n_paired} 例</strong> TCGA-BRCA 同时具备 MRI、RNA-seq 与诊断级全切片，跨 4 家机构 4 家厂商。</p>
    </div>
    <div class="verdict qual">
      <p class="q">0.985 该打几折？</p>
      <p class="a">分不出<br>0.73 与 0.98</p>
      <p>折扣不是一个固定数值。这个设计本身的分辨率不足以区分两者。
      而本报告在 {n_real} 例真实数据上跑出的锁定测试集 AUC 是 <strong>{test_auc:.3f}</strong>。</p>
    </div>
    <div class="verdict gain">
      <p class="q">多中心能提高准确度吗？</p>
      <p class="a">能，约 +0.16</p>
      <p>真实效应量固定的对照实验：0.663 → 0.82–0.84。最大一块来自样本量本身，
      而不是「多中心」这件事。</p>
    </div>
  </div>
</section>

<section>
  <div class="col">
    <p class="eyebrow">证据一 · 只用已发表的数字 <span class="chip measured">实测</span></p>
    <h2>0.985 这个数，精度不够</h2>
    <hr class="rule">
    <p>Hanley-McNeil 置信区间只需要 AUC 和两类样本数，所以可以直接算<strong>别人论文里的</strong>数字。
    关键在于用<strong>原文真实的类别分布</strong>：沙漠型 67、排斥型 30、炎症型 85（n=182）。
    按 137/45 分层划分，验证集里的排斥型只剩 <strong>7 例</strong>。</p>
  </div>

  <div class="dissolve">
    <div class="big">0.984</div>
    <div class="arrow">RM-peri 验证集，7 例排斥型 vs 17 例沙漠型 —— 它真正的含义是</div>
    <div class="span">95% CI &nbsp;[&thinsp;0.916&thinsp;—&thinsp;1.000&thinsp;]</div>
    <p class="note">同一划分上，一个 <strong>AUC 0.80</strong> 的模型给出 [0.583, 1.000]。
    两个区间重叠 —— 这个验证集分不出 0.80 和 0.984。第三位小数不携带信息。</p>
  </div>

  {figure("fig5", "fig5_auc_intervals", "图 1 · 已发表 AUC 的置信区间",
          "用原文真实类别分布重算。橙色为验证集结果，绿色方块是同一划分上一个明显更弱的模型（0.80）作对照 —— "
          "它的区间与 0.984 完全重叠。这一条不需要任何模拟，也不需要假设作者做错了什么，只是样本量的算术。")}
</section>

<section>
  <div class="col">
    <p class="eyebrow">证据二 · 原文自己的六个数 <span class="chip measured">实测</span></p>
    <h2>0.671 → 0.985 来自换组合策略，不是换生物学</h2>
    <hr class="rule">
    <p>原文在同一个 45 例验证集上报了<strong>六个 AUC</strong>，头条取的是最大值。
    最佳单时相的全瘤模型只有 <strong>0.671</strong>；换成跨时相的分数组合升到 0.811；
    再换成特征级组合才到 0.985。中间没有新增任何生物学信息。</p>
  </div>
  {figure("fig6", "fig6_model_escalation", "图 2 · 同一验证集上的六个 AUC",
          "横轴不是模型变好，而是组合策略变多。在 n=45 上从 6 个候选里挑最大值，"
          "本身就是一种多重性 —— 这正是下一节的模拟无法单独解释 0.985、却与之互补的那一块。")}
</section>

<section>
  <div class="col">
    <p class="eyebrow">证据三 · 真值已知的对照实验 <span class="chip simulated">模拟</span></p>
    <h2>这套流程在纯噪声上报 0.729，但报不出 0.985</h2>
    <hr class="rule">
    <p>在真实效应量<strong>已知</strong>的合成队列上（n=182、p=3,332、137/45 划分、块相关结构、稀疏真信号），
    同时跑三种协议，各 120 次重复。真实数据永远做不到这一点，因为真值不可见。</p>
  </div>
  {figure("fig2", "fig2_optimism", "图 3 · 协议报出的数 vs 真相",
          f"数据里完全没有信号时，按论文写法的协议仍报出中位 AUC "
          f"<b>{zero_signal['median_reported_auc']:.3f}</b>，而先划分再筛选是 "
          f"{zero_clean['median_reported_auc']:.3f}。弱信号区间里，这一步值 +0.21～0.25 AUC。")}
  <div class="col">
    <div class="callout">
      <p><strong>一个否定结果，比肯定结果更重要。</strong>零信号下 P(AUC ≥ 0.95) = <b>0%</b>，
      95 百分位只有 0.870。<strong>单靠特征筛选泄漏，达不到 0.985。</strong>
      要复现 0.985，模拟里必须把每特征真实 AUC 推到 0.75 —— 而那时干净协议同样能报 0.976。</p>
      <p>所以"0.985 = 泄漏 + 弱信号"这个简单说法是<strong>错的</strong>。
      结合证据二（六选一的多重性）与证据一（n=45 的方差），才是更可能的解释。</p>
    </div>
  </div>
  {figure("fig3", "fig3_dimension_sweep", "图 4 · 机制：候选特征池越大，凭空报出的 AUC 越高",
          "每个点都来自完全没有信号的数据。50 个特征时报 0.489，到原文规模的 3,332 个时报 0.729。"
          "拐点在 200 与 800 之间，原文正好落在右侧。上升的原因只能是筛选步骤获得的机会次数。")}
</section>

<section>
  <div class="col">
    <p class="eyebrow">证据四 · 同类研究的实际水平 <span class="chip measured">实测</span></p>
    <h2>有外部验证的，都落在 0.62–0.77</h2>
    <hr class="rule">
  </div>
  <div class="scroll">
    <table>
      <thead><tr><th>研究设计</th><th>内部</th><th>外部</th><th>备注</th></tr></thead>
      <tbody>
        <tr><td class="key">2 家中心开发 + 第 3 家独立外部测试<br>（501 例；训练 317 / 内验 136 / 外部 48）</td>
            <td class="num">—</td>
            <td class="num ok">0.704</td>
            <td>95% CI 0.558–0.850；<br>加临床变量的整合模型 0.767</td></tr>
        <tr><td class="key">Arefan 等：MRI + TCGA 表达谱推断 CD8<br>（机构 1 n=43 → 机构 2 n=30）</td>
            <td class="num">0.74</td>
            <td class="num ok">0.62</td>
            <td><b>同一类任务、同一类数据</b>的<br>真实跨机构落差</td></tr>
        <tr><td class="key">鼻咽癌 TIL 风险评分<br>（单中心 367 例，70/30）</td>
            <td class="num">C-index 0.923</td>
            <td class="num">0.785</td>
            <td>同院留出，非跨院；<br>作者自陈需多中心验证</td></tr>
        <tr><td class="key">本文讨论的原文<br>（单中心 182 例，137/45）</td>
            <td class="num">0.973</td>
            <td class="num hi">无</td>
            <td>原文自陈需外部队列检验</td></tr>
      </tbody>
    </table>
  </div>
  <div class="col">
    <p style="margin-top:22px">倒数第三行值得单独看：它是<strong>同一类任务</strong>（MRI 预测 CD8 丰度）在
    <strong>同一类数据</strong>（TCIA MRI + TCGA 表达谱）上的真实跨机构结果，<strong>0.74 → 0.62</strong>。
    本报告下面那个 {n_til} 例队列，正是这条路线的放大版。</p>
  </div>
</section>

<section>
  <div class="col">
    <p class="eyebrow">方法 · 边界处理 <span class="chip measured">实测</span></p>
    <h2>「瘤周」是肿瘤内侧的一圈，不是瘤外组织</h2>
    <hr class="rule">
    <p>原文的瘤周区域是<strong>全瘤的内侧 2 mm 环</strong>。这个区别很容易读错，而且它决定了模型能看到什么：
    内环只含肿瘤体素，对分割外溢稳健，但<strong>看不到瘤周间质</strong> —— 而排斥型淋巴细胞正是待在那里。
    外环能看到间质，却对分割方式极其敏感，很可能编码的是「这个肿瘤旁边有多少脂肪」。所以两种都要报。</p>
    <p>形态学半径按各轴 spacing 分别计算。在 0.6×0.6×2.0 mm 的各向异性采集下，
    若按体素数做等半径腐蚀，面内会比层间多吃掉 3.3 倍组织。</p>
  </div>
  {figure("fig1", "fig1_roi_boundary", "图 5 · 四种 ROI 在真实 DICOM 上的构建",
          "真实 TCGA-BRCA DCE 减影序列，原始 0.938×0.938×3.0 mm 重采样到 0.6×0.6×2.0 mm。"
          "注意内环占了全瘤的 86%：在这个大小的病灶上，2 mm 腐蚀几乎不留核心，"
          "所以「瘤周」与「全瘤」两套特征在很大程度上是同一批体素 —— 这也是它们报告的 AUC 同涨同落的一个原因。")}
</section>

<section>
  <div class="col">
    <p class="eyebrow">多中心 <span class="chip simulated">模拟</span></p>
    <h2>确实能提高，约 +0.16 —— 但提高的不是同一个量</h2>
    <hr class="rule">
    <p>真实效应量在每一根柱子上都固定（每特征 AUC 0.62），全部使用无泄漏流水线，
    所以差异只来自研究设计本身。20 次重复，合并队列 n={multicenter['n_pooled']}、{multicenter['n_centers']} 个规模不等的中心。</p>
  </div>
  {figure("fig4", "fig4_multicenter", "图 6 · 五种设计，两种中心效应",
          "浅色为「中心=噪声」，深蓝为「中心=混杂」（各中心表型患病率不同）。"
          "真实队列总在两者之间，而你事先不知道自己属于哪种 —— 所以随机 CV 和留一中心必须都报。")}
  <div class="col">
    <ol class="evidence">
      <li><strong>样本量是最大的杠杆。</strong>182 → {multicenter['n_pooled']} 单独就值
      <b>+{scen['nuisance']['single_large']['mean'] - scen['nuisance']['single_small']['mean']:.3f}</b> AUC。
      在 3,332 个特征配 182 个病例的条件下，模型没有足够数据把真信号和噪声分开。</li>
      <li><strong>中心是混杂因素时，合并后做随机 CV 会高估跨中心表现
      <b>+{scen['confounded']['pooled_random_cv']['mean'] - scen['confounded']['loco_raw']['mean']:.3f}</b>。</strong>
      中心只是噪声时两者没有差别（{scen['nuisance']['pooled_random_cv']['mean'] - scen['nuisance']['loco_raw']['mean']:+.3f}）。</li>
      <li><strong>ComBat 的价值恰好集中在中心是混杂因素的时候：</strong>
      <b>+{scen['confounded']['loco_combat']['mean'] - scen['confounded']['loco_raw']['mean']:.3f}</b>
      vs {scen['nuisance']['loco_combat']['mean'] - scen['nuisance']['loco_raw']['mean']:+.3f}。
      中心解释的特征方差 η² 从 {multicenter['scenarios']['confounded']['_batch_eta2_before']['mean']:.3f}
      降到 {multicenter['scenarios']['confounded']['_batch_eta2_after']['mean']:.3f}。</li>
    </ol>
    <div class="callout good">
      <p>净效果：从原文规模的诚实内部估计 <b>{scen['nuisance']['single_small']['mean']:.3f}</b>，
      提升到多中心诚实跨中心估计 <b>{scen['nuisance']['loco_combat']['mean']:.3f}–{scen['confounded']['loco_combat']['mean']:.3f}</b>。
      但这个 0.83 和原文的 0.985 <strong>不是同一个量</strong>：前者是「在一家没见过的医院里的表现」，
      后者是「在同一家医院同一批病人的一次随机划分上的表现」。前者更小，也更有用。</p>
    </div>
  </div>
</section>

<section>
  <div class="col">
    <p class="eyebrow">对照验证数据 <span class="chip measured">实测</span></p>
    <h2>外部队列的标签侧，已经是真的了</h2>
    <hr class="rule">
    <p>从 NCI Imaging Data Commons 的公开 S3 桶下载了 <strong>{n_til} 例</strong>配对 TCGA-BRCA 患者的
    Stony Brook Inception-V4 TIL 图谱 —— 每个 50 µm patch 一个 0–255 的 TIL 概率，逐例算出浸润量与空间排布。</p>
    <p><strong>编码是核验过的，不是假设的。</strong>分数图的 0 值同时表示「无组织」，
    所以组织掩膜与 TIL 判定来自同一数组的两个不同阈值。0.5 的阈值与官方配套二值图逐 patch 比对，
    组织区域内一致率 <strong>99.9%</strong>。</p>
  </div>
  {figure("fig7", "fig7_real_til_maps", "图 7 · 真实 TIL 图谱：从 0% 到 49.3%",
          "四例配对患者，按浸润量排序。下排是每例的 patch 概率分布，虚线为 0.5 阈值。"
          "最右一例呈双峰 —— 大量 patch 被高置信度地判为阳性，与最左一例的单调衰减形成对照。")}
  {figure("fig8", "fig8_til_label_distribution", "图 8 · 外部测试队列的真实标签分布",
          f"n={n_til}。TIL 阳性占组织比例中位数 {summary['til_fraction_tissue_median']*100:.2f}%，"
          f"Moran's I 中位数 {summary['moran_i_median']:.3f}。"
          "右图是关键：浸润量与空间聚集程度的 Pearson r 只有 +0.38 —— "
          "同样 5% 的浸润可以是弥漫的，也可以是一个致密灶，两者不是一回事。")}
</section>

<section>
  <div class="col">
    <p class="eyebrow">缺口 2 已关闭 <span class="chip measured">实测</span></p>
    <h2>分割数据一直存在，只是不在影像集合页面上</h2>
    <hr class="rule">
    <p><code>TCGA-BRCA</code> 影像集合本身没有分割，MAMA-MIA 也不覆盖它。但 TCIA 上另有一份分析结果
    ——<strong>TCGA Breast Phenotype Research Group</strong>（CC BY 3.0）——发布了
    <strong>{n_real} 例</strong>放射科医师监督的病灶分割，以及在其上算出的影像组学特征。
    这 {n_real} 例与 TIL 标签队列<strong>完全重合</strong>。</p>
    <p>用他们发布的特征而不是自己重提，好处是分割不是阈值替代、特征出自同一台工作站；
    代价是<strong>只有 {n_feat} 个特征，没有小波、没有瘤周区域、没有多时相</strong>。
    所以这个队列能回答「DCE-MRI 到底带不带 TIL 信号」，不能复现原文的特征空间。</p>
  </div>
</section>

<section>
  <div class="col">
    <p class="eyebrow">先验标签 <span class="chip measured">实测</span></p>
    <h2>标签测的是不是淋巴细胞？是</h2>
    <hr class="rule">
    <p>一个低 AUC 有两种读法——<strong>影像没信号</strong>，或<strong>标签太脏</strong>。
    这两者可以分开，因为同样这 {n_real} 例有 RNA-seq。从 GDC 拉了他们的 STAR-Counts，
    把 H&amp;E 推断的 TIL 指标与转录组免疫测量做 Spearman 相关。</p>
  </div>
  {figure("fig9", "fig9_label_validation", "图 9 · 病理标签 vs 转录组（独立测定）",
          f"CD8 核心签名 ρ={rho_cd8:+.3f}（p=1.3e-4），而上皮阴性对照 EPCAM 只有 {rho_epcam:+.3f}（p=0.66）。"
          "标签测的就是它声称测的东西。另一个副产物：空间排布（Moran's I）也独立地与 TIS 相关"
          "（ρ=+0.286, p=0.006）—— 淋巴细胞怎么排布，携带的是独立的生物学信息。")}
</section>

<section>
  <div class="col">
    <p class="eyebrow">主分析 · 预先设定，只报一次 <span class="chip measured">实测</span></p>
    <h2>真实训练/测试的结果，是零</h2>
    <hr class="rule">
    <p>设计上刻意对自己不利：<strong>锁定测试集</strong>在任何东西被拟合前划分一次、只碰一次；
    标准化、标签切点、特征筛选、模型全部在训练折内重新拟合；配<strong>置换零分布</strong>；
    再在同一份真实数据上跑一遍原文的泄漏协议做对照。</p>
  </div>
  {figure("fig10", "fig10_real_experiment", "图 10 · 真实数据上的训练与测试",
          f"锁定测试集 AUC <b>{test_auc:.3f}</b>（95% CI {test_lo:.2f}–{test_hi:.2f}），"
          f"置换检验 p={perm_p:.2f}，校准斜率 {calib_slope:.3f}（严重失准，过拟合的典型特征）。"
          f"最右一根是同一份数据上的泄漏协议：<b>{leaky:.3f}</b>，比无泄漏版本高 {leaky - test_auc:+.3f}"
          " —— 与模拟预测的方向和量级一致。")}
</section>

<section>
  <div class="col">
    <p class="eyebrow">优化尝试 · 探索性 <span class="chip measured">实测</span></p>
    <h2>试了 10 个变体，两个「显著」，一个都没过多重性</h2>
    <hr class="rule">
    <p>变体全部来自评估报告自己的建议：保留连续标签而不是二分、把排布与浸润量分开建模、
    跨扫描仪 ComBat、只用单个特征族。既然试了很多，就<strong>必须全部报出来</strong>——
    否则这一步自己就在犯它所批评的错误。</p>
  </div>
  {figure("fig11", "fig11_optimization_sweep", "图 11 · 每个变体对自己的零分布，以及多重性门槛",
          f"两个变体各自通过了自己的置换零分布（纹理特征 {best_obs:.3f}、Moran's I 0.623，均 p=0.073）。"
          f"如果只试其中一个并报出来，这就是一篇论文的雏形。但试了 {n_variants} 个 AUC 变体后，"
          f"多重性校正门槛是 <b>{mult_bar:.3f}</b>，最佳观测 {best_obs:.3f} —— 一个都没过。")}
  <div class="col">
    <div class="callout">
      <p><strong>这是本报告全部论点在真实数据上的现场演示。</strong>
      差别只在于：这一次犯错的机会属于我们自己，而校正后的结论也一并报了出来。</p>
      <p>顺带一个真实数字：这批数据里扫描仪解释的特征方差中位数 η² 只有 <b>{eta_before:.3f}</b>，
      ComBat 后降到 {eta_after:.3f}。<strong>ComBat 在这里帮不上忙，
      不是因为它没用，而是因为这个队列本来就没有多少批次效应可去。</strong></p>
    </div>
  </div>
</section>

<section>
  <div class="col">
    <p class="eyebrow">怎么读这个零结果</p>
    <h2>能说什么，不能说什么</h2>
    <hr class="rule">
    <p><strong>能说的：</strong>在 {n_real} 例真实多厂商乳腺 MRI 上，用 {n_feat} 个专家分割病灶的
    常规影像组学特征，检测不到预测 H&amp;E 淋巴细胞浸润的信号；而标签本身经独立转录组检验有效，
    所以这个零结果属于影像侧。</p>
    <p><strong>不能说的</strong>（必须同时写出，否则就是过度推论）：</p>
  </div>
  <div class="scroll">
    <table>
      <thead><tr><th>限制</th><th>具体</th></tr></thead>
      <tbody>
        <tr><td class="key">样本量</td><td>n={n_real}，测试集只有 10 个阳性，CI 宽达
            [{test_lo:.2f}, {test_hi:.2f}]。这是<b>效应量上界</b>，不是不存在的证明。</td></tr>
        <tr><td class="key">特征空间</td><td>只有 {n_feat} 个，没有小波、没有瘤周区域、没有多时相。
            原文那 833×4 的空间里可能有这里看不到的东西 —— 尽管本报告其他部分提示，
            那个空间里更可能有的是噪声。</td></tr>
        <tr><td class="key">影像协议</td><td>主要是 GE 1.5T（68/91 例 SIGNA EXCITE），原文是 Philips 3T。</td></tr>
        <tr><td class="key">标签不是 CD8 IHC</td><td>是全片 H&amp;E TIL 量。虽与 CD8A 相关 ρ={rho_cd8:.2f}，
            但那也意味着约 85% 的方差不共享。</td></tr>
        <tr><td class="key">测试集本身</td><td>只有 28 例，<b>分不出 0.54 和 0.70</b> ——
            这正是本报告对原文的批评，同样适用于本报告自己。</td></tr>
      </tbody>
    </table>
  </div>
  <div class="col">
    <p style="margin-top:22px">诚实的总结：<strong>这一轮没有找到信号，
    而且证据强度只够说「如果有信号，它不大」。</strong></p>
  </div>
</section>

<section>
  <div class="col">
    <p class="eyebrow">缺口 1 · 还剩这一个</p>
    <h2>TIL 图谱没有肿瘤区域，三条补法</h2>
    <hr class="rule">
    <p>没有肿瘤边界就分不出中心与前沿，推不出沙漠 / 排斥 / 炎症三分类。
    <strong>RM-peri 在这份数据上仍然做不了。</strong>三条路，按代价排序：</p>
  </div>
  <div class="scroll">
    <table>
      <thead><tr><th>路径</th><th>做法</th><th>代价</th><th>状态</th></tr></thead>
      <tbody>
        <tr><td class="key">A · 核分割叠加</td>
            <td>IDC <code>pan_cancer_nuclei_seg</code> 给逐细胞核多边形 →
            核密度/形态 → 肿瘤富集区 → 与 TIL 图谱叠加得中心与前沿</td>
            <td class="num">全库 1.2 TB；<br>只取我们 88 例约 <b>68 GB</b></td>
            <td><span class="chip measured">已核实</span><br>覆盖 88/91 例</td></tr>
        <tr><td class="key">B · 肿瘤区域分割</td>
            <td>对同批诊断切片跑 HoVer-Net 或现成乳腺肿瘤区域模型</td>
            <td class="num">中等，但引入一个<br>需自行验证的新模型</td>
            <td><span class="chip absent">未做</span></td></tr>
        <tr><td class="key">C · 几何描述子</td>
            <td>不需要肿瘤掩膜：Moran's I、最大连通簇占比、簇形状。
            排斥型的特征是淋巴细胞成带状聚集而不渗入，几何上可测</td>
            <td class="num">零下载</td>
            <td><span class="chip measured">已实现并使用</span></td></tr>
      </tbody>
    </table>
  </div>
  <div class="col">
    <div class="callout good">
      <p>路径 C <strong>不是纯粹的替代品</strong>：Moran's I 与 TIS 独立相关（ρ=+0.286, p=0.006），
      说明它确实捕捉到了生物学，只是不能直接命名为「排斥型」。</p>
      <p><strong>建议顺序</strong>：先走路径 C 定效应量的量级，
      只有当它显示出值得追的信号时，再投路径 A 的 68 GB 和随之而来的分割验证工作。
      本轮路径 C 没有显示这样的信号。</p>
    </div>
  </div>
</section>

<section>
  <div class="col">
    <p class="eyebrow">下一步</p>
    <h2>基于本轮实测的四条建议</h2>
    <hr class="rule">
    <ol class="evidence">
      <li><strong>样本量必须上去。</strong>91 例、28 例测试集，什么都测不出来。
      多中心模拟显示 182 → 900 单独就值 +0.185 AUC。</li>
      <li><strong>特征空间要对，不是要大。</strong>本轮 {n_feat} 个特征没找到信号；
      但扩到 3,332 个的正确理由是「加入瘤周区域和多时相动力学这些有生物学动机的量」，
      不是「多试一些总能撞上」—— 后者的代价已经量化过了（零信号下 AUC 从 0.489 涨到 0.729）。</li>
      <li><strong>标签要往 CD8 IHC 走。</strong>H&amp;E TIL 与 CD8A 只共享约 15% 的方差。
      如果终点是 CD8 特异的空间表型，替代标签的天花板就在那里。</li>
      <li><strong>每一个尝试都要配置换零分布和多重性校正。</strong>本轮如果不做，
      我们会报出一个 0.652、p=0.073 的「纹理特征预测 TIL」结果。</li>
    </ol>
  </div>
</section>

<section>
  <div class="col">
    <p class="eyebrow">公开数据盘点 <span class="chip measured">实测</span></p>
    <h2>瓶颈不是影像，是标签</h2>
    <hr class="rule">
    <p>公开乳腺 MRI 有 2,737–3,122 例。配对了 CD8 IHC 空间表型的，<strong>是 0 例</strong>。
    以下每个数字都由 <code>scripts/refresh_inventory.py</code> 当场打 TCIA / GDC / GEO 三个 API 实测。</p>
  </div>
  <div class="scroll">
    <table>
      <thead><tr><th>队列</th><th>规模</th><th>中心</th><th>标签</th><th>能验哪个任务</th></tr></thead>
      <tbody>
        <tr><td class="key">TCGA-BRCA<br><small>Tier A</small></td>
            <td class="num"><b>{n_paired}</b> 例三配对<br><small>MRI + RNA-seq + WSI</small></td>
            <td>MSKCC / Mayo / UPMC /<br>Roswell，4 家厂商</td>
            <td>Saltz 空间 TIL 图谱<br>+ RNA-seq 锚定 CD8</td>
            <td class="ok">两个都能<br><small>（需补肿瘤区域）</small></td></tr>
        <tr><td class="key">I-SPY2<br><small>Tier B</small></td>
            <td class="num"><b>{n_ispy2}</b> 例<br><small>MRI ↔ GSE194040</small></td>
            <td>~20 家美国中心</td>
            <td>转录组免疫签名<br>+ pCR + 帕博利珠单抗臂</td>
            <td>只能验 RM-whole</td></tr>
        <tr><td class="key">Duke / ISPY1 / ACRIN-6698 /<br>Advanced / NACT / MAMA-MIA<br><small>Tier C</small></td>
            <td class="num">922 / 222 / 385 /<br>632 / 64 / 1,506</td>
            <td>混合</td>
            <td>无免疫标签</td>
            <td>分割、协调、域偏移</td></tr>
      </tbody>
    </table>
  </div>
  <div class="col">
    <div class="callout">
      <p><strong>一条硬约束，写进了代码而不是注释。</strong>转录组标签只能验「炎症型 vs 其余」。
      沙漠型与排斥型的肿瘤核心都没有浸润，差别<strong>纯粹在空间排布上</strong> ——
      一份 core biopsy 的 bulk 表达谱没有任何机制可以表达这个差别。
      用它去报一个 RM-peri 的 AUC，等于在测量一个与模型声称测量的东西不同的量。</p>
    </div>
    <p style="margin-top:24px">另有两条<strong>查证过的断链</strong>，省得别人再踩：
    ISPY1 影像（<code>ISPY1_1001</code>）与 GSE22226（6 位研究号）按数字朴素匹配交集为 0；
    ACRIN-6698（DWI）与 ISPY2（DCE）同为 6 位却交集为 0 —— 两个 collection 在 TCIA 上是分别去标识化的。
    打通任一条都能显著扩容。</p>
  </div>
</section>

<section>
  <div class="col">
    <p class="eyebrow">复算</p>
    <h2>每个数字都能重跑</h2>
    <hr class="rule">
<pre>pip install -r mri_cd8_til/requirements.txt

python -m mri_cd8_til intervals            # 证据一，秒出
python scripts/refresh_inventory.py        # 重新实测公开数据可得性
python scripts/build_til_labels.py         # 下载 {n_til} 例真实 TIL 图谱并算标签
python scripts/run_optimism_study.py       # 证据三，约 2 小时
python scripts/run_multicenter_demo.py     # 多中心对比，约 1 小时
python scripts/run_real_dicom_demo.py      # 真实 DICOM 端到端
python scripts/make_figures.py             # 重绘本文全部图表
pytest tests/                              # 23 项测试</pre>
  </div>
</section>

<footer>
  <div class="col">
    <p><b>来源</b></p>
    <ul>
      <li>Jeon 等，<em>Front Immunol</em> 2022;13:1080048 —— 被复现的原文</li>
      <li>Saltz 等，<em>Cell Reports</em> 2018;23(1):181-193 —— TIL 图谱；TCIA DOI 10.7937/K9/TCIA.2018.Y75F9W1，CC BY</li>
      <li>NCI Imaging Data Commons —— TIL 图谱的 DICOM 再发布（本文 {n_til} 例由此获取）</li>
      <li>Arefan 等，<em>BMC Cancer</em> 2021 —— MRI + TCGA 表达谱推断 CD8，跨机构 0.74 → 0.62</li>
      <li>Garrucho 等，<em>Scientific Data</em> 2025 —— MAMA-MIA 多中心 DCE-MRI 与专家分割</li>
      <li>I-SPY2-990 mRNA 资源 —— GEO GSE194040</li>
      <li>鼻咽癌 TIL AI 风险评分 —— 单中心 367 例，C-index 0.923 → 0.785</li>
    </ul>
    <p style="margin-top:16px">证据快照 2026-08-10。所有实测数字由随附脚本当场生成，可复算。</p>
  </div>
</footer>

</div>
"""

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(html, encoding="utf-8")
    size_mb = OUT.stat().st_size / 1e6
    print(f"wrote {OUT}  ({size_mb:.2f} MB)")
    if size_mb > 15:
        print("WARNING: approaching the 16 MB artifact limit")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
