#!/usr/bin/env python3
"""Chinese decision memo for the RFP feasibility review. Numbers from results/."""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import results_path  # noqa: E402


def jsn(n):
    p = results_path(n)
    return json.load(open(p)) if os.path.exists(p) else {}


def esc(s):
    return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


CSS = """
:root{
  --paper:#FCFCFD;--ground:#F1F5F9;--surface:#FFF;--line:#D9E1EA;
  --ink:#151E28;--ink-2:#3A4A5B;--muted:#5E7185;
  --accent:#1F6F8B;--accent-soft:#E6F1F5;
  --no:#A93B32;--no-soft:#FBEDEC;--risk:#9E6516;--risk-soft:#FAF2E5;
  --yes:#2C7A5C;--yes-soft:#E9F4EF;
  --serif:"Noto Serif SC",Songti SC,SimSun,Georgia,serif;
  --sans:"Noto Sans SC",-apple-system,"PingFang SC","Hiragino Sans GB","Microsoft YaHei",sans-serif;
  --mono:"IBM Plex Mono",ui-monospace,Menlo,monospace;
}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){
  --paper:#0E141A;--ground:#161F27;--surface:#18222B;--line:#2A3742;
  --ink:#E3EAF1;--ink-2:#B4C4D1;--muted:#8598AA;
  --accent:#63A8C4;--accent-soft:#152730;
  --no:#DE8178;--no-soft:#2A1918;--risk:#D39E55;--risk-soft:#291F14;
  --yes:#6DB795;--yes-soft:#14251E;}}
:root[data-theme="dark"]{
  --paper:#0E141A;--ground:#161F27;--surface:#18222B;--line:#2A3742;
  --ink:#E3EAF1;--ink-2:#B4C4D1;--muted:#8598AA;
  --accent:#63A8C4;--accent-soft:#152730;
  --no:#DE8178;--no-soft:#2A1918;--risk:#D39E55;--risk-soft:#291F14;
  --yes:#6DB795;--yes-soft:#14251E;}
*{box-sizing:border-box}
body{background:var(--paper);color:var(--ink);font-family:var(--sans);
  font-size:16px;line-height:1.85;margin:0}
.wrap{max-width:58rem;margin:0 auto;padding:0 20px;padding-block:0 4.5rem}
.col{max-width:36rem}
h1,h2,h3{font-family:var(--serif);font-weight:700;text-wrap:balance;line-height:1.35}
h1{font-size:clamp(1.7rem,4.2vw,2.4rem);margin:0 0 .5rem}
h2{font-size:clamp(1.25rem,2.6vw,1.5rem);margin:0 0 .9rem}
h3{font-size:1.04rem;margin:1.8rem 0 .5rem}
p{margin:0 0 1rem}
code{font-family:var(--mono);font-size:.86em;background:var(--ground);
  padding:.1em .35em;border-radius:3px}
.kicker{font-size:.72rem;font-weight:600;letter-spacing:.18em;text-transform:uppercase;
  color:var(--accent);margin:0 0 .8rem}
header{border-bottom:2px solid var(--ink);padding-block:3rem 1.5rem;margin-bottom:2rem}
.lede{font-size:1.1rem;color:var(--ink-2)}
section{padding-block:2.2rem 0}
hr.rule{border:0;border-top:1px solid var(--line);margin:2.6rem 0 0}
.buckets{display:grid;grid-template-columns:repeat(auto-fit,minmax(170px,1fr));
  gap:1px;background:var(--line);border:1px solid var(--line);margin:1.5rem 0}
.bk{background:var(--surface);padding:1.05rem 1.1rem;min-width:0}
.bk .n{font-family:var(--serif);font-size:2rem;font-weight:700;line-height:1}
.bk .l{font-size:.74rem;font-weight:600;letter-spacing:.09em;text-transform:uppercase;
  color:var(--muted);margin:.45rem 0 .3rem}
.bk .s{font-size:.83rem;color:var(--ink-2);line-height:1.6}
.note{border-left:3px solid var(--accent);background:var(--accent-soft);
  padding:1rem 1.15rem;margin:1.4rem 0;font-size:.95rem;line-height:1.78}
.note.no{border-color:var(--no);background:var(--no-soft)}
.note.risk{border-color:var(--risk);background:var(--risk-soft)}
.note.yes{border-color:var(--yes);background:var(--yes-soft)}
.note p:last-child{margin-bottom:0}
.note .t{font-weight:600;font-family:var(--serif)}
.tw{overflow-x:auto;margin:1.2rem 0 1.5rem;border:1px solid var(--line)}
table{border-collapse:collapse;width:100%;font-size:.85rem;background:var(--surface)}
th,td{text-align:left;padding:.54rem .7rem;border-bottom:1px solid var(--line);
  vertical-align:top}
th{background:var(--ground);font-weight:600;font-size:.73rem;letter-spacing:.05em;
  color:var(--muted);white-space:nowrap}
tbody tr:last-child td{border-bottom:0}
.num td:nth-child(n+2),.num th:nth-child(n+2){text-align:right;
  font-variant-numeric:tabular-nums;font-family:var(--mono);font-size:.82rem}
td.l,th.l{text-align:left!important;font-family:var(--sans)!important}
.tag{display:inline-block;font-size:.7rem;font-weight:600;padding:.08rem .45rem;
  border-radius:2px;white-space:nowrap}
.tag.no{background:var(--no-soft);color:var(--no)}
.tag.risk{background:var(--risk-soft);color:var(--risk)}
.tag.yes{background:var(--yes-soft);color:var(--yes)}
.item{border-top:1px solid var(--line);padding-block:1.4rem;
  display:grid;grid-template-columns:5.4rem 1fr;gap:0 1.2rem}
.item:last-of-type{border-bottom:1px solid var(--line)}
.item>div{min-width:0}
.item h3{margin:0 0 .3rem;font-size:1rem}
.item p{margin:0 0 .5rem;font-size:.92rem}
.item .can{color:var(--ink-2)}
.item .cant{color:var(--muted);font-size:.88rem}
@media(max-width:600px){.item{grid-template-columns:1fr;gap:.4rem}}
ul{margin:0 0 1rem;padding-left:1.15rem}li{margin-bottom:.4rem}
footer{margin-top:3rem;padding-top:1.3rem;border-top:1px solid var(--line);
  font-size:.8rem;color:var(--muted);line-height:1.7}
@media (prefers-reduced-motion:reduce){*{animation:none!important}}
"""

FONTS = ('<link rel="preconnect" href="https://fonts.googleapis.com">'
         '<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>'
         '<link rel="stylesheet" href="https://fonts.googleapis.com/css2?'
         'family=IBM+Plex+Mono:wght@400;500&family=Noto+Sans+SC:wght@400;500;600&'
         'family=Noto+Serif+SC:wght@600;700&display=swap">')


# Chinese rendering of the C7 matrix. C7 itself stays English because it is the
# machine-readable record; only the wording is localised here, keyed on the
# stable ids C7 emits, so the numbers still come from the JSON.
ITEM_ZH = {
 "step1_il2": ("第 1 步 · MHC-II 依赖的 IL-2 读数（DO11.10 / I-A\u1d48）",
   "复现了 cOVA323-339 在 I-A\u1d48 上的文献 register（核心 {core}，%Rank {rank}），"
   "确认杂交瘤体系的限制性元件行为与文献一致。",
   "IL-2 分泌是细胞层面的功能读数。没有任何序列模型能预测细胞因子输出。"
   "结合是必要条件，不是充分条件。"),
 "step1_panel": ("第 1 步 · 五条肽的组合（8-mer、9-mer+4D、17-mer、两条 30-mer）",
   "六个成员里评分了 {scored} 个。17-mer uniblock 与 30-mer diblock 在同一扫描窗下"
   "返回相同的核心和 %Rank\u2014\u2014即接上自组装嵌段这件事，序列模型完全看不见。",
   "{floor} 条低于类 II 接口的 11 残基下限，根本无法提交；"
   "{dee} 条含 D-氨基酸，没有任何类 II 预测器能表示立体化学——"
   "D-取代肽会被当成它的全-L 孪生体来打分，那是另一个分子。"
   "而第 1 步的核心比较正是全-L 对照 vs D-取代试验，所以这个比较计算做不了。"),
 "step2_binding": ("第 2 步 · 无细胞 peptide-MHC 结合 / 稳定性筛选",
   "针对这两个等位基因本身做了实测，对手是 IEDB 上它们的全部有标签记录："
   "{auc}。这把预筛从一个主张变成了一条校准过的规则：{skip}。",
   "预测的是结合，不是复合物稳定性或解离速率。RFP 要的是稳定性——"
   "那才是决定四聚体能否扛过染色的性质，而没有序列模型会报半衰期。"),
 "step2_tetramer": ("第 2 步 · 定制 2 个等位基因 + 8 个 uniblock 复合物",
   "对全部 {poss} 个可能的（肽，等位基因）复合物排了序，判定 {build} 个值得做。"
   "不折叠的单体永远不会变成四聚体，而这件事通常在钱花完之后才发现。",
   "无法确认重折叠、生物素化或多聚化产率。预测会结合的肽，仍然可能折不起来。"),
 "step2_precursor": ("第 2 步 · IL-7/IL-15 扩增后用 FACS 测基线前体频率",
   "什么也做不了。序列层面没有任何替身。",
   "前体频率是某一个供体 T 细胞库的属性，由他本人的胸腺选择和暴露史决定。"
   "它只能被测量，不能被预测。这是整份 RFP 里最不可替代的一条。"),
 "step3_mlr": ("第 3 步 · MLR（增殖、细胞因子、ICS）",
   "什么也做不了。",
   "增殖和细胞因子输出是加工、呈递、共刺激和 TCR 库共同决定的细胞层面结果。"
   "序列模型只看得见第一步。"),
 "donors": ("供体 · 携带指定 HLA 抗原的人 PBMC",
   "美欧加权携带率：{a} {aw}，{b} {bw}，有其一 {ew}，两个都有 {bothw}。"
   "要凑 10 名供体：有其一约需分型 {n_either} 人，两个都要约需 {n_both} 人。",
   "无法告诉你某个商业 HLA 分型 PBMC 库里现成有什么；"
   "两位点按独立相乘也违反连锁不平衡的事实。"),
}
EV_ZH = {"measured (C6)": "实测（C6）", "measured (C1, C3)": "实测（C1、C3）",
         "measured (C4)": "实测（C4）", "measured (C5)": "实测（C5）",
         "reasoning, not measurement": "推理，非实测"}

ZH = {"REPLACE": "可替代", "DE-RISK": "可降风险", "NO": "不可替代"}
CLS = {"REPLACE": "yes", "DE-RISK": "risk", "NO": "no"}


def main():
    feas = jsn("c7_feasibility.json")
    cal = jsn("c3_calibration.json")
    pre = jsn("c4_prescreen.json")
    donor = jsn("c5_donor_feasibility.json")
    cls2 = jsn("c6_nanopeptide_classII.json")
    bench = jsn("c1_benchmark_summary.json")

    g = feas["group_counts"]
    A, B = "HLA-A*32:01", "HLA-B*57:01"
    ca, cb = cal["per_allele"][A], cal["per_allele"][B]
    w = donor["weighted_us_eu"]
    bur = donor["donors_to_screen"]
    v = cls2["verdict"]
    pc = cls2["positive_control"]["result"]

    H = []
    A_ = H.append

    A_(f'''<header><div class="wrap">
<p class="kicker">RFP 可行性审阅 · 体外/体内验证 vs 计算方法 · RUO</p>
<h1>这份 RFP 里，哪些能用计算做，哪些不能</h1>
<p class="lede col">把 RFP 拆成 7 个条目，逐条判定。结论先给：
<strong>没有任何一条能被计算替代</strong>；有 {g['DE-RISK']} 条可以在下单之前被显著降风险，
而且其中一条可能会改变研究设计本身。</p>
</div></header><main class="wrap">''')

    # ---- buckets
    A_(f'''<section>
<div class="buckets">
<div class="bk"><div class="n" style="color:var(--yes)">{g['REPLACE']}</div>
  <div class="l">可替代</div>
  <div class="s">计算直接给出答案、实验不必做。<strong>这一栏是空的</strong>&mdash;&mdash;
  这是本审阅最有用的一句话</div></div>
<div class="bk"><div class="n" style="color:var(--risk)">{g['DE-RISK']}</div>
  <div class="l">可降风险</div>
  <div class="s">答不了，但能改变你下单的内容、数量，或者判断这件事根本做不做得成</div></div>
<div class="bk"><div class="n" style="color:var(--no)">{g['NO']}</div>
  <div class="l">不可替代</div>
  <div class="s">物理测量，序列层面没有替身；或者分子本身超出模型能表示的范围</div></div>
</div>
<p class="col">下面每一条都带证据出处。凡是没有实测支撑的判断，标注为"推理"而非"实测"。</p>
</section><hr class="rule">''')

    # ---- item by item
    A_('<section><h2>一、逐条判定</h2>')
    fmt = {
        "core": pc["core"], "rank": pc["rank"],
        "scored": v["scored"], "floor": v["blocked_by_length_floor"],
        "dee": v["blocked_by_d_amino_acids"],
        "auc": f"{A} AUC {ca['EL']['auc']:.3f}；{B} AUC {cb['EL']['auc']:.3f}",
        "skip": (f"{A} 在 %Rank ≥ {ca['recommended_skip_cut']} 跳过"
                 f"（NPV {ca['recommended_skip_npv']:.2f}）；"
                 f"{B} 在 %Rank ≥ {cb['recommended_skip_cut']} 跳过"
                 f"（NPV {cb['recommended_skip_npv']:.2f}）"),
        "poss": pre["complexes_possible"], "build": pre["n_build"],
        "a": A, "b": B, "aw": f"{w[A]*100:.2f}%", "bw": f"{w[B]*100:.2f}%",
        "ew": f"{w['either']*100:.2f}%", "bothw": f"{w['both_independent']*100:.2f}%",
        "n_either": f"{bur['either']['10']:.0f}", "n_both": f"{bur['both_independent']['10']:.0f}",
    }
    for it in feas["items"]:
        cls = CLS[it["group"]]
        title, can, cant = ITEM_ZH[it["id"]]
        A_(f'''<div class="item">
<div><span class="tag {cls}">{ZH[it['group']]}</span>
  <div style="font-size:.72rem;color:var(--muted);margin-top:.4rem;font-family:var(--mono)">
  {EV_ZH.get(it['evidence_type'], esc(it['evidence_type']))}</div></div>
<div><h3>{title}</h3>
<p class="can"><strong>计算做到了：</strong>{can.format(**fmt)}</p>
<p class="cant"><strong>做不到：</strong>{cant.format(**fmt)}</p></div></div>''')
    A_('</section><hr class="rule">')

    # ---- the three de-risk items in detail
    A_('<section><h2>二、三件真能降风险的事</h2>')

    A_('<h3>1 · 供体可行性&mdash;&mdash;这条最可能改变研究设计</h3>')
    A_(f'''<p class="col">RFP 要求"携带指定 HLA 抗原的人 PBMC"，但没说是<strong>两个等位基因都要</strong>
还是<strong>有其一即可</strong>。这个区别是数量级的：</p>''')
    A_(f'''<div class="tw"><table class="num">
<thead><tr><th class="l">要求</th><th>美欧加权携带率</th>
<th>要 5 人需分型</th><th>要 10 人</th><th>要 20 人</th></tr></thead><tbody>
<tr><td class="l">{A}</td><td>{w[A]*100:.2f}%</td>
  <td>{bur[A]['5']:.0f}</td><td>{bur[A]['10']:.0f}</td><td>{bur[A]['20']:.0f}</td></tr>
<tr><td class="l">{B}</td><td>{w[B]*100:.2f}%</td>
  <td>{bur[B]['5']:.0f}</td><td>{bur[B]['10']:.0f}</td><td>{bur[B]['20']:.0f}</td></tr>
<tr><td class="l">有其一即可</td><td>{w['either']*100:.2f}%</td>
  <td>{bur['either']['5']:.0f}</td><td>{bur['either']['10']:.0f}</td>
  <td>{bur['either']['20']:.0f}</td></tr>
<tr><td class="l"><strong>两个都要</strong></td><td>{w['both_independent']*100:.2f}%</td>
  <td>{bur['both_independent']['5']:.0f}</td>
  <td><strong>{bur['both_independent']['10']:.0f}</strong></td>
  <td>{bur['both_independent']['20']:.0f}</td></tr>
</tbody></table></div>''')
    A_(f'''<div class="note risk col"><p><span class="t">先把这个问题问清楚，再谈报价。</span>
要凑齐 10 名供体：<strong>有其一即可 → 约筛 {bur['either']['10']:.0f} 人</strong>；
<strong>两个都要 → 约筛 {bur['both_independent']['10']:.0f} 人</strong>。
后者基本不是一个能在 30 天内完成的研究。</p>
<p style="font-size:.88rem;color:var(--muted);margin-top:.6rem">
两位点按独立相乘是近似。HLA-A 与 -B 相距约 1.3&nbsp;Mb、存在强连锁不平衡，
而欧洲人群的 B*57:01 主要搭载在 57.1 祖先单倍型上（A*01:01-B*57:01-C*06:02-DRB1*07:01），
带的是 A*01:01 而不是 A*32:01&mdash;&mdash;所以对<em>这一对</em>而言，
独立估计更可能是<strong>高估</strong>，真实情况比上表还要难。正式报价前应查单倍型表。</p></div>''')

    A_('<h3>2 · 四聚体该做哪几个&mdash;&mdash;钱花在这里</h3>')
    A_(f'''<p class="col">RFP 要定制 2 个 HLA 等位基因 + <strong>8 个 uniblock 复合物</strong>。
不折叠的单体不会变成四聚体，而这件事通常是在钱和时间都花掉之后才发现的。
5 条候选肽 × 2 个等位基因 = {pre['complexes_possible']} 个可能复合物，预筛结果：
<strong>{pre['n_build']} 个值得做，{pre['n_skip']} 个建议跳过</strong>。</p>''')

    rows = pre["rows"]
    A_('<div class="tw"><table class="num"><thead><tr>'
       '<th class="l">候选</th><th class="l">肽</th><th class="l">等位基因</th>'
       '<th>EL %Rank</th><th>IC50 nM</th><th class="l">判定</th>'
       '<th class="l">已知真值</th></tr></thead><tbody>')
    for r in rows:
        truth = {1: '<span class="tag yes">结合</span>', 0: '<span class="tag no">不结合</span>',
                 None: '<span style="color:var(--muted)">&mdash;</span>'}[r["known_label"]]
        call = ('<span class="tag risk">做</span>' if r["call"] == "BUILD"
                else '<span class="tag" style="background:var(--ground);color:var(--muted)">跳过</span>')
        A_(f'<tr><td class="l">{esc(r["candidate"])}</td>'
           f'<td class="l"><code>{esc(r["peptide"])}</code></td>'
           f'<td class="l">{esc(r["allele"])}</td>'
           f'<td>{r["el_rank"]:.2f}</td><td>{r["ic50_nM"]:.0f}</td>'
           f'<td class="l">{call}</td><td class="l">{truth}</td></tr>')
    A_('</tbody></table></div>')
    ag = pre["agreement_with_known_labels"]
    A_(f'''<p class="col">这 5 条候选肽是<strong>从校准集里留出的基准肽</strong>（真值已知、
且没有参与阈值设定），所以这张表是在演示预筛做真实判断而不是自我验证。
在 {ag['n_checkable']} 个真值可查的复合物上，判定与 IEDB 标签一致 {ag['n_agree']} 个。
唯一一处不一致（<code>RPKPDYSAM</code> × {B}）是规则<em>按设计</em>发生的：
漏掉一个真表位无法挽回，多做一个只是钱，所以规则偏向多做。</p>''')
    A_(f'''<div class="note col"><p><span class="t">顺带一个对 RFP 第 2 步第一段有直接影响的观察。</span>
经典的 <strong>IC50 &lt; 500 nM</strong> 亲和力标准在这 {pre['complexes_possible']} 个复合物上只通过
{pre['n_passing_ic50_standard']} 个，比 EL 判定的少。更要紧的是：
<code>FSPEVIPMF</code>&mdash;&mdash;一条有文献记载的 B*57:01 表位&mdash;&mdash;
EL %Rank 0.47（强结合）但 BA IC50 <strong>3017 nM</strong>，比 500 nM 标准差了六倍。
如果第一步的无细胞筛选用 IC50&lt;500 nM 做门槛，它会把真表位筛掉。两个指标都报，不要让亲和力标准单独否决。</p></div>''')

    A_('<h3>3 · 预筛本身可不可信&mdash;&mdash;在这两个等位基因上实测</h3>')
    A_(f'''<p class="col">类 I 预测"平均而言"很好，而平均值是被 A*02:01 那类研究透彻的等位基因撑起来的。
A*32:01 不属于那一类。所以在允许预筛说"别做这个"之前，先把 IEDB 上这两个等位基因的
全部有标签记录（{bench['raw_records']:,} 条原始记录 &rarr; {bench['labelled_pairs']:,} 对）拿来实测：</p>''')
    A_(f'''<div class="tw"><table class="num">
<thead><tr><th class="l">等位基因</th><th>校准 n</th><th>阳性占比</th><th>ROC AUC</th>
<th class="l">95% CI</th><th class="l">建议跳过阈值</th><th>能跳过</th></tr></thead><tbody>
<tr><td class="l">{A}</td><td>{ca['n_calibration']:,}</td>
  <td>{ca['positive_fraction']*100:.0f}%</td><td>{ca['EL']['auc']:.3f}</td>
  <td class="l">[{ca['EL']['auc_ci95'][0]:.3f}, {ca['EL']['auc_ci95'][1]:.3f}]</td>
  <td class="l">%Rank &ge; {ca['recommended_skip_cut']}</td>
  <td>{ca['expected_skip_fraction_at_assumed_prev']*100:.0f}%</td></tr>
<tr><td class="l">{B}</td><td>{cb['n_calibration']:,}</td>
  <td>{cb['positive_fraction']*100:.0f}%</td><td>{cb['EL']['auc']:.3f}</td>
  <td class="l">[{cb['EL']['auc_ci95'][0]:.3f}, {cb['EL']['auc_ci95'][1]:.3f}]</td>
  <td class="l">%Rank &ge; {cb['recommended_skip_cut']}</td>
  <td>{cb['expected_skip_fraction_at_assumed_prev']*100:.0f}%</td></tr>
</tbody></table></div>''')
    A_(f'''<p class="col">两个等位基因的判别力都很好（AUC 约 0.91&ndash;0.92），
包括数据量少得多的 A*32:01&mdash;&mdash;这是预筛可以被信任的实测依据。
但两者的<strong>实用价值差别很大</strong>：在"保留 95% 真结合肽"这个安全线下，
A*32:01 可以安全跳过约 {ca['expected_skip_fraction_at_assumed_prev']*100:.0f}% 的候选，
而 B*57:01 只能跳过约 {cb['expected_skip_fraction_at_assumed_prev']*100:.0f}%。
对 B*57:01，预筛的省钱能力有限，该做的还是得做。</p>''')
    A_(f'''<div class="note col"><p><span class="t">一个必须说明的统计陷阱。</span>
IEDB 的类 I 数据绝大部分是洗脱配体&mdash;&mdash;按构造就是阳性，
所以 A*32:01 回来是 <strong>{ca['positive_fraction']*100:.0f}% 阳性</strong>。
在这个比例下算出来的 PPV/NPV 描述的是 IEDB 的收录政策，不是你的决策。
上表的 PPV/NPV 一律按<strong>声明的先验</strong>（{cal['assumed_prevalence']:.0%}，
即"为某个等位基因设计的候选肽真的结合"的概率）重算，并做了敏感性扫描。
另外 A*32:01 的 AUC 只靠 {ca['n_negative']} 条阴性撑着，置信区间明显更宽&mdash;&mdash;
不要把数据贫乏的等位基因和数据充足的等位基因当成同样可信。</p></div>''')
    A_('</section><hr class="rule">')

    # ---- step 1
    A_('<section><h2>三、第 1 步（纳米肽）：计算基本插不上手</h2>')
    A_(f'''<p class="col">先说好消息：<strong>体系本身的对照跑通了</strong>。
DO11.10 杂交瘤的限制性元件是小鼠 I-A<sup>d</sup>，参考表位是 cOVA<sub>323-339</sub>。
把它送进 NetMHCIIpan，返回的结合核心是 <code>{esc(pc['core'])}</code>、
%Rank {pc['rank']}&mdash;&mdash;<strong>与文献 register 一致</strong>。
体系设定没有问题，这一条本身值得单独发一封邮件确认。</p>''')
    A_(f'''<p class="col">但五条肽里，能被评分的只有 <strong>{v['scored']} 条</strong>：</p>''')
    A_('<div class="tw"><table class="num"><thead><tr>'
       '<th class="l">肽</th><th>长度</th><th class="l">结果</th></tr></thead><tbody>')
    for it in cls2["panel"]:
        if it.get("result"):
            res = (f'核心 <code>{esc(it["result"]["core"])}</code>，'
                   f'%Rank {it["result"]["rank"]}')
        elif it.get("api_refusal"):
            res = '<span class="tag no">API 拒收</span> 低于类 II 的 11 残基下限'
        else:
            res = '<span class="tag no">排除</span> 含 D-氨基酸'
        A_(f'<tr><td class="l"><code>{esc(it["name"])}</code></td>'
           f'<td>{it["length"]}</td><td class="l">{res}</td></tr>')
    A_('</tbody></table></div>')
    A_(f'''<div class="note no col"><p><span class="t">D-氨基酸不是"注意事项"，是拒绝评分的理由。</span>
在用的每一个 MHC 预测器都只在 L-氨基酸肽上训练，输入里根本没有立体化学这一维。
把含 4 个 D-残基的肽送进去，返回的是<strong>它的全-L 孪生体</strong>的分数，那是另一个分子。
而第 1 步的核心比较恰恰是"全-L 对照 vs D-取代试验"&mdash;&mdash;
<strong>这个比较计算做不了</strong>，只能表征 L 那一臂。</p></div>''')
    sw = cls2["padding_sweep"]
    A_(f'''<div class="note col"><p><span class="t">顺带堵掉一个看起来很聪明的捷径。</span>
8-mer 和裸 9-mer 核心低于类 II 接口的 11 残基下限，根本提交不了。
有人会想"那就补几个甘氨酸凑长度"。实测：把 <code>{esc(sw['core'])}</code>
补到 11/13/15/17，核心始终不变，但 %Rank 从 {sw['rank_range'][0]} 变到 {sw['rank_range'][1]}
&mdash;&mdash;<strong>{sw['fold_range']} 倍的摆动，完全来自补了多少个甘氨酸</strong>。
这不是绕过限制，是把一个自己造的假象当数据。类 II 的槽是两端开放的，
侧翼残基本来就参与结合&mdash;&mdash;裸核心和有上下文的同一核心是两个分子，
这大概也正是 RFP 把 8-mer 放进来当对照的原因。</p></div>''')
    A_(f'''<p class="col">还有一点：17-mer uniblock 和 30-mer diblock 在同一个 15-mer 扫描窗下
返回<strong>完全相同的核心和 %Rank</strong>。也就是说，接上自组装嵌段这件事，
序列模型完全看不见&mdash;&mdash;而那恰恰是 diblock 存在的理由。</p>''')
    A_('</section><hr class="rule">')

    # ---- recommendation
    A_('<section><h2>四、建议</h2>')
    A_(f'''<div class="note yes col"><p><span class="t">这项研究该做，但可以做得更小、更准。</span>
计算替代不了其中任何一条，但在下单之前跑一遍上面这套（几小时，无试剂成本），
可以：把四聚体合成清单从 {pre['complexes_possible']} 个收敛到 {pre['n_build']} 个、
确认体系对照成立、并且把供体可行性这个可能致命的问题提前暴露出来。</p></div>''')
    A_('''<h3>下单之前需要 RFP 方回答的两个问题</h3>
<ol class="col">
<li><strong>供体是要同时携带两个等位基因，还是有其一即可？</strong>
这决定了筛查量是几十人还是几千人，也决定了 30 天时限是否现实。
这是本次审阅里价值最高的一个问题&mdash;&mdash;它比任何一个预测都重要。</li>
<li><strong>第 1 步的 D-取代对照，打算怎么在计算上处理？</strong>
建议的答案是：不处理，明确写进假设，由湿实验独占。
任何声称能预测 D-取代肽 MHC 结合的方案都应当被要求出示其立体化学表示方式。</li>
</ol>''')
    A_('''<h3>可以直接照搬的两点方法建议</h3>
<ul class="col">
<li>第 2 步第一段的无细胞筛选，<strong>不要只用 IC50&lt;500 nM 做门槛</strong>，
洗脱配体评分和亲和力两个指标都记录。上面已给出一个反例。</li>
<li>四聚体清单按"<strong>不确定就做</strong>"的原则定：漏掉一个真表位没人会回头补，
多做一个只是钱。上面的阈值就是按这个不对称性选的。</li>
</ul>''')
    A_('</section>')

    A_(f'''<footer class="col">
<p>本文所有数字由 <code>results/</code> 下的输出生成。第 2 步的候选肽为<strong>构造的替身</strong>
（从校准集留出的基准肽），真实序列由 RFP 方提供后原样替换即可，流程不变。
第 1 步的 17-mer/30-mer 同为构造替身，以 cOVA<sub>323-339</sub> 为基础，已在文中标明。</p>
<p>预测工具：NetMHCpan-4.1（EL + BA）、NetMHCIIpan，均经 IEDB API；
有标签数据来自 IEDB query API；等位基因频率来自 IEDB population-coverage 表。</p>
<p><strong>仅供研究使用（RUO）。本文是一份下单前的可行性审阅，不是免疫原性结论。</strong></p>
</footer></main>''')

    html = (f"<title>RFP 可行性审阅</title>\n{FONTS}\n<style>{CSS}</style>\n"
            + "\n".join(H))
    out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "rfp_feasibility_zh.html")
    with open(out, "w") as f:
        f.write(html)
    print(f"wrote {out} ({len(html)//1024} KB)")


if __name__ == "__main__":
    main()
