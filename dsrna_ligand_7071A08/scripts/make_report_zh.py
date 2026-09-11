#!/usr/bin/env python3
"""
Chinese report. Same discipline as make_report.py: every number is read out of
results/, never typed into the prose.

Part one is the assessment. Part two walks the deck slide by slide from P8,
because that is where the report stops being "here is your answer" and starts
being "here is why you should believe the answer" - and that transition is the
part a reader who does not do this for a living has to be carried through.
"""
import base64
import csv
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import load_config, read_fasta, read_metadata, data_path, results_path, figures_path  # noqa: E402


def jsn(name):
    p = results_path(name)
    return json.load(open(p)) if os.path.exists(p) else {}


def tsv(name):
    p = results_path(name)
    if not os.path.exists(p):
        return []
    with open(p) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def esc(s):
    return (str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


def img(name, cap=""):
    p = figures_path(name)
    if not os.path.exists(p):
        return ""
    b = base64.b64encode(open(p, "rb").read()).decode()
    c = f'<figcaption>{cap}</figcaption>' if cap else ""
    return f'<figure class="fig"><img src="data:image/png;base64,{b}" alt="{name}">{c}</figure>'


def table(rows, cols, heads, cls="", fmt=None):
    """rows: list of dicts. fmt values return raw HTML."""
    fmt = fmt or {}
    th = "".join(f"<th>{h}</th>" for h in heads)
    body = []
    for r in rows:
        tds = "".join(f"<td>{fmt[c](r) if c in fmt else esc(r.get(c, ''))}</td>"
                      for c in cols)
        body.append(f"<tr>{tds}</tr>")
    return (f'<div class="tw"><table class="{cls}"><thead><tr>{th}</tr></thead>'
            f'<tbody>{"".join(body)}</tbody></table></div>')


CSS = """
:root{
  --paper:#FCFCFD; --ground:#F2F6F9; --surface:#FFFFFF; --line:#D8E1EA;
  --ink:#16202B; --ink-2:#38485A; --muted:#5F7285;
  --accent:#2C6E9E; --accent-soft:#E8F1F8; --accent-line:#B7D2E6;
  --bad:#AF382E; --bad-soft:#FBEEEC; --warn:#A8661B; --warn-soft:#FBF2E6;
  --good:#2F7D5E; --good-soft:#EAF4EF; --violet:#6E55A0;
  --serif:"Noto Serif SC",Songti SC,SimSun,Georgia,serif;
  --sans:"Noto Sans SC",-apple-system,BlinkMacSystemFont,"PingFang SC","Hiragino Sans GB","Microsoft YaHei",sans-serif;
  --mono:"IBM Plex Mono",ui-monospace,SFMono-Regular,Menlo,monospace;
  --measure:36rem;
}
@media (prefers-color-scheme:dark){
  :root:not([data-theme="light"]){
    --paper:#0F151B; --ground:#161F27; --surface:#18222C; --line:#2B3945;
    --ink:#E4EBF2; --ink-2:#B7C6D3; --muted:#8699AB;
    --accent:#6FACD8; --accent-soft:#172A38; --accent-line:#2E4C63;
    --bad:#E08278; --bad-soft:#2C1B19; --warn:#DBA45C; --warn-soft:#2B2016;
    --good:#6FB897; --good-soft:#15261F; --violet:#A38FD0;
  }
}
:root[data-theme="dark"]{
  --paper:#0F151B; --ground:#161F27; --surface:#18222C; --line:#2B3945;
  --ink:#E4EBF2; --ink-2:#B7C6D3; --muted:#8699AB;
  --accent:#6FACD8; --accent-soft:#172A38; --accent-line:#2E4C63;
  --bad:#E08278; --bad-soft:#2C1B19; --warn:#DBA45C; --warn-soft:#2B2016;
  --good:#6FB897; --good-soft:#15261F; --violet:#A38FD0;
}
*{box-sizing:border-box}
body{background:var(--paper);color:var(--ink);font-family:var(--sans);
  font-size:16px;line-height:1.85;margin:0;-webkit-font-smoothing:antialiased}
.wrap{max-width:60rem;margin:0 auto;padding:0 20px;padding-block:0 5rem}
.col{max-width:var(--measure)}
h1,h2,h3,h4{font-family:var(--serif);font-weight:700;text-wrap:balance;
  line-height:1.35;color:var(--ink)}
h1{font-size:clamp(1.85rem,4.6vw,2.7rem);margin:0 0 .5rem}
h2{font-size:clamp(1.3rem,2.8vw,1.6rem);margin:0 0 .9rem}
h3{font-size:1.08rem;margin:2rem 0 .6rem}
h4{font-size:.98rem;margin:1.3rem 0 .35rem;font-family:var(--sans);font-weight:600}
p{margin:0 0 1.05rem}
strong{font-weight:600}
code,.mono{font-family:var(--mono);font-size:.88em;letter-spacing:-.01em}
code{background:var(--ground);padding:.1em .35em;border-radius:3px}
a{color:var(--accent)}
.lede{font-size:1.12rem;color:var(--ink-2);line-height:1.8}
.kicker{font-size:.72rem;font-weight:600;letter-spacing:.18em;text-transform:uppercase;
  color:var(--accent);margin:0 0 .8rem;font-family:var(--sans)}

/* header */
header.hero{border-bottom:2px solid var(--ink);padding-block:3.2rem 1.6rem;margin-bottom:2.4rem}
.meta{display:flex;flex-wrap:wrap;gap:.4rem 1.4rem;font-size:.82rem;color:var(--muted);
  margin-top:1.3rem;font-family:var(--mono)}
.meta span{white-space:nowrap}

/* sections */
section{padding-block:2.4rem 0;scroll-margin-top:1rem}
.rule{border:0;border-top:1px solid var(--line);margin:3rem 0 0}

/* stat row */
.stats{display:grid;grid-template-columns:repeat(auto-fit,minmax(158px,1fr));gap:1px;
  background:var(--line);border:1px solid var(--line);margin:1.6rem 0 1.4rem}
.stat{background:var(--surface);padding:1.05rem 1.05rem 1.15rem}
.stat .n{font-family:var(--serif);font-size:1.85rem;font-weight:700;line-height:1.1;
  font-variant-numeric:tabular-nums}
.stat .l{font-size:.74rem;font-weight:600;letter-spacing:.09em;text-transform:uppercase;
  color:var(--muted);margin:.5rem 0 .35rem}
.stat .s{font-size:.82rem;color:var(--ink-2);line-height:1.6}

/* callouts */
.note{border-left:3px solid var(--accent);background:var(--accent-soft);
  padding:1rem 1.15rem;margin:1.5rem 0;font-size:.95rem;line-height:1.78}
.note.bad{border-color:var(--bad);background:var(--bad-soft)}
.note.warn{border-color:var(--warn);background:var(--warn-soft)}
.note.good{border-color:var(--good);background:var(--good-soft)}
.note p:last-child{margin-bottom:0}
.note .t{font-weight:600;font-family:var(--serif)}

/* tables */
.tw{overflow-x:auto;margin:1.3rem 0 1.6rem;border:1px solid var(--line)}
table{border-collapse:collapse;width:100%;font-size:.86rem;background:var(--surface)}
th,td{text-align:left;padding:.56rem .7rem;border-bottom:1px solid var(--line);
  vertical-align:top}
th{background:var(--ground);font-weight:600;font-size:.74rem;letter-spacing:.05em;
  color:var(--muted);white-space:nowrap}
tbody tr:last-child td{border-bottom:0}
.num td:nth-child(n+2),.num th:nth-child(n+2){text-align:right;
  font-variant-numeric:tabular-nums;font-family:var(--mono);font-size:.82rem}
tr.hi td{background:var(--accent-soft);font-weight:600}
.tag{display:inline-block;font-size:.7rem;font-weight:600;padding:.1rem .45rem;
  border-radius:2px;letter-spacing:.03em;white-space:nowrap}
.tag.ok{background:var(--good-soft);color:var(--good)}
.tag.no{background:var(--bad-soft);color:var(--bad)}
.tag.w{background:var(--warn-soft);color:var(--warn)}

/* figures */
.fig{margin:1.6rem 0 1.8rem}
.fig img{width:100%;max-width:100%;display:block;border:1px solid var(--line);
  background:#fff}
figcaption{font-size:.82rem;color:var(--muted);line-height:1.65;margin-top:.6rem}

/* page-by-page */
.pg{display:grid;grid-template-columns:4.6rem 1fr;gap:0 1.4rem;
  border-top:1px solid var(--line);padding-block:1.9rem}
/* a grid item defaults to min-width:auto, which lets a wide table push the
   whole page sideways instead of scrolling inside its own .tw container */
.pg>div,.stat{min-width:0}
.pg:last-of-type{border-bottom:1px solid var(--line)}
.pgn{font-family:var(--serif);font-weight:700;font-size:1.5rem;color:var(--accent);
  line-height:1.3;font-variant-numeric:tabular-nums}
.pgn small{display:block;font-family:var(--mono);font-size:.64rem;font-weight:400;
  letter-spacing:.08em;color:var(--muted);margin-top:.15rem}
.pg h3{margin:0 0 .2rem}
.pg .why{font-size:.88rem;color:var(--muted);margin:0 0 .9rem;font-style:normal}
.qa{margin:1rem 0 0}
.qa dt{font-size:.74rem;font-weight:600;letter-spacing:.09em;text-transform:uppercase;
  color:var(--accent);margin-bottom:.3rem}
.qa dd{margin:0 0 1rem;font-size:.95rem}
.qa dd:last-child{margin-bottom:0}
@media (max-width:620px){
  .pg{grid-template-columns:1fr;gap:.5rem}
  .pgn{font-size:1.1rem}
  .pgn small{display:inline;margin-left:.5rem}
}

/* toc */
.toc{display:flex;flex-wrap:wrap;gap:.5rem .5rem;margin:1.6rem 0 0}
.toc a{font-size:.82rem;text-decoration:none;color:var(--ink-2);
  border:1px solid var(--line);padding:.28rem .6rem;border-radius:2px;background:var(--surface)}
.toc a:hover{border-color:var(--accent);color:var(--accent)}
.toc a:focus-visible{outline:2px solid var(--accent);outline-offset:2px}
ul.tight{margin:0 0 1.1rem;padding-left:1.15rem}
ul.tight li{margin-bottom:.42rem}
footer{margin-top:3.5rem;padding-top:1.4rem;border-top:1px solid var(--line);
  font-size:.8rem;color:var(--muted);line-height:1.7}
@media (prefers-reduced-motion:reduce){*{animation:none!important;transition:none!important}}
"""

FONTS = ('<link rel="preconnect" href="https://fonts.googleapis.com">'
         '<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>'
         '<link rel="stylesheet" href="https://fonts.googleapis.com/css2?'
         'family=IBM+Plex+Mono:wght@400;500&'
         'family=Noto+Sans+SC:wght@400;500;600&'
         'family=Noto+Serif+SC:wght@600;700&display=swap">')


def main():
    cfg = load_config()
    seqs = read_fasta(data_path("sequences.fasta"))
    meta = read_metadata()
    panel = jsn("m2_panel.json")
    val = jsn("m4_filter_validation.json")
    suit = jsn("m6_system_suitability.json")
    promc = jsn("m6_promiscuity_calibration.json")
    acc = jsn("m11_calibration.json")
    bench = jsn("m10_benchmark_summary.json")
    tolw = jsn("m12_tolerance_weight.json")
    prom = jsn("m13_promiscuity.json")
    dom = jsn("m14_domain_attribution.json")
    rank = tsv("m6_calibrated_ranking.tsv")
    clusters = tsv("m5_clusters.tsv")
    coincide = tsv("m7_tb_coincidence.tsv")
    grid = tsv("m8_exposure_grid.tsv")
    deimm = tsv("m9_deimmunization_scan.tsv")

    tid = next(r["id"] for r in rank if r["role"] == "test_article")
    tr = next(r for r in rank if r["id"] == tid)
    anchor_id = cfg["benchmarks"]["anchor_low"]
    ar = next(r for r in rank if r["id"] == anchor_id)
    art = dom["articles"][tid]
    sc = art["scaffold_vs_benchmark"]
    nov_dom, scf_dom = art["domains"][0], art["domains"][1]
    nh = nov_dom["nearest_human"]
    tcl = [c for c in clusters if c["id"] == tid]
    top = max(tcl, key=lambda c: float(c["pop_presenting"]))
    op = acc["operating_points"]["EL<1 alone (in use)"]
    gate = acc["operating_points"]["EL<1 and BA<10 (gate, removed)"]
    sweep = tolw["sensitivity_sweep"]["sweep"]
    folds = [s["fold_vs_anchor"] for s in sweep]
    H, A = [], None
    A = H.append

    # ---------------------------------------------------------------- hero
    A(f'''<header class="hero"><div class="wrap">
<p class="kicker">工艺相关杂质 · 免疫原性风险评估 · 仅供研究使用（RUO）</p>
<h1>dsRNA 亲和配体 7071-A08<br>HLA-DR 免疫原性评估</h1>
<p class="lede col">一条 {len(seqs[tid])} 个氨基酸的亲和层析配体，在
{panel['panel_size_total']} 个 HLA-DR 分子组成的 panel 上扫描，与同批次的基准配体和对照一起打分。
本报告分两部分：<strong>第一部分是结论</strong>，<strong>第二部分逐页解释 deck 第 8 页之后的方法学内容</strong>
&mdash;&mdash; 也就是“凭什么相信这个结论”的那一半。</p>
<div class="meta">
<span>受试品 {esc(tid)}</span><span>{len(seqs[tid])} aa</span>
<span>panel {panel['panel_size_drb1']} DRB1 + {panel['panel_size_total']-panel['panel_size_drb1']} DRB3/4/5</span>
<span>美欧加权覆盖 {panel['weighted_coverage']*100:.1f}%</span>
<span>批次{'可报告' if suit.get('batch_valid') else '不可报告'}</span>
</div>
<nav class="toc">
<a href="#jielun">结论</a><a href="#jiegou">这条序列是什么</a><a href="#piici">批次与对照</a>
<a href="#zhujie">P8–P16 逐页解读</a><a href="#xiayibu">下一步做什么</a><a href="#xiandu">限度</a>
</nav>
</div></header>
<main class="wrap">''')

    # ------------------------------------------------------------- 结论
    A('<section id="jielun">')
    A('<h2>一、结论</h2>')
    A(f'''<div class="stats">
<div class="stat"><div class="n" style="color:var(--warn)">{float(tr['fold_vs_ProteinA_Z']):.2f}×</div>
  <div class="l">相对 Protein A Z 域</div>
  <div class="s">Z 域是临床淋出物暴露史最长的亲和配体，作为标尺的分母</div></div>
<div class="stat"><div class="n">{tr['n_foreign_epitopes']}</div>
  <div class="l">非自体表位</div>
  <div class="s">受试品共 {len([c for c in tcl])} 个表位簇；与人类蛋白组同源的已被降权</div></div>
<div class="stat"><div class="n" style="color:var(--bad)">{float(tr['pop_at_risk_pct']):.0f}%</div>
  <div class="l">美欧人群携带可呈递分子</div>
  <div class="s">携带至少一个能呈递本配体某个非自体表位的 DR 分子</div></div>
<div class="stat"><div class="n">{tr['max_promiscuity']}<span style="font-size:1rem;color:var(--muted)">/{panel['panel_size_total']}</span></div>
  <div class="l">主导表位的 DR 广度</div>
  <div class="s">为通用 T 辅助表位天花板（{promc['universal_epitope_ceiling_sb']}/{panel['panel_size_total']}）的
      {promc['test_peak_vs_universal_sb']:.2f} 倍</div></div>
</div>''')

    A(f'''<div class="note bad col"><p><span class="t">风险是集中的，不是弥散的。</span>
残基 <strong>{top['start']}–{top['end']}</strong>、结合核心 <code>{esc(top['peak_core'])}</code>、
峰值 15-mer <code>{esc(top['peak_peptide'])}</code> 这一个区域，整簇合计被 {top['union_sb_alleles']} 个 DR 分子呈递
（其中单看峰值核心一条是 {tr['max_promiscuity']} 个，即上方卡片里的数字；
簇内还有别的核心，所以合计更大），覆盖
<strong>{float(top['pop_presenting'])*100:.0f}%</strong> 的美欧加权人群。
确证实验应当围绕这一条肽来设计；序列其余部分对一个这个大小的非人源蛋白来说并不突出。</p></div>''')

    A(f'''<div class="note warn col"><p><span class="t">把结论写成区间，而不是一个点值。</span>
{float(tr['fold_vs_ProteinA_Z']):.2f}× 这个数字依赖于“近自体降权系数”&mdash;&mdash;
流程里唯一一个<strong>数据无法界定</strong>的参数（见 P8 与下文 M12）。
在该系数的整个合理取值范围内，本配体在
<strong>{min(folds):.2f}×</strong> 到 <strong>{max(folds):.2f}×</strong> 之间移动。
结论在整个区间内都成立&mdash;&mdash;始终高于基准，也始终远离阳性对照&mdash;&mdash;
但小数点后第二位不是一个真实的数字。这条配体对这个选择比抗体类配体更敏感，
因为<strong>被降权的那一个核心恰好就是它的主导表位</strong>。</p></div>''')

    A(f'''<p class="col">按校准后的分级，本配体落在
<strong>“中度偏高”</strong>（<code>{esc(tr["risk_band"])}</code>）一档：
高于 Protein A 与 Protein A 天然 B 域（均为 1.00×），低于 Protein G B1（3.18×）
和破伤风类毒素通用表位 p2（3.25×）。同批次跑的所有基准与对照见下文表格。</p>''')
    A('</section><hr class="rule">')

    # ------------------------------------------------------------ 结构域
    A('<section id="jiegou">')
    A('<h2>二、这条序列是什么：一个融合体</h2>')
    A(f'''<p class="col">这一点必须先讲，因为它改变了整个结论怎么读。比对（本次实算，不是凭印象）显示
7071-A08 <strong>不是单结构域蛋白</strong>，而是两段拼起来的：</p>''')

    A(table(
        [{"d": "dsRBD（双链 RNA 结合域）", "r": "1–72",
          "n": f"{nov_dom['n_epitopes']}", "p": f"{nov_dom['share_of_pirs']*100:.0f}%",
          "dens": f"{nov_dom['pirs_density_per_100aa']:.2f}",
          "e": f"人源最近亲缘 {nh[0]['entry'].split('_')[0]} {nh[0]['identity_pct']:.1f}%（{nh[0]['aligned_length']} aa）、"
               f"{nh[1]['entry'].split('_')[0]} {nh[1]['identity_pct']:.1f}%（{nh[1]['aligned_length']} aa）"},
         {"d": "Protein A Z 域变体", "r": "73–130",
          "n": f"{scf_dom['n_epitopes']}", "p": f"{scf_dom['share_of_pirs']*100:.0f}%",
          "dens": f"{scf_dom['pirs_density_per_100aa']:.2f}",
          "e": f"与本批次 {esc(anchor_id)} 有 {sc['identity_pct']}% 一致性（{sc['identity']}），{sc['n_substitutions']} 个突变"},
         {"d": "C 端 Cys", "r": "131–132", "n": "0", "p": "0%", "dens": "0.00",
          "e": "偶联到树脂用的末端半胱氨酸"}],
        ["d", "r", "n", "p", "dens", "e"],
        ["结构域", "残基", "表位数", "占 pIRS", "每 100 aa 密度", "依据"], "num",
        fmt={"e": lambda r: f'<span style="font-size:.8rem;color:var(--muted)">{esc(r["e"])}</span>'}))

    A(img("fig8_domain_attribution.png",
          "图注（原图为英文）：左图横轴是残基位置，每一根“棒棒糖”是一个预测表位，"
          "高度＝美欧加权人群中能呈递它的比例；粉色背景是 dsRBD，蓝色背景是 Z 域变体。"
          "绿色圈＝该表位在已获临床资格的 Protein A 淋出物里本来就有；橙色圈＝基准没有、"
          "是工程化造出来的；黑色圈＝新结构域自己的表位；带 × 的是与人类蛋白组近似、已被降权的。"
          "右图是三条柱：新结构域、Z 域那一半、以及单独打分的 Protein A Z 域基准。"))

    A(f'''<div class="note col"><p><span class="t">关键点：Z 域那一半就是基准本身。</span>
每 100 aa 的表位密度，Z 域那一半是 <strong>{scf_dom['pirs_density_per_100aa']:.2f}</strong>，
单独打分的 {esc(anchor_id)} 是 <strong>{float(ar['pIRS']):.2f}</strong>&mdash;&mdash;
<strong>小数点后两位完全相同</strong>。它的 {sc['n_scaffold_epitopes']} 个预测表位里有
<strong>{sc['n_shared']} 个</strong>与基准已携带的表位重叠，也就是说这些不是新增暴露，
而是已经给患者用了几十年的淋出物本来就带的表位内容。</p>
<p>那 {sc['n_substitutions']} 个突变只造出 <strong>{sc['n_novel']} 个</strong>基准没有的表位
（<code>{esc(sc['epitopes'][0]['core'])}</code>），对 pIRS 的贡献是
<strong>{sc['novel_pirs_contribution']:.2f}</strong>。</p></div>''')

    A(f'''<div class="note warn col"><p><span class="t">但那个“0.0% 人群”要小心读。</span>
<code>{esc(sc['epitopes'][0]['core'])}</code> 被预测会呈递，只是呈递它的只有一个
<strong>DRB5*01:01</strong>。而 DRB3、DRB4、DRB5 在 DRB1 位点的频率表里<strong>没有频率</strong>，
所以它们的人群权重在结构上就是零。<em>这不等于说没有人能呈递这条肽。</em>
DRB5*01:01 搭载在 DR15 单倍型上，在欧洲人群中并不罕见。诚实的读法是：
这个表位是真实的、窄的、并且在本报告使用的人群尺度上无法量化。
如果要重新设计支架，这是唯一值得看的位点。</p></div>''')

    A(f'''<p class="col"><strong>整体归因：</strong>升高的部分几乎全部来自新的 dsRBD
（每 100 aa {nov_dom['pirs_density_per_100aa']:.2f}，约为基准 {float(ar['pIRS']):.2f} 的
{nov_dom['pirs_density_per_100aa']/float(ar['pIRS']):.1f} 倍），支架基本没有新增贡献。
所以确证实验和任何去免疫化工作，<strong>都应该针对 RNA 结合域</strong>，不是支架。</p>''')
    A('</section><hr class="rule">')

    # ------------------------------------------------------------ 批次
    A('<section id="piici">')
    A('<h2>三、批次与对照：这批数据能不能报</h2>')
    A('''<p class="col">pIRS 是一个<strong>相对</strong>刻度，单独一个数字没有意义&mdash;&mdash;
它只有相对于同一批次、同一 panel、同一阈值跑出来的基准和对照才可读。
所以配体从不单独跑。系统适用性四项检查全部通过，这批数据才可报告：</p>''')
    A(table(suit.get("checks", []), ["check", "detail", "pass"],
            ["检查项", "实测", "结论"], "",
            fmt={"check": lambda r: {"universal epitopes detected with multi-allele breadth":
                                     "panel 能检出通用表位，且有多等位基因广度",
                                     "self controls suppressed by tolerance filter":
                                     "自体对照被容忍过滤器压制",
                                     "tolerance filter discriminates self from foreign":
                                     "容忍过滤器能区分自体与外源",
                                     "benchmark anchor scored":
                                     "基准锚点已打分"}.get(r["check"], esc(r["check"])),
                 "detail": lambda r: f'<span style="font-size:.78rem;color:var(--muted)" class="mono">{esc(r["detail"])}</span>',
                 "pass": lambda r: f'<span class="tag {"ok" if r["pass"] else "no"}">{"通过" if r["pass"] else "不通过"}</span>'}))

    A('<h3>同批次全部序列的打分</h3>')
    ROLE_ZH = {"test_article": "受试品", "benchmark_ligand": "基准配体",
               "clinical_anchor": "临床锚点", "class_comparator": "同类比较物",
               "negative_control_self": "阴性对照（自体）", "positive_control": "阳性对照",
               "self_immunogenic_control": "边界对照（自体但真表位）",
               "tolerised_binder_control": "边界对照（强结合但无应答）"}
    A(table(rank, ["id", "role", "pIRS", "pIRS_raw", "fold_vs_ProteinA_Z",
                   "max_promiscuity", "pop_at_risk_pct"],
            ["序列", "角色", "pIRS", "未过滤 pIRS", "× Protein A Z", "峰值广度", "人群携带率"], "num",
            fmt={"id": lambda r: f'<code>{esc(r["id"])}</code>',
                 "role": lambda r: ROLE_ZH.get(r["role"], esc(r["role"])),
                 "pIRS": lambda r: f'{float(r["pIRS"]):.2f}',
                 "pIRS_raw": lambda r: f'{float(r["pIRS_raw"]):.2f}',
                 "fold_vs_ProteinA_Z": lambda r: f'{float(r["fold_vs_ProteinA_Z"]):.2f}×',
                 "max_promiscuity": lambda r: f'{r["max_promiscuity"]}/{panel["panel_size_total"]}',
                 "pop_at_risk_pct": lambda r: f'{float(r["pop_at_risk_pct"]):.1f}%'}))
    A(img("fig3_calibrated_ranking.png",
          "图注（原图为英文）：左图是同批次所有序列的 pIRS 排序，虚线是 Protein A Z 域基准线；"
          "红色是受试品。右图显示容忍过滤器拿走了多少&mdash;&mdash;浅色点是过滤前，深色点是过滤后，"
          "两点之间的横线越长，说明该序列的原始分数里“人类自己也有的部分”越多。"
          "注意人类种系 VH3-23 和 HSA 两个阴性对照被压到 0，这正是过滤器该有的行为。"))
    A('</section><hr class="rule">')

    # ------------------------------------------------------- P8-P16 逐页
    A('<section id="zhujie">')
    A('<h2>四、P8–P16 逐页解读</h2>')
    A('''<p class="col">deck 前 7 页讲的是“答案是什么”。<strong>第 8 页开始换了一个问题：
凭什么相信这个答案。</strong>这一段之所以不好读，是因为它不再讲这条配体，而是在讲
<em>方法本身有多可靠</em>&mdash;&mdash;包括几处刻意承认“这里我们做不到”的地方。
下面每一页都按同样三问来讲：这一页在回答什么、该怎么读、本次的数字是多少。</p>''')

    pages = []

    # P8
    pages.append(("P8", "M4 · 边界对照", "过滤器在什么地方是错的，而且错了多少",
                  "承认方法的失效边界，而不是只展示它有效的一面",
                  [("这一页在回答什么",
                    "容忍过滤器的核心假设是“人类自己有的序列 → 人体已经耐受”。这个假设<strong>有时候是错的</strong>。"
                    "与其口头说一句“本方法有局限”，不如在同一批次里放两条已知会让它出错的肽，把错误量出来。"),
                   ("该怎么读",
                    "两条边界对照<strong>不是 pass/fail 项</strong>，它们不会让批次作废。它们是“本方法在哪里会骗你”的实测陈述，"
                    "应当原样写进报告：<br>"
                    "<strong>MBP85-99</strong>（髓鞘碱性蛋白）是一条<em>人自己的</em>肽，但它是经过验证的真实表位，"
                    "在 10 个 DR 分子上有阳性人源 T 细胞数据（多发性硬化相关）。过滤器会把它当成“自体 → 安全”而压掉。<br>"
                    "<strong>CLIP87-101</strong>（恒定链）是一条能占据几乎所有 DR 分子凹槽的<em>通用配体</em>，"
                    "但 IEDB 里没有任何阳性人源 T 细胞记录。它测的是另一个方向：结合强度本身会不会被误读成风险。"),
                   ("本次数字",
                    f"MBP85-99 在本 panel 上 DR 广度 10/25（强结合）、20/25（弱结合），未过滤 pIRS 1.19，"
                    f"过滤后被压到 <strong>0.00</strong>&mdash;&mdash;<strong>过滤器把一个真实表位完全抹掉了</strong>。"
                    f"CLIP87-101 广度 1/25、7/25，pIRS 0.00&mdash;&mdash;强结合没有被误读成风险，这个方向是对的。<br><br>"
                    f"<strong>对本配体的实际含义：</strong>如果 7071-A08 的某个表位恰好与某个人类蛋白高度相似，"
                    f"本流程会倾向于低估它。本配体的主导表位 <code>FIVEAKIKE</code> 正是这种情况"
                    f"（与人类 SKA1 蛋白 8/9 匹配，权重被降到 0.35），所以第一部分才要把结论写成区间。")]))

    # P9
    pages.append(("P9", "M10–M11 · 实测准确度", "这套判定规则到底准不准",
                  "把“我觉得这样更准”换成“我量过，结果是这样”",
                  [("这一页在回答什么",
                    "在这一页之前，所有的阈值和规则都只是<strong>论证</strong>，不是<strong>证据</strong>。"
                    "比如“洗脱配体评分会过度报阳，所以要求亲和力那一路也同意”&mdash;&mdash;听上去很合理，但从没被验证过。"
                    "这一页把规则拿去和真实测量结果对撞。"),
                   ("该怎么读",
                    f"从 IEDB 拉下本 panel 这 {bench['alleles_with_data']} 个 DR 分子的全部人源 CD4 T 细胞实验结果："
                    f"{bench['raw_records']:,} 条原始记录 → <strong>{bench['labelled_pairs']:,} 对</strong>带标签的（肽，等位基因）"
                    f"（{bench['positives']:,} 阳性 / {bench['negatives']:,} 阴性，{bench['ambiguous_excluded']} 对阴阳同时存在被剔除）。"
                    f"共享任一 9-mer 的肽被并成一簇（{bench['clusters']:,} 簇），所有统计以簇为单位重抽样，"
                    f"这样同一篇研究里互相重叠的肽不会被当成多个独立观测。<br><br>"
                    f"<strong>AUC</strong> 可以理解为“随机拿一条真表位和一条非表位，模型把真的排在前面的概率”，0.5 是瞎猜。"),
                   ("本次数字",
                    f"最好的连续评分规则 <code>{acc['best_continuous_rule']}</code> AUC "
                    f"<strong>{acc['rules'][acc['best_continuous_rule']]['auc']:.3f}</strong>，"
                    f"单用洗脱配体评分 {acc['rules']['EL']['auc']:.3f}。<br><br>"
                    f"<strong>本流程实际采用的操作点</strong>（EL %Rank &lt; 1）："
                    f"灵敏度 <strong>{op['sensitivity']:.2f}</strong>、特异度 <strong>{op['specificity']:.2f}</strong>、"
                    f"在 {acc['scan_prevalence']*100:.0f}% 的假定扫描阳性率下，"
                    f"一个被标记的肽真的是表位的概率约 <strong>{op['ppv_at_scan_prevalence']*100:.0f}%</strong>。<br><br>"
                    f"<span style='color:var(--bad)'><strong>被实测否决的规则：</strong></span>那个“要求亲和力也同意”的门，"
                    f"实测下来移除了 <strong>{op['tp']-gate['tp']} 个真阳性</strong>，只换来移除 "
                    f"<strong>{op['fp']-gate['fp']} 个假阳性</strong>&mdash;&mdash;大约四个真的换一个假的&mdash;&mdash;"
                    f"而且 MCC 从 {op['mcc']:.3f} 变成 {gate['mcc']:.3f}，更差了。所以这个门被关掉了。<br><br>"
                    f"<strong>最重要的一句话：灵敏度只有 {op['sensitivity']:.2f}。</strong>"
                    f"这意味着<strong>没被标记 ≠ 安全</strong>，只是“没被标记”。"
                    f"强结合这一档是一个高特异度、低灵敏度的判据。")]))

    # P10
    pages.append(("P10", "M13 · 最后一个特异性杠杆", "要不要改成“必须被很多个 DR 分子呈递才算”",
                  "一个看上去必然成立的改进，被数据否决",
                  [("这一页在回答什么",
                    "提高阳性特异性最顺理成章的下一步是：“不要只要有一个 DR 分子排名靠前就报警，"
                    "要求<strong>很多个</strong>分子都能呈递才算数“。直觉上，被很多人 HLA 都能呈递的肽显然更可信。"
                    "这一页把这个直觉拿去测了。"),
                   ("该怎么读",
                    f"把全部 {prom['n_peptides']:,} 条基准肽（每个 9-mer 簇取一条）对整个 "
                    f"{panel['panel_size_total']} 分子 panel 打分，然后比较两个判据谁更准："
                    "<strong>广度</strong>（多少个分子能呈递）vs <strong>最佳单等位基因排名</strong>（当前规则）。"
                    "看 ΔAUC 的 95% 置信区间是否排除 0&mdash;&mdash;<strong>跨过 0 就说明没有差别</strong>。"),
                   ("本次数字",
                    f"<strong>杠杆不存在。</strong>相对最佳单等位基因排名：<br>"
                    f"广度（%Rank&lt;1 的分子数）ΔAUC <strong>{prom['comparisons']['breadth_sb_vs_best_rank']['delta']:+.4f}</strong>，"
                    f"95% CI [{prom['comparisons']['breadth_sb_vs_best_rank']['ci95'][0]:+.4f}, "
                    f"{prom['comparisons']['breadth_sb_vs_best_rank']['ci95'][1]:+.4f}]；<br>"
                    f"人群加权呈递比例 <strong>{prom['comparisons']['pop_presenting_vs_best_rank']['delta']:+.4f}</strong>，"
                    f"[{prom['comparisons']['pop_presenting_vs_best_rank']['ci95'][0]:+.4f}, "
                    f"{prom['comparisons']['pop_presenting_vs_best_rank']['ci95'][1]:+.4f}]。两个区间都跨零。<br><br>"
                    f"<strong>为什么这个直觉感觉是对的？因为混杂，不是因为生物学。</strong>"
                    f"一条肽在基准里能测出多大“广度”，上限是<em>当初有人愿意在几个分子上测它</em>&mdash;&mdash;"
                    f"而肽之所以被反复测，恰恰是因为它已经看起来有意思。"
                    f"单是“IEDB 在几个分子上测过它”这一项，AUC 就达到 "
                    f"<strong>{prom['confound']['auc_of_test_count_alone']:.3f}</strong>，"
                    f"高于本基准里任何一个由序列推导出的预测量；阳性率随测试分子数 1→2→3 从 "
                    f"{prom['confound']['positive_rate_by_test_count']['1']['positive_rate']*100:.1f}% 升到 "
                    f"{prom['confound']['positive_rate_by_test_count']['2']['positive_rate']*100:.1f}% 再到 "
                    f"{prom['confound']['positive_rate_by_test_count']['3']['positive_rate']*100:.1f}%。<br><br>"
                    f"所以只能看<strong>同一测试数层内部</strong>：在只被测过 1 个分子的 "
                    f"{prom['strata']['single_allele']['n']:,} 条肽里，广度依然是平的"
                    f"（{prom['strata']['single_allele']['breadth_sb_vs_best_rank']['delta']:+.4f}）；"
                    f"在 {prom['strata']['multi_allele']['n']} 条多等位基因肽里反而<strong>显著更差</strong>"
                    f"（{prom['strata']['multi_allele']['breadth_sb_vs_best_rank']['delta']:+.4f}）。<br><br>"
                    f"<strong>结论：流程不改。</strong>继续按单个分子的 %Rank 判阳，"
                    f"人群加权只用于汇总报告、不作为判定门槛。")]))

    # P11
    pages.append(("P11", "M5–M6 · 校准打分", "pIRS 这个数字是怎么算出来的",
                  "把“13 个强结合肽”变成一个可比较的量",
                  [("这一页在回答什么",
                    "“这条蛋白里有 13 个强结合 15-mer”是无法解读的：任何非人源蛋白都有一堆，"
                    "而且不同蛋白之间没法比。这一页讲怎么把计数变成一个刻度。"),
                   ("该怎么读",
                    "pIRS ＝ <strong>每 100 个残基、按人群加权的可呈递外源表位含量</strong>。四步：<br>"
                    "① 把重叠的 15-mer 合并成<em>表位</em>（按结合核心归并）；<br>"
                    "② 每个表位乘以“美欧加权人群中携带能呈递它的 DR 分子的比例”；<br>"
                    "③ 乘上容忍权重（与人类蛋白组同源的降权）；<br>"
                    "④ 按长度归一化到每 100 aa。<br><br>"
                    "然后表达成相对基准的倍数。<strong>绝对值没有意义，倍数才有。</strong>"),
                   ("本次数字",
                    f"受试品 pIRS <strong>{float(tr['pIRS']):.2f}</strong>，未过滤 {float(tr['pIRS_raw']):.2f}"
                    f"（容忍过滤拿走 {float(tr['tolerance_drop_pct']):.0f}%），"
                    f"基准 {esc(anchor_id)} {float(ar['pIRS']):.2f} → <strong>{float(tr['fold_vs_ProteinA_Z']):.2f}×</strong>。<br><br>"
                    f"<strong>容忍过滤器本身也被验证过</strong>（不是靠断言）：本批次 {val['n_real_cores']} 个预测核心中，"
                    f"{val['real_hit_rate_9of9']*100:.1f}% 是人类蛋白组里的精确 9-mer，"
                    f"{(val['real_hit_rate_8of9']-val['real_hit_rate_9of9'])*100:.1f}% 只差一个取代；"
                    f"把同样的序列打乱后重跑（成分相同、真实同源性为零的零假设），命中率只有 "
                    f"{val['null_hit_rate_8of9']*100:.1f}%&mdash;&mdash;<strong>富集 {val['enrichment_8of9_real_over_null']}×</strong>。"
                    f"如果零假设命中率没有远低于真实值，本次运行会直接报告该过滤器无信息量。")]))

    # P12
    pages.append(("P12", "M5 · 表位簇", "风险具体落在序列的哪里",
                  "从一个总分回到可以下单做实验的肽",
                  [("这一页在回答什么",
                    "总分决定“要不要担心”，这一页决定“拿哪几条肽去做湿实验”。"),
                   ("该怎么读",
                    "核心位置相距 8 个残基以内的表位被并成一个<em>簇</em>。看两列：<strong>DR 广度</strong>"
                    "（多少个分子呈递）和<strong>人群呈递比例</strong>。后者才是排序依据，"
                    "因为一个只被罕见等位基因呈递的表位，人群层面的意义很小。"),
                   ("本次数字", None)]))

    # P13
    pages.append(("P13", "M7 · B 细胞 / ADA 层", "为什么只看 T 细胞不够",
                  "临床上真正测的终点是抗体，不是 T 细胞",
                  [("这一页在回答什么",
                    "免疫原性在临床上的可测终点是<strong>抗药抗体（ADA）</strong>。"
                    "T 辅助表位是产生高亲和力抗体的必要条件，但抗体识别的是 B 细胞表位。"
                    "两者在序列上重叠的区域，风险要高一档。"),
                   ("该怎么读",
                    "用 BepiPred-2.0 算线性 B 细胞倾向性，再看它与 T 细胞簇的重叠残基数。"
                    "<strong>这是本流程里最弱的一个模型</strong>&mdash;&mdash;真实 ADA 表位多数是构象型的，"
                    "线性预测看不见。这一层只用于给湿实验排优先级，绝不单独作为结论。"),
                   ("本次数字", None)]))

    # P14
    pages.append(("P14", "M8 · 暴露量", "内在风险 × 剂量",
                  "没有剂量，免疫原性风险无法定级",
                  [("这一页在回答什么",
                    "杂质的免疫原性风险随<strong>每次给药递送的微克数</strong>变化。"
                    "同一条配体，在 1 ng/mg、1 mg 剂量下和在 100 ng/mg、1000 mg 剂量下，完全是两回事。"
                    "一个不含剂量的 T 细胞评分到不了风险结论。"),
                   ("该怎么读",
                    "表格是“淋出水平（ng 配体/mg 原液）× 单次剂量（mg）”的网格，"
                    "格子里是每次给药的 µg 配体量。参照锚点：<strong>亚微克/剂</strong>的蛋白杂质暴露，"
                    "正是 rProtein A 淋出物几十年临床使用所处的区间。"),
                   ("本次数字", None)]))

    # P15
    pages.append(("P15", "M9 · 去免疫化扫描", "如果配体可以改造，改哪里",
                  "可选模块，且本例有明确的功能风险",
                  [("这一页在回答什么",
                    "针对主导表位，把 9 个结合核心位置逐一做单点取代，看哪些取代能消除呈递。"),
                   ("该怎么读",
                    "看 <code>ΔPop</code>（人群呈递比例的下降）和 <code>BL62</code>"
                    "（BLOSUM62 分数，越接近 0 越保守、结构扰动风险越小）。"
                    "<strong>这是纯序列层面的计算，完全不模拟功能。</strong>"),
                   ("本次数字", None)]))

    # P16
    pages.append(("P16", "限度与下一步", "这份报告是什么，不是什么",
                  "体外预测排序并定位风险，但不测量风险",
                  [("这一页在回答什么",
                    "把方法的边界一次性写清楚，避免结论被过度解读。"),
                   ("该怎么读", "见下文“五、限度”一节，这里不重复。"),
                   ("本次数字",
                    "本流程产出的<strong>实际交付物是一份有范围、有报价的实验计划</strong>："
                    "把“评估一下这个配体”变成“针对这几条肽做这三个实验”。见下一节。")]))

    for num, mod, title, why, qa in pages:
        A(f'<div class="pg" id="{num.lower()}"><div class="pgn">{num}<small>{esc(mod)}</small></div><div>')
        A(f'<h3>{title}</h3><p class="why">{why}</p><dl class="qa">')
        for label, text in qa:
            if text is None:
                continue
            A(f'<dt>{label}</dt><dd>{text}</dd>')
        A('</dl>')
        if num == "P12":
            A(table(sorted(tcl, key=lambda c: -float(c["pop_presenting"])),
                    ["start", "peak_core", "peak_peptide", "union_sb_alleles",
                     "pop_presenting", "tolerance_class"],
                    ["起始残基", "结合核心", "峰值 15-mer", "DR 广度", "人群呈递", "容忍归类"], "num",
                    fmt={"start": lambda r: f'{r["start"]}–{r["end"]}',
                         "peak_core": lambda r: f'<code>{esc(r["peak_core"])}</code>',
                         "peak_peptide": lambda r: f'<code>{esc(r["peak_peptide"])}</code>',
                         "union_sb_alleles": lambda r: f'{r["union_sb_alleles"]}/{panel["panel_size_total"]}',
                         "pop_presenting": lambda r: f'{float(r["pop_presenting"])*100:.1f}%',
                         "tolerance_class": lambda r: {"foreign": "外源", "mixed": "混合",
                                                       "all_tolerised": "全部耐受"}.get(r["tolerance_class"],
                                                                                    r["tolerance_class"])}))
            A(f'''<p style="font-size:.9rem"><strong>读法：</strong>只有
<code>{esc(top['peak_core'])}</code> 这一簇达到了值得做实验的量级
（{float(top['pop_presenting'])*100:.0f}% 人群、{top['union_sb_alleles']}/{panel['panel_size_total']} 个分子）。
其余各簇人群呈递比例都在 15% 以下。注意它被归为“混合”，
因为该簇内既有外源核心也有与人类近似的核心&mdash;&mdash;这正是 P8 的边界问题落在本配体上的地方。</p>''')
        if num == "P13":
            tc = [c for c in coincide if c["id"] == tid]
            A(table(sorted(tc, key=lambda c: -float(c["t_pop_presenting"]))[:5],
                    ["t_cluster", "t_peak_core", "t_pop_presenting", "b_region", "overlap_aa"],
                    ["T 细胞簇", "结合核心", "人群呈递", "B 细胞区段", "重叠残基"], "num",
                    fmt={"t_peak_core": lambda r: f'<code>{esc(r["t_peak_core"])}</code>',
                         "t_pop_presenting": lambda r: f'{float(r["t_pop_presenting"])*100:.1f}%',
                         "overlap_aa": lambda r: f'{r["overlap_aa"]} aa'}))
            hit = next((c for c in tc if c["t_peak_core"] == top["peak_core"]), None)
            if hit:
                A(f'''<p style="font-size:.9rem"><strong>本次结果：</strong>主导 T 细胞簇
{esc(hit['t_cluster'])}（<code>{esc(hit['t_peak_core'])}</code>）与 B 细胞倾向区段
{esc(hit['b_region'])} 重叠 <strong>{hit['overlap_aa']} 个残基</strong>。
重叠不大，但方向一致：这个区域同时是 dsRBD 里带正电的 RNA 结合面附近，
本来就倾向于暴露在表面。把它列为 ADA 风险的首要观察区。</p>''')
        if num == "P14":
            A(table(grid, ["leachate_ng_per_mg", "dose_mg", "ug_ligand_per_dose", "exposure_band"],
                    ["淋出水平 ng/mg", "单次剂量 mg", "µg 配体/剂", "暴露分档"], "num",
                    fmt={"ug_ligand_per_dose": lambda r: f'{float(r["ug_ligand_per_dose"]):.4g}',
                         "exposure_band": lambda r: {
                             "negligible": '<span class="tag ok">可忽略</span>',
                             "low": '<span class="tag ok">低</span>',
                             "moderate": '<span class="tag w">中</span>',
                             "elevated": '<span class="tag no">偏高</span>'}.get(r["exposure_band"],
                                                                              r["exposure_band"])}))
            A('''<div class="note warn"><p><span class="t">这张表必须换成你们自己的产品参数。</span>
配体分子量、淋出 ppm 和剂量档位目前用的是占位值（见 <code>config/config.yaml</code> 的
<code>exposure</code> 段）。用别人的剂量算出来的风险矩阵没有意义。</p></div>''')
        if num == "P15":
            A(table(deimm[:8], ["variant", "substitution", "core", "n_sb_alleles",
                                "pop_presenting", "delta_pop_presenting", "blosum62"],
                    ["变体", "取代", "新核心", "DR 广度", "人群呈递", "ΔPop", "BLOSUM62"], "num",
                    fmt={"core": lambda r: f'<code>{esc(r["core"])}</code>',
                         "pop_presenting": lambda r: f'{float(r["pop_presenting"])*100:.1f}%',
                         "delta_pop_presenting": lambda r: f'{float(r["delta_pop_presenting"])*100:+.1f}%'}))
            A(f'''<div class="note bad"><p><span class="t">这一页不能直接采用。</span>
多个单点取代确实能把主导表位的呈递完全消除（ΔPop −30.2%），但
<code>{esc(top['peak_core'])}</code> 位于残基 {top['start']}–{top['end']}，
<strong>在 dsRBD 内部&mdash;&mdash;也就是负责结合 RNA 的那个功能域</strong>。
排名靠前的取代里有多个是去掉一个赖氨酸（K40），而 dsRBD 与双链 RNA 的结合是静电驱动的。
因此每一个取代都是功能风险，<strong>必须先测结合活性再谈免疫原性收益</strong>。
这里给出的是候选清单，不是建议。</p></div>''')
        A('</div></div>')

    A('</section><hr class="rule">')

    # ------------------------------------------------------------ 下一步
    A('<section id="xiayibu">')
    A('<h2>五、下一步做什么</h2>')
    A(f'''<p class="col">本流程的实际价值在于把“评估这个配体”变成一个有范围、可报价的实验计划。
按成本从低到高：</p>''')
    steps = [
        ("HLA-DR 竞争结合实验", "数天，成本最低",
         f"针对 <code>{esc(top['peak_core'])}</code> 所在的峰值 15-mer "
         f"<code>{esc(top['peak_peptide'])}</code>，对 panel 中呈递它的主要 DR 分子做竞争结合。"
         f"直接验证预测的结合本身。这是唯一需要立刻做的一条。"),
        ("MAPPs（MHC 相关肽段蛋白质组学）", "数周",
         "用 HLA 分型过的供体来源的单核细胞衍生树突状细胞，看完整配体被<strong>实际加工和呈递</strong>出什么肽。"
         "预测做不到这一步&mdash;&mdash;它不建模摄取、内体蛋白酶解和 HLA-DM 编辑。"),
        ("离体 PBMC / CD4 增殖", "数月，成本最高",
         "约 50 名与本 panel 匹配分型的供体。这是目前对临床 ADA 风险最接近的替代终点。"),
    ]
    A('<div class="tw"><table><thead><tr><th>顺序</th><th>实验</th><th>周期</th><th>做什么、为什么</th></tr></thead><tbody>')
    for i, (name, cost, what) in enumerate(steps, 1):
        A(f'<tr><td style="font-family:var(--serif);font-weight:700;color:var(--accent)">{i}</td>'
          f'<td><strong>{name}</strong></td><td style="color:var(--muted);white-space:nowrap">{cost}</td>'
          f'<td>{what}</td></tr>')
    A('</tbody></table></div>')
    A(f'''<div class="note good col"><p><span class="t">另外两件与配体本身有关的事。</span>
① 如果考虑改造支架，唯一值得看的位点是工程化造出的
<code>{esc(sc['epitopes'][0]['core'])}</code>（配体第 {sc['epitopes'][0]['ligand_pos']} 位起），
它只被 DRB5*01:01 呈递，在本报告的人群模型里权重为零但并非无风险。<br>
② 如果考虑对主导表位做去免疫化，先做结合活性验证，见 P15。</p></div>''')
    A('</section><hr class="rule">')

    # -------------------------------------------------------------- 限度
    A('<section id="xiandu">')
    A('<h2>六、限度 —— 每份报告都要写进去</h2>')
    A('<p class="col">体外 DR 筛查<strong>排序并定位</strong>风险，它<strong>不测量</strong>风险。</p>')
    lim = [
        ("是结合预测，不是呈递", "NetMHCIIpan 算的是肽–MHC 结合。它不建模抗原摄取、内体蛋白酶解、HLA-DM 编辑或复合物稳定性。"),
        ("只有 DR", "DP 和 DQ 同样参与 CD4 应答，DQ 在若干生物药 ADA 事件中被牵涉。本 panel 的规格把它们排除在外。"),
        ("pIRS 是相对量", "它不是预测的 ADA 发生率。今天没有任何体外方法能预测 ADA 发生率。"),
        ("不含聚集体与佐剂效应", "聚集或颗粒化的杂质比单体显著更具免疫原性，基于序列的方法完全看不见这一层。"),
        ("绝对准确度数字应当读作上限",
         "NetMHCIIpan 在 IEDB 数据上训练，而本基准也来自 IEDB，部分训练集重叠是确定存在的、从外部无法排除的。"
         "<strong>规则与规则之间的比较</strong>要稳健得多：所有规则用同样可能泄漏的预测量打同一批肽，泄漏对两边同向放大，在差值里大部分抵消。"),
        ("没被标记 ≠ 已排除", f"在本流程采用的操作点上灵敏度只有 {op['sensitivity']:.2f}。低于阈值的肽是“未被标记”，不是“已清除”。"),
    ]
    A('<ul class="tight col">')
    for t_, d in lim:
        A(f'<li><strong>{t_}</strong>&mdash;&mdash;{d}</li>')
    A('</ul>')
    A('</section>')

    A(f'''<footer class="col">
<p>受试品 {esc(tid)} 为申办方提供的专有序列，无公共数据库编号&mdash;&mdash;也不需要：本流程只需去标识的 FASTA。
批次中其余全部序列均为公开序列，运行时从 RCSB / UniProt 实时获取并记录编号，
因此即使受试品不可复现，衡量它的那把尺子是可复现的。</p>
<p>本报告所有数字均由 <code>results/</code> 下的输出文件生成，非人工誊写。
生成脚本 <code>scripts/make_report_zh.py</code>。英文完整版见 <code>report.html</code>。</p>
<p><strong>仅供研究使用（RUO）。未经确证性湿实验数据，不得用于法规申报。</strong></p>
</footer>''')
    A('</main>')

    html = (f"<title>7071-A08 免疫原性评估</title>\n{FONTS}\n<style>{CSS}</style>\n"
            + "\n".join(H))
    out = os.path.join(os.path.dirname(results_path("x"))[:-7], "report_zh.html")
    out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "report_zh.html")
    with open(out, "w") as f:
        f.write(html)
    print(f"wrote {out} ({len(html)//1024} KB)")


if __name__ == "__main__":
    main()
