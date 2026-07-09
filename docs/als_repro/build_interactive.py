#!/usr/bin/env python3
"""Build the formal INTERACTIVE ALS campaign report (single self-contained HTML):
- all campaign data/tables, 9 figures (base64)
- 2 live 3Dmol.js viewers (nanobody-SOD1 complex ipTM 0.886; small-molecule-SOD1)
- top sequences (VHH + ASO), SMILES, score tables
3Dmol.js + CIF structures inlined -> works fully offline.
"""
import base64, json, pathlib, html

D = pathlib.Path("docs/als_repro")
FIG = D/"figures"; STR = D/"structures"; DATA = D/"data"

def b64img(p):
    return "data:image/png;base64,"+base64.b64encode(pathlib.Path(p).read_bytes()).decode()
def cif(p): return pathlib.Path(p).read_text()
def tsv_rows(p, n=None):
    ls = pathlib.Path(p).read_text().splitlines()
    hdr = ls[0].split("\t"); rows=[l.split("\t") for l in ls[1:]]
    return hdr, (rows[:n] if n else rows)

threed = (STR/"3Dmol-min.js").read_text()
nb_cif = cif(STR/"nanobody_top.cif")
sm_cif = cif(STR/"smallmol_top.cif")
seqs = json.loads((pathlib.Path("scratchpad/als_campaign/design_sequences.json")).read_text())
nb_seq = seqs["nanobody_top"]["NANO1"]
state = json.loads((DATA/"state.json").read_text())

# score tables
nb_hdr, nb_rows = tsv_rows(DATA/"sod1_nanobody_design_top.tsv", 15)   # ipTM,struct,plddt
sm_hdr, sm_rows = tsv_rows(DATA/"sod1_smallmol_design_top.tsv", 10)   # bc,iptm,opt,solub,perm,logp,smiles
aso_hdr, aso_rows = tsv_rows(DATA/"aso_SOD1.tsv", 5)

FIGS = {n: b64img(FIG/f"fig_als_{n}.png") for n in
        ["market","axis","biomarker","synbiomarker","network","matrix","designs","ranking","delivery"]}
RENDER_NB = b64img(STR/"render_nanobody.png")
RENDER_SM = b64img(STR/"render_smallmol.png")

def tbl(hdr, rows, cls=""):
    h="".join(f"<th>{html.escape(c)}</th>" for c in hdr)
    b="".join("<tr>"+"".join(f"<td>{html.escape(str(c))}</td>" for c in r)+"</tr>" for r in rows)
    return f'<table class="{cls}"><thead><tr>{h}</tr></thead><tbody>{b}</tbody></table>'

# nanobody CDR3 highlight in sequence
cdr3="AASPSGGTKFTS"
nb_disp = nb_seq.replace(cdr3, f'<span class="cdr">{cdr3}</span>')

CSS = """
*{box-sizing:border-box} html{-webkit-print-color-adjust:exact}
body{font-family:"Noto Sans","PingFang SC","Microsoft YaHei","WenQuanYi Zen Hei",system-ui,sans-serif;
 color:#1a2230;line-height:1.62;margin:0;background:#f4f6f9}
.wrap{max-width:960px;margin:0 auto;padding:40px 44px 90px;background:#fff}
h1{font-size:25px;border-bottom:3px solid #0b6b5e;padding-bottom:12px;margin:0 0 4px}
h2{font-size:19px;color:#0e2a45;border-left:4px solid #0b6b5e;padding-left:11px;margin:34px 0 10px}
h3{font-size:15.5px;margin:20px 0 6px;color:#22303f}
p{margin:8px 0} a{color:#0d47a1;text-decoration:none}
.lead{background:#f5f9f8;border-left:4px solid #0b6b5e;padding:12px 16px;border-radius:0 6px 6px 0;font-size:14px;color:#33414f}
table{border-collapse:collapse;width:100%;margin:12px 0;font-size:12.5px}
th,td{border:1px solid #e3e8ef;padding:6px 9px;text-align:left;vertical-align:top}
thead th{background:#0e2a45;color:#fff;font-weight:600}
tbody tr:nth-child(even){background:#f7fafc}
img{display:block;max-width:100%;height:auto;margin:12px auto 4px;border:1px solid #e3e8ef;border-radius:8px}
.cap{font-size:12px;color:#667;text-align:center;font-style:italic;margin:0 0 14px}
code,.seq{font-family:ui-monospace,Menlo,Consolas,monospace}
.seqbox{background:#0e2233;color:#d5e6ff;padding:12px 14px;border-radius:8px;font-size:12px;
 word-break:break-all;line-height:1.9;margin:8px 0}
.cdr{background:#e8b800;color:#000;padding:1px 3px;border-radius:3px;font-weight:700}
.viewer{width:100%;height:440px;position:relative;border:1px solid #cfd8e3;border-radius:10px;
 margin:10px 0;background:#0b0f14}
.vwrap{margin:14px 0}
.badge{display:inline-block;padding:2px 9px;border-radius:11px;font-size:11.5px;font-weight:700;color:#fff}
.b-real{background:#0b6b5e}.b-hyp{background:#b8860b}
.grid2{display:grid;grid-template-columns:1fr 1fr;gap:18px}
@media(max-width:760px){.grid2{grid-template-columns:1fr}}
.kpi{display:flex;gap:14px;flex-wrap:wrap;margin:12px 0}
.kpi div{flex:1;min-width:120px;background:#eef3f9;border-radius:10px;padding:12px 14px;text-align:center}
.kpi b{display:block;font-size:22px;color:#0b6b5e}
.controls{font-size:12px;margin:6px 0}
.controls button{font:inherit;border:1px solid #0b6b5e;background:#fff;color:#0b6b5e;border-radius:6px;
 padding:3px 9px;margin-right:5px;cursor:pointer}
.controls button:hover{background:#0b6b5e;color:#fff}
.toc{columns:2;font-size:13px;background:#f7fafc;border:1px solid #e3e8ef;border-radius:8px;padding:12px 18px}
.note{font-size:11px;color:#889;margin-top:22px}
@media print{ .viewer,.controls{display:none!important} .wrap{padding:0} body{background:#fff}
 h2,h3{page-break-after:avoid} table,img{page-break-inside:avoid} }
@page{size:A4;margin:14mm 12mm}
"""

def viewer(vid, cif_data, kind, config_js):
    esc = cif_data.replace("\\","\\\\").replace("`","\\`").replace("$","\\$")
    return f'''
<div class="vwrap">
 <div class="controls">
   <button onclick="vresetsurf_{vid}()">卡通</button>
   <button onclick="vsurf_{vid}()">表面</button>
   <button onclick="v_{vid}.spin(!spin_{vid});spin_{vid}=!spin_{vid}">旋转开关</button>
   <span style="color:#667">拖动旋转 · 滚轮缩放</span>
 </div>
 <div id="{vid}" class="viewer"></div>
</div>
<script>
var v_{vid}, spin_{vid}=false;
(function(){{
  var cifdata=`{esc}`;
  v_{vid}=$3Dmol.createViewer(document.getElementById("{vid}"),{{backgroundColor:"0x0b0f14"}});
  v_{vid}.addModel(cifdata,"cif");
  {config_js}
  v_{vid}.zoomTo(); v_{vid}.render();
}})();
function vresetsurf_{vid}(){{ v_{vid}.removeAllSurfaces(); v_{vid}.render(); }}
function vsurf_{vid}(){{ v_{vid}.addSurface($3Dmol.SurfaceType.VDW,{{opacity:0.65,color:"white"}},{{chain:"A"}}); v_{vid}.render(); }}
</script>'''

# nanobody viewer: SOD1 chain A grey cartoon, NANO1 teal cartoon, interface sticks
nb_cfg = '''
  v_VID.setStyle({chain:"A"},{cartoon:{color:"0x9aa6b2"}});
  v_VID.setStyle({chain:"NANO1"},{cartoon:{color:"0x0b6b5e"}});
  v_VID.addStyle({chain:"NANO1",resi:"98-110"},{stick:{radius:0.18,colorscheme:"yellowCarbon"}});
'''.replace("VID","nbview")
# smallmol viewer: SOD1 cartoon + ligand sticks
sm_cfg = '''
  v_VID.setStyle({chain:"A"},{cartoon:{color:"0x9aa6b2"}});
  v_VID.setStyle({resn:"LIG"},{stick:{radius:0.25,colorscheme:"cyanCarbon"}});
  v_VID.addStyle({resn:"LIG"},{sphere:{scale:0.28,colorscheme:"cyanCarbon"}});
'''.replace("VID","smview")

BODY = f"""
<h1>ALS · TDP-43 蛋白稳态轴 —— 从头多模态药物发现<br>正式全流程报告(含交互式 3D 结构)</h1>
<p class="cap">denovo-drug-campaign skill · 11 阶段 · 全流程真实数据 + GPU 设计 · 预算 ≤100 候选/模态</p>

<div class="lead">
<b>一句话结论。</b>从 8 个候选适应症里数据驱动选定 <b>ALS</b>(最高未满足 + 最强遗传学 + 真多模态 + FDA 认可的 NfL 生物标志物),锁定跨层级 <b>TDP-43 蛋白稳态轴</b>,三方法交叉验证 biomarker(Jaccard 0.74),真跑 5 模态设计:<b>de-novo 纳米抗体 top ipTM 0.886</b>、ebselen 重利用 ipTM 0.824、3×ASO×100,而 de-novo 小分子仅 0.31 —— <b>数据自证 SOD1 上 ASO/抗体 优于小分子</b>。全程 <span class="badge b-real">真数据</span> vs <span class="badge b-hyp">假设/示意</span> 分开标;ipTM 只验 pose,不等于亲和力。
</div>

<div class="kpi">
 <div><b>6,292</b>OT 关联靶点(ALS)</div>
 <div><b>0.886</b>纳米抗体 top ipTM</div>
 <div><b>0.74</b>biomarker Jaccard</div>
 <div><b>δ=1</b>网络双稳(NfL 预警)</div>
 <div><b>$7.6</b>GPU 花费</div>
</div>

<h2>目录</h2>
<div class="toc">
1. 市场分析与适应症选择<br>2. 跨层级靶点轴<br>3. Biomarker 三方法交叉验证<br>4. 合成 biomarker<br>
5. 网络不可约 + 临界转变<br>6. 靶点×模态可行性矩阵<br>7. 多模态 GPU 设计(真实分数)<br>
8. <b>交互式 3D 结构 + 序列 + 对接</b><br>9. 安全/可开发性排序<br>10. CNS 递送 formulation<br>11. 诚实边界
</div>

<h2>1 · 市场分析与适应症选择 <span class="badge b-real">真数据</span></h2>
<p>8 个候选按 Open Targets 关联靶点数 + ClinicalTrials 在招募数(均 live)作硬锚,叠加未满足/多模态/biomarker/遗传/拥挤(专家维度)综合打分。AD 名义第一为<b>研究体量伪影</b>(常见病被研究多→关联/试验数虚高,且拥挤+单模态);<b>ALS 的分数来自对的维度</b>。</p>
<img src="{FIGS['market']}"><p class="cap">图1 适应症综合打分:ALS 在决策相关维度胜出。</p>

<h2>2 · 跨层级靶点轴 <span class="badge b-real">真数据</span></h2>
<p>Open Targets 真实分值:SOD1 0.887 / TARDBP 0.844 / FUS 0.840 / TBK1 0.828 / C9orf72 0.658。轴:<code>SOD1/C9orf72/FUS(基因)→ TDP-43(蛋白枢纽,&gt;95% ALS)→ STMN2/UNC13A(隐蔽外显子)→ NfL</code>。</p>
<img src="{FIGS['axis']}"><p class="cap">图2 ALS TDP-43 蛋白稳态跨层级轴。</p>

<h2>3 · Biomarker 三方法交叉验证 <span class="badge b-real">真数据</span></h2>
<p>pathway-biomarker-triangulation skill 真跑:STRING 38 节点,最小核心 20/40,RWR 重叠 17,<b>Jaccard 0.74</b>;共享核心 = TARDBP/FUS/SOD1/TBK1/OPTN/VCP/C9orf72/ATXN2/<b>NEFH</b>(呼应 NfL)。STMN2/UNC13A 孤儿救回(TDP-43 隐蔽外显子生物学)。</p>
<img src="{FIGS['biomarker']}"><p class="cap">图3 biomarker 三方法(方法一归因 + 方法二 RWR)。</p>
<p><b>临床落地</b>:入组=基因型分层+基线 NfL;PD=NfL 动态(tofersen 已验证)+CSF pTDP-43;效果=ALSFRS-R+生存。</p>

<h2>4 · 合成 biomarker <span class="badge b-hyp">机制示意</span></h2>
<p>Syn-TDP43-STMN2sensor:STMN2 隐蔽外显子剪接传感器,把"TDP-43 核功能丢失"转成可测剪接开关。⚠️ 概念设计(Evo2/Proto GPU 部署 WIP,未做序列自然性打分)。</p>
<img src="{FIGS['synbiomarker']}"><p class="cap">图4 合成 biomarker 盒架构。</p>

<h2>5 · 网络不可约 + 临界转变 <span class="badge b-real">真数据(kernel)</span></h2>
<p>kernel 真跑:3 纤维({{SOD1,C9orf72,FUS}}/{{TDP-43}}/{{STMN2,UNC13A}})、<b>CRNT δ=1</b>(可溶↔聚集双稳)、方差 5.15×、AR1 0.989 → <b>NfL 方差↑/AR1↑ = 早期预警 biomarker</b>。⚠️ 临床态映射与 Schlögl 常数为假设。</p>
<img src="{FIGS['network']}"><p class="cap">图5 TDP-43 蛋白稳态双稳 + NfL 临界慢化早期预警。</p>

<h2>6 · 靶点 × 模态可行性矩阵 <span class="badge b-real">真数据</span></h2>
<p>核心洞见:靶点生物学决定最优模态 —— RNA 层(SOD1/C9orf72/STMN2)→ ASO/siRNA 最优。</p>
<img src="{FIGS['matrix']}"><p class="cap">图6 ALS 靶点×模态可行性(★=选中先导)。</p>

<h2>7 · 多模态 GPU 设计(真实分数)<span class="badge b-real">真数据</span></h2>
<p>全部 Boltz-2(Modal GPU,CD Bioparticles)真跑,每模态 ≤100 候选。</p>
<img src="{FIGS['designs']}"><p class="cap">图7 五模态真实产出:抗体(0.89)&gt;ebselen(0.82)&gt;&gt;小分子(0.31)。</p>

<h2>8 · 交互式 3D 结构 · 序列 · 对接效果 <span class="badge b-real">真数据(可交互)</span></h2>
<p class="lead">以下为 Boltz-2 真实预测结构(可用鼠标旋转/缩放/切换表面)。ipTM/pLDDT 反映界面几何置信度,<b>不等于</b>结合亲和力。</p>

<h3>8.1 冠军先导 L4:SOD1 × de-novo 纳米抗体复合物 —— ipTM 0.886(interaction PAE 2.05 Å)</h3>
<div class="grid2">
<div>
<img src="{RENDER_NB}"><p class="cap">静态渲染(真实坐标,任何环境可见)。灰=SOD1 β-桶 · 青绿=de-novo VHH 结合于一面。</p>
{viewer("nbview", nb_cif, "protein", nb_cfg)}
<p class="cap">↑ 可交互查看器(真实浏览器中拖动/缩放)。灰=SOD1 · 青绿=VHH · 黄=CDR3。</p>
</div>
<div>
<h3 style="margin-top:0">纳米抗体序列(VHH,119 aa)</h3>
<div class="seqbox">{nb_disp}</div>
<p style="font-size:12px">CDR3(黄标)= <code>{cdr3}</code>。框架为 VHH(camelid)骨架。<b>ipTM 0.886 · struct_conf 0.703 · PAE 2.05 Å</b> = 高置信界面几何。</p>
<h3>Top-15 设计评分(100 条中检出)</h3>
{tbl(["ipTM","结构置信","界面 PAE (Å)"], [[f"{float(r[0]):.3f}",f"{float(r[1]):.3f}",f"{float(r[2]):.2f}"] for r in nb_rows])}
</div>
</div>

<h3>8.2 对照先导 L6:SOD1 × de-novo 小分子复合物 —— binding_confidence 0.31(ipTM 0.674)</h3>
<div class="grid2">
<div>
<img src="{RENDER_SM}"><p class="cap">静态渲染(真实坐标)。灰=SOD1 · 红=de-novo 小分子(26 原子)结合于浅表面。</p>
{viewer("smview", sm_cif, "ligand", sm_cfg)}
<p class="cap">↑ 可交互查看器。灰=SOD1 · 青=小分子(球+棒)。</p>
</div>
<div>
<h3 style="margin-top:0">Top 小分子 SMILES + ADME</h3>
<div class="seqbox" style="color:#a8e6cf">{html.escape(sm_rows[0][6])}</div>
<p style="font-size:12px"><b>binding_confidence 0.31 · ipTM 0.67</b>。整体 100 条 0.15–0.31 偏低,坐实 <b>SOD1 是小分子难靶</b>;2 条 ADME 良好(高溶解、cLogP≈2.5)。</p>
<h3>Top-10 设计评分 + ADME(100 条中检出)</h3>
{tbl(["bindConf","ipTM","optScore","溶解","cLogP"], [[f"{float(r[0]):.3f}",f"{float(r[1]):.3f}",f"{float(r[2]):.3f}",r[3],f"{float(r[5]):.2f}"] for r in sm_rows])}
</div>
</div>

<h3>8.3 基因模态 L1:SOD1 gapmer ASO(top 序列,机制=切断/敲低,tofersen 类)</h3>
{tbl(["排名","位点","反义 gapmer 5'→3'","GC%","Tm","CpG"], aso_rows)}
<p style="font-size:12px">另设计 STMN2(剪接开关,恢复全长)+ C9orf72(敲低重复)各 100 条 —— 见附带 TSV。</p>

<h3>8.4 重利用 L5:ebselen × SOD1 co-fold</h3>
<p>Ebselen(SPI-1005,共价 Cys111,Phase 3)co-fold SOD1 → <b>ipTM 0.824 · 结合置信度 0.660</b>(co-fold 不建模共价键,故非上限)。SMILES:<code>O=c1c2ccccc2[se]n1-c1ccccc1</code>。</p>

<h2>9 · 安全 / 可开发性排序 <span class="badge b-real">真数据(kernel)</span></h2>
<p>kernel rank_leads 加权。<b>基因/ASO 模态包揽前三</b>,与靶点 RNA-层生物学一致。</p>
<img src="{FIGS['ranking']}"><p class="cap">图8 6 条先导综合排序:L1 SOD1-ASO 4.54 居首。</p>

<h2>10 · CNS 靶向递送 formulation <span class="badge b-hyp">机制推荐</span></h2>
<p>ASO→鞘内注射(tofersen/nusinersen 类,绕 BBB,最成熟);纳米抗体→AAV 脑内表达 / TfR 穿梭;小分子→脑穿透优化+CNS-LNP。⚠️ 基于已批准先例的机制推荐,非湿实验。</p>
<img src="{FIGS['delivery']}"><p class="cap">图9 CNS 递送:据模态定制。</p>

<h2>11 · 诚实边界</h2>
<ul style="font-size:13px">
<li><span class="badge b-real">真数据/严谨</span> OT/CT/STRING/Ensembl 数据、Boltz co-fold/设计分数(纳米抗体 0.886、小分子 0.31、ebselen 0.824)、kernel 网络层(纤维/δ/DNB 几何)。</li>
<li><span class="badge b-hyp">假设/示意</span> "高吸引子=临床发病态"映射、Schlögl 速率常数、合成 biomarker 序列、CNS 递送方案。</li>
<li><b>ipTM/binding_confidence 只验 pose/界面,不等于亲和力/中和/敲低效力</b> —— 均需 SPR/细胞/动物湿实验确证。</li>
<li>设计结果为 100 条中检出的 top(纳米抗体 n=15、小分子 n=10 已检视);全部 100 条在 Boltz 平台。</li>
</ul>

<p class="note">数据来源:Open Targets Platform GraphQL(ALS MONDO_0004976)· ClinicalTrials.gov v2 · STRING v12 · Ensembl REST · ChEMBL v34(ebselen CHEMBL51085)· Boltz-2.1(Modal GPU,CD Bioparticles)· pathway-biomarker-triangulation + denovo-drug-campaign skill。3D 结构为 Boltz-2.1 真实预测(.cif 内联,3Dmol.js 渲染)。本文为研究性内容,非临床或投资建议。</p>
"""

doc = f'''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>ALS TDP-43 从头多模态药物发现 · 正式交互报告</title>
<style>{CSS}</style>
<script>{threed}</script>
</head><body><div class="wrap">{BODY}</div></body></html>'''

out = D.parent/"als_report_interactive.html"
out.write_text(doc, encoding="utf-8")
print(f"wrote {out} : {out.stat().st_size} bytes ({out.stat().st_size//1024} KB)")
print("viewers: nanobody(SOD1+VHH ipTM0.886), smallmol(SOD1+LIG)")
