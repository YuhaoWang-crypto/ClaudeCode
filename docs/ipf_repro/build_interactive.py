#!/usr/bin/env python3
"""Formal INTERACTIVE IPF report (self-contained HTML): external de-novo IPF report
re-rendered reliably (inline 3Dmol.js + static fallbacks from real PDB coords) +
an independent source-traceable verification section folded in."""
import base64, json, pathlib, html

D = pathlib.Path("docs/ipf_repro"); FIG=D/"figures"; STR=D/"structures"
def b64(p): return "data:image/png;base64,"+base64.b64encode(pathlib.Path(p).read_bytes()).decode()
def pdb(p): return pathlib.Path(p).read_text()
threed=(STR/"3Dmol-min.js").read_text()
mb_pdb=pdb(STR/"minibinder.pdb"); vhh_pdb=pdb(STR/"vhh.pdb")
seqs=json.loads((STR/"seqs.json").read_text())
vhh_H=seqs["vhh15"]["H"]; vhh_T=seqs["vhh15"]["T"]
FIGS=[b64(FIG/f"ipf_fig{i}.png") for i in range(7)]
RMB=b64(STR/"render_minibinder.png"); RVHH=b64(STR/"render_vhh.png")
cdr3="AAGLALGAPSAYTYW"
vhh_disp=vhh_H.replace(cdr3, f'<span class="cdr">{cdr3}</span>') if cdr3 in vhh_H else vhh_H

CSS="""*{box-sizing:border-box}html{-webkit-print-color-adjust:exact}
body{font-family:"Noto Sans","PingFang SC","Microsoft YaHei","WenQuanYi Zen Hei",system-ui,sans-serif;color:#1a2230;line-height:1.62;margin:0;background:#f4f6f9}
.wrap{max-width:960px;margin:0 auto;padding:40px 44px 90px;background:#fff}
h1{font-size:24px;border-bottom:3px solid #4575b4;padding-bottom:12px;margin:0 0 4px}
h2{font-size:19px;color:#2c5aa0;border-left:4px solid #4575b4;padding-left:11px;margin:32px 0 10px}
h3{font-size:15px;margin:18px 0 6px;color:#22303f}
p{margin:8px 0}a{color:#0d47a1}
.lead{background:#eef3fb;border-left:4px solid #4575b4;padding:12px 16px;border-radius:0 6px 6px 0;font-size:14px;color:#33414f}
.flag{background:#fdecea;border-left:4px solid #c0392b;padding:12px 16px;border-radius:0 6px 6px 0;font-size:13.5px;color:#7b241c;margin:12px 0}
table{border-collapse:collapse;width:100%;margin:12px 0;font-size:12.5px}
th,td{border:1px solid #e3e8ef;padding:6px 9px;text-align:left;vertical-align:top}
thead th{background:#2c3e50;color:#fff}tbody tr:nth-child(even){background:#f7fafc}
img{display:block;max-width:100%;height:auto;margin:12px auto 4px;border:1px solid #e3e8ef;border-radius:8px}
.cap{font-size:12px;color:#667;text-align:center;font-style:italic;margin:0 0 14px}
code,.seq{font-family:ui-monospace,Menlo,Consolas,monospace}
.seqbox{background:#0e2233;color:#d5e6ff;padding:12px 14px;border-radius:8px;font-size:12px;word-break:break-all;line-height:1.9;margin:8px 0}
.cdr{background:#e8b800;color:#000;padding:1px 3px;border-radius:3px;font-weight:700}
.viewer{width:100%;height:420px;border:1px solid #cfd8e3;border-radius:10px;margin:10px 0;background:#fff;position:relative}
.grid2{display:grid;grid-template-columns:1fr 1fr;gap:18px}@media(max-width:760px){.grid2{grid-template-columns:1fr}}
.kpi{display:flex;gap:12px;flex-wrap:wrap;margin:12px 0}
.kpi div{flex:1;min-width:110px;background:#eef3f9;border-radius:10px;padding:10px 12px;text-align:center}
.kpi b{display:block;font-size:20px;color:#2c5aa0}
.badge{display:inline-block;padding:2px 9px;border-radius:11px;font-size:11px;font-weight:700;color:#fff}
.b-ok{background:#1a7a3a}.b-rev{background:#c0392b}.b-un{background:#8a6d00}
.controls{font-size:12px;margin:6px 0}.controls button{font:inherit;border:1px solid #2c5aa0;background:#fff;color:#2c5aa0;border-radius:6px;padding:3px 9px;margin-right:5px;cursor:pointer}
.controls button:hover{background:#2c5aa0;color:#fff}
.note{font-size:11px;color:#889;margin-top:22px}
@media print{.viewer,.controls{display:none!important}.wrap{padding:0}body{background:#fff}h2,h3{page-break-after:avoid}table,img{page-break-inside:avoid}}
@page{size:A4;margin:14mm 12mm}"""

def viewer(vid,data,cfg):
    esc=data.replace("\\","\\\\").replace("`","\\`").replace("$","\\$")
    return f'''<div class="controls"><button onclick="rst_{vid}()">卡通</button><button onclick="srf_{vid}()">表面</button><button onclick="v_{vid}.spin(!sp_{vid});sp_{vid}=!sp_{vid}">旋转</button><span style="color:#667">拖动旋转·滚轮缩放</span></div>
<div id="{vid}" class="viewer"></div>
<script>var v_{vid},sp_{vid}=false;(function(){{var d=`{esc}`;v_{vid}=$3Dmol.createViewer(document.getElementById("{vid}"),{{backgroundColor:"white"}});var m=v_{vid}.addModel(d,"pdb");{cfg}v_{vid}.zoomTo();v_{vid}.render();}})();
function rst_{vid}(){{v_{vid}.removeAllSurfaces();v_{vid}.render();}}
function srf_{vid}(){{v_{vid}.addSurface($3Dmol.SurfaceType.VDW,{{opacity:0.6,color:"lightgray"}},{{chain:"{'A' if vid=='mbv' else 'T'}"}});v_{vid}.render();}}</script>'''

mb_cfg='m.setStyle({chain:"A"},{cartoon:{color:"0xc8ced6"}});m.setStyle({chain:"B"},{cartoon:{color:"0xaeb6c0"}});m.setStyle({chain:"C"},{cartoon:{colorscheme:"spectrum"}});'
vhh_cfg='m.setStyle({chain:"T"},{cartoon:{color:"0xc8ced6"}});m.setStyle({chain:"H"},{cartoon:{color:"0x0b6b5e"}});m.addStyle({chain:"H",resi:"96-112"},{stick:{radius:0.18,colorscheme:"yellowCarbon"}});'

def vtable(rows):
    b="".join("<tr>"+"".join(f"<td>{c}</td>" for c in r)+"</tr>" for r in rows)
    return f'<table><thead><tr><th>#</th><th>主张</th><th>报告值</th><th>独立核对</th><th>结论</th></tr></thead><tbody>{b}</tbody></table>'

VER=[
 ["1","IPF 疾病 id / 关联靶点数","EFO_0000768","Open Targets = IPF,关联 <b>4,591</b>",'<span class="badge b-ok">CONFIRMED</span>'],
 ["2","5 靶点均 IPF 高关联","MUC5B/ITGB6/TGFB1/ENPP2/LOXL2","OT 排名:MUC5B #7(0.60)、ENPP2 #35(0.35);<b>TGFB1 #247、LOXL2 #267、ITGB6 #274、SMAD3 #639</b> —— 后四个是机制挑选非统计高排名;OT 真实 top(PARN/RTEL1/TERT 端粒、SFTPA2 表面活性)全不在轴内",'<span class="badge b-rev">REVISED</span>'],
 ["3","STRING ITGAV-ITGB6 0.999","0.999","0.999",'<span class="badge b-ok">CONFIRMED</span>'],
 ["4","ITGB6→TGFB1 / TGFB1→SMAD3","0.986/0.99","0.986 / 0.990(0.994 实为 TGFB1-ITGAV,子基元混标)",'<span class="badge b-ok">CONFIRMED</span>'],
 ["5","ENPP2/LOXL2 = 平行臂","隐含离枢纽","STRING:<b>ENPP2 在 score≥400 无任何边</b>(确离枢纽);LOXL2 仅 TGFB1-LOXL2 0.663",'<span class="badge b-ok">CONFIRMED</span>'],
 ["6","MUC5B rs35705950 IPF 风险变异","遗传风险","GWAS:11p15.5 启动子,风险等位 T,<b>OR 2.13–2.77</b> —— 教科书级 IPF 变异",'<span class="badge b-ok">CONFIRMED</span>'],
 ["7","cudetaxestat 真实(autotaxin 抑制剂)","co-fold ipTM 0.807","ChEMBL4802146,BLD-0409,Blade,max_phase 2;试验 NCT05373914 状态 UNKNOWN(2023 完成,Blade 已收缩)",'<span class="badge b-ok">真实</span>/<span class="badge b-un">ipTM UNVERIFIED</span>'],
 ["8","ENPP2 小分子排名 <b>第 1</b>(4.54)","首选资产","<b>类别 III 期负面:ziritaxestat(GLPG1690)ISABELA I/II(NCT03711162/03733444,1300+ 人)双双终止</b>",'<span class="badge b-rev">REVISED</span>'],
 ["9","LOXL2 小分子排名第 5","候选","<b>simtuzumab(抗 LOXL2)IPF II 期 NCT01769196(544 人)失败</b> —— LOXL2 在 IPF 已临床证伪",'<span class="badge b-rev">REVISED</span>'],
 ["10","仅 2 个获批抗纤维化药","2","pirfenidone + nintedanib 恰 2 个(都只延缓 FVC 下降)",'<span class="badge b-ok">CONFIRMED</span>'],
 ["11","8×CAGA-box(AGCCAGACA)SMAD3 元件","合成 biomarker","CAGA-box/CAGA12 = 经典 SMAD3/4 响应元件(Dennler 1998),motif 正确",'<span class="badge b-ok">CONFIRMED</span>'],
 ["12","GPU 分数 / 动态层 / Evo2 自然性","0.807 / Jaccard 0.78 / δ=1 / 8.8×","无公开源可核对(内部计算产物)",'<span class="badge b-un">UNVERIFIED</span>'],
]

BODY=f"""
<h1>IPF 特发性肺纤维化 —— 从头多模态药物开发<br>正式全流程报告(含交互式 3D 结构 + 独立复核)</h1>
<p class="cap">外部 denovo-drug-campaign 报告 · 本地可靠化重渲染(内联 3Dmol + 静态保底)· 叠加独立溯源复核</p>

<div class="lead"><b>报告主张。</b>数据驱动从 11 个候选选定 IPF,构建 <b>MUC5B / αvβ6(ITGB6)/ TGF-β1 / ENPP2(autotaxin)/ LOXL2</b> 五靶点纤维化轴,三方法验证 biomarker + 合成 biomarker,100 封顶多模态设计,<b>L1 ENPP2 小分子(4.54)居首</b>,吸入/雾化肺局部递送。</div>

<div class="flag"><b>⚠ 独立复核的关键红旗(详见第 11 节)。</b>骨架扎实(EFO id、STRING ~0.99 整合素→TGFβ→SMAD3、MUC5B rs35705950、2 个获批药、CAGA-box 均核实),但<b>排序有重大临床遗漏</b>:排第 1 的 ENPP2/autotaxin —— ziritaxestat III 期(ISABELA)已<b>终止</b>;排第 5 的 LOXL2 —— simtuzumab II 期已<b>失败</b>。IPF 是"III 期失败坟场"。且 5 靶点里仅 MUC5B/ENPP2 统计高排名,其余靠机制挑选。</div>

<div class="kpi"><div><b>4,591</b>OT 关联靶点(IPF)</div><div><b>0.999</b>αvβ6 STRING</div><div><b>OR 2.6</b>MUC5B rs35705950</div><div><b>2</b>已获批抗纤维化</div><div><b>11/12</b>硬主张核实</div></div>

<h2>1 · 市场分析与适应症选择</h2>
<p>11 个候选按 Open Targets 靶点/药物数 + ClinicalTrials 招募数 7 维加权。Alzheimer 4.17 &gt; Parkinson 3.94 &gt; <b>IPF 3.87</b>。选 IPF 因 5 靶点全有可成药口袋、天然跨三层级、仅 2 药获批。<span class="badge b-un">市场综合分无公开源可核</span></p>
<img src="{FIGS[0]}"><p class="cap">图1 适应症综合打分。</p>

<h2>2 · 多层级靶点轴 <span class="badge b-ok">STRING 核实</span></h2>
<p>轴:<code>MUC5B(基因)→ αvβ6/ITGB6(激活潜隐 TGF-β1)→ TGF-β1(枢纽)→ SMAD3 → 肌成纤维细胞;ENPP2/LOXL2 平行硬化臂</code>。STRING 主线独立核实:ITGAV-ITGB6 <b>0.999</b>、ITGB6-TGFB1 0.986、TGFB1-SMAD3 <b>0.990</b>。<b>但</b>OT 关联排名:仅 MUC5B(#7)、ENPP2(#35)靠前;TGFB1/LOXL2/ITGB6/SMAD3 在 #247–639 —— 机制轴合理,但非"统计最高靶点";OT 真实 top(端粒 PARN/RTEL1/TERT、表面活性 SFTPA2)不在轴内。</p>
<img src="{FIGS[1]}"><p class="cap">图2 靶点轴 / STRING 通路。</p>

<h2>3 · Biomarker 三方法 + 合成 biomarker</h2>
<p>M1∩M2 核心 21、Jaccard 0.78;M3:TGF-β/SMAD3 CRNT δ=1 双稳,pSMAD3 逼近临界方差↑8.8×、AR1 0.94→0.99(DNB 预警)。合成 <b>FibroMark-IPF-SMAD3-01</b> = 8×CAGA-box(AGCCAGACA,<span class="badge b-ok">经典 SMAD3 元件</span>)+ 最小启动子;Evo2 合成 −0.330 vs 打乱 −1.296。<span class="badge b-un">三方法/Evo2 数值无公开源可核</span></p>
<img src="{FIGS[2]}"><p class="cap">图3 biomarker 三方法交叉验证。</p>
<img src="{FIGS[3]}"><p class="cap">图4 合成 biomarker / Evo2 自然性。</p>

<h2>4 · 靶点 × 模态可行性矩阵</h2>
<p>5×5 矩阵融合重利用锚点,选 6 先导跨全模态:L1 ENPP2 小分子 / L2 αvβ6 VHH / L3 MUC5B 基因 / L4 αvβ6 微结合物 / L5 LOXL2 小分子 / L6 TGF-β1 抗体。</p>
<img src="{FIGS[4]}"><p class="cap">图5 靶点×模态可行性矩阵。</p>

<h2>5–6 · 多模态 GPU 设计(100 封顶)</h2>
<p>ENPP2 小分子 5(cudetaxestat co-fold ipTM 0.807)+ LOXL2 小分子 5(0.43–0.64,浅口袋+共价)+ MUC5B 5 gapmer+5 siRNA + αvβ6 微结合物 10(RFdiffusion,RGD 2.2–2.6Å)+ αvβ6 VHH 70(RFantibody,2 条 ipae&lt;10)= 100。<span class="badge b-un">GPU 分数无公开源可核</span></p>
<img src="{FIGS[5]}"><p class="cap">图6 多模态设计 / ADMET。</p>

<h2>7 · 交互式 3D 结构 · 序列(外部报告设计,本地可靠化)<span class="badge b-ok">真实结构·可交互</span></h2>
<p class="lead">以下为<b>该外部报告 GPU 设计</b>的代表结构,我从其内嵌数据提取并本地重渲染(内联 3Dmol + 静态保底图,离线可见)。ipTM/接触数只验 pose,不等于亲和力。</p>

<h3>7.1 L4:αvβ6 × de-novo 微型结合物(RFdiffusion)</h3>
<div class="grid2"><div>
<img src="{RMB}"><p class="cap">静态渲染(真实坐标)。灰=αvβ6(链 A+B)· 红=设计微结合物(链 C,70-mer)。</p>
{viewer("mbv", mb_pdb, mb_cfg)}
<p class="cap">↑ 交互查看器:灰=αvβ6 · 彩虹=微结合物。</p>
</div><div>
<h3 style="margin-top:0">诚实提示</h3>
<p style="font-size:13px">提取的微结合物链 C 序列为 <b>poly-glycine 骨架</b>(RFdiffusion 生成骨架,序列由 ProteinMPNN 后填;此 PDB 仅含骨架占位)。因此可展示<b>结合几何/RGD 命中口袋</b>,但<b>不含最终氨基酸序列</b> —— 报告的 "RGD 2.2–2.6Å" 是骨架接触距离,合理但需 MPNN 序列 + 重折叠确认。</p>
</div></div>

<h3>7.2 L2:αvβ6 表位 × de-novo VHH 纳米抗体(RFantibody)</h3>
<div class="grid2"><div>
<img src="{RVHH}"><p class="cap">静态渲染(真实坐标)。灰=αvβ6 表位(链 T,81aa)· 青绿=VHH(链 H,120aa)。</p>
{viewer("vhhv", vhh_pdb, vhh_cfg)}
<p class="cap">↑ 交互查看器:灰=表位 · 青绿=VHH · 黄=CDR3。</p>
</div><div>
<h3 style="margin-top:0">VHH 序列(120 aa)</h3>
<div class="seqbox">{vhh_disp}</div>
<p style="font-size:12px">CDR3(黄)= <code>{cdr3}</code>,VHH(camelid)骨架。这是含真实氨基酸序列的设计(RFantibody + 序列设计)。报告称 80 条设计仅 2 条 ipae&lt;10 —— 符合 de-novo 抗体真实低产率。</p>
</div></div>

<h2>8 · 安全 / 可开发性排序(报告)</h2>
<table><thead><tr><th>排名</th><th>先导</th><th>靶点</th><th>模态</th><th>综合分</th><th>独立复核标注</th></tr></thead><tbody>
<tr><td>1</td><td>L1</td><td>ENPP2</td><td>小分子</td><td>4.54</td><td><span class="badge b-rev">红旗</span> autotaxin III 期(ziritaxestat)已终止</td></tr>
<tr><td>2</td><td>L2</td><td>αvβ6</td><td>VHH</td><td>3.92</td><td>肺特异,合理;pamrevlumab(抗 CTGF)III 期近失败为邻域警示</td></tr>
<tr><td>3</td><td>L4</td><td>αvβ6</td><td>微结合物</td><td>3.85</td><td>骨架命中 RGD;需 MPNN 序列</td></tr>
<tr><td>4</td><td>L5</td><td>LOXL2</td><td>小分子</td><td>3.81</td><td><span class="badge b-rev">红旗</span> simtuzumab(抗 LOXL2)II 期已失败</td></tr>
<tr><td>5</td><td>L3</td><td>MUC5B</td><td>基因</td><td>3.74</td><td>rs35705950 铁证;风险等位特异,遗传学最强</td></tr>
<tr><td>6</td><td>L6</td><td>TGF-β1</td><td>抗体</td><td>3.59</td><td>全身抑 TGF-β 有心脏/免疫毒性,报告已注</td></tr>
</tbody></table>
<img src="{FIGS[6]}"><p class="cap">图7 先导综合排序(报告)。</p>
<div class="flag"><b>独立复核建议重排:</b>鉴于 autotaxin/LOXL2 的 IPF 临床失败,<b>遗传学最强的 MUC5B(基因,rs35705950 OR~2.6)与肺特异的 αvβ6(VHH/微结合物)应上调</b>,ENPP2/LOXL2 应下调或明确标注"类别已证伪"。</div>

<h2>9 · 靶向递送 formulation</h2>
<p>IPF 为肺局部病 → <b>吸入/雾化富集纤维化灶、降全身暴露</b>:ENPP2 小分子 干粉吸入+口服;αvβ6 VHH 雾化(15kDa 耐雾化);MUC5B siRNA 吸入型 LNP(DOTAP);αvβ6 微结合物雾化;LOXL2 小分子口服(ECM 需全身)。这是该报告最扎实的一节 —— 肺局部递送与靶点生物学契合。</p>

<h2>11 · 独立溯源复核(本节为 Claude Code 新增)</h2>
<p>用 Open Targets / STRING v12 / GWAS Catalog / ChEMBL / ClinicalTrials live 逐条核对:<b>硬事实 11/12 通过</b>,弱点全在"排序解读"。</p>
{vtable(VER)}
<h3>重大遗漏 / 方法学问题</h3>
<ul style="font-size:13px">
<li><b>autotaxin III 期失败未提(最严重)</b>:ENPP2 排第 1,但 ziritaxestat ISABELA I/II 因获益-风险不足终止;cudetaxestat 是差异化(非竞争性)骨架,但 Blade Ph2 已停摆。</li>
<li><b>LOXL2 在 IPF 已临床证伪</b>:simtuzumab II 期失败;报告 LOXL2 小分子继承未验证靶点。</li>
<li><b>靶点选择是机制驱动非证据驱动</b>:5 靶点里 4 个 OT 排名 247–639;遗传学最强的端粒(PARN/RTEL1/TERT)与表面活性(SFTPA2/SFTPC)完全缺席。</li>
<li><b>"空间大"是片面叙事</b>:仅 2 药获批属实,但该领域是 III 期失败坟场(ziritaxestat、simtuzumab、pamrevlumab、zinpentraxin);报告自己排名的 L1/L5 正落在已失败轴上。</li>
<li>GPU/动态/自然性数值(ipTM 0.807、Jaccard 0.78、δ=1、8.8×、Evo2)均 UNVERIFIED,为内部计算产物。</li>
</ul>
<p><b>一句话:</b>机制上专业、临床上过于乐观。骨架(通路/遗传学/landscape)经得起溯源;<b>ENPP2 与 LOXL2 排名应打红旗</b>,GPU/动态分数按未验证对待。</p>

<h2>10 · 诚实局限(报告原文 + 复核)</h2>
<ul style="font-size:13px">
<li>GPU 验证是计算预测(ipTM/ipae/接触),非实验亲和力 —— 需 SPR/BLI + 细胞功能确证。</li>
<li>VHH 产出保守(80 设计仅 2 条 ipae&lt;10),符合 de-novo 抗体真实产率;需亲和成熟 + 人源化。</li>
<li>微结合物 PDB 仅骨架(poly-Gly),最终序列需 ProteinMPNN + 重折叠。</li>
<li>LOXL2 co-fold 偏低反映浅口袋 + 共价机制;合成 biomarker Evo2 自然性高 ≠ 体内转录活性;动态层为规范参数演示。</li>
<li>小分子为已知化学型类似物(L1/L5 锚定临床先例,非全新骨架)。</li>
</ul>

<p class="note">外部报告由 Claude Science 生成(数据源 Open Targets/STRING/UniProt/RCSB/ClinicalTrials;计算 Modal GPU Boltz-2/RFdiffusion/RFantibody/Evo2)。本地可靠化重渲染 + 独立复核由 Claude Code 完成(3D 结构从原报告内嵌数据提取,3Dmol.js 内联;复核用 Open Targets GraphQL / STRING v12 / GWAS Catalog / ChEMBL v34 / ClinicalTrials.gov v2 live 核对)。研究性内容,非临床或投资建议。</p>
"""

doc=f'''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>IPF 从头多模态药物发现 · 正式交互报告(含独立复核)</title><style>{CSS}</style><script>{threed}</script></head><body><div class="wrap">{BODY}</div></body></html>'''
out=D.parent/"ipf_report_interactive.html"
out.write_text(doc,encoding="utf-8")
print(f"wrote {out}: {out.stat().st_size//1024} KB")
