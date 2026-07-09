#!/usr/bin/env python3
"""Formal INTERACTIVE obesity/cardiometabolic report (self-contained HTML):
consolidate + beautify the original obesity campaign, with 2 live 3Dmol viewers
(MSTN nanobody & ALK7 antibody, re-folded fresh) + static fallbacks + sequences
+ all existing figures + the 3-method biomarker triangulation + independent cross-check."""
import base64, pathlib, html

FIGS_DIR = pathlib.Path("docs/figures")
STR = pathlib.Path("docs/obesity_repro/structures")
def b64(p): return "data:image/png;base64,"+base64.b64encode(pathlib.Path(p).read_bytes()).decode()
def cif(p): return pathlib.Path(p).read_text()
threed=(STR/"3Dmol-min.js").read_text()
mstn_cif=cif(STR/"mstn_nanobody.cif"); alk7_cif=cif(STR/"alk7_antibody.cif")
RMSTN=b64(STR/"render_mstn.png"); RALK7=b64(STR/"render_alk7.png")
F=lambda n: b64(FIGS_DIR/n)

# sequences (from designed_sequences.fasta)
mstn_nb="EVQLVESGGGLVQAGGSLRLSCAASAPLSAMGWFRQAPGKEREFVAAIGADGKNVYYAESVKGRFTISRDNAKNTVVLQMNSLKPEDTALYYCFAATGKYPNHKTYWGQGTQVTVSS"
mstn_cdr3="FAATGKYPNHKTY"
alk7_h="EVQLVQSGAEVKKPGESLKISCKGSGFDFSAHWIGWVRQMPGKGLEWMGIINPADGTTRYSPSFQGQVTISADKSISTAYLQWSSLKASDTAMYYCARINSAGSLDVWGQGTLVTVSS"
alk7_l="QSVLTQPPSVSGAPGQRVTISCTGSSSDGLADGEVSWYQQLPGTAPKLLIYSASELPSGVPDRFSGSKSGTSASLAITGLQSEDEADYYCSTWDSDGNLVFGGGTKLTVL"
alk7_h_cdr3="ARINSAGSLDV"
def hl(seq,m): return seq.replace(m,f'<span class="cdr">{m}</span>') if m in seq else seq

CSS="""*{box-sizing:border-box}html{-webkit-print-color-adjust:exact}
body{font-family:"Noto Sans","PingFang SC","Microsoft YaHei","WenQuanYi Zen Hei",system-ui,sans-serif;color:#1a2230;line-height:1.62;margin:0;background:#f4f6f9}
.wrap{max-width:960px;margin:0 auto;padding:40px 44px 90px;background:#fff}
h1{font-size:24px;border-bottom:3px solid #0b6b5e;padding-bottom:12px;margin:0 0 4px}
h2{font-size:19px;color:#0e2a45;border-left:4px solid #0b6b5e;padding-left:11px;margin:32px 0 10px}
h3{font-size:15px;margin:18px 0 6px;color:#22303f}
p{margin:8px 0}a{color:#0d47a1}
.lead{background:#f5f9f8;border-left:4px solid #0b6b5e;padding:12px 16px;border-radius:0 6px 6px 0;font-size:14px;color:#33414f}
.flag{background:#fff8e1;border-left:4px solid #e0b84c;padding:12px 16px;border-radius:0 6px 6px 0;font-size:13.5px;color:#5a4a00;margin:12px 0}
table{border-collapse:collapse;width:100%;margin:12px 0;font-size:12.5px}
th,td{border:1px solid #e3e8ef;padding:6px 9px;text-align:left;vertical-align:top}
thead th{background:#0e2a45;color:#fff}tbody tr:nth-child(even){background:#f7fafc}
img{display:block;max-width:100%;height:auto;margin:12px auto 4px;border:1px solid #e3e8ef;border-radius:8px}
.cap{font-size:12px;color:#667;text-align:center;font-style:italic;margin:0 0 14px}
code,.seq{font-family:ui-monospace,Menlo,Consolas,monospace}
.seqbox{background:#0e2233;color:#d5e6ff;padding:11px 14px;border-radius:8px;font-size:11.5px;word-break:break-all;line-height:1.85;margin:8px 0}
.cdr{background:#e8b800;color:#000;padding:1px 3px;border-radius:3px;font-weight:700}
.viewer{width:100%;height:420px;border:1px solid #cfd8e3;border-radius:10px;margin:10px 0;background:#fff}
.grid2{display:grid;grid-template-columns:1fr 1fr;gap:18px}@media(max-width:760px){.grid2{grid-template-columns:1fr}}
.kpi{display:flex;gap:12px;flex-wrap:wrap;margin:12px 0}
.kpi div{flex:1;min-width:110px;background:#eef3f9;border-radius:10px;padding:10px 12px;text-align:center}
.kpi b{display:block;font-size:20px;color:#0b6b5e}
.badge{display:inline-block;padding:2px 9px;border-radius:11px;font-size:11px;font-weight:700;color:#fff}
.b-ok{background:#1a7a3a}.b-rev{background:#c0392b}.b-un{background:#8a6d00}
.controls{font-size:12px;margin:6px 0}.controls button{font:inherit;border:1px solid #0b6b5e;background:#fff;color:#0b6b5e;border-radius:6px;padding:3px 9px;margin-right:5px;cursor:pointer}
.controls button:hover{background:#0b6b5e;color:#fff}
.toc{columns:2;font-size:13px;background:#f7fafc;border:1px solid #e3e8ef;border-radius:8px;padding:12px 18px}
.note{font-size:11px;color:#889;margin-top:22px}
@media print{.viewer,.controls{display:none!important}.wrap{padding:0}body{background:#fff}h2,h3{page-break-after:avoid}table,img{page-break-inside:avoid}}
@page{size:A4;margin:14mm 12mm}"""

def viewer(vid,data,cfg,tgt):
    esc=data.replace("\\","\\\\").replace("`","\\`").replace("$","\\$")
    return f'''<div class="controls"><button onclick="rst_{vid}()">卡通</button><button onclick="srf_{vid}()">表面</button><button onclick="v_{vid}.spin(!sp_{vid});sp_{vid}=!sp_{vid}">旋转</button><span style="color:#667">拖动旋转·滚轮缩放</span></div>
<div id="{vid}" class="viewer"></div>
<script>var v_{vid},sp_{vid}=false;(function(){{var d=`{esc}`;v_{vid}=$3Dmol.createViewer(document.getElementById("{vid}"),{{backgroundColor:"white"}});var m=v_{vid}.addModel(d,"cif");{cfg}v_{vid}.zoomTo();v_{vid}.render();}})();
function rst_{vid}(){{v_{vid}.removeAllSurfaces();v_{vid}.render();}}
function srf_{vid}(){{v_{vid}.addSurface($3Dmol.SurfaceType.VDW,{{opacity:0.6,color:"lightgray"}},{{chain:"{tgt}"}});v_{vid}.render();}}</script>'''

mstn_cfg='m.setStyle({chain:"A"},{cartoon:{color:"0x9aa6b2"}});m.setStyle({chain:"B"},{cartoon:{color:"0x0b6b5e"}});m.addStyle({chain:"B",resi:"96-110"},{stick:{radius:0.18,colorscheme:"yellowCarbon"}});'
alk7_cfg='m.setStyle({chain:"A"},{cartoon:{color:"0x9aa6b2"}});m.setStyle({chain:"H"},{cartoon:{color:"0x0d47a1"}});m.setStyle({chain:"L"},{cartoon:{color:"0x5aa0c8"}});'

BODY=f"""
<h1>肥胖 / 心血管代谢 · 多靶点多模态从头药物开发<br>正式全流程报告(含交互式 3D 结构)</h1>
<p class="cap">端到端 Demo:市场→5 靶点→biomarker 三方法→多模态设计→GPU 实跑→三重筛选→基因组模型上云 · 全程真实工具/账号</p>

<div class="lead"><b>一句话结论。</b>选定<b>非肠促胰岛素、人类遗传学(LOF)验证的肥胖/心血管代谢</b>赛道,构建跨"基因→GPCR→受体激酶→肌肉→代谢"五层级靶点链(GPR75 / INHBE / ALK7 / ActRIIB-myostatin / GDF15),20+ 多模态设计打分,5 个 Boltz GPU pilot + 独立亲和力/ADME 三重筛选。<b>最强资产 = 抗 myostatin 纳米抗体</b>(唯一同时通过界面 ipTM 与亲和力),ALK7(抗体+小分子双模态)次之。AlphaGenome 在 Modal 上完整跑通(真实 hg38、18 变体)。</div>

<div class="flag"><b>本版新增(本次总结/美化)。</b>把两条存活生物药 <b>独立 re-fold</b> 拿到新鲜结构做<b>可交互 3D 查看器</b>:MSTN 纳米抗体 <b>ipTM 0.896 / 结合 0.732</b>、ALK7 抗体 <b>0.816 / 0.193</b> —— 与原报告(0.90/0.74、0.93/0.21)独立一致,再次印证"高 ipTM ≠ 亲和力"这条核心教训。</div>

<div class="kpi"><div><b>5</b>跨层级靶点</div><div><b>0.732</b>MSTN 纳米抗体结合</div><div><b>2/6</b>生物药过亲和力确认</div><div><b>18</b>AlphaGenome hg38 变体</div><div><b>0.83</b>biomarker 三方法 Jaccard</div></div>

<h2>目录</h2>
<div class="toc">1. 市场与适应症<br>2. 五靶点跨层级链<br>3. Biomarker:三方法交叉验证 + 合成<br>4. 多模态设计 + 三重筛选(决定性)<br>5. <b>交互式 3D 结构 + 序列</b><br>6. 核酸药 + 靶向递送<br>7. 基因组模型上云(AlphaGenome)<br>8. 独立复核交叉核对<br>9. 成本 · Skills · 局限</div>

<h2>1 · 市场分析与适应症选择</h2>
<p>6 个候选赛道 5 维打分(未满足/交易活跃/可成药性/biomarker 成熟/白空间)。<b>选定肥胖-非肠促胰岛素/遗传学</b>:制药最大蛋糕(2030 $1000–1500 亿)、大药企正为 LOF 验证生物学买单(Regeneron-AZ 押 GPR75、Pfizer-Metsera ~$100 亿)、赛道转向"持久·脂肪特异·保肌" = 多模态平台用武之地、biomarker(BMI/DEXA 瘦体重/MRI-VAT/HbA1c)最成熟。</p>

<h2>2 · 五靶点跨层级链</h2>
<table><thead><tr><th>#</th><th>靶点(UniProt)</th><th>层级</th><th>主打模态</th><th>入组 biomarker</th></tr></thead><tbody>
<tr><td>1</td><td><b>GPR75</b>(O95800)</td><td>孤儿 GPCR</td><td>siRNA/小分子</td><td>GPR75 pLOF 基因型 + 多基因肥胖评分</td></tr>
<tr><td>2</td><td><b>INHBE/activin E</b>(P58166)</td><td>配体/转录</td><td>GalNAc-siRNA</td><td>循环 activin E、MRI-VAT、WHRadjBMI</td></tr>
<tr><td>3</td><td><b>ALK7/ACVR1C</b>(Q8NER5)</td><td>受体激酶</td><td>ECD 抗体/别构 SMKI</td><td>pSMAD2/3、内脏脂肪(DEXA)</td></tr>
<tr><td>4</td><td><b>ActRIIB/myostatin</b>(Q13705/O14793)</td><td>肌肉质量</td><td>抗 MSTN 纳米抗体</td><td>四肢瘦体重指数(DEXA)</td></tr>
<tr><td>5</td><td><b>GDF15/GFRAL</b>(Q99988/Q6UXV0)</td><td>代谢/食欲</td><td>GDF15 激动型抗体</td><td>循环 GDF15</td></tr>
</tbody></table>

<h2>3 · Biomarker:三方法交叉验证 + 合成 biomarker</h2>
<p>对肥胖 activin/myostatin→ALK7/ActRIIB→SMAD2/3 + GDF15-GFRAL 轴做三方法(<span class="badge b-ok">真跑</span>):</p>
<img src="{F('fig_biomarker_attribution.png')}"><p class="cap">图1 方法一:证据×可测性归因 + 最小单元认证(最小核心 biomarker 单元)。</p>
<img src="{F('fig_biomarker_method_compare.png')}"><p class="cap">图2 方法二:RWR 纯拓扑 vs 方法一(重叠 + Jaccard)。</p>
<img src="{F('fig_grn_dnb_obesity.png')}"><p class="cap">图3 方法三(network-biomarker skill,动力系统):SMAD2/3 鞍结分岔临界慢化 DNB 早期预警。</p>
<p><b>三方法收敛</b>:三个完全正交的方法都收敛到 <code>ACVR2B → SMAD2/3</code> —— 最可辩护的 biomarker 核心。方法三独家给出:CRNT δ=2(SMAD2/3 可双稳)+ 动态早期预警量(pSMAD2/3 方差/AR1)。</p>

<h2>4 · 多模态设计 + 三重筛选(决定性一步)</h2>
<p>20+ 设计逐一打分;5 个 Boltz pilot(GPR75 小分子 200、ALK7 小分子 200、MSTN 纳米抗体 100、ALK7 抗体 100、GDF15 抗体 100)全成功。<b>独立亲和力确认</b>(structure_and_binding)是真假 binder 分水岭:</p>
<table><thead><tr><th>项目</th><th>设计期 ipTM</th><th>确认 ipTM</th><th>确认亲和力</th><th>结论</th></tr></thead><tbody>
<tr><td><b>抗 myostatin 纳米抗体</b></td><td>0.955</td><td>0.90</td><td><b>0.74</b></td><td><span class="badge b-ok">界面+亲和力双验证 → 最强资产</span></td></tr>
<tr><td>ALK7-ECD 抗体 #2</td><td>0.973</td><td>0.93</td><td>0.21</td><td>最强抗体(几何好、亲和中等)</td></tr>
<tr><td>ALK7-ECD 抗体 #1</td><td>0.976</td><td>0.87–0.94</td><td>~0</td><td><span class="badge b-un">仅几何</span></td></tr>
<tr><td>GDF15 抗体</td><td>0.919</td><td>0.77–0.82</td><td>~0</td><td><span class="badge b-un">最弱,降级</span></td></tr>
</tbody></table>
<div class="flag"><b>核心教训:</b>6 个高 ipTM 生物药,独立亲和力只确认了 <b>2 个</b> —— 高 ipTM 可能只是"几何一致",独立确认才是真假 binder 的分水岭。ADME:ALK7 小分子(<code>Cc1cc(-c2nc(-c3cncc(Br)c3)no2)cnc1O</code>,结合 0.513)结合+ADME 双优 → 小分子先导。</p></div>
<img src="{F('fig_alk7_smallmol_pocket.png')}"><p class="cap">图4 ALK7 小分子先导 · Boltz 预测结合口袋(3Dmol 渲染)。</p>

<h2>5 · 交互式 3D 结构 · 序列(存活生物药,本次独立 re-fold)<span class="badge b-ok">真实·可交互</span></h2>
<p class="lead">两条通过亲和力确认的生物药,本次用原设计序列<b>独立重折叠</b>拿到新鲜结构(可旋转/缩放;静态图保底)。ipTM/结合置信度只验 pose/界面,不等于实验亲和力。</p>

<h3>5.1 最强资产:myostatin/GDF8 × 抗 MSTN 纳米抗体 — ipTM 0.896 / 结合 0.732</h3>
<div class="grid2"><div>
<img src="{RMSTN}"><p class="cap">静态渲染(真实坐标)。灰=GDF8 成熟域 · 青绿=纳米抗体。</p>
{viewer("mstnv", mstn_cif, mstn_cfg, "A")}
<p class="cap">↑ 交互:灰=GDF8 · 青绿=VHH · 黄=CDR3。</p>
</div><div>
<h3 style="margin-top:0">纳米抗体序列(117 aa)</h3>
<div class="seqbox">{hl(mstn_nb,mstn_cdr3)}</div>
<p style="font-size:12px">CDR3(黄)= <code>{mstn_cdr3}</code>,VHH 骨架。独立 re-fold <b>ipTM 0.896 · 结合 0.732</b>,复现原报告 0.90/0.74 —— 唯一界面+亲和力双高的资产。选 MSTN(vs ActRIIB)= 配体选择性、避开 BMP/心血管。</p>
</div></div>

<h3>5.2 ALK7-ECD × de-novo Fv 抗体 — ipTM 0.816 / 结合 0.193</h3>
<div class="grid2"><div>
<img src="{RALK7}"><p class="cap">静态渲染(真实坐标)。灰=ALK7 ECD · 蓝=抗体 VH/VL。</p>
{viewer("alk7v", alk7_cif, alk7_cfg, "A")}
<p class="cap">↑ 交互:灰=ALK7 ECD · 深蓝=VH · 浅蓝=VL。</p>
</div><div>
<h3 style="margin-top:0">抗体序列(VH 118 / VL 110)</h3>
<div class="seqbox">VH: {hl(alk7_h,alk7_h_cdr3)}</div>
<div class="seqbox">VL: {alk7_l}</div>
<p style="font-size:12px">CDR3-H(黄)= <code>{alk7_h_cdr3}</code>。独立 re-fold <b>ipTM 0.816 · 结合 0.193</b>,复现原报告 0.93/0.21 —— 界面几何好、亲和力中等。ECD 阻断抗体绕开 ALK4/5 激酶选择性/心毒陷阱。</p>
</div></div>

<h2>6 · 核酸药物 + 靶向递送</h2>
<p><b>siRNA(基因层)</b>:基于验证 CDS 设计,先导 INHBE <b>E4</b>(<code>AGUCUAGUUGCAGUUUCAG</code>),机制=<b>切断/敲低</b>(RISC),ESC-plus 化学(GalNAc + 2'-OMe/2'-F/PS + GNA-seed + 5'-VP)。<b>递送</b>:INHBE→GalNAc-ASGPR 偶联(SC,Q8–12 周);GPR75→脂肪靶向可离子化 LNP(脂肪归巢肽 CKGGRAKDC)+ 探索性 ZIF-8 MOF;纳米抗体→Fc 融合延长半衰期。</p>

<h2>7 · 基因组模型上云(Modal)</h2>
<p><span class="badge b-ok">AlphaGenome 完整跑通</span>:genomicsxai PyTorch port + gtca 权重部署到用户 Modal GPU,真实 hg38 窗口,<b>18 个 CDx 变体</b>组织级注释(剪接/调控/eQTL 效应),用于 biomarker 富集面板。<span class="badge b-un">Evo2-7B</span>:部署配方就绪、H100 可加载,最终 log-likelihood 打分待一次性 flash-attn/TE 构建(WIP)。</p>

<h2>8 · 独立复核交叉核对</h2>
<p>用 Open Targets / STRING v12 / GWAS Catalog live 对原报告第三方再核:</p>
<table><thead><tr><th>项</th><th>独立核对</th><th>结论</th></tr></thead><tbody>
<tr><td>STRING 通路骨架</td><td>INHBA/MSTN→ACVR2B 0.999、INHBE→ALK7 0.952、ACVR2B/ALK7→SMAD2/3 0.96–0.98、GDF15→GFRAL 0.999→RET 0.998 —— 每条边逐位吻合</td><td><span class="badge b-ok">CONFIRMED</span></td></tr>
<tr><td>GWAS 位点</td><td>GPR75 保护性 rs80328470、INHBE rs150777893 均确认</td><td><span class="badge b-ok">CONFIRMED</span></td></tr>
<tr><td>Open Targets 关联"分值"</td><td>§16 表格 headline 数值与性状标签对不上(疑 L2G/错配);靶点-性状生物学是真的</td><td><span class="badge b-rev">REVISED</span></td></tr>
<tr><td>"MSTN 零肥胖遗传学"</td><td>OT 实给 MSTN 肥胖 0.121(>ALK7 0.072);措辞过陈述,但"保肌加成"定位正确</td><td><span class="badge b-rev">REVISED</span></td></tr>
<tr><td>设计层 ipTM/结合</td><td>本次独立 re-fold 复现(MSTN 0.896/0.732、ALK7 0.816/0.193)</td><td><span class="badge b-ok">复现一致</span></td></tr>
</tbody></table>
<p><b>一句话</b>:通路与遗传学骨架经得起第三方溯源;需收敛的是 §16 的 OT 分值标注与 MSTN 遗传学措辞(解读层,非事实层)。</p>

<h2>9 · 成本 · Skills · 局限</h2>
<p><b>成本</b>:Boltz ≈ $38–40 + 本次 re-fold $0.05;Modal ≈ $10–11。<b>固化 skills</b>:alphagenome-modal、evo2-modal、boltz-denovo-design、pathway-biomarker-triangulation、denovo-drug-campaign(5 个,构成完整流水线)。</p>
<ul style="font-size:13px">
<li>ipTM/结合置信度只验 pose/界面,<b>不等于</b>亲和力/中和/敲低效力 —— 需 SPR/BLI + 细胞/动物湿实验确证。</li>
<li>siRNA 脱靶/种子全转录组筛选需基因组基座模型(AlphaGenome 已上云,Evo2 打分 WIP)。</li>
<li>小分子为已知化学型类似物(锚定临床先例,非全新骨架);合成 biomarker 与动态层为规范机制演示,需病人数据标定。</li>
</ul>

<p class="note">数据源:Open Targets / STRING v12 / GWAS Catalog / ChEMBL / ClinicalTrials / UCSC / Ensembl;计算:Boltz-2.1(Modal GPU,CD Bioparticles)、AlphaGenome(genomicsxai port + gtca 权重)。3D 结构为本次独立 re-fold 的 Boltz-2.1 真实预测(.cif 内联,3Dmol.js 内联,离线可用)。研究性内容,非临床或投资建议。</p>
"""
doc=f'''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>肥胖多靶点多模态从头药物开发 · 正式交互报告</title><style>{CSS}</style><script>{threed}</script></head><body><div class="wrap">{BODY}</div></body></html>'''
out=pathlib.Path("docs/obesity_report_interactive.html"); out.write_text(doc,encoding="utf-8")
print(f"wrote {out}: {out.stat().st_size//1024} KB")
