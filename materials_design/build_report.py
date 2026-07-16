# -*- coding: utf-8 -*-
"""生成自包含 HTML 报告(内嵌所有图为 base64)。"""
import base64, os
D = "/home/user/ClaudeCode/materials_design"
def img(fn, w="100%"):
    p = os.path.join(D, fn)
    if not os.path.exists(p): return f"<p><i>[缺图 {fn}]</i></p>"
    b64 = base64.b64encode(open(p,"rb").read()).decode()
    return f'<img src="data:image/png;base64,{b64}" style="width:{w};max-width:1100px;border:1px solid #e2e4e8;border-radius:8px;margin:10px 0"/>'

SECTIONS = f"""
<h1>计算材料设计与验证 — 项目报告</h1>
<p class="sub">FairChem UMA · MatterGen · GPAW · Quantum ESPRESSO，全部在 Modal GPU 上运行。
每条结论都标注 ✅已验证 / ⚠️限制。</p>

<div class="toc"><b>目录</b>
<ol>
<li>三个独立设计项目（量子点 / 上转换 / 超导）</li>
<li>模拟 vs 真实实验数据</li>
<li>同族替换 &amp; 掺杂的性质模拟</li>
<li>真实 DFT 电子结构（n 型能级）</li>
<li>高通量筛选（描述符能行也会失效）</li>
<li>生成式逆向设计 + 多描述符 + DFPT 电声耦合</li>
<li>方法论与诚实边界</li>
</ol></div>

<h2>1 · 三个独立设计项目</h2>
<p>三个各自独立的设计任务（<b>不是</b>同一材料）。每个给出材料/结构/机理/验证。</p>
<table>
<tr><th>项目</th><th>代表材料</th><th>结构</th><th>为什么能实现</th></tr>
<tr><td>量子点</td><td>CsPbBr₃ / CdSe / PbS / Si</td><td>Pm-3m / P6₃mc / Fm-3m / Fd-3m</td><td>尺寸&lt;激子玻尔半径→三维限域→分立可调能级(∝1/R²)</td></tr>
<tr><td>上转换</td><td>β-NaYF₄:Yb,Er</td><td>P6₃/m 氟化物</td><td>Yb敏化→Er 4f阶梯能级ETU;低声子基质抑制多声子猝灭</td></tr>
<tr><td>超导</td><td>MgB₂ / Nb₃Sn / Nb</td><td>P6/mmm / A15 / bcc</td><td>离域电子经电声耦合配库珀对;MgB₂靠σ带空穴+E₂g声子</td></tr>
</table>

<h2>2 · 模拟 vs 真实实验数据（同图对比）</h2>
<p>线/柱=模拟，点=实验。(A)量子点尺寸-发光(Brus vs 实验尺寸曲线)；(B)上转换基质声子氟化物&lt;氧化物；
(C)超导 MgB₂ E₂g；(D)13 材料晶格 parity，<b>MAE=0.057 Å (1.2%)</b>。</p>
{img("sim_vs_experiment.png")}

<h2>3 · 同族替换 &amp; 掺杂的性质模拟</h2>
<p><b>Demo A｜同族替换 Ca→Sr→Ba</b>：晶格↑、M–F键↑、声子↓、体模量↓，UMA vs 实验高度吻合。</p>
{img("demoA_substitution.png")}
<p><b>Demo B｜掺杂 Nb:SrTiO₃</b>：Vegard 定律，晶格随 Nb 线性膨胀（Nb⁵⁺&gt;Ti⁴⁺），电子学后果是 n 型载流子。</p>
{img("demoB_doping.png")}
<p><b>①同族替换扩展</b>：CsPbBr₃→CsSnBr₃(Pb→Sn) 与 Nb₃Sn→V₃Si(A15)。</p>
{img("step2_substitution.png")}

<h2>4 · 真实 DFT 电子结构（n 型能级）</h2>
<p>GPAW/PBE：纯 SrTiO₃ 绝缘体(带隙 1.68 eV，E_F 在带隙)；Nb:SrTiO₃ E_F 进导带→<b>简并 n 型</b>。
这是 UMA 给不了的电子结构。</p>
{img("step3_dft_dos.png")}

<h2>5 · 高通量筛选：描述符能行也会失效</h2>
<p>枚举替换→算描述符→排序。(A)上转换基质按声子排序，<b>描述符有效</b>，氟化物全优于氧化物；
(B)(C)硼化物按声子排序，<b>MgB₂(唯一39K超导)排倒数第3，声子 vs Tc 无相关→单一描述符失效</b>。</p>
{img("screening.png")}

<h2>6 · 逆向设计 + 多描述符 + DFPT</h2>
<p><b>#1 生成式逆向设计</b>：只给 MatterGen "体模量=300 GPa"，生成的全是已知超硬家族(重金属硼化物+锇族合金)，
UMA 独立验证 291–589 GPa，全部≥目标。</p>
{img("mattergen_inverse_design.png")}
<p><b>#2 多描述符筛选</b>：DFT 的 B-p(σ)占比把范围缩小到 sp 硼化物 σ 家族，但 YB₂ 也高却不超导→仍需 λ。</p>
{img("screening_2descriptor.png")}
<p><b>#3 QE DFPT 电声耦合</b>：MgB₂ 声子 E₂g=485 cm⁻¹(驱动超导的模，谐波LDA)；
Allen-Dynes 用 DFPT ω_log(627K)，39K 对应 λ≈0.9-1.0。⚠️ λ 未完全从头收敛(用文献 λ 演示公式)。</p>
{img("step3b_qe_epc.png")}

<h2>7 · 方法论与诚实边界</h2>
<p class="method"><b>流程</b>：生成式逆向设计(按性质造) → 廉价描述符粗筛(声子/弹性,UMA) →
电子结构描述符缩小范围(N(E_F)/σ占比,DFT) → DFPT 精算 Tc(电声耦合) → 实验。
<b>性质越依赖电子结构(如 Tc)，越要往流程后端走。</b></p>
<ul class="honest">
<li>✅ 已验证：结构是 UMA 势能面稳定极小点，晶格接近实验(MAE 1.2%)；氟化物&lt;氧化物声子；
MgB₂ E₂g 声子；量子点尺寸-带隙趋势；MatterGen 生成结构 UMA 体模量≥目标；SrTiO₃ n 型 DOS。</li>
<li>⚠️ 限制：UMA 给能量/力/声子/弹性，不直接给带隙/Tc；Brus &lt;3nm 高估；PBE 低估带隙；
QE 的 λ 未完全收敛(需密网格双δ+非谐)；LiYF₄ 初始结构不合格。</li>
</ul>
<p class="foot">工具：FairChem UMA <code>uma-s-1p1</code> · MatterGen · GPAW/PBE · Quantum ESPRESSO(LDA DFPT) ·
pymatgen/ASE · Modal(A10G/CPU)。全部脚本与数据见 <code>materials_design/</code>。</p>
"""

HTML = f"""<!doctype html><html lang="zh"><head><meta charset="utf-8">
<title>计算材料设计与验证 — 项目报告</title>
<style>
body{{font-family:-apple-system,"Segoe UI","Noto Sans CJK SC","WenQuanYi Zen Hei",sans-serif;
 max-width:1180px;margin:0 auto;padding:32px 40px;color:#1b1f24;line-height:1.65;background:#fff}}
h1{{font-size:26px;border-bottom:3px solid #0072B2;padding-bottom:10px}}
h2{{font-size:20px;color:#0a4d73;margin-top:34px;border-left:4px solid #56B4E9;padding-left:10px}}
.sub{{color:#555;font-size:14px}} .toc{{background:#f5f8fb;border:1px solid #dce3ea;border-radius:8px;padding:12px 20px;font-size:14px}}
table{{border-collapse:collapse;width:100%;margin:12px 0;font-size:14px}}
th,td{{border:1px solid #dce3ea;padding:7px 10px;text-align:left}} th{{background:#eef3f8}}
code{{background:#f0f2f5;padding:1px 5px;border-radius:4px;font-size:13px}}
.method{{background:#E9F7EF;border:1px solid #009E73;border-radius:8px;padding:12px 16px}}
.honest li{{margin:6px 0}} .foot{{color:#666;font-size:13px;margin-top:24px;border-top:1px solid #eee;padding-top:12px}}
img{{display:block}}
</style></head><body>{SECTIONS}</body></html>"""

open(os.path.join(D,"report.html"),"w").write(HTML)
print("wrote report.html,", len(HTML), "bytes")
