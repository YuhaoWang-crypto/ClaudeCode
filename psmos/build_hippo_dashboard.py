"""
build_hippo_dashboard.py — regenerate the Hippo–YAP regeneration PSMOS dashboard
from computed data (Evo2 live, Compara live, AlphaGenome adapter-ready).

Mirrors the PSMOS regeneration dashboard: five role cards, a fidelity×utility
scatter, a G/D/N/R/E/X heatmap (with computed Evo2-constraint and Compara-
redundancy columns), and weight-profile tabs — but the numbers are now driven by
computed overlays where available, with provenance banners.

Output: psmos/hippo_dashboard_live.html
"""

import json
import os

from .scoring_psmos import compute_hippo

OUT = os.path.join(os.path.dirname(__file__), "hippo_dashboard_live.html")


def build():
    D = compute_hippo()
    html = TEMPLATE.replace("__DATA__", json.dumps(D, ensure_ascii=False))
    with open(OUT, "w") as fh:
        fh.write(html)
    print(f"wrote {OUT}")
    print(f"  Evo2 live: {D['evo2_live']}")
    print(f"  AlphaGenome R: {D['alphagenome_status']}")
    return OUT


TEMPLATE = r"""<!DOCTYPE html><html lang="zh"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>PSMOS · Hippo–YAP 再生 模式生物组合 (Evo2 live)</title>
<style>
:root{--s1:#fcfcfb;--s2:#f4f3f0;--bd:#e2e0da;--tp:#0b0b0b;--ts:#52514e;--tm:#8a887f;
 --c-mammal:#2a78d6;--c-aqua:#1baf7a;--c-invert:#eb6834;--c-ref:#8a887f;--evo:#7a3ff2;--ok:#1baf7a;
 --seq1:#cde2fb;--seq2:#86b6ef;--seq3:#2a78d6;--seq4:#184f95;}
@media(prefers-color-scheme:dark){:root{--s1:#1a1a19;--s2:#242422;--bd:#38382f;--tp:#fff;--ts:#c3c2b7;--tm:#8f8e82;
 --c-mammal:#3987e5;--c-aqua:#199e70;--c-invert:#d95926;--c-ref:#8f8e82;--evo:#a889f7;--ok:#199e70;
 --seq1:#104281;--seq2:#184f95;--seq3:#3987e5;--seq4:#86b6ef;}}
*{box-sizing:border-box}
body{margin:0;font-family:-apple-system,BlinkMacSystemFont,"Segoe UI","PingFang SC","Microsoft YaHei",sans-serif;
 background:var(--s1);color:var(--tp);line-height:1.5;padding:20px;max-width:1000px;margin:0 auto}
h1{font-size:20px;margin:0 0 4px}h2{font-size:15px;margin:24px 0 6px}
.sub{color:var(--ts);font-size:13px;margin-bottom:6px}
.banner{border:1px solid var(--bd);border-radius:6px;padding:9px 12px;font-size:12.5px;color:var(--ts);margin:8px 0;background:var(--s2)}
.banner.evo{border-left:3px solid var(--evo)}.banner.ag{border-left:3px solid var(--c-invert)}
.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(165px,1fr));gap:10px;margin-top:12px}
.card{border:1px solid var(--bd);border-radius:10px;padding:11px 12px;background:var(--s2)}
.card .role{font-size:11px;color:var(--tm)}.card .org{font-size:14px;font-weight:600;margin:3px 0}.card .why{font-size:11px;color:var(--ts)}
.legend{font-size:11px;color:var(--tm);margin:8px 0}
.dot{display:inline-block;width:9px;height:9px;border-radius:50%;margin-right:3px;vertical-align:middle}
.chip{display:inline-block;border-radius:3px;padding:1px 5px;margin-right:4px;font-size:10px}
svg{width:100%;height:auto;display:block;margin-top:4px}
.tick{fill:var(--tm);font-size:10px}.axlab{fill:var(--ts);font-size:11px}.pt-label{fill:var(--tp);font-size:10.5px;font-weight:500}
table{border-collapse:collapse;width:100%;font-size:11.5px;margin-top:6px}
th,td{border:1px solid var(--bd);padding:4px 6px;text-align:center}th{background:var(--s2)}
td.l,th.l{text-align:left}.hm{display:block;padding:4px 0;border-radius:3px}
.tabs{display:flex;gap:6px;flex-wrap:wrap;margin:8px 0}
.tab{padding:6px 10px;border:1px solid var(--bd);border-radius:8px;background:var(--s2);cursor:pointer;font-size:12px;color:var(--ts)}
.tab.on{background:var(--c-mammal);color:#fff;border-color:var(--c-mammal)}
.row{display:flex;align-items:center;gap:8px;margin:4px 0;font-size:12.5px}
.rank{width:16px;color:var(--tm);text-align:right}.name{width:220px;white-space:nowrap;overflow:hidden;text-overflow:ellipsis}
.barwrap{flex:1;background:var(--s2);border-radius:4px;height:18px;min-width:60px}
.bar{display:block;height:100%;border-radius:4px;background:var(--seq3);min-width:3px}
.val{width:34px;text-align:right;color:var(--ts);font-variant-numeric:tabular-nums}
.foot{font-size:11px;color:var(--tm);margin-top:22px;border-top:1px solid var(--bd);padding-top:12px}
b{color:var(--tp)}code{background:var(--s2);padding:1px 4px;border-radius:3px}
</style></head><body>
<h1>PSMOS · Hippo–YAP/TAZ–TEAD 再生 模式生物「组合」推荐</h1>
<div class="sub">研究问题 <b>regeneration</b> · 参考 <b>Homo sapiens</b>。不给唯一答案,给<b>五类分工模型</b> + 六层打分(G/D/N/R/E/X)。</div>
<div class="banner evo" id="evobanner"></div>
<div class="banner ag" id="agbanner"></div>

<h2>① 五类分工模型(核心输出)</h2>
<div class="cards" id="cards"></div>

<h2>② 保真度 × 可实施性(为什么没有"唯一最佳")</h2>
<div class="legend">
 <span class="dot" style="background:var(--c-mammal)"></span>哺乳
 <span class="dot" style="background:var(--c-aqua)"></span>鱼/两栖
 <span class="dot" style="background:var(--c-invert)"></span>无脊椎
 <span class="dot" style="background:var(--c-ref)"></span>参考(人) · 右上=理想;左上=最像最难做;右下=最好做最不像。</div>
<svg id="scatter" viewBox="0 0 700 440"></svg>

<h2>③ 六层打分热图 (G/D/N/R/E/X) + 实测覆盖</h2>
<div class="legend" id="provlegend"></div>
<div id="heat"></div>

<h2>④ 权重随研究问题变化 → 排序变化 (PSMOS)</h2>
<div class="tabs" id="tabs"></div>
<div id="rank"></div>

<div class="foot" id="foot"></div>
<script>
const D=__DATA__;
const CATCOL={mammal:'var(--c-mammal)',aqua:'var(--c-aqua)',invert:'var(--c-invert)',ref:'var(--c-ref)'};

document.getElementById('evobanner').innerHTML =
 `<b>Evo2 状态:已接入(live)。</b> evo2_7b 在 Modal H100 上为 <b>${D.evo2_live.length}</b> 个物种(${D.evo2_live.join(', ')})的 YAP/TAZ + TEAD 真实 CDS 计算了对数似然;`+
 `G 层("核心基因正交 + 序列保守")对这些物种已是<b>curated ⊕ Evo2 实测</b>混合。旁系同源冗余由 Ensembl Compara <b>实测</b>。`;
document.getElementById('agbanner').innerHTML =
 `<b>AlphaGenome(R 层)状态:${D.alphagenome_status}。</b> R 层(启动子/增强子/TF motif/剪接)只能覆盖<b>人与小鼠</b>(AlphaGenome 模型范围);`+
 `适配器已就绪、人鼠基因组区间已解析,设置 <code>ALPHAGENOME_API_KEY</code> 即用人鼠调控轨道实测值替换该层。`;

// role cards
const roles={mechanistic:['机制解析 Mechanistic','拆解通路因果、上位性'],
 imaging:['动态成像 Imaging','活体观察激活时程/命运'],
 genetic:['遗传筛选 Genetic','高通量扰动/修饰子筛选'],
 translational:['人体转化 Translational','向人外推、临床前验证'],
 comparative:['因果/阴性对照 Comparative','天然重连→验证充分必要性']};
const rowByKey={};D.rows.forEach(r=>rowByKey[r.key]=r);
let ch='';
for(const k of ['mechanistic','imaging','genetic','translational','comparative']){
 const w=D.winners[k];const r=rowByKey[w.key];
 ch+=`<div class="card"><div class="role">${roles[k][0]}</div>
  <div class="org" style="color:${CATCOL[r.clade]}">${w.name.split(' (')[0]}</div>
  <div class="why">${roles[k][1]}<br><span style="color:var(--tm)">拟合度 ${w.score.toFixed(0)}</span></div></div>`;}
document.getElementById('cards').innerHTML=ch;

// scatter
(function(){const W=700,H=440,ml=52,mr=20,mt=16,mb=44;
 const X=v=>ml+(v)*(W-ml-mr),Y=v=>H-mb-(v)*(H-mt-mb);
 let s='';
 for(let g=0;g<=1;g+=0.25){s+=`<line x1="${X(g)}" y1="${mt}" x2="${X(g)}" y2="${H-mb}" stroke="var(--bd)"/>`;
  s+=`<line x1="${ml}" y1="${Y(g)}" x2="${W-mr}" y2="${Y(g)}" stroke="var(--bd)"/>`;
  s+=`<text class="tick" x="${X(g)}" y="${H-mb+14}" text-anchor="middle">${(g*100).toFixed(0)}</text>`;
  s+=`<text class="tick" x="${ml-6}" y="${Y(g)+3}" text-anchor="end">${(g*100).toFixed(0)}</text>`;}
 s+=`<text class="axlab" x="${(ml+W-mr)/2}" y="${H-6}" text-anchor="middle">实验可实施性 Experimental Utility →</text>`;
 s+=`<text class="axlab" transform="translate(14,${(mt+H-mb)/2}) rotate(-90)" text-anchor="middle">生物学保真度 Biological Fidelity →</text>`;
 D.rows.forEach(r=>{const col=CATCOL[r.clade];const x=X(r.utility),y=Y(r.fidelity);
  s+=`<circle cx="${x}" cy="${y}" r="6.5" fill="${col}" stroke="var(--s1)" stroke-width="2"/>`;
  const lab=r.name.split(' (')[1]? r.name.split('(')[1].replace(')',''):r.name;
  const anchor=x>W-140?'end':'start';const dx=x>W-140?-10:10;
  s+=`<text class="pt-label" x="${x+dx}" y="${y+3}" text-anchor="${anchor}">${lab}</text>`;});
 document.getElementById('scatter').innerHTML=s;})();

// heatmap
document.getElementById('provlegend').innerHTML=
 `G基因·D结构域·N网络·R调控·E表达·X实验 + 两列实测覆盖。来源 `+
 `<span class="chip" style="background:var(--evo);color:#fff">computed:evo2</span>`+
 `<span class="chip" style="background:var(--c-mammal);color:#fff">computed:compara</span>`+
 `<span class="chip" style="background:var(--s2);color:var(--ts);border:1px solid var(--bd)">curated</span>。`;
function ramp(v,base){const st=base||[[205,226,251],[134,182,239],[42,120,214],[24,79,149]];
 const t=Math.max(0,Math.min(1,v))*(st.length-1);const i=Math.floor(t),f=t-i;
 const a=st[i],b=st[Math.min(i+1,st.length-1)];return`rgb(${a.map((x,k)=>Math.round(x+(b[k]-x)*f)).join(',')})`;}
const EVORAMP=[[237,231,251],[199,180,244],[150,110,235],[104,52,205]];
const L=[['G','G 基因'],['D','D 结构域'],['N','N 网络'],['R','R 调控'],['E','E 表达'],['X','X 实验']];
let h='<table><tr><th class="l">物种</th>'+L.map(l=>`<th>${l[1]}</th>`).join('')+
 '<th style="background:var(--evo);color:#fff">Evo2约束</th><th style="background:var(--c-mammal);color:#fff">冗余(Compara)</th></tr>';
D.rows.forEach(r=>{
 h+=`<tr><td class="l" style="border-left:3px solid ${CATCOL[r.clade]}">${r.name}</td>`+
  L.map(l=>{const v=r.layers[l[0]];const ink=v>.55?'#fff':'var(--tp)';
   const star=(r.prov[l[0]]||'').startsWith('computed')?'*':'';
   return`<td title="${r.prov[l[0]]}"><span class="hm" style="background:${ramp(v)};color:${ink}">${(v*100).toFixed(0)}${star}</span></td>`;}).join('');
 if(r.evo2_constraint!=null){const c=ramp(r.evo2_constraint,EVORAMP);const ink=r.evo2_constraint>.5?'#fff':'var(--tp)';
  h+=`<td title="meanLL ${r.evo2_meanLL}"><span class="hm" style="background:${c};color:${ink}">${(r.evo2_constraint*100).toFixed(0)}</span></td>`;}
 else h+=`<td><span class="hm" style="background:var(--s2);color:var(--tm)">—</span></td>`;
 if(r.compara_simplicity!=null){const c=ramp(r.compara_simplicity);const ink=r.compara_simplicity>.55?'#fff':'var(--tp)';
  h+=`<td><span class="hm" style="background:${c};color:${ink}">${(r.compara_simplicity*100).toFixed(0)}</span></td>`;}
 else h+=`<td><span class="hm" style="background:var(--s2);color:var(--tm)">—</span></td>`;
 h+='</tr>';});
h+='</table><div class="legend">* = 该层含实测覆盖(Evo2 混入 G / Compara 提供冗余)。</div>';
document.getElementById('heat').innerHTML=h;

// weight-profile tabs
const profs=Object.keys(D.profiles),tabs=document.getElementById('tabs'),rank=document.getElementById('rank');
function draw(p){const rows=D.profiles[p];const mx=Math.max(...rows.map(r=>r.psmos))||1;
 rank.innerHTML=rows.map((r,i)=>`<div class="row"><span class="rank">${i+1}</span>
  <span class="name" style="color:${CATCOL[r.clade]}">${r.name.split(' (')[0]}</span>
  <span class="barwrap"><span class="bar" style="width:${Math.max(3,r.psmos/mx*100)}%;background:${CATCOL[r.clade]}"></span></span>
  <span class="val">${r.psmos.toFixed(0)}</span></div>`).join('');}
profs.forEach((p,i)=>{const b=document.createElement('div');b.className='tab'+(i==0?' on':'');b.textContent=p;
 b.onclick=()=>{document.querySelectorAll('.tab').forEach(t=>t.classList.remove('on'));b.classList.add('on');draw(p);};tabs.appendChild(b);});
draw(profs[0]);

document.getElementById('foot').innerHTML=
 `<b>方法学:PSMOS 六层 v1(Evo2 live)。</b> `+
 `<b>实测层:</b>① G 的序列保守 ← Evo2-7B log-likelihood(Modal H100,真实 CDS,${D.evo2_live.length} 物种);② 冗余 ← Ensembl Compara 旁系同源计数;`+
 `③ R(人/鼠)← AlphaGenome(适配器就绪,待 API key)。<b>策划层:</b>N 网络拓扑、E 时空表达、X 工具/通量为 Hippo 比较生物学策划先验。`+
 `<br><b>核心结论:</b>同一通路在同一问题下,机制/成像/遗传/转化/对照应选<b>不同</b>物种——果蝇 G 层因 Evo2 实测偏低仍是机制/遗传首选(单拷贝、遗传学无敌);`+
 `斑马鱼夺成像/再生;涡虫(天然重连)是最佳阴性/对照。"唯一最相似物种"是伪命题。`;
</script></body></html>"""


if __name__ == "__main__":
    build()
