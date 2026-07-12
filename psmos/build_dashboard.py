"""
build_dashboard.py — regenerate the Notch model-organism dashboard from COMPUTED
data, with Evo2 wired in live and per-axis provenance shown.

Differences from the v0 pilot dashboard (which was fully curated):
  * Evo2 status banner is LIVE when psmos/data/evo2/<pw>_scores.json exists —
    it names the species whose sequence-constraint layer is now model-measured.
  * A computed "Evo2 constraint" column joins the heatmap; a provenance legend
    marks computed:uniprot / computed:ensembl / computed:evo2 / curated.
  * The hard gate is driven by the real UniProt/Ensembl ortholog search.
  * Footer reports Pearson r between Evo2 constraint and the curated fidelity
    proxy — corroboration, not laundering.

Output: psmos/notch_dashboard_live.html
"""

import json
import os

from .scoring import compute, evo2_vs_curated_correlation, GOALS

OUT = os.path.join(os.path.dirname(__file__), "notch_dashboard_live.html")

AX = [("completeness", "完整"), ("fidelity", "保真(代理)"), ("simplicity", "低冗余"),
      ("arch_match", "架构"), ("tract", "可操作"), ("thru", "通量")]


def build(pathway_key="Notch"):
    pw, scores = compute(pathway_key)
    r, n_live = evo2_vs_curated_correlation(scores)
    live_species = [s.key for s in scores if s.evo2]

    # rankings per goal
    rankings = {}
    for gname in GOALS:
        rows = sorted(scores, key=lambda s: -s.goal_scores[gname])
        rankings[gname] = [
            dict(rank=i + 1, key=s.key, name=s.name, score=s.goal_scores[gname],
                 gate=s.gate_ok,
                 evo2=(s.evo2.get("constraint") if s.evo2 else None),
                 meanLL=(s.evo2.get("mean_LL") if s.evo2 else None))
            for i, s in enumerate(rows)
        ]

    heat = []
    for s in scores:
        heat.append(dict(
            key=s.key, name=s.name, gate=s.gate_ok,
            gate_prov=s.gate_provenance,
            axes=s.axes,
            evo2_constraint=(s.evo2.get("constraint") if s.evo2 else None),
            evo2_meanLL=(s.evo2.get("mean_LL") if s.evo2 else None),
            evo2_genes=(s.evo2.get("per_gene") if s.evo2 else None),
        ))

    D = dict(
        pathway=pw.label, reference="Homo sapiens",
        human_mean_paralogs=round(pw.human_mean_paralogs, 2),
        live_species=live_species,
        evo2_live=bool(live_species),
        corr=(None if r is None else round(r, 2)), corr_n=n_live,
        goals={g: GOALS[g]["note"] for g in GOALS},
        rankings=rankings, heat=heat,
    )

    html = TEMPLATE.replace("__DATA__", json.dumps(D, ensure_ascii=False))
    with open(OUT, "w") as fh:
        fh.write(html)
    print(f"wrote {OUT}")
    print(f"  Evo2 live: {D['evo2_live']}  species={live_species}")
    if r is not None:
        print(f"  Evo2 constraint vs curated fidelity: r={r:+.2f} (n={n_live})")
    return OUT


TEMPLATE = r"""<!DOCTYPE html><html lang="zh"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>通路感知模式生物推荐 · Notch (Evo2 live)</title>
<style>
:root{--s1:#fcfcfb;--s2:#f4f3f0;--bd:#e2e0da;--tp:#0b0b0b;--ts:#52514e;--tm:#8a887f;
 --series:#2a78d6;--seq1:#cde2fb;--seq2:#86b6ef;--seq3:#2a78d6;--seq4:#184f95;--reject:#b8b6ad;
 --evo:#7a3ff2;--ok:#1baf7a;}
@media(prefers-color-scheme:dark){:root{--s1:#1a1a19;--s2:#242422;--bd:#38382f;--tp:#fff;--ts:#c3c2b7;--tm:#8f8e82;
 --series:#3987e5;--seq1:#104281;--seq2:#184f95;--seq3:#3987e5;--seq4:#86b6ef;--reject:#4a4a42;--evo:#a889f7;--ok:#199e70;}}
*{box-sizing:border-box}
body{margin:0;font-family:-apple-system,BlinkMacSystemFont,"Segoe UI","PingFang SC","Microsoft YaHei",sans-serif;
 background:var(--s1);color:var(--tp);line-height:1.5;padding:20px;max-width:1000px;margin:0 auto}
h1{font-size:20px;margin:0 0 4px}h2{font-size:15px;margin:24px 0 8px}
.sub{color:var(--ts);font-size:13px;margin-bottom:12px}
.banner{border:1px solid var(--bd);border-left:3px solid var(--ok);border-radius:6px;padding:9px 12px;
 font-size:12.5px;color:var(--ts);margin:10px 0;background:var(--s2)}
.banner.evo{border-left-color:var(--evo)}
.tabs{display:flex;gap:6px;flex-wrap:wrap;margin:10px 0}
.tab{padding:8px 12px;border:1px solid var(--bd);border-radius:8px;background:var(--s2);cursor:pointer;font-size:13px;color:var(--ts)}
.tab.on{background:var(--series);color:#fff;border-color:var(--series)}
.panel{display:none}.panel.on{display:block}
.goalnote{font-size:12px;color:var(--tm);margin:0 0 10px}
.row{display:flex;align-items:center;gap:8px;margin:5px 0;font-size:13px}
.rank{width:20px;color:var(--tm);text-align:right;flex:none}
.name{width:220px;flex:none;white-space:nowrap;overflow:hidden;text-overflow:ellipsis}
.barwrap{flex:1;background:var(--s2);border-radius:4px;height:20px;position:relative;min-width:60px}
.bar{display:block;height:100%;border-radius:4px;background:var(--seq3);min-width:3px}
.val{width:38px;flex:none;text-align:right;font-variant-numeric:tabular-nums;color:var(--ts)}
.rej .bar{background:var(--reject)}.rej .name{color:var(--tm)}
.tag{font-size:10px;color:#e34948;margin-left:4px}
.elive{font-size:9.5px;color:#fff;background:var(--evo);border-radius:3px;padding:1px 4px;margin-left:5px;vertical-align:middle}
.tt{position:fixed;pointer-events:none;background:var(--tp);color:var(--s1);padding:8px 10px;border-radius:6px;
 font-size:12px;max-width:280px;opacity:0;transition:opacity .1s;z-index:9}
.legend{font-size:11px;color:var(--tm);margin:6px 0 10px}
.chip{display:inline-block;border-radius:3px;padding:1px 5px;margin-right:4px;font-size:10px}
table{border-collapse:collapse;width:100%;font-size:11.5px;margin-top:6px}
th,td{border:1px solid var(--bd);padding:4px 6px;text-align:center}th{background:var(--s2)}
td.l,th.l{text-align:left}.hm{display:inline-block;width:100%;padding:4px 0;border-radius:3px}
.foot{font-size:11px;color:var(--tm);margin-top:22px;border-top:1px solid var(--bd);padding-top:12px}
code{background:var(--s2);padding:1px 4px;border-radius:3px}
</style></head><body>
<h1>通路感知模式生物推荐 · Notch 信号通路</h1>
<div class="sub">参考物种 = 人。同一物种面板,在三种研究目标下重新排序。分数 0–100。
硬门控:缺受体或 CSL 转录因子 → 判定通路缺失、直接淘汰。</div>
<div class="banner evo" id="evobanner"></div>
<div class="banner" id="gatebanner"></div>

<div class="tabs" id="tabs"></div>
<div id="panels"></div>

<h2>各物种打分维度 + Evo2 实测约束(热图)</h2>
<div class="legend" id="provlegend"></div>
<div id="heat"></div>

<div class="foot" id="foot"></div>
<div class="tt" id="tt"></div>
<script>
const D=__DATA__;
const goals=Object.keys(D.rankings);
const tt=document.getElementById('tt');
function show(e,html){tt.innerHTML=html;tt.style.opacity=1;let x=e.clientX+14,y=e.clientY+14;
 if(x+290>innerWidth)x=e.clientX-290;tt.style.left=x+'px';tt.style.top=y+'px';}
function hide(){tt.style.opacity=0;}

document.getElementById('evobanner').innerHTML = D.evo2_live
 ? `<b>Evo2 状态:已接入(live)。</b> 通过 Modal H100 对 <b>${D.live_species.length}</b> 个物种的 Notch 核心基因编码序列(真实 CDS,来自 Ensembl)计算了 Evo2 对数似然(约束/自然度)。这条“序列保真”维度已从<b>分歧时间代理</b>升级为<b>模型实测</b>。`
 : `<b>Evo2 状态:适配器就绪,尚未运行。</b> 运行 <code>python3 -m psmos.run_evo2_scoring Notch</code> 即可为下列物种填入实测约束分。`;

document.getElementById('gatebanner').innerHTML =
 `<b>硬门控 = 实测。</b> 受体 + CSL 两个门控家族的存在与否,由 UniProt/Ensembl 直系同源检索实证得到:` +
 `酵母、拟南芥两个物种对两个门控基因均<b>零命中</b> → Notch 通路确实缺失 → 淘汰(不是断言,是检索结果)。`;

document.getElementById('provlegend').innerHTML =
 `列:6 维打分 + <b>Evo2 约束</b>。来源标注 ` +
 `<span class="chip" style="background:var(--evo);color:#fff">computed:evo2</span>` +
 `<span class="chip" style="background:var(--series);color:#fff">computed:uniprot/ensembl</span>` +
 `<span class="chip" style="background:var(--s2);color:var(--ts);border:1px solid var(--bd)">curated 策划</span>。灰行=通路缺失。`;

// tabs + panels
const tabs=document.getElementById('tabs'),panels=document.getElementById('panels');
goals.forEach((g,i)=>{
 const b=document.createElement('div');b.className='tab'+(i==0?' on':'');b.textContent=g;
 b.onclick=()=>{document.querySelectorAll('.tab').forEach(t=>t.classList.remove('on'));
  document.querySelectorAll('.panel').forEach(t=>t.classList.remove('on'));
  b.classList.add('on');document.getElementById('p'+i).classList.add('on');};
 tabs.appendChild(b);
 const p=document.createElement('div');p.className='panel'+(i==0?' on':'');p.id='p'+i;
 const nd=document.createElement('div');nd.className='goalnote';nd.textContent='目标:'+D.goals[g];p.appendChild(nd);
 const rows=D.rankings[g];const max=Math.max(...rows.map(r=>r.score))||1;
 rows.forEach(r=>{
  const row=document.createElement('div');row.className='row'+(r.gate?'':' rej');
  const evtag=(r.evo2!=null)?`<span class="elive">Evo2 ${(r.evo2*100).toFixed(0)}</span>`:'';
  const detail=`<b>${r.name}</b><br>`+
   (r.evo2!=null?`Evo2 约束(实测): ${(r.evo2*100).toFixed(0)}/100 · meanLL ${r.meanLL}<br>`:`Evo2: 无 live 数据<br>`)+
   `门控: ${r.gate?'通过':'淘汰(通路缺失)'}`;
  row.innerHTML=`<span class="rank">${r.rank}</span>`+
   `<span class="name" title="${r.name}">${r.name}</span>`+
   `<span class="barwrap"><span class="bar" style="width:${Math.max(2,r.score/max*100)}%"></span></span>`+
   `<span class="val">${r.score.toFixed(0)}</span>`+(r.gate?'':'<span class="tag">淘汰</span>')+evtag;
  row.onmousemove=e=>show(e,detail);row.onmouseleave=hide;row.ontouchstart=e=>show(e.touches[0],detail);
  p.appendChild(row);
 });
 panels.appendChild(p);
});

// heatmap
const AX=[["completeness","完整"],["fidelity","保真代理"],["simplicity","低冗余"],
 ["arch_match","架构"],["tract","可操作"],["thru","通量"]];
function ramp(v,base){const st=base||[[205,226,251],[134,182,239],[42,120,214],[24,79,149]];
 const t=Math.max(0,Math.min(1,v))*(st.length-1);const i=Math.floor(t),f=t-i;
 const a=st[i],b=st[Math.min(i+1,st.length-1)];return`rgb(${a.map((x,k)=>Math.round(x+(b[k]-x)*f)).join(',')})`;}
const EVORAMP=[[237,231,251],[199,180,244],[150,110,235],[104,52,205]];
let h='<table><tr><th class="l">物种</th>'+AX.map(a=>`<th>${a[1]}</th>`).join('')+
 '<th style="background:var(--evo);color:#fff">Evo2约束</th></tr>';
D.heat.forEach(r=>{
 const rej=!r.gate;
 h+=`<tr><td class="l">${r.name}</td>`+AX.map(x=>{
  const v=r.axes[x[0]];const col=rej?'var(--reject)':ramp(v);
  const ink=(!rej&&v>0.55)?'#fff':'var(--tp)';
  return`<td><span class="hm" style="background:${col};color:${ink}">${(v*100).toFixed(0)}</span></td>`;
 }).join('');
 if(r.evo2_constraint!=null){
  const col=ramp(r.evo2_constraint,EVORAMP);const ink=r.evo2_constraint>0.5?'#fff':'var(--tp)';
  const genes=r.evo2_genes?Object.entries(r.evo2_genes).map(([k,v])=>`${k}:${v.toFixed(2)}`).join(' · '):'';
  h+=`<td title="meanLL ${r.evo2_meanLL} | ${genes}"><span class="hm" style="background:${col};color:${ink}">`+
     `${(r.evo2_constraint*100).toFixed(0)}</span></td>`;
 }else{
  h+=`<td><span class="hm" style="background:var(--s2);color:var(--tm)">—</span></td>`;
 }
 h+='</tr>';
});
h+='</table>';document.getElementById('heat').innerHTML=h;

let corrTxt = (D.corr!=null)
 ? `Evo2 实测约束 vs 策划保真代理:Pearson r = <b>${D.corr>=0?'+':''}${D.corr}</b>(n=${D.corr_n} 个 live 物种)——二者${Math.abs(D.corr)>0.5?'同向印证':'相关性弱,说明 Evo2 捕捉到分歧时间代理之外的信号'}。`
 : `Evo2 live 物种不足以估计相关性。`;
document.getElementById('foot').innerHTML=
 `<b>方法学:PSMOS v1(Evo2 live)。</b> `+
 `<b>实测层:</b>① 硬门控 + 门控家族存在 ← UniProt/Ensembl 直系同源检索;② 序列约束 ← Evo2-7B log-likelihood(Modal H100,真实 CDS)。`+
 `<b>策划层(待接入):</b>架构相似度、低冗余(旁系同源计数)、可操作性/通量为 Notch 比较基因组学策划先验;`+
 `跨物种调控语法(R 层)需 AlphaGenome(人/鼠)。`+
 `<br>${corrTxt} 人均旁系同源/家族 = ${D.human_mean_paralogs}。`+
 `<br><b>诚实边界:</b>Evo2 log-likelihood 衡量的是“序列在生命之树上的约束/自然度”,<b>不是</b>“与人的等价性”;因此它作为独立实测维度呈现,而非直接覆盖策划的保真轴。`;
</script></body></html>"""


if __name__ == "__main__":
    build("Notch")
