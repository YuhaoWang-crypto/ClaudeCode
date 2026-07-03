#!/usr/bin/env python3
"""Build the self-contained full-humanized-IgG anti-Tirzepatide report."""
import json

P = json.load(open("report_payload.json"))
ORDER = ["ab2-mat1","ab2-mat2","ab2-wt","ab1-wt","ab2-mat3"]  # display order (by tier then name)
rankById = {r["ab"]: r for r in P["rank"]}

# ---- LEADS array for the ranking table (JS) ----
leads = []
for ab in ORDER:
    r = rankById[ab]
    leads.append([ab, P["mut"][ab], r["bind"], r["iptm_mean"], r["iptm_best"], r["iptm_sd"],
                  r["sconf"], r["core_retained"], r["n_epi"], r["scfv_iptm"], r["scfv_bind"],
                  P["tier"][ab], P["note"][ab]])

DATA = {
    "LEADS": leads,
    "SEQ": P["seqs"], "CDRH3": P["cdrh3"], "MUT": P["mut"], "EPI": P["epi"],
    "CIF": P["cif"], "TIER": P["tier"], "CONSTS": P["consts"], "TIRZ": P["tirz"],
}
DATA_JSON = json.dumps(DATA, separators=(",", ":"))

HTML = r"""<!DOCTYPE html><html lang="en"><head>
<meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Anti-Tirzepatide Humanized IgG Report</title>
<script>(function(){var s=["https://cdn.jsdelivr.net/npm/3dmol@2.1.0/build/3Dmol-min.js",
"https://unpkg.com/3dmol@2.1.0/build/3Dmol-min.js","https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.1.0/3Dmol-min.js"];
var i=0;function n(){if(i>=s.length){window._3dfail=1;return;}var e=document.createElement('script');
e.src=s[i++];e.onload=function(){window._3dok=1};e.onerror=n;document.head.appendChild(e);}n();})();</script>
<style>
:root{--bg:#0f1420;--panel:#171e2e;--ink:#e8edf6;--mut:#93a1bd;--acc:#4C9AFF;--acc2:#E4572E;--ok:#54A24B;--line:#26304a}
*{box-sizing:border-box}body{margin:0;font:15px/1.55 -apple-system,Segoe UI,Roboto,Helvetica,Arial,"PingFang SC","Microsoft YaHei",sans-serif;background:var(--bg);color:var(--ink)}
header{padding:22px 26px;border-bottom:1px solid var(--line);background:linear-gradient(180deg,#141b2b,#0f1420)}
h1{margin:0 0 4px;font-size:21px}.sub{color:var(--mut);font-size:13px;max-width:1000px}
nav{display:flex;gap:6px;padding:10px 20px;position:sticky;top:0;background:#0f1420ee;backdrop-filter:blur(6px);border-bottom:1px solid var(--line);z-index:5;flex-wrap:wrap}
nav button{background:transparent;color:var(--mut);border:1px solid transparent;padding:7px 13px;border-radius:8px;cursor:pointer;font-size:14px}
nav button.active{color:var(--ink);background:var(--panel);border-color:var(--line)}
main{padding:20px 26px;max-width:1180px;margin:0 auto}section{display:none}section.on{display:block}
.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(160px,1fr));gap:12px;margin:14px 0}
.card{background:var(--panel);border:1px solid var(--line);border-radius:12px;padding:14px}
.card .n{font-size:24px;font-weight:700}.card .l{color:var(--mut);font-size:12px;margin-top:2px}
table{border-collapse:collapse;width:100%;font-size:13.5px;margin-top:8px}
th,td{padding:8px 9px;border-bottom:1px solid var(--line);text-align:left;vertical-align:top}
th{color:var(--mut);cursor:pointer;user-select:none}tr:hover td{background:#141c2e}
.seq{font-family:ui-monospace,Menlo,monospace;font-size:12px;letter-spacing:.4px;color:#bfe0ff;word-break:break-all}
.mono{font-family:ui-monospace,Menlo,monospace}
.copy{cursor:pointer;border:1px solid var(--line);background:#1d2740;color:var(--mut);border-radius:6px;padding:2px 7px;font-size:11px;margin-left:6px}
.copy:hover{color:var(--ink)}
.viewwrap{display:grid;grid-template-columns:330px 1fr;gap:16px}@media(max-width:860px){.viewwrap{grid-template-columns:1fr}}
#viewer{position:relative;width:100%;height:560px;background:#0a0e17;border:1px solid var(--line);border-radius:12px;overflow:hidden}
select,.ctl{width:100%;background:#1d2740;color:var(--ink);border:1px solid var(--line);border-radius:8px;padding:9px;font-size:14px}
.legend span{display:inline-flex;align-items:center;gap:6px;margin-right:12px;font-size:12px}.dot{width:12px;height:12px;border-radius:3px;display:inline-block}
.kv{display:flex;justify-content:space-between;border-bottom:1px dashed var(--line);padding:5px 0;font-size:13px}
.note{color:var(--mut);font-size:12.5px}a{color:var(--acc)}.pill{font-size:11px;color:var(--mut);border:1px solid var(--line);border-radius:20px;padding:1px 8px}
.warn{background:#2a1e14;border:1px solid #5b3b1e;border-radius:10px;padding:12px 14px;color:#f0c7a0;font-size:13px}
.good{background:#12251a;border:1px solid #1f5b37;border-radius:10px;padding:12px 14px;color:#a7e8c0;font-size:13px}
h2{font-size:17px;border-left:3px solid var(--acc);padding-left:9px;margin-top:26px}code{background:#1d2740;padding:1px 5px;border-radius:4px}
.badge{display:inline-block;padding:2px 8px;border-radius:20px;font-size:11px;font-weight:600}
.b-A{background:#1c3a2a;color:#7ee2a8}.b-B{background:#28324e;color:#8fb6ff}.b-C{background:#3a2f1c;color:#e2c07e}
.bar{height:7px;border-radius:4px;background:#1d2740;overflow:hidden;margin-top:3px}.bar>i{display:block;height:100%;background:var(--acc)}
</style></head><body>
<header><h1>Anti-Tirzepatide Antibodies — Full Humanized IgG1: Re-dock, Affinity &amp; Ranking</h1>
<div class="sub">All previously designed binders, re-assembled as <b>native VH/VL paired into complete humanized IgG1</b> (VH–CH1–hinge–CH2–CH3 heavy · VL–Cλ light), then re-docked against <b>modified Tirzepatide</b> (Aib2/Aib13, K20 lipid). Boltz-2.1 co-fold of the antigen-binding <b>Fab</b> + peptide, 5 samples each. This supersedes the scFv / scFv-Fc round.</div></header>
<nav>
 <button data-t="rec" class="active">Recommendation</button>
 <button data-t="rank">Ranking &amp; Scores</button>
 <button data-t="fmt">scFv → IgG (before/after)</button>
 <button data-t="dock">Docking (3D)</button>
 <button data-t="seq">IgG Sequences</button>
 <button data-t="meth">Methods</button>
</nav><main>

<section id="rec" class="on">
 <div class="cards">
  <div class="card"><div class="n">ab2-mat1</div><div class="l">primary lead · humanized IgG1 · H3:A8Y</div></div>
  <div class="card"><div class="n">0.453</div><div class="l">binding-head score (highest of 5)</div></div>
  <div class="card"><div class="n">5 / 5</div><div class="l">constructs retain F22·V23·L26·I27 epitope</div></div>
  <div class="card"><div class="n">0.79–0.84</div><div class="l">Fab interface ipTM range (all confident)</div></div>
 </div>
 <div class="good"><b>Recommendation.</b> Advance <b>ab2-mat1 (H3:A8Y)</b> as the primary lead <b>as a full human IgG1</b> (VH–CH1–hinge–CH2–CH3 / VL–Cλ, bivalent → avidity). It is the top binding-head scorer in <b>both</b> the earlier scFv round (0.731) and this full-IgG round (0.453), keeps the clean 5/5 canonical epitope, and folds cleanly with the constant domains attached. Carry <b>ab2-mat2 (A9Y)</b> and <b>ab1-wt</b> (orthogonal lineage) as backups. <b>Do not</b> advance the triple mutant ab2-mat3. Confirm empirically by SPR/BLI + DSF/SEC.</div>
 <h2>What changed vs the scFv round — 4 findings</h2>
 <ul>
  <li><b>Reformatting is safe.</b> Pairing VH/VL natively and bolting on the human CH1/Cλ + IgG1 Fc did <b>not</b> move the paratope: all five constructs still dock onto the C-terminal amphipathic face (<b>F22·V23·Q24·L26·I27</b>) — the same epitope found for the scFvs.</li>
  <li><b>Interface confidence softens a little.</b> Adding the constant domains lowers atomistic interface ipTM from the scFv range (0.86–0.90) to the Fab range (0.79–0.84). All stay above 0.79 (confident); this is expected geometry dilution, not loss of binding.</li>
  <li><b>The lead is robust.</b> ab2-mat1 (A8Y) is the #1 binding-head scorer in both formats — the ranking survives the format change, which is the decision that matters.</li>
  <li><b>The A9Y edge does not survive.</b> In scFv, ab2-mat2 (A9Y) slightly beat WT on the binder head; in the full-IgG context that gap collapses (A9Y 0.343 ≈ WT 0.372). A8Y, not A9Y, is the durable improvement.</li>
 </ul>
 <p class="note">Target epitope = Tirzepatide C-terminal amphipathic face <b>F22·V23·L26·I27</b> (+ flanking Q19/K20/Q24). The K20 fatty-diacid sits on the opposite face and is distal to every paratope (kept as bare Lys here for a controlled cross-construct comparison).</p>
</section>

<section id="rank">
 <h2>Ranking &amp; scores <span class="pill">click a header to sort</span></h2>
 <table id="tblLead"><thead><tr>
  <th>Tier</th><th>IgG</th><th>CDR-H3 mut</th><th data-num="1">binding-head</th><th data-num="1">ipTM mean</th>
  <th data-num="1">ipTM best</th><th data-num="1">±sd</th><th data-num="1">fold conf.</th><th data-num="1">epitope core</th><th data-num="1"># contacts</th><th>note</th>
 </tr></thead><tbody></tbody></table>
 <p class="note">binding-head = Boltz protein–protein binding_confidence (affinity proxy, 0–1, higher=stronger); ipTM = atomistic Fab↔peptide protein_ipTM averaged over 5 samples (interface geometry confidence); fold conf. = structure_confidence; epitope core = how many of the 5 canonical epitope residues (F22 V23 Q24 L26 I27) the Fab contacts; # contacts = total peptide residues within 4.5 Å of the antibody.</p>
 <div class="warn"><b>How to read this.</b> Two "affinity" signals can diverge: <b>binding-head</b> tracks predicted affinity magnitude, <b>ipTM</b> tracks how confident the model is in the interface <i>shape</i>. ab1-wt has the best ipTM/fold but the weakest binding-head — a rigid, well-defined, but modest-affinity grip. We rank primarily by the binding head. No single co-fold score predicts developability — pair with wet-lab SPR/BLI + liability screening.</div>
</section>

<section id="fmt">
 <h2>scFv → full humanized IgG — does adding the constant region / Fc change binding?</h2>
 <div class="good"><b>Short answer: it preserves the binder, softens the interface score, and keeps the lead on top.</b> The Fc/CH1/Cλ are separate modules that never touch the paratope; per-arm binding is governed by the Fv, and the bivalent IgG adds avidity. The only glycosylation introduced is the canonical, distal <b>Fc N297</b> glycan — not a paratope liability.</div>
 <table id="tblFmt"><thead><tr><th>IgG</th><th data-num="1">scFv ipTM (prev)</th><th data-num="1">Fab ipTM (now, μ)</th><th data-num="1">Δ ipTM</th><th data-num="1">scFv bind (prev)</th><th data-num="1">Fab bind (now)</th><th>interpretation</th></tr></thead><tbody></tbody></table>
 <p class="note">Δ ipTM is the change from the isolated scFv interface to the Fab interface (constant domains present). A small negative Δ is expected — the interface is a smaller fraction of a larger complex — and does not indicate weaker binding.</p>
</section>

<section id="dock">
 <h2>Docking viewer <span class="pill">drag=rotate · scroll=zoom</span></h2>
 <div class="viewwrap"><div>
   <select id="pick"></select>
   <div class="legend" style="margin-top:10px">
    <span><i class="dot" style="background:#4C9AFF"></i>heavy (VH-CH1)</span>
    <span><i class="dot" style="background:#8FB6FF"></i>light (VL-Cλ)</span>
    <span><i class="dot" style="background:#9aa7bd"></i>Tirzepatide</span>
    <span><i class="dot" style="background:#E4572E"></i>epitope F22/V23/L26/I27</span></div>
   <div id="mx" style="margin-top:14px"></div>
 </div><div id="viewer"></div></div>
 <p class="note">Atomistic Boltz-2.1 co-fold of the <b>Fab</b> (heavy VH-CH1 + light VL-Cλ) with the modified <b>Tirzepatide</b> peptide (best of 5 samples). The Fc (CH2-CH3) is omitted from the co-fold — it is distal to the paratope and, when co-folded whole, misfolds into an artifact (shown in the prior round). Every construct docks the peptide at the CDR paratope on the canonical epitope.</p>
</section>

<section id="seq"><h2>Full humanized IgG1 sequences</h2>
<p class="note">Wet-lab-ready. Heavy = VH–CH1–hinge–CH2–CH3 (human IgG1). Light = VL–Cλ. One click to copy. Machine-readable: <code>results/antibody/humanized_igg_constructs.fasta</code>.</p>
<div id="seqcards"></div>
<div class="card" style="margin:14px 0"><b>Shared human constant regions</b>
 <div class="kv"><span>CH1</span><span class="seq" style="flex:1" id="cCH1"></span></div>
 <div class="kv"><span>hinge–CH2–CH3 (Fc)</span><span class="seq" style="flex:1" id="cFC"></span></div>
 <div class="kv"><span>Cλ</span><span class="seq" style="flex:1" id="cCL"></span></div>
</div>
</section>

<section id="meth">
 <h2>Methods</h2>
 <p class="note"><b>Constructs.</b> Each previously designed binder was assembled into a complete human IgG1: heavy chain = VH + human IgG1 CH1–hinge–CH2–CH3; light chain = VL + human Cλ. VH/VL are the native cognate pairs (no scFv linker). Two lineages: ab2 (matured CDR-H3, four variants incl. parent) and ab1 (orthogonal). <b>Docking / affinity.</b> Boltz-2.1 (AF3-class) structure_and_binding co-fold of the antigen-binding <b>Fab</b> (heavy VH-CH1 + light VL-Cλ, single-sequence mode) with the modified Tirzepatide peptide; 5 samples/construct. The Fc is excluded from the co-fold because it does not contact antigen and co-folding the whole IgG dilutes/artefacts the interface (established last round). <b>Metrics.</b> protein_ipTM (atomistic Fab↔peptide interface, mean & best over 5 samples), binding_confidence (Boltz binder head, affinity proxy), structure_confidence (fold). <b>Epitope.</b> peptide residues with any heavy atom within 4.5 Å of the antibody. <b>Target.</b> Tirzepatide 39-mer with Aib2/Aib13 as CCD AIB; K20 fatty-diacid on the opposite face kept as bare Lys for a controlled cross-construct comparison (validated non-paratope last round). <b>Not modelled:</b> Fc N-glycan, disulfides, MD refinement — this is a co-fold screen, to be confirmed by SPR/BLI + DSF/SEC.</p>
</section>

</main>
<script>
var D=__DATA__;
// nav
document.querySelectorAll('nav button').forEach(b=>b.onclick=()=>{
 document.querySelectorAll('nav button').forEach(x=>x.classList.remove('active'));
 document.querySelectorAll('section').forEach(x=>x.classList.remove('on'));
 b.classList.add('active');document.getElementById(b.dataset.t).classList.add('on');
 if(b.dataset.t==='dock') setTimeout(()=>{if(!window._vi)initV();else glv.resize();},60);});
document.addEventListener('click',e=>{const c=e.target.closest('.copy');if(!c)return;
 navigator.clipboard.writeText(c.dataset.s);const t=c.textContent;c.textContent='✓';setTimeout(()=>c.textContent=t,900);});
function f(v,d){return v==null?'–':(typeof v==='number'?v.toFixed(d==null?2:d):v);}

// ---- leads table ----
(function(){const tb=document.querySelector('#tblLead tbody');let sk=2,dir=1;
 function paint(rows){tb.innerHTML='';rows.forEach(r=>{const tr=document.createElement('tr');
  var bindbar=Math.round(r[2]*100/0.5);
  tr.innerHTML=`<td><span class="badge b-${r[11]}">${r[11]}</span></td><td class="mono"><b>${r[0]}</b></td><td class="mono note">${r[1]}</td>
   <td><b>${f(r[2],3)}</b><div class="bar"><i style="width:${bindbar}%"></i></div></td><td>${f(r[3],3)}</td><td>${f(r[4],3)}</td><td class="note">${f(r[5],3)}</td><td>${f(r[6],3)}</td><td>${r[7]}/5</td><td>${r[8]}</td><td class="note">${r[12]}</td>`;
  tb.appendChild(tr);});}
 function srt(k){const th=document.querySelector(`#tblLead th:nth-child(${k+1})`);const num=th.dataset.num;
  if(k===sk)dir=-dir;else{sk=k;dir=1;}
  const rows=[...D.LEADS].sort((a,b)=>{let x=a[k],y=b[k];if(x==null)x=-1;if(y==null)y=-1;
   return (num?(x-y):(''+x).localeCompare(y))*dir;});paint(rows);}
 document.querySelectorAll('#tblLead th').forEach((th,i)=>th.onclick=()=>srt(i));
 // default: sort by binding-head desc
 paint([...D.LEADS].sort((a,b)=>b[2]-a[2]));
})();

// ---- format before/after table ----
(function(){const tb=document.querySelector('#tblFmt tbody');
 D.LEADS.forEach(r=>{const scfv=r[9],fab=r[3],d=(scfv==null?null:(fab-scfv));
  var interp = r[0]=='ab2-mat1'?'lead preserved; top binder in both formats':
               r[0]=='ab2-mat2'?'A9Y edge over WT lost in IgG context':
               r[0]=='ab2-wt'?'parent; binder score now on par with A9Y':
               r[0]=='ab1-wt'?'no scFv atomistic ref; strong fold, weak binder':
               'regresses in both formats — drop';
  const tr=document.createElement('tr');
  tr.innerHTML=`<td class="mono"><b>${r[0]}</b></td><td>${f(scfv,3)}</td><td>${f(fab,3)}</td>
   <td class="${d!=null&&d<0?'':''}">${d==null?'–':(d>0?'+':'')+d.toFixed(3)}</td>
   <td>${f(r[10],3)}</td><td><b>${f(r[2],3)}</b></td><td class="note">${interp}</td>`;
  tb.appendChild(tr);});
})();

// ---- sequences ----
(function(){const box=document.getElementById('seqcards');
 const ord=['ab2-mat1','ab2-mat2','ab2-wt','ab1-wt','ab2-mat3'];
 ord.forEach(ab=>{const s=D.SEQ[ab];
  const c=document.createElement('div');c.className='card';c.style.margin='12px 0';
  c.innerHTML=`<div style="display:flex;justify-content:space-between;flex-wrap:wrap;gap:6px">
    <b style="font-size:15px">${ab} <span class="badge b-${D.TIER[ab]}" style="margin-left:6px">Tier ${D.TIER[ab]}</span></b>
    <span class="pill">CDR-H3: <span class="mono">${D.CDRH3[ab]}</span> · ${D.MUT[ab]}</span></div>
   <div class="kv"><span style="min-width:230px">Heavy chain IgG1 <span class="note">(${s.lenHC} aa)</span></span>
     <span class="seq" style="flex:1">${s.HC}<button class="copy" data-s="${s.HC}">copy</button></span></div>
   <div class="kv"><span style="min-width:230px">Light chain (Cλ) <span class="note">(${s.lenLC} aa)</span></span>
     <span class="seq" style="flex:1">${s.LC}<button class="copy" data-s="${s.LC}">copy</button></span></div>
   <div class="kv"><span style="min-width:230px">epitope contacts</span><span class="note" style="flex:1">${D.EPI[ab].join('  ')}</span></div>`;
  box.appendChild(c);});
 document.getElementById('cCH1').innerHTML=D.CONSTS.CH1+`<button class="copy" data-s="${D.CONSTS.CH1}">copy</button>`;
 document.getElementById('cFC').innerHTML=D.CONSTS.HINGE_CH2_CH3+`<button class="copy" data-s="${D.CONSTS.HINGE_CH2_CH3}">copy</button>`;
 document.getElementById('cCL').innerHTML=D.CONSTS.CL+`<button class="copy" data-s="${D.CONSTS.CL}">copy</button>`;
})();

// ---- 3D viewer ----
var glv;
function b64toStr(b){return atob(b);}
function initV(){
 if(window._3dfail||!window.$3Dmol){document.getElementById('viewer').innerHTML=
   '<div style="padding:20px;color:#93a1bd">3D library blocked offline. Structures are embedded (results/antibody/cif/*.cif) — open in PyMOL/ChimeraX.</div>';return;}
 window._vi=1;
 glv=$3Dmol.createViewer(document.getElementById('viewer'),{backgroundColor:'#0a0e17'});
 const sel=document.getElementById('pick');
 const ord=['ab2-mat1','ab2-mat2','ab2-wt','ab1-wt','ab2-mat3'];
 ord.forEach(ab=>{if(D.CIF[ab]){const o=document.createElement('option');o.value=ab;
   o.textContent=ab+'  ('+D.MUT[ab]+')';sel.appendChild(o);}});
 sel.onchange=()=>show(sel.value);
 show(ord[0]);
}
function show(ab){
 glv.clear();
 const cif=b64toStr(D.CIF[ab]);
 glv.addModel(cif,'cif');
 const m=glv.getModel();
 // chains: H (heavy VH-CH1), L (light VL-CL), P (peptide)
 m.setStyle({chain:'H'},{cartoon:{color:'#4C9AFF',opacity:.9}});
 m.setStyle({chain:'L'},{cartoon:{color:'#8FB6FF',opacity:.9}});
 m.setStyle({chain:'P'},{cartoon:{color:'#9aa7bd'},stick:{colorscheme:'grayCarbon',radius:.18}});
 // epitope core highlight on peptide (resi 22,23,24,26,27)
 [22,23,24,26,27].forEach(ri=>{
   m.setStyle({chain:'P',resi:ri},{stick:{color:'#E4572E',radius:.28}});
 });
 glv.zoomTo({chain:'P'});glv.zoom(0.75);glv.render();
 // metrics box
 var r=D.LEADS.find(x=>x[0]===ab);
 document.getElementById('mx').innerHTML=
   `<div class="kv"><span>tier</span><b>${r[11]}</b></div>
    <div class="kv"><span>binding-head</span><b>${r[2].toFixed(3)}</b></div>
    <div class="kv"><span>ipTM (mean / best)</span><b>${r[3].toFixed(3)} / ${r[4].toFixed(3)}</b></div>
    <div class="kv"><span>fold confidence</span><b>${r[6].toFixed(3)}</b></div>
    <div class="kv"><span>epitope core</span><b>${r[7]}/5</b></div>
    <div class="note" style="margin-top:8px">${r[12]}</div>`;
}
</script>
</body></html>"""

HTML = HTML.replace("__DATA__", DATA_JSON)
open("antibody_report_igg.html", "w").write(HTML)
print("wrote antibody_report_igg.html", len(HTML), "bytes")
