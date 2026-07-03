#!/usr/bin/env python3
"""Standalone interactive ANTIBODY report:
   recommendation · ranking · scores · docking (3D viewer) · Fc before/after · sequences.
Self-contained single HTML; only 3Dmol.js is fetched from a CDN at view time.
"""
import base64, json, os, re

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
def b64(path):
    with open(os.path.join(ROOT, path), "rb") as fh: return base64.b64encode(fh.read()).decode()

# ---- 3D structures to embed ----
STRUCTS = [
 ("ab2-mat1 (A8Y) · scFv + modified peptide (atomistic)", "results/antibody/atomistic/ab2mat1_A8Y_best.cif", "pep"),
 ("ab2-mat2 (A9Y) · scFv + modified peptide (atomistic)", "results/antibody/atomistic/ab2mat2_A9Y_best.cif", "pep"),
 ("ab2 WT · scFv + modified peptide (atomistic)",         "results/antibody/atomistic/ab2WT_sample1_best.cif", "pep"),
 ("ab2-mat1 (A8Y) · scFv + WHOLE-drug ligand (dock)",     "results/antibody/ligand_dock/A8Y_wholedrug.cif", "lig"),
 ("ab2-mat2 (A9Y) · scFv + WHOLE-drug ligand (dock)",     "results/antibody/ligand_dock/A9Y_wholedrug.cif", "lig"),
]
structs_js = json.dumps([{"label": l, "cif": b64(p), "kind": k} for l, p, k in STRUCTS])

# ---- ranking / scores (curated, verified) ----
# [lead, mut, unforced_bind, atom_iptm, md5_dg, md20_dg, md20_rmsd, md20_ret, ligand_cdr_contacts, tier, note]
LEADS = [
 ["ab2-mat1", "H3:A8Y", 0.731, 0.885, -47.9, -53.1, 0.48, 0.93, 14, "A", "Primary lead — strongest MM/GBSA (5 & 20 ns), stable, no CDR glycosite"],
 ["ab2-mat2", "H3:A9Y", 0.720, 0.895, -37.5, -41.6, 0.66, 1.02, 10, "A", "Co-lead — best atomistic interface; robust"],
 ["ab2 WT",   "(parent)", 0.635, 0.855, -38.9, None, None, None, None, "B", "Baseline / parent"],
 ["ab2-mat3", "H3:A2Y+A4W+A8Y", None, 0.800, None, None, None, None, None, "C", "Triple stack REGRESSES (epistasis) — do not advance"],
 ["ab1-wt",   "(orthogonal)", 0.694, None, None, None, None, None, None, "B", "Backup, different lineage (λ)"],
]
leads_js = json.dumps(LEADS)

# ---- Fc before/after ----
FC = {
 "bare": {"WT": 0.638, "A8Y": 0.664, "A9Y": 0.621},
 "mono": {"WT": (0.137, 0.12), "A8Y": (0.025, 0.008), "A9Y": (0.010, 0.002)},
 "dimer": {"WT": 1.7e-6, "A8Y": 5.6e-7, "A9Y": 5.3e-6},
}
fc_js = json.dumps(FC)

# ---- sequences ----
_constructs = json.load(open(os.path.join(ROOT, "results/antibody/matured_constructs.json")))
by_lead = {}
for r in _constructs:
    lead, part = r["name"].split("__", 1); by_lead.setdefault(lead, {})[part] = r["seq"]
LEAD_ORDER = [("ab2-mat1_H3-A8Y","ab2-mat1","H3:A8Y"),("ab2-mat2_H3-A9Y","ab2-mat2","H3:A9Y"),
              ("ab2-mat3_H3-A2Y-A4W-A8Y","ab2-mat3","H3:A2Y+A4W+A8Y"),("ab1-wt","ab1-wt","wildtype")]
PART_ORDER=[("__VH","VH"),("VL_ab2_lambda","VL (λ)"),("scFv_VH-G4Sx3-VL","scFv  (VH–(G4S)₃–VL)"),
            ("scFv-Fc","scFv-Fc  (+ human IgG1 hinge-CH2-CH3)"),("IgG1_HeavyChain","IgG1 Heavy chain"),
            ("IgG1_LightChain","IgG1 Light chain (Cλ)")]
def h3(vh):
    m=re.search(r"Y[YF]C[A-Z]{2}([A-Z]+?)WG.G", vh or ""); return m.group(1) if m else ""
def seqrow(label,n,seq):
    return (f'<div class="kv" style="align-items:flex-start"><span style="min-width:210px">{label} '
            f'<span class="note">({n} aa)</span></span><span class="seq" style="flex:1">{seq}'
            f'<button class="copy" data-s="{seq}">copy</button></span></div>')
seq_html=['<p class="note">Every wet-lab-ready sequence, one click to copy. Light chain = λ. '
          'Machine-readable: <code>results/antibody/matured_constructs.fasta</code>.</p>']
for key,short,mut in LEAD_ORDER:
    parts=by_lead.get(key,{}); vh=next((v for k,v in parts.items() if k.endswith("VH")),None)
    seq_html.append(f'<div class="card" style="margin:12px 0"><div style="display:flex;justify-content:space-between;flex-wrap:wrap;gap:6px">'
                    f'<b style="font-size:15px">{short}</b><span class="pill">CDR-H3: <span class="mono">{h3(vh)}</span> · {mut}</span></div>')
    for sub,lab in PART_ORDER:
        k=next((k for k in parts if sub in k),None)
        if k: seq_html.append(seqrow(lab,len(parts[k]),parts[k]))
    seq_html.append('</div>')
const=[]; hdr=None
for line in open(os.path.join(ROOT,"results/antibody/human_constant_regions.fasta")):
    line=line.rstrip("\n")
    if line.startswith(">"): hdr=line[1:]
    elif line and hdr: const.append((hdr,len(line),line)); hdr=None
seq_html.append('<h2>Human constant-region cassettes</h2><div class="card" style="margin:10px 0">'
                +"".join(seqrow(h,n,s) for h,n,s in const)+'</div>')
SEQ_HTML="\n".join(seq_html)

html = r"""<!DOCTYPE html><html lang="en"><head>
<meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Anti-Tirzepatide Antibody Report</title>
<script>(function(){var s=["https://cdn.jsdelivr.net/npm/3dmol@2.1.0/build/3Dmol-min.js",
"https://unpkg.com/3dmol@2.1.0/build/3Dmol-min.js","https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.1.0/3Dmol-min.js"];
var i=0;function n(){if(i>=s.length){window._3dfail=1;return;}var e=document.createElement('script');
e.src=s[i++];e.onload=function(){window._3dok=1};e.onerror=n;document.head.appendChild(e);}n();})();</script>
<style>
:root{--bg:#0f1420;--panel:#171e2e;--ink:#e8edf6;--mut:#93a1bd;--acc:#4C9AFF;--acc2:#E4572E;--ok:#54A24B;--line:#26304a}
*{box-sizing:border-box}body{margin:0;font:15px/1.55 -apple-system,Segoe UI,Roboto,Helvetica,Arial,"PingFang SC","Microsoft YaHei",sans-serif;background:var(--bg);color:var(--ink)}
header{padding:22px 26px;border-bottom:1px solid var(--line);background:linear-gradient(180deg,#141b2b,#0f1420)}
h1{margin:0 0 4px;font-size:21px}.sub{color:var(--mut);font-size:13px}
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
</style></head><body>
<header><h1>Anti-Tirzepatide Antibody — Design, Scoring &amp; Docking Report</h1>
<div class="sub">Target = <b>modified Tirzepatide</b> (Aib2/Aib13 + K20 C20-diacid-γGlu-(AEEA)₂). Boltz-2.1 (AF3-class) co-folding + OpenMM MD/MM-GBSA. λ light chain.</div></header>
<nav>
 <button data-t="rec" class="active">Recommendation</button>
 <button data-t="rank">Ranking &amp; Scores</button>
 <button data-t="dock">Docking (3D)</button>
 <button data-t="fc">Fc before/after</button>
 <button data-t="seq">Sequences</button>
 <button data-t="meth">Methods</button>
</nav><main>

<section id="rec" class="on">
 <div class="cards">
  <div class="card"><div class="n">ab2-mat1</div><div class="l">primary lead · CDR-H3 A8Y</div></div>
  <div class="card"><div class="n">−53.1</div><div class="l">A8Y MM/GBSA ΔG, 20 ns (kcal/mol)</div></div>
  <div class="card"><div class="n">0.885/0.895</div><div class="l">A8Y/A9Y atomistic protein_ipTM</div></div>
  <div class="card"><div class="n">0</div><div class="l">CDR glycosylation sequons (leads)</div></div>
 </div>
 <div class="good"><b>Recommendation.</b> Advance <b>ab2-mat1 (H3:A8Y)</b> as the primary lead and
 <b>ab2-mat2 (H3:A9Y)</b> as co-lead, formatted as <b>human IgG1 or scFv-Fc</b> (bivalent → avidity), λ light.
 Both bind the modified drug at the CDR paratope (matured CDR-H3 + VL CDRs), form complexes stable over 20 ns MD,
 and carry no CDR glycosylation liability. <b>ab1-wt</b> is an orthogonal-lineage backup. Do <b>not</b> advance the
 triple mutant (epistasis regression). Confirm empirically by SPR/BLI + DSF/SEC.</div>
 <h2>Why these leads (5 independent, consistent lines of evidence)</h2>
 <ul>
  <li><b>Affinity maturation</b> of scFv <b>ab2</b>: CDR-H3 Ala→Tyr at pos 8/9 raised un-forced binding_confidence from 0.635 (WT) to <b>0.731 (A8Y) / 0.720 (A9Y)</b>.</li>
  <li><b>Atomistic interface</b> (Boltz structure_and_binding): protein_ipTM <b>0.885 / 0.895</b> (WT 0.855).</li>
  <li><b>MD MM/GBSA</b>: A8Y ΔG −47.9 (5 ns) → <b>−53.1 (20 ns)</b>; A9Y −37.5 → −41.6; complexes stable (binder RMSD &lt;0.7 nm, contact retention ≥0.93).</li>
  <li><b>Orthogonal whole-drug docking</b>: the entire drug (as one ligand) docks at the <b>CDR paratope</b> (CDR-H3 + VL CDRs; A8Y 14 CDR contacts).</li>
  <li><b>Developability</b>: no N-glyc sequon in the Fv/CDRs; no free Cys.</li>
 </ul>
 <p class="note">Epitope = Tirzepatide C-terminal amphipathic face <b>F22·V23·L26·I27</b> (+ E3, PPPS). K20 lipid sits on the opposite face.</p>
</section>

<section id="rank">
 <h2>Ranking &amp; scores <span class="pill">click a header to sort</span></h2>
 <table id="tblLead"><thead><tr>
  <th>Tier</th><th>Lead</th><th>CDR-H3 mut</th><th data-num="1">bind (un-forced)</th><th data-num="1">atomistic ipTM</th>
  <th data-num="1">ΔG 5 ns</th><th data-num="1">ΔG 20 ns</th><th data-num="1">binder RMSD 20 ns</th><th data-num="1">contact ret.</th><th data-num="1">ligand-dock CDR contacts</th><th>note</th>
 </tr></thead><tbody></tbody></table>
 <p class="note">bind = Boltz binding_confidence (un-forced, +lipid); ipTM = atomistic protein_ipTM (scFv+peptide); ΔG = chain-split MM/GBSA (kcal/mol, ranking signal within scaffold); RMSD/retention = MD pose stability.</p>
 <div class="warn"><b>Honest note on ranking.</b> MM/GBSA is compared <b>within the ab2 scaffold</b> (A8Y &gt; A9Y ≈ WT). Do not rank across different scaffolds by raw ΔG (interface-size bias). No single affinity score predicts developability failures — pair with liability screening.</div>
</section>

<section id="dock">
 <h2>Docking viewer <span class="pill">drag=rotate · scroll=zoom</span></h2>
 <div class="viewwrap"><div>
   <select id="pick"></select>
   <div class="legend" style="margin-top:10px">
    <span><i class="dot" style="background:#4C9AFF"></i>antibody (scFv)</span>
    <span><i class="dot" style="background:#9aa7bd"></i>Tirzepatide peptide</span>
    <span><i class="dot" style="background:#E4572E"></i>epitope F22/V23/L26/I27</span>
    <span><i class="dot" style="background:#F2C14E"></i>whole-drug ligand</span></div>
   <div id="mx" style="margin-top:14px"></div>
 </div><div id="viewer"></div></div>
 <p class="note">Two orthogonal representations: (1) atomistic scFv + modified <b>peptide</b> (protein_ipTM 0.86–0.90); (2) scFv + the entire drug as one <b>small-molecule ligand</b> — it docks at the CDR paratope, confirming the binding site. Both agree.</p>
</section>

<section id="fc">
 <h2>Fc before / after — does adding the humanized IgG1 Fc change binding?</h2>
 <div class="good"><b>No.</b> The Fc is a separate module across a flexible hinge and never contacts the CDRs, so per-arm binding = the bare-scFv result (unchanged); the bivalent Fc adds <b>avidity</b> (a gain). The only glycosylation site introduced is the canonical, distal <b>Fc N297</b> glycan — not a paratope liability.</div>
 <table id="tblFc"><thead><tr><th>construct</th><th data-num="1">WT</th><th data-num="1">A8Y</th><th data-num="1">A9Y</th><th>interpretation</th></tr></thead><tbody></tbody></table>
 <div class="warn"><b>Why the Fc rows collapse.</b> Whole-antibody Boltz co-folding cannot score the Fc format: a lone scFv-Fc chain misfolds (CH3 is an obligate homodimer → structure_confidence ≈ 0), and the ~980-aa flexible dimer gives diluted global scores (all constructs → ~1e-6, uniformly). These are <b>modeling artifacts, not loss of binding</b>. Quantify the Fc format by MD (dimer+peptide) or experiment (SPR/BLI), using the bare-scFv values as the valid per-arm readout.</div>
</section>

<section id="seq"><h2>Antibody sequences</h2>%SEQ_HTML%</section>

<section id="meth">
 <h2>Methods</h2>
 <ul>
  <li><b>Co-folding / docking:</b> Boltz-2.1 (AlphaFold3-class; no separate AF3 endpoint available) — protein_screen (epitope-forced &amp; un-forced, +K20 lipid), structure_and_binding (atomistic protein_ipTM + CIF), and whole-drug-as-ligand co-fold.</li>
  <li><b>MD:</b> OpenMM ff14SB/TIP3P, 4 fs HMR, 0.2 ns NVT + 0.2 ns NPT + 5/20 ns production, on Modal GPUs (<code>mdscreen</code>). Chain-split protein–peptide MM/GBSA (GBn2, PBC-unwrapped) + pose stability.</li>
  <li><b>Target:</b> peptide chain (Aib2/Aib13 as CCD AIB) + K20 acyl chain (from user SMILES) as co-folded ligand; the whole modified drug also modelled as a single SMILES ligand for the orthogonal dock.</li>
 </ul>
 <div class="kv"><span>binding_confidence</span><span>Boltz binder head; un-forced = honest.</span></div>
 <div class="kv"><span>protein_ipTM</span><span>atomistic scFv↔peptide interface, 0–1; &gt;0.8 confident. Cleanest metric.</span></div>
 <div class="kv"><span>MM/GBSA ΔG</span><span>dynamics-averaged, entropy omitted → ranking within a scaffold, not absolute affinity.</span></div>
 <p class="note">Full write-ups: <code>results/MASTER_SUMMARY.md</code>, <code>MATURATION_REPORT.md</code>, <code>MD_CROSSVALIDATION.md</code>, <code>antibody/fc/FC_EVALUATION.md</code>, <code>antibody/ligand_dock/LIGAND_DOCK.md</code>.</p>
</section>
</main>
<script>
const STRUCTS=%STRUCTS%, LEADS=%LEADS%, FC=%FC%;
document.querySelectorAll('nav button').forEach(b=>b.onclick=()=>{
 document.querySelectorAll('nav button').forEach(x=>x.classList.remove('active'));
 document.querySelectorAll('section').forEach(x=>x.classList.remove('on'));
 b.classList.add('active');document.getElementById(b.dataset.t).classList.add('on');
 if(b.dataset.t==='dock') setTimeout(()=>{if(!window._vi)initV();else glv.resize();},60);});
document.addEventListener('click',e=>{const c=e.target.closest('.copy');if(!c)return;
 navigator.clipboard.writeText(c.dataset.s);const t=c.textContent;c.textContent='✓';setTimeout(()=>c.textContent=t,900);});
function f(v,d){return v==null?'–':(typeof v==='number'?v.toFixed(d==null?2:d):v);}
// leads table
(function(){const tb=document.querySelector('#tblLead tbody');let sk=6,dir=1;
 function paint(rows){tb.innerHTML='';rows.forEach(r=>{const tr=document.createElement('tr');
  tr.innerHTML=`<td><span class="badge b-${r[9]}">${r[9]}</span></td><td class="mono"><b>${r[0]}</b></td><td class="mono note">${r[1]}</td>
   <td><b>${f(r[2])}</b></td><td>${f(r[3],3)}</td><td>${f(r[4],1)}</td><td><b>${f(r[5],1)}</b></td><td>${f(r[6])}</td><td>${f(r[7])}</td><td>${r[8]==null?'–':r[8]}</td><td class="note">${r[10]}</td>`;
  tb.appendChild(tr);});}
 function srt(k){const num=document.querySelector(`#tblLead th:nth-child(${k+1})`).dataset.num;
  if(k===sk)dir=-dir;else{sk=k;dir=num?1:1;}
  const rows=[...LEADS].sort((a,b)=>{let x=a[k],y=b[k];if(x==null)x=num?-1e9:'';if(y==null)y=num?-1e9:'';
   // for ΔG (more negative=better) sort ascending by default handled by user clicks
   return num?(x-y)*dir:(''+x).localeCompare(''+y)*dir;});paint(rows);}
 document.querySelectorAll('#tblLead th').forEach((th,k)=>th.onclick=()=>srt(k));
 paint(LEADS);})();
// Fc table
(function(){const tb=document.querySelector('#tblFc tbody');
 const rows=[['bare scFv (paratope) — forced+lipid',FC.bare.WT,FC.bare.A8Y,FC.bare.A9Y,'valid readout: healthy binding'],
  ['scFv-Fc monomer — binding_conf',FC.mono.WT[0],FC.mono.A8Y[0],FC.mono.A9Y[0],'lone Fc misfolds (CH3 homodimer) → artifact'],
  ['scFv-Fc monomer — structure_conf',FC.mono.WT[1],FC.mono.A8Y[1],FC.mono.A9Y[1],'≈0 fold confidence (artifact)'],
  ['scFv-Fc homodimer — binding_conf',FC.dimer.WT,FC.dimer.A8Y,FC.dimer.A9Y,'~1e-6 uniformly → non-discriminative artifact']];
 rows.forEach(r=>{const tr=document.createElement('tr');
  tr.innerHTML=`<td>${r[0]}</td><td>${(+r[1]).toExponential? (r[1]<0.01?(+r[1]).toExponential(1):(+r[1]).toFixed(3)):r[1]}</td>
   <td>${r[2]<0.01?(+r[2]).toExponential(1):(+r[2]).toFixed(3)}</td><td>${r[3]<0.01?(+r[3]).toExponential(1):(+r[3]).toFixed(3)}</td><td class="note">${r[4]}</td>`;
  tb.appendChild(tr);});})();
// viewer
let glv;
function initV(){ if(typeof $3Dmol==='undefined'){ if(window._3dfail){document.getElementById('viewer').innerHTML=
  '<div style="padding:24px;color:#93a1bd">3Dmol.js (CDN) blocked here. CIFs are in results/antibody/atomistic/ &amp; ligand_dock/ — open in PyMOL/ChimeraX.</div>';return;} setTimeout(initV,400);return;}
 window._vi=1;glv=$3Dmol.createViewer(document.getElementById('viewer'),{backgroundColor:'#0a0e17'});
 const sel=document.getElementById('pick');STRUCTS.forEach((s,i)=>{const o=document.createElement('option');o.value=i;o.textContent=s.label;sel.appendChild(o);});
 sel.onchange=()=>load(+sel.value);load(0);}
function load(i){glv.clear();glv.addModel(atob(STRUCTS[i].cif),'cif');const k=STRUCTS[i].kind;
 glv.setStyle({chain:'A'},{cartoon:{color:'#4C9AFF'}});
 if(k==='pep'){glv.setStyle({chain:'T'},{cartoon:{color:'#9aa7bd'},stick:{radius:0.1,color:'#9aa7bd'}});
   glv.addStyle({chain:'T',resi:[3,18,19,22,23,26,27]},{stick:{radius:0.3,color:'#E4572E'}});}
 else{glv.setStyle({chain:'L'},{stick:{radius:0.18,color:'#F2C14E'}});}
 glv.zoomTo();glv.render();
 const mx=document.getElementById('mx');const L=LEADS.find(x=>STRUCTS[i].label.indexOf(x[0])>=0);
 mx.innerHTML=L?`<div class="kv"><span>lead</span><b>${L[0]} (${L[1]})</b></div>
  <div class="kv"><span>atomistic ipTM</span><span>${f(L[3],3)}</span></div>
  <div class="kv"><span>ΔG 20 ns</span><span>${f(L[5],1)} kcal/mol</span></div>
  <div class="note" style="margin-top:6px">${L[10]}</div>`:'';}
</script></body></html>"""

html = (html.replace("%SEQ_HTML%", SEQ_HTML).replace("%STRUCTS%", structs_js)
            .replace("%LEADS%", leads_js).replace("%FC%", fc_js))
out = os.path.join(ROOT, "report", "antibody_report.html")
open(out, "w").write(html)
print("wrote", out, round(os.path.getsize(out)/1024), "KB")
