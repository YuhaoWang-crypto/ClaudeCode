#!/usr/bin/env python
"""Build the interactive Raygun alignment viewer HTML from the results JSON."""
import json, html

with open("/home/user/ClaudeCode/scratchpad/raygun_8_8M_results.json") as f:
    data = json.load(f)

# Human-friendly variant descriptions
DESC = {
    "v1_same_lownoise":  "Same length, low noise — faithful reconstruction",
    "v2_same_highnoise": "Same length, high noise — substitution scan",
    "v3_shortened":      "Shortened — deletion design",
    "v4_extended":       "Extended — insertion design",
}
for v in data["variants"]:
    v["desc"] = DESC.get(v["label"], v["label"])
    v["delta"] = v["len"] - data["template"]["len"]

payload = json.dumps(data)

HTML = """<title>Raygun × mCherry — Redesign Alignment Explorer</title>
<style>
:root{
  --bg:#faf7fb; --panel:#ffffff; --ink:#1c1620; --muted:#6b6270;
  --line:#e7dfe9; --line-soft:#f0eaf2;
  --accent:#d81b60; --accent-soft:#fce4ec;
  --hydro:#2a9d8f; --polar:#e9c46a; --acidic:#e63946; --basic:#4361ee;
  --special:#9b5de5; --other:#9aa0a6;
  --match:#cbb; --sub:#f4a259; --gap:#e63946;
  --shadow:0 1px 2px rgba(20,10,30,.05),0 8px 24px rgba(20,10,30,.06);
  --mono:ui-monospace,"SF Mono","JetBrains Mono",Menlo,Consolas,monospace;
  --sans:system-ui,-apple-system,"Segoe UI",Roboto,sans-serif;
}
@media (prefers-color-scheme:dark){
  :root{
    --bg:#14101a; --panel:#1e1826; --ink:#f2ecf4; --muted:#a89db2;
    --line:#2f2739; --line-soft:#241d2e;
    --accent:#ff5c8a; --accent-soft:#3a1424;
    --shadow:0 1px 2px rgba(0,0,0,.3),0 10px 30px rgba(0,0,0,.35);
  }
}
:root[data-theme="dark"]{
  --bg:#14101a; --panel:#1e1826; --ink:#f2ecf4; --muted:#a89db2;
  --line:#2f2739; --line-soft:#241d2e;
  --accent:#ff5c8a; --accent-soft:#3a1424;
  --shadow:0 1px 2px rgba(0,0,0,.3),0 10px 30px rgba(0,0,0,.35);
}
:root[data-theme="light"]{
  --bg:#faf7fb; --panel:#ffffff; --ink:#1c1620; --muted:#6b6270;
  --line:#e7dfe9; --line-soft:#f0eaf2;
  --accent:#d81b60; --accent-soft:#fce4ec;
  --shadow:0 1px 2px rgba(20,10,30,.05),0 8px 24px rgba(20,10,30,.06);
}
*{box-sizing:border-box}
body{background:var(--bg);color:var(--ink);font-family:var(--sans);
  line-height:1.55;margin:0;-webkit-font-smoothing:antialiased}
.wrap{max-width:1080px;margin:0 auto;padding:clamp(20px,4vw,52px) clamp(16px,3vw,32px) 72px}
.eyebrow{font-size:12px;letter-spacing:.16em;text-transform:uppercase;
  color:var(--accent);font-weight:700;margin:0 0 10px}
h1{font-size:clamp(28px,5vw,44px);line-height:1.08;margin:0 0 14px;
  text-wrap:balance;font-weight:750;letter-spacing:-.02em}
.lede{color:var(--muted);font-size:17px;max-width:64ch;margin:0}
.badges{display:flex;flex-wrap:wrap;gap:8px;margin:22px 0 0}
.badge{font-family:var(--mono);font-size:12px;padding:5px 11px;border-radius:999px;
  border:1px solid var(--line);background:var(--panel);color:var(--muted)}
.badge b{color:var(--ink);font-weight:600}
.card{background:var(--panel);border:1px solid var(--line);border-radius:16px;
  box-shadow:var(--shadow);padding:clamp(18px,3vw,28px);margin-top:26px}
.section-label{font-size:12px;letter-spacing:.14em;text-transform:uppercase;
  color:var(--muted);font-weight:700;margin:0 0 16px}
/* pill selector */
.pills{display:flex;flex-wrap:wrap;gap:10px}
.pill{cursor:pointer;border:1px solid var(--line);background:var(--bg);
  border-radius:12px;padding:12px 15px;text-align:left;min-width:0;flex:1 1 190px;
  transition:border-color .15s,transform .05s;font-family:inherit;color:inherit}
.pill:hover{border-color:var(--accent)}
.pill:focus-visible{outline:2px solid var(--accent);outline-offset:2px}
.pill[aria-selected="true"]{border-color:var(--accent);
  background:var(--accent-soft);box-shadow:inset 0 0 0 1px var(--accent)}
.pill .pt{font-weight:650;font-size:14px;margin-bottom:3px}
.pill .ps{font-size:12px;color:var(--muted);font-family:var(--mono)}
/* stats */
.stats{display:grid;grid-template-columns:repeat(auto-fit,minmax(120px,1fr));
  gap:1px;background:var(--line);border:1px solid var(--line);
  border-radius:12px;overflow:hidden;margin-top:20px}
.stat{background:var(--panel);padding:14px 16px}
.stat .k{font-size:11px;letter-spacing:.08em;text-transform:uppercase;color:var(--muted)}
.stat .v{font-size:24px;font-weight:700;font-variant-numeric:tabular-nums;
  margin-top:4px;font-family:var(--mono)}
.stat .v small{font-size:13px;color:var(--muted);font-weight:500}
.delta-pos{color:var(--hydro)} .delta-neg{color:var(--accent)}
/* length bars */
.lenbars{margin-top:20px;display:flex;flex-direction:column;gap:9px}
.lenrow{display:grid;grid-template-columns:78px 1fr auto;gap:12px;align-items:center;font-size:13px}
.lenrow .nm{color:var(--muted);font-family:var(--mono)}
.track{height:12px;background:var(--line-soft);border-radius:6px;overflow:hidden}
.fill{height:100%;border-radius:6px;background:var(--accent)}
.fill.tmpl{background:var(--muted)}
.lenrow .amt{font-family:var(--mono);font-variant-numeric:tabular-nums;font-size:12px;color:var(--muted)}
/* alignment */
.aln-head{display:flex;justify-content:space-between;align-items:center;
  flex-wrap:wrap;gap:12px;margin-bottom:14px}
.toggle{display:inline-flex;border:1px solid var(--line);border-radius:10px;overflow:hidden}
.toggle button{border:0;background:var(--panel);color:var(--muted);cursor:pointer;
  padding:8px 14px;font-family:var(--sans);font-size:13px;font-weight:600}
.toggle button[aria-pressed="true"]{background:var(--accent);color:#fff}
.aln{overflow-x:auto;border-radius:12px;border:1px solid var(--line);background:var(--bg)}
.blk{font-family:var(--mono);font-size:13.5px;line-height:1.15;padding:14px 16px;
  white-space:pre;min-width:max-content}
.blk .lab{color:var(--muted);user-select:none}
.blk .rlr{color:var(--muted);opacity:.7}
.res{border-radius:2px}
.mode-diff .m-sub{background:var(--sub);color:#1c1206;font-weight:600}
.mode-diff .m-gap{background:var(--gap);color:#fff}
.mode-diff .m-match{color:var(--muted)}
.midline{color:var(--muted)}
/* legend */
.legend{display:flex;flex-wrap:wrap;gap:8px 16px;margin-top:16px;font-size:12.5px;color:var(--muted)}
.legend span{display:inline-flex;align-items:center;gap:6px}
.sw{width:13px;height:13px;border-radius:3px;display:inline-block}
.foot{color:var(--muted);font-size:12.5px;margin-top:30px;line-height:1.7}
.foot b{color:var(--ink)}
a{color:var(--accent)}
</style>

<div class="wrap">
  <p class="eyebrow">Template-based protein redesign</p>
  <h1>Raygun redesigns of mCherry</h1>
  <p class="lede">The Raygun autoencoder compresses a template protein into a
  fixed-length latent from ESM-2 650M embeddings, then decodes new sequences at
  a chosen length &mdash; performing controlled insertions, deletions and
  substitutions. Each variant below is aligned back to the mCherry template
  (Needleman&ndash;Wunsch, BLOSUM62) to show exactly what changed.</p>
  <div class="badges" id="badges"></div>

  <div class="card">
    <p class="section-label">Choose a variant</p>
    <div class="pills" id="pills" role="tablist"></div>
    <div class="stats" id="stats"></div>
    <div class="lenbars" id="lenbars"></div>
  </div>

  <div class="card">
    <div class="aln-head">
      <p class="section-label" style="margin:0">Sequence alignment</p>
      <div class="toggle" role="group" aria-label="Color mode">
        <button id="mDiff" aria-pressed="true">Difference</button>
        <button id="mProp" aria-pressed="false">Property</button>
      </div>
    </div>
    <div class="aln" id="aln"></div>
    <div class="legend" id="legend"></div>
  </div>

  <p class="foot" id="foot"></p>
</div>

<script>
const DATA = __PAYLOAD__;
const T = DATA.template;
let cur = 0, mode = "diff";

// amino-acid property groups -> css var
const PROP = {};
"AVLIMFWPGC".split("").forEach(a=>PROP[a]="hydro");
"STYNQ".split("").forEach(a=>PROP[a]="polar");
"DE".split("").forEach(a=>PROP[a]="acidic");
"KRH".split("").forEach(a=>PROP[a]="basic");
const PROPCOL={hydro:"--hydro",polar:"--polar",acidic:"--acidic",basic:"--basic",other:"--other"};

function esc(s){return s.replace(/&/g,"&amp;").replace(/</g,"&lt;")}

// badges
document.getElementById("badges").innerHTML =
  `<span class="badge">model <b>${esc(DATA.model)}</b></span>`+
  `<span class="badge">template <b>${esc(T.name)}</b></span>`+
  `<span class="badge">length <b>${T.len} aa</b></span>`+
  `<span class="badge">variants <b>${DATA.variants.length}</b></span>`+
  `<span class="badge">encoder <b>ESM-2 650M</b></span>`;

// pills
const pills = document.getElementById("pills");
DATA.variants.forEach((v,i)=>{
  const b=document.createElement("button");
  b.className="pill"; b.role="tab"; b.setAttribute("aria-selected", i===0);
  b.innerHTML=`<div class="pt">${esc(v.desc)}</div>`+
    `<div class="ps">len ${v.len} · noise ${v.noise} · id ${v.pct_identity}%</div>`;
  b.onclick=()=>{cur=i;render()};
  pills.appendChild(b);
});

function statCard(k,val,cls){return `<div class="stat"><div class="k">${k}</div>`+
  `<div class="v ${cls||""}">${val}</div></div>`}

function renderStats(){
  const v=DATA.variants[cur];
  const dsign = v.delta>0?`+${v.delta}`:`${v.delta}`;
  const dcls = v.delta>0?"delta-pos":(v.delta<0?"delta-neg":"");
  document.getElementById("stats").innerHTML =
    statCard("Length", `${v.len}<small> aa</small>`)+
    statCard("Length &Delta;", `${dsign}<small> aa</small>`, dcls)+
    statCard("Identity", `${v.pct_identity}<small>%</small>`)+
    statCard("Similarity", `${v.pct_similarity}<small>%</small>`)+
    statCard("Substitutions", v.mismatches)+
    statCard("Gap columns", v.gap_cols)+
    statCard("Noise", v.noise)+
    statCard("BLOSUM score", v.score);
}

function renderLenBars(){
  const v=DATA.variants[cur];
  const mx=Math.max(T.len, ...DATA.variants.map(x=>x.len));
  const bar=(nm,len,cls)=>`<div class="lenrow"><span class="nm">${nm}</span>`+
    `<span class="track"><span class="fill ${cls}" style="width:${100*len/mx}%"></span></span>`+
    `<span class="amt">${len} aa</span></div>`;
  document.getElementById("lenbars").innerHTML =
    bar("template",T.len,"tmpl")+bar("variant",v.len,"");
}

function colorRes(ch, other){
  if(mode==="prop"){
    if(ch==="-") return `<span class="res m-gap" style="background:var(--gap);color:#fff">-</span>`;
    const g=PROP[ch]||"other";
    return `<span class="res" style="color:var(${PROPCOL[g]});font-weight:600">${ch}</span>`;
  }
  // difference mode
  if(ch==="-") return `<span class="res m-gap">-</span>`;
  if(other==="-") return `<span class="res m-sub">${ch}</span>`; // aligned to a gap = indel side
  if(ch===other) return `<span class="res m-match">${ch}</span>`;
  return `<span class="res m-sub">${ch}</span>`;
}

function renderAln(){
  const v=DATA.variants[cur];
  const a=v.ref_aln, b=v.var_aln, W=60;
  let out="";
  for(let i=0;i<a.length;i+=W){
    const ra=a.slice(i,i+W), rb=b.slice(i,i+W);
    let mid="";
    for(let j=0;j<ra.length;j++){
      const x=ra[j], y=rb[j];
      if(x==="-"||y==="-") mid+=" ";
      else if(x===y) mid+="|";
      else mid+="."; // could be similar; kept simple
    }
    let la="", lb="";
    for(let j=0;j<ra.length;j++){la+=colorRes(ra[j],rb[j]); lb+=colorRes(rb[j],ra[j]);}
    out+=`<span class="rlr">${String(i+1).padStart(4)} </span>\n`;
    out+=`<span class="lab">mCherry  </span>${la}\n`;
    out+=`<span class="lab">         </span><span class="midline">${mid}</span>\n`;
    out+=`<span class="lab">variant  </span>${lb}\n\n`;
  }
  const el=document.getElementById("aln");
  el.innerHTML=`<div class="blk mode-${mode}">${out}</div>`;
}

function renderLegend(){
  const L=document.getElementById("legend");
  if(mode==="diff"){
    L.innerHTML =
      `<span><i class="sw" style="background:var(--sub)"></i>substitution / indel residue</span>`+
      `<span><i class="sw" style="background:var(--gap)"></i>gap</span>`+
      `<span><i class="sw" style="background:var(--muted)"></i>identical (dimmed)</span>`+
      `<span>| identical &nbsp; . substitution</span>`;
  } else {
    L.innerHTML =
      `<span><i class="sw" style="background:var(--hydro)"></i>hydrophobic (AVLIMFWPGC)</span>`+
      `<span><i class="sw" style="background:var(--polar)"></i>polar (STYNQ)</span>`+
      `<span><i class="sw" style="background:var(--acidic)"></i>acidic (DE)</span>`+
      `<span><i class="sw" style="background:var(--basic)"></i>basic (KRH)</span>`+
      `<span><i class="sw" style="background:var(--gap)"></i>gap</span>`;
  }
}

function render(){
  [...pills.children].forEach((p,i)=>p.setAttribute("aria-selected", i===cur));
  renderStats(); renderLenBars(); renderAln(); renderLegend();
}

document.getElementById("mDiff").onclick=()=>{mode="diff";
  mDiff.setAttribute("aria-pressed",true); mProp.setAttribute("aria-pressed",false); renderAln(); renderLegend();};
document.getElementById("mProp").onclick=()=>{mode="prop";
  mProp.setAttribute("aria-pressed",true); mDiff.setAttribute("aria-pressed",false); renderAln(); renderLegend();};

document.getElementById("foot").innerHTML =
  `<b>Method.</b> mCherry (236 aa) embedded with ESM-2 650M (layer 33); Raygun `+
  `<b>${esc(DATA.model)}</b> decodes each variant at the target length on CPU. `+
  `Alignments: global Needleman&ndash;Wunsch, BLOSUM62 (gap open &minus;10, extend &minus;0.5). `+
  `Identity = identical columns / aligned non-gap columns; similarity adds positive-BLOSUM pairs. `+
  `These are <b>sequence</b>-level differences (primary structure); 3D-fold validation would require a `+
  `structure predictor and is not part of this demo.`;

render();
</script>
"""

out = HTML.replace("__PAYLOAD__", payload)
with open("/home/user/ClaudeCode/scratchpad/raygun_viewer.html", "w") as f:
    f.write(out)
print("wrote viewer html,", len(out), "bytes")
