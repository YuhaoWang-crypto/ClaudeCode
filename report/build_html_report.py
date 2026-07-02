#!/usr/bin/env python3
"""Build a single self-contained interactive HTML report:
   - sortable sequence/metrics tables
   - embedded 3D docking viewer (3Dmol.js, CIF embedded as base64)
   - maturation figure embedded as base64
No server needed; only 3Dmol.js is fetched from CDN at view time.
"""
import base64, json, os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# ---- structures to embed in the viewer (label -> path, chains present) ----
STRUCTS = [
 ("TZP-S1 · 3-helix bundle (75 aa) + modified target", "results/cif_stable/TZP-S1_bundle.cif", True),
 ("TZP-S2 · mixed α/β (49 aa) + modified target",       "results/cif_stable/TZP-S2_mixedab.cif", True),
 ("TZP-S3 · α-hairpin (47 aa) + modified target",       "results/cif_stable/TZP-S3_hairpin.cif", True),
 ("TZP-S4 · stapled helix (34 aa) + modified target",   "results/cif_stable/TZP-S4_stapled.cif", True),
 ("TZP-P3 · lead peptide (34 aa) + target",             "results/cif_leads/R2_035.cif", False),
 ("TZP-P2 · lead peptide (26 aa) + target",             "results/cif_leads/R2_002.cif", False),
 ("TZP-B1 · spec_13 mini-protein + modified target",    "results/cif_modified/TZP-B1_spec13.cif", True),
 ("spec_13 · existing binder + target (lipid-free)",    "results/cif/design_spec_13.cif", False),
 ("ab2-mat2 (H3:A9Y) scFv + peptide · atomistic",       "results/antibody/atomistic/ab2mat2_A9Y_best.cif", False),
 ("ab2-mat1 (H3:A8Y) scFv + peptide · atomistic",       "results/antibody/atomistic/ab2mat1_A8Y_best.cif", False),
 ("ab2 WT scFv + peptide · atomistic (baseline)",       "results/antibody/atomistic/ab2WT_sample1_best.cif", False),
 ("ab2-mat3 (triple) scFv + peptide · atomistic",       "results/antibody/atomistic/ab2mat3_triple_best.cif", False),
]

def b64(path):
    with open(os.path.join(ROOT, path), "rb") as fh:
        return base64.b64encode(fh.read()).decode()

structs_js = [{"label": lbl, "cif": b64(p), "lipid": lip} for lbl, p, lip in STRUCTS]

fig_b64 = b64("results/figure_maturation.png")

# ---- Antibody / maturation data (loaded from result JSONs) ----
_mat  = json.load(open(os.path.join(ROOT,"results/maturation_ranked.json")))
_val  = {x["id"]:x for x in json.load(open(os.path.join(ROOT,"results/validation_ranked.json")))}
_atom = json.load(open(os.path.join(ROOT,"results/antibody/atomistic/metrics.json")))
_pf   = {x["id"]:x for x in json.load(open(os.path.join(ROOT,"results/prior_forced_ranked.json")))}
_pu   = {x["id"]:x for x in json.load(open(os.path.join(ROOT,"results/prior_unforced_ranked.json")))}
_plip = {x["id"]:x for x in json.load(open(os.path.join(ROOT,"results/prior_lipid_sensitivity.json")))}
_prior= {o["id"]:o for o in json.load(open(os.path.join(ROOT,"design/prior_designs.json")))}

# curated matured-lead summary: [lead, mutation, forced, unforced, atom_piptm, atom_ppbind, verdict]
AB_LEADS = [
 ["ab2-mat2","H3:A9Y",0.687,_val["V_A9Y"]["bind_unforced"],_atom["ab2mat2_A9Y"]["protein_iptm_best"],_atom["ab2mat2_A9Y"]["pp_binding_confidence"],"Best atomistic interface (ipTM 0.895, PAE 0.84 Å)"],
 ["ab2-mat1","H3:A8Y",0.709,_val["V_A8Y"]["bind_unforced"],_atom["ab2mat1_A8Y"]["protein_iptm_best"],_atom["ab2mat1_A8Y"]["pp_binding_confidence"],"Top honest score; highest pp-binding 0.533"],
 ["ab2-WT","(parent)",0.558,_val["V_ab2WT"]["bind_unforced"],_atom["ab2WT"]["protein_iptm_best"],_atom["ab2WT"]["pp_binding_confidence"],"Baseline lead"],
 ["ab2-mat3","H3:A2Y+A4W+A8Y",0.693,None,_atom["ab2mat3_triple"]["protein_iptm_best"],_atom["ab2mat3_triple"]["pp_binding_confidence"],"Triple stack REGRESSES — epistasis"],
]
# round-1 single-point scan (top 12), [lib_id, mut, bind_forced, iptm]
AB_SINGLES = [[r["lib_id"],r["mut"],round(r["bind_conf"],3),round(r["iptm"],2)]
              for r in _mat if r["lib_id"] not in ("M00","M01")][:12]
# prior designs failure table, [their_rank, id, forced, unforced+lipid, naked, dlip, their_iptm, liab, glycosite]
_gly = {"design_spec_7":"CDR-H2 (NGS@55)","design_spec_0":"CDR-H2 ×2","design_spec_4":"CDR-H1 (NGT@30)",
        "design_spec_9":"CDR-H2 (NGT@54)","design_spec_5":"CDR-H2 (NGS@54)"}
_pord = ["design_spec_4","design_spec_9","design_spec_8","design_spec_7","design_spec_1",
         "design_spec_2","design_spec_0","design_spec_3","design_spec_6","design_spec_5"]
AB_PRIOR = []
for did in _pord:
    ext = next(k for k in _pf if k.startswith(did))
    AB_PRIOR.append([_prior[did]["final_rank"], did, round(_pf[ext]["bind"],3), round(_pu[ext]["bind"],3),
                     round(_plip[did]["bind_nolipid"],3), round(_prior[did]["their_iptm"],3),
                     _prior[did]["liability"], _gly.get(did,"—")])

# MD + MM/GBSA cross-validation (5 ns explicit-solvent, chain-split MM/GBSA)
_mdnote = {"ab2WT":"ab2 family · no glycosite","ab2mat1_A8Y":"ab2 family · no glycosite",
           "ab2mat2_A9Y":"ab2 family · no glycosite",
           "design_spec_4":"failed WL · CDR-H1 glycosite","design_spec_7":"failed WL · CDR-H2 glycosite"}
_md = json.load(open(os.path.join(ROOT,"results/antibody/md/md_mmgbsa_ranked.json")))
AB_MD = [[m["system"], round(m["dg_bind_kcal"],1), round(m["dg_std_kcal"],1),
          round(m["binder_rmsd_mean_nm"],2), round(m["contact_retention"],2),
          round(m["interface_contacts_mean"]), round(m["dg_per_contact"],3),
          _mdnote.get(m["system"],"")] for m in _md]

# ---- candidate table data ----
CANDS = [
 # cls, id, src, len, topology, seq, metric_label, metric, fold, note
 ["Stable","TZP-S1.1","S1o_short",66,"3-helix bundle (opt)","SLEELSKLLEEISKLIEEFSKLGSPGEELLKKLEEALKKLEEHLKKLGNSGSEEIRKLAEEIKKLG","bind(screen)",0.75,0.87,"Best of optimization; shortest confident bundle"],
 ["Stable","TZP-S1","LS3_bundle2",75,"3-helix bundle","SLEELSKLLEEISKLIEEWSKLGSPGSEELLKKLEEALKKLEEHLKKLGSDGTSEEIRKLAEEIKKLAEELKKLG","protein_iptm",0.96,0.96,"Atomistic-confirmed; most stable fold"],
 ["Stable","TZP-S2","LS5_c",49,"mixed α/β (β-hairpin+helix)","KTVTLTVEGKTYTVTVSDGSKLEGSPSLEELSKLLEEISKLIEEFSKLG","protein_iptm",0.97,0.96,"Has real β-sheet (33%) + helix"],
 ["Stable","TZP-S3","LS4_hairpinF",47,"α-hairpin (2 helices)","SLEELSKLLEEISKLIEEFSKLGPNGKELKKLAEELKKLAEELKKLA","protein_iptm",0.97,0.96,"Compact two-helix fold"],
 ["Stable","TZP-S4","LS2_b",34,"disulfide-stapled helix","SAEELLKKLCELSKLLEEICKLIEEWLKKLEELG","protein_iptm",0.96,0.97,"Stapled; protease-resistant"],
 ["Peptide","TZP-P3","R2_035",34,"amphipathic α-helix","SLSTLENELSTLENEISTIENEWSTGLENEISTG","bind(un-forced)",0.41,None,"Highest peptide binding signal"],
 ["Peptide","TZP-P1","R2_008",26,"amphipathic α-helix","SLSTLENELSTLENEISTIENEFSTG","bind(un-forced)",0.37,None,"Best-balanced 26-mer; F-anchor"],
 ["Peptide","TZP-P2","R2_002",26,"amphipathic α-helix","SLSTLENELSTLENEISTIENEWSTG","bind(un-forced)",0.37,None,"Cleanest; W-anchor; no liabilities"],
 ["Peptide","TZP-P4","R2_010",26,"amphipathic α-helix","SLSTLENELSTLENEISTIENEYSTG","bind(un-forced)",0.35,None,"Tyr-anchor"],
 ["Peptide","TZP-P5","R2_007",26,"amphipathic α-helix","SLSTWENELSTLENEISTIENEWSTG","bind(un-forced)",0.34,None,"di-Trp"],
 ["Peptide","TZP-P6","R2_034",21,"amphipathic α-helix","SLSTLENELSTLENEISTIEG","bind(un-forced)",0.33,None,"Minimal 21-mer; cheapest SPPS"],
 ["Existing","TZP-B1","design_spec_13",96,"mixed α/β mini-protein","MVSKTFTVTLPTGATLTVTVTYDYETNVLTVKVTLPQTGKSDEATFKVAAGVTVETILPGVGLVVKAHYDAEKNTITITLECPALGATFTFTLNLS","bind(un-forced)",0.20,None,"Best existing binder; orthogonal scaffold"],
]

IGG1=("EPKSCDKTHTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSRDELTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK")
LINK="GGGGSGGGGSGGGGS"

data_js = json.dumps(CANDS)
structs_json = json.dumps(structs_js)
ableads_js = json.dumps(AB_LEADS)
absingles_js = json.dumps(AB_SINGLES)
abprior_js = json.dumps(AB_PRIOR)
abmd_js = json.dumps(AB_MD)

html = """<!DOCTYPE html>
<html lang="en"><head>
<meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Tirzepatide Binder Report — interactive</title>
<script>
// load 3Dmol from whichever CDN the viewer's network allows
(function(){var srcs=["https://cdn.jsdelivr.net/npm/3dmol@2.1.0/build/3Dmol-min.js",
 "https://unpkg.com/3dmol@2.1.0/build/3Dmol-min.js",
 "https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.1.0/3Dmol-min.js"];
 var i=0;function nx(){if(i>=srcs.length){window._3dfail=true;return;}
  var s=document.createElement('script');s.src=srcs[i++];s.onload=function(){window._3dok=true;};s.onerror=nx;document.head.appendChild(s);} nx();})();
</script>
<style>
:root{--bg:#0f1420;--panel:#171e2e;--ink:#e8edf6;--mut:#93a1bd;--acc:#4C9AFF;--acc2:#E4572E;--ok:#54A24B;--line:#26304a}
*{box-sizing:border-box}
body{margin:0;font:15px/1.55 -apple-system,Segoe UI,Roboto,Helvetica,Arial,"PingFang SC","Microsoft YaHei",sans-serif;background:var(--bg);color:var(--ink)}
header{padding:22px 26px;border-bottom:1px solid var(--line);background:linear-gradient(180deg,#141b2b,#0f1420)}
h1{margin:0 0 4px;font-size:21px}
.sub{color:var(--mut);font-size:13px}
nav{display:flex;gap:6px;padding:10px 20px;position:sticky;top:0;background:#0f1420ee;backdrop-filter:blur(6px);border-bottom:1px solid var(--line);z-index:5;flex-wrap:wrap}
nav button{background:transparent;color:var(--mut);border:1px solid transparent;padding:7px 13px;border-radius:8px;cursor:pointer;font-size:14px}
nav button.active{color:var(--ink);background:var(--panel);border-color:var(--line)}
main{padding:20px 26px;max-width:1180px;margin:0 auto}
section{display:none} section.on{display:block}
.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:12px;margin:14px 0}
.card{background:var(--panel);border:1px solid var(--line);border-radius:12px;padding:14px}
.card .n{font-size:26px;font-weight:700}
.card .l{color:var(--mut);font-size:12px;margin-top:2px}
.badge{display:inline-block;padding:2px 8px;border-radius:20px;font-size:11px;font-weight:600}
.b-Stable{background:#1c3a2a;color:#7ee2a8}.b-Peptide{background:#28324e;color:#8fb6ff}.b-Existing{background:#3a2f1c;color:#e2c07e}
table{border-collapse:collapse;width:100%;font-size:13.5px;margin-top:8px}
th,td{padding:8px 9px;border-bottom:1px solid var(--line);text-align:left;vertical-align:top}
th{color:var(--mut);cursor:pointer;user-select:none;position:sticky;top:52px;background:#131a29}
th:hover{color:var(--ink)}
tr:hover td{background:#141c2e}
.seq{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:12px;letter-spacing:.5px;color:#bfe0ff;word-break:break-all}
.mono{font-family:ui-monospace,Menlo,monospace}
.copy{cursor:pointer;border:1px solid var(--line);background:#1d2740;color:var(--mut);border-radius:6px;padding:2px 7px;font-size:11px;margin-left:6px}
.copy:hover{color:var(--ink)}
.viewwrap{display:grid;grid-template-columns:320px 1fr;gap:16px}
@media(max-width:860px){.viewwrap{grid-template-columns:1fr}}
#viewer{position:relative;width:100%;height:560px;background:#0a0e17;border:1px solid var(--line);border-radius:12px;overflow:hidden}
select,.ctl{width:100%;background:#1d2740;color:var(--ink);border:1px solid var(--line);border-radius:8px;padding:9px;font-size:14px}
.ctlrow{display:flex;align-items:center;gap:8px;margin:8px 0;color:var(--mut);font-size:13px}
.legend span{display:inline-flex;align-items:center;gap:6px;margin-right:12px;font-size:12px}
.dot{width:12px;height:12px;border-radius:3px;display:inline-block}
.kv{display:flex;justify-content:space-between;border-bottom:1px dashed var(--line);padding:5px 0;font-size:13px}
.kv b{color:var(--acc)}
.note{color:var(--mut);font-size:12.5px}
a{color:var(--acc)} .pill{font-size:11px;color:var(--mut);border:1px solid var(--line);border-radius:20px;padding:1px 8px}
.warn{background:#2a1e14;border:1px solid #5b3b1e;border-radius:10px;padding:12px 14px;color:#f0c7a0;font-size:13px}
h2{font-size:17px;border-left:3px solid var(--acc);padding-left:9px;margin-top:26px}
code{background:#1d2740;padding:1px 5px;border-radius:4px}
</style></head>
<body>
<header>
  <h1>Tirzepatide Binder Program — Interactive Report</h1>
  <div class="sub">Boltz-2.1 evaluation &amp; de novo design · target = Tirzepatide (39-aa GIP/GLP-1 agonist, Aib2/Aib13 + K20 C20-diacid-γGlu-(AEEA)₂) · generated 2026-07-01</div>
</header>
<nav>
  <button data-t="ov" class="active">Overview</button>
  <button data-t="seq">Sequences</button>
  <button data-t="ab">Antibody &amp; Maturation</button>
  <button data-t="view">3D Docking Viewer</button>
  <button data-t="meth">Methods &amp; Caveats</button>
</nav>
<main>

<section id="ov" class="on">
  <div class="cards">
    <div class="card"><div class="n">90→</div><div class="l">existing binders re-evaluated</div></div>
    <div class="card"><div class="n">135</div><div class="l">new candidates designed &amp; screened</div></div>
    <div class="card"><div class="n">0.96–0.97</div><div class="l">stable-lead fold pTM &amp; interface iptm</div></div>
    <div class="card"><div class="n">11</div><div class="l">wet-lab-ready leads (3 classes)</div></div>
  </div>
  <h2>Key findings</h2>
  <ul>
    <li>The 5 yellow-highlighted picks were <b>exactly the top-5</b> by an independent composite; <code>spec_13</code> (TZP-B1) is the strongest existing binder on 3 independent signals.</li>
    <li><b>Consensus epitope</b> = Tirzepatide C-terminal amphipathic face <b>F22 · V23 · L26 · I27</b> (+ E3 and the PPPS tail). All new designs target it.</li>
    <li>De novo peptide <b>NB094 → affinity-matured</b> ~2× to the TZP-P series (un-forced binding_confidence 0.33–0.41).</li>
    <li>Longer <b>autonomously-folding</b> designs (TZP-S1..S4) confirmed atomistically: monomer pTM 0.96–0.97 <i>and</i> clean interface (protein_iptm 0.96–0.97) against the <b>fully-modified</b> drug (K20 lipid included).</li>
    <li>Human <b>IgG1 / IgG4-S228P Fc-fusions</b> preserve binding (interface iptm 0.88–0.95).</li>
  </ul>
  <h2>Affinity-maturation &amp; robustness</h2>
  <img src="data:image/png;base64,%FIG%" style="max-width:100%;border:1px solid var(--line);border-radius:10px">
  <p class="warn">Scores are <b>relative triage signals</b>, not affinity guarantees. Tirzepatide is a small, heavily-modified peptide. With the K20 lipid present, screen <code>binding_confidence</code> is an upper bound (binders may also touch the greasy chain); the cleaner number is the atomistic <code>protein_iptm</code>. Confirm empirically by SPR/BLI and CD.</p>
</section>

<section id="seq">
  <h2>Recommended candidate panel <span class="pill">click a column header to sort</span></h2>
  <div class="note" style="margin:6px 0">Metric = the most trustworthy score per class: <b>protein_iptm</b> (atomistic interface) for stable leads, <b>binding_confidence (un-forced)</b> for peptides/existing. Fold = monomer/structure confidence.</div>
  <table id="tbl"><thead><tr>
    <th data-k="0">Class</th><th data-k="1">ID</th><th data-k="2">Design</th><th data-k="3" data-num="1">AA</th>
    <th data-k="4">Topology</th><th data-k="6">Metric</th><th data-k="7" data-num="1">Score</th><th data-k="8" data-num="1">Fold</th>
    <th data-k="5">Sequence</th><th data-k="9">Notes</th>
  </tr></thead><tbody></tbody></table>
  <p class="note">Full sequences incl. IgG1/IgG4-Fc fusions: <code>results/wetlab_constructs.fasta</code>. Fc cassette appended to any lead as <code>[binder]-(G4S)₃-[human IgG1 hinge-CH2-CH3]</code>.</p>
</section>

<section id="ab">
  <h2>Antibody scFv — affinity maturation of lead <span class="mono">ab2</span></h2>
  <p class="note">CDR point-mutation of the anti-Tirzepatide scFv <b>ab2</b>, co-folded against the <b>fully-modified</b> drug
  (Aib2/Aib13 + K20 acyl chain), epitope-directed. 81 co-folds across 3 rounds + atomistic validation.
  Primary signal = <code>binding_confidence</code>; leads confirmed by un-forced re-score and atomistic <code>protein_iptm</code>.</p>
  <div class="cards">
    <div class="card"><div class="n">48/50</div><div class="l">single CDR mutations that improved WT</div></div>
    <div class="card"><div class="n">+0.10</div><div class="l">un-forced binding gain (WT 0.63→0.73)</div></div>
    <div class="card"><div class="n">0.895</div><div class="l">best matured protein_ipTM (A9Y, atomistic)</div></div>
    <div class="card"><div class="n">0.84 Å</div><div class="l">best interface PAE (A9Y)</div></div>
  </div>

  <h2>Recommended matured leads <span class="pill">click a header to sort</span></h2>
  <div class="note" style="margin:6px 0">Advance the two single-point leads (liability-clean, no CDR glycosite). The triple stack regresses — a real epistasis signal caught by validation. View complexes in the <b>3D Docking Viewer</b> tab (the <code>ab2-mat…·atomistic</code> entries).</div>
  <table id="tblLead" class="sortable"><thead><tr>
    <th>Lead</th><th>CDR-H3 mutation</th><th data-num="1">bind (forced)</th><th data-num="1">bind (un-forced)</th>
    <th data-num="1">protein_ipTM</th><th data-num="1">pp-binding</th><th>verdict</th>
  </tr></thead><tbody></tbody></table>

  <h2>Round 1 — single-point CDR-H3 scan (top 12)</h2>
  <table id="tblSing" class="sortable"><thead><tr>
    <th>Variant</th><th>mutation</th><th data-num="1">binding_confidence (forced)</th><th data-num="1">ipTM</th>
  </tr></thead><tbody></tbody></table>
  <p class="note">Aromatic (W/Y/F) substitutions dominate — consistent with engaging Tirzepatide's hydrophobic C-terminal face. Full narrative: <code>results/MATURATION_REPORT.md</code> · atomistic: <code>results/MATURATION_ATOMISTIC.md</code>.</p>

  <h2>Why the prior (unmodified-target) scFv failed in wet-lab</h2>
  <p class="warn"><b>Honest headline:</b> Boltz <code>binding_confidence</code> alone does <b>not</b> reproduce the failure — several priors score fine vs the modified drug. The real, sequence-level causes it doesn't model are below. <b>5/10 designs (incl. the preferred <span class="mono">design_spec_7</span>) carry an N-glycosylation sequon inside a CDR</b> → paratope glycan blocks binding in mammalian expression.</p>
  <table id="tblPrior" class="sortable"><thead><tr>
    <th data-num="1">their rank</th><th>design</th><th data-num="1">forced</th><th data-num="1">un-forced (+lipid)</th>
    <th data-num="1">naked peptide</th><th data-num="1">their ipTM (unmod)</th><th data-num="1">liability</th><th>CDR glycosite</th>
  </tr></thead><tbody></tbody></table>
  <p class="note"><b>naked peptide</b> = binding_confidence with the lipid removed; most priors collapse to ≈0, i.e. weak intrinsic peptide grip. Full analysis: <code>results/PRIOR_DESIGNS_EVALUATION.md</code>. Matured constructs (scFv/scFv-Fc/IgG1): <code>results/antibody/matured_constructs.fasta</code>.</p>

  <h2>MD + MM/GBSA cross-validation <span class="pill">5 ns explicit-solvent · OpenMM/Modal</span></h2>
  <p class="note" style="margin:6px 0">Independent physics check of the Boltz static scores: 5 ns MD per complex, then chain-split end-state <b>MM/GBSA</b> ΔG + pose stability (binder RMSD after target superposition; inter-chain contact retention). Trajectories PBC-unwrapped + target min-imaged before scoring.</p>
  <table id="tblMD" class="sortable"><thead><tr>
    <th data-num="1">MM/GBSA ΔG</th><th>system</th><th data-num="1">± </th><th data-num="1">binder RMSD (nm)</th>
    <th data-num="1">contact retention</th><th data-num="1">⟨contacts⟩</th><th data-num="1">ΔG / contact</th><th>note</th>
  </tr></thead><tbody></tbody></table>
  <p class="warn"><b>Read this correctly.</b> All 5 Boltz poses are <b>dynamically stable</b> (binder RMSD 0.26–0.50 nm, contact retention ≈0.9–1.3) — the docking modes are physically plausible. Within the ab2 scaffold, <b>A8Y improves MM/GBSA ~+9 kcal/mol over WT</b> (a 3rd independent vote for ab2-mat1). BUT raw MM/GBSA ranks the <b>failed</b> spec_7/spec_4 highest — a well-known <b>interface-size bias</b> (spec_7 buries the most atoms). On <b>ΔG-per-contact</b> the ab2 family is the most efficient, and both Boltz <i>and</i> MM/GBSA miss spec_7's CDR-glycosylation failure: affinity scores can't predict a developability failure. Use MM/GBSA only <b>within a scaffold</b>, paired with liability screening. Detail: <code>results/MD_CROSSVALIDATION.md</code>.</p>
  <div class="cards">
    <div class="card"><div class="n">0</div><div class="l">N-glyc sequons in ab2 WT/A8Y/A9Y (vs spec_7: CDR-H2)</div></div>
    <div class="card"><div class="n">5/5</div><div class="l">complexes stable over 5 ns MD</div></div>
    <div class="card"><div class="n">+9</div><div class="l">kcal/mol A8Y MM/GBSA gain vs WT (same scaffold)</div></div>
  </div>
</section>

<section id="view">
  <h2>3D docking viewer <span class="pill">drag=rotate · scroll=zoom · Boltz-2 predicted complex</span></h2>
  <div class="viewwrap">
    <div>
      <select id="pick"></select>
      <div class="ctlrow"><label><input type="checkbox" id="spin"> spin</label>
        <label><input type="checkbox" id="epi" checked> highlight epitope</label>
        <label><input type="checkbox" id="lip" checked> show K20 lipid</label></div>
      <select id="style" class="ctl">
        <option value="cartoon">Cartoon</option><option value="stick">Sticks</option><option value="surface">Binder surface</option>
      </select>
      <div class="legend" style="margin-top:10px">
        <span><i class="dot" style="background:#9aa7bd"></i>Tirzepatide (target)</span>
        <span><i class="dot" style="background:#4C9AFF"></i>Binder</span>
        <span><i class="dot" style="background:#E4572E"></i>epitope F22/V23/L26/I27</span>
        <span><i class="dot" style="background:#F2C14E"></i>K20 lipid</span>
      </div>
      <div id="mx" style="margin-top:14px"></div>
    </div>
    <div id="viewer"></div>
  </div>
</section>

<section id="meth">
  <h2>Methods</h2>
  <ul>
    <li><b>Model:</b> Boltz-2.1 (co-folding + binding module) via Boltz API; local biophysical/developability analysis.</li>
    <li><b>Target:</b> peptide chain with CCD <code>AIB</code> at 2 &amp; 13; K20 acyl (C20-diacid-γGlu-(AEEA)₂) extracted verbatim from the supplied SMILES, co-folded as ligand chain L.</li>
    <li><b>Epitope mapping:</b> heavy-atom contacts (≤4.5 Å) across top complexes → consensus F22/V23/L26/I27 + E3 + PPPS.</li>
    <li><b>Screens:</b> protein_screen (epitope-directed &amp; un-forced); structure_and_binding for atomistic protein_iptm + CIF. ~200 Boltz-2 predictions total.</li>
  </ul>
  <h2>Interpreting the scores</h2>
  <div class="kv"><span>protein_iptm</span><span>binder↔peptide interface, 0-1, &gt;0.8 confident. <b>Cleanest metric.</b></span></div>
  <div class="kv"><span>binding_confidence</span><span>Boltz binder head, 0-1. Un-forced = honest; lipid-present = upper bound.</span></div>
  <div class="kv"><span>pTM / structure_confidence</span><span>monomer autonomous-fold confidence, 0-1.</span></div>
  <div class="kv"><span>ipTM</span><span>whole-complex interface; inflated under epitope-forcing → not used as primary.</span></div>
  <h2>Exact covalent model</h2>
  <p class="note">Hosted API bonds only to CCD ligands, so the lipid is co-folded (not covalently bonded) and the C-terminal amide is omitted. For the exact Lys20-Nζ covalent model, run <code>design/tirzepatide_fully_modified.boltz.yaml</code> in a local Boltz-2 install.</p>
</section>
</main>

<script>
const CANDS=%DATA%, STRUCTS=%STRUCTS%;
const AB_LEADS=%ABLEADS%, AB_SINGLES=%ABSINGLES%, AB_PRIOR=%ABPRIOR%, AB_MD=%ABMD%;
// tabs
document.querySelectorAll('nav button').forEach(b=>b.onclick=()=>{
  document.querySelectorAll('nav button').forEach(x=>x.classList.remove('active'));
  document.querySelectorAll('section').forEach(x=>x.classList.remove('on'));
  b.classList.add('active');document.getElementById(b.dataset.t).classList.add('on');
  if(b.dataset.t==='view') setTimeout(()=>{if(!window._vinit)initViewer();else glv.resize();},60);
});
// table
const tb=document.querySelector('#tbl tbody');
function render(rows){tb.innerHTML='';rows.forEach(r=>{
  const tr=document.createElement('tr');
  const sc=(r[7]==null?'–':r[7].toFixed(2)), fo=(r[8]==null?'–':r[8].toFixed(2));
  tr.innerHTML=`<td><span class="badge b-${r[0]}">${r[0]}</span></td>
   <td class="mono"><b>${r[1]}</b></td><td class="mono note">${r[2]}</td><td>${r[3]}</td>
   <td>${r[4]}</td><td class="note">${r[6]}</td><td><b>${sc}</b></td><td>${fo}</td>
   <td><span class="seq">${r[5]}</span><button class="copy" data-s="${r[5]}">copy</button></td>
   <td class="note">${r[9]}</td>`;
  tb.appendChild(tr);});
  document.querySelectorAll('.copy').forEach(c=>c.onclick=()=>{navigator.clipboard.writeText(c.dataset.s);c.textContent='✓';setTimeout(()=>c.textContent='copy',900);});
}
let sortk=7,dir=-1;
function sortby(k){const num=document.querySelector(`th[data-k="${k}"]`).dataset.num;
  if(k===sortk)dir=-dir;else{sortk=k;dir=(k===3||k===7||k===8)?-1:1;}
  const rows=[...CANDS].sort((a,b)=>{let x=a[k],y=b[k];if(x==null)x=-1;if(y==null)y=-1;
    if(num){return (x-y)*dir;} return (''+x).localeCompare(''+y)*dir;});
  render(rows);}
document.querySelectorAll('#tbl th').forEach(th=>th.onclick=()=>sortby(+th.dataset.k));
sortby(7);

// ---- antibody / maturation tables (generic sortable) ----
function fmt(v){return v==null?'–':(typeof v==='number'?(Math.abs(v)<10?v.toFixed(3):v):v);}
function cellHTML(v,ci,kind){
  if(kind==='lead'){ if(ci===0)return `<b class="mono">${v}</b>`; if(ci===1)return `<span class="mono">${v}</span>`;
    if(ci>=2&&ci<=5)return `<b>${fmt(v)}</b>`; return `<span class="note">${v}</span>`; }
  if(kind==='sing'){ if(ci===0)return `<span class="mono">${v}</span>`; if(ci===1)return `<span class="mono note">${v}</span>`; return fmt(v); }
  if(kind==='prior'){ if(ci===1)return `<span class="mono">${v}</span>`;
    if(ci===6)return `<b style="color:${v>=180?'#E4572E':'var(--ink)'}">${v}</b>`;
    if(ci===7)return v==='—'?'<span class="note">—</span>':`<b style="color:#E4572E">${v}</b>`;
    if(ci===4)return `<span style="color:${v<0.1?'#E4572E':'var(--ink)'}">${fmt(v)}</span>`; return fmt(v); }
  if(kind==='md'){ if(ci===0)return `<b>${fmt(v)}</b>`; if(ci===1)return `<span class="mono">${v}</span>`;
    if(ci===7)return `<span class="note" style="color:${/failed/.test(v)?'#E4572E':'#7ee2a8'}">${v}</span>`;
    return fmt(v); }
  return fmt(v);
}
function mkTable(tid,data,kind,defcol){
  const tb=document.querySelector('#'+tid+' tbody');
  const ths=[...document.querySelectorAll('#'+tid+' th')];
  let sk=defcol,sd=-1;
  function paint(rows){tb.innerHTML='';rows.forEach(r=>{const tr=document.createElement('tr');
    tr.innerHTML=r.map((v,ci)=>`<td>${cellHTML(v,ci,kind)}</td>`).join('');tb.appendChild(tr);});}
  function srt(k){const num=ths[k].dataset.num; if(k===sk)sd=-sd;else{sk=k;sd=num?-1:1;}
    const rows=[...data].sort((a,b)=>{let x=a[k],y=b[k];if(x==null)x=num?-1:'';if(y==null)y=num?-1:'';
      return num?(x-y)*sd:(''+x).localeCompare(''+y)*sd;}); paint(rows);}
  ths.forEach((th,k)=>th.onclick=()=>srt(k)); srt(defcol);
}
mkTable('tblLead',AB_LEADS,'lead',4);   // sort by protein_ipTM
mkTable('tblSing',AB_SINGLES,'sing',2); // sort by forced bind
mkTable('tblPrior',AB_PRIOR,'prior',6); // sort by liability
mkTable('tblMD',AB_MD,'md',0);          // sort by MM/GBSA dG
// viewer
let glv;
function initViewer(){
  if(typeof $3Dmol==='undefined'){
    if(window._3dfail){document.getElementById('viewer').innerHTML=
      '<div style="padding:24px;color:#93a1bd;font-size:14px">3D library (3Dmol.js) could not load — this browser/network blocked the CDN.<br><br>Open this file on a machine with internet, or download 3Dmol-min.js and place a local &lt;script&gt; tag. The docked CIF files are also in <code>results/cif_stable/</code> and can be opened in PyMOL/ChimeraX.</div>';return;}
    setTimeout(initViewer,400);return; // still loading
  }
  window._vinit=true;
  glv=$3Dmol.createViewer(document.getElementById('viewer'),{backgroundColor:'#0a0e17'});
  const sel=document.getElementById('pick');
  STRUCTS.forEach((s,i)=>{const o=document.createElement('option');o.value=i;o.textContent=s.label;sel.appendChild(o);});
  sel.onchange=()=>load(+sel.value);
  ['spin','epi','lip'].forEach(id=>document.getElementById(id).onchange=()=>style());
  document.getElementById('style').onchange=style;
  load(0);
}
function cur(){return STRUCTS[+document.getElementById('pick').value];}
function load(i){
  glv.clear();
  const cif=atob(STRUCTS[i].cif);
  glv.addModel(cif,'cif');
  style(); glv.zoomTo(); glv.render();
  // metrics panel
  const lbl=STRUCTS[i].label.split(' · ')[0].replace('spec_13','TZP-B1');
  const c=CANDS.find(x=>x[1]===lbl)||CANDS.find(x=>lbl.indexOf(x[1])>=0);
  const mx=document.getElementById('mx');
  if(c){mx.innerHTML=`<div class="kv"><span>ID</span><b>${c[1]}</b></div>
    <div class="kv"><span>length</span><span>${c[3]} aa</span></div>
    <div class="kv"><span>topology</span><span>${c[4]}</span></div>
    <div class="kv"><span>${c[6]}</span><b>${c[7]==null?'–':c[7].toFixed(2)}</b></div>
    <div class="kv"><span>fold</span><span>${c[8]==null?'–':c[8].toFixed(2)}</span></div>
    <div class="note" style="margin-top:6px">${c[9]}</div>`;}
  else mx.innerHTML='';
}
function style(){
  const s=document.getElementById('style').value;
  const showEpi=document.getElementById('epi').checked;
  const showLip=document.getElementById('lip').checked;
  glv.setStyle({},{});
  // target peptide (chain T)
  glv.setStyle({chain:'T'},{cartoon:{color:'#9aa7bd'},stick:{radius:0.12,color:'#9aa7bd'}});
  // binder (chain A)
  if(s==='surface'){glv.setStyle({chain:'A'},{cartoon:{color:'#4C9AFF'}});
     glv.addSurface($3Dmol.SurfaceType.VDW,{opacity:0.75,color:'#4C9AFF'},{chain:'A'});}
  else if(s==='stick'){glv.setStyle({chain:'A'},{stick:{colorscheme:'cyanCarbon',radius:0.18}});}
  else{glv.setStyle({chain:'A'},{cartoon:{color:'#4C9AFF'}});}
  // lipid ligand (chain L)
  if(showLip) glv.setStyle({chain:'L'},{stick:{radius:0.18,color:'#F2C14E'}});
  else glv.setStyle({chain:'L'},{});
  // epitope highlight on target
  if(showEpi) glv.addStyle({chain:'T',resi:[3,18,19,22,23,26,27]},{stick:{radius:0.28,color:'#E4572E'}});
  glv.spin(document.getElementById('spin').checked?'y':false);
  glv.render();
}
</script>
</body></html>
"""

html = (html.replace("%DATA%", data_js).replace("%STRUCTS%", structs_json).replace("%FIG%", fig_b64)
            .replace("%ABLEADS%", ableads_js).replace("%ABSINGLES%", absingles_js).replace("%ABPRIOR%", abprior_js)
            .replace("%ABMD%", abmd_js))
out = os.path.join(ROOT, "report", "tzp_report.html")
with open(out, "w") as fh:
    fh.write(html)
print("wrote", out, round(os.path.getsize(out)/1024), "KB")
