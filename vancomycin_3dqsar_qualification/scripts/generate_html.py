# -*- coding: utf-8 -*-
"""Render the full regulatory report as self-contained HTML (base64 images) for PDF via Chromium."""
import json, base64
D="/tmp/claude-0/-home-user-ClaudeCode/0e1c2284-3690-5d68-a96a-480ebe908c1c/scratchpad"
R=json.load(open(f"{D}/analysis_results.json"))
MCS=json.load(open(f"{D}/mcs_core_alignment.json"))
Q=json.load(open(f"{D}/qsar_results.json"))
QA=json.load(open(f"{D}/qsar_augment.json"))
M7=json.load(open(f"{D}/m7_results.json"))
def img64(f):
    return "data:image/png;base64,"+base64.b64encode(open(f"{D}/{f}","rb").read()).decode()
EPS=[r["endpoint"] for r in Q["cv_metrics"]]
aucs=[r["CV_AUC"] for r in Q["cv_metrics"]]
meanauc=Q["mean_cv_auc"]

def fig(f,cap,w="100%"):
    return f'<figure><img src="{img64(f)}" style="width:{w}"><figcaption>{cap}</figcaption></figure>'

# QSAR table rows
qrows=""
for ep in EPS:
    m=[r for r in Q["cv_metrics"] if r["endpoint"]==ep][0]
    meas=QA["measured_vancomycin"][ep]; hp=QA["heldout_parent_pred"][ep]; p75=QA["heldout_pred_rrt075"][ep]; d=p75-hp
    mc='style="color:#c00;font-weight:700"' if meas=="active" else ""
    qrows+=f"<tr><td>{ep} <span class=sub>({m['description']})</span></td><td class=c>{m['CV_AUC']:.3f}</td><td class=c {mc}>{meas}</td><td class=c>{hp:.3f}</td><td class=c>{p75:.3f} ({'+' if d>=0 else ''}{d:.3f})</td></tr>"

# M7 table rows
m7rows=""
for c,r in M7["per_compound"].items():
    ec="#2e7d32" if not r["expert_rulebase_alerts"] else "#c00"
    sc="#2e7d32" if r["stat_Ames_call"]=="negative" else "#c00"
    m7rows+=f'<tr><td>{c}</td><td class=c style="color:{ec}">{r["expert_call"]}</td><td class=c style="color:{sc}">{r["stat_Ames_call"]} (p={r["stat_Ames_prob"]:.2f}, AD {r["stat_AD_maxTanimoto"]:.2f})</td><td class=c style="color:#2e7d32;font-weight:700">5</td></tr>'

n_alerts=M7["methodology_1"].split(";")[-1].strip().split()[0]

HTML=f"""<!doctype html><html><head><meta charset="utf-8"><style>
@page {{ size:A4; margin:16mm 15mm; }}
* {{ box-sizing:border-box; }}
body {{ font-family:'Helvetica Neue',Arial,sans-serif; font-size:10.3pt; line-height:1.5; color:#1a1a1a; margin:0; }}
h1 {{ font-size:15pt; color:#1F3864; border-bottom:2px solid #1F3864; padding-bottom:3px; margin:20px 0 8px; page-break-after:avoid; }}
h2 {{ font-size:12pt; color:#2E74B5; margin:14px 0 5px; page-break-after:avoid; }}
h3 {{ font-size:10.6pt; color:#2E74B5; margin:10px 0 4px; }}
p {{ text-align:justify; margin:6px 0; }}
.title {{ text-align:center; margin-top:40px; }}
.title h1 {{ border:0; color:#1F3864; font-size:22pt; margin:2px; }}
.title .sub {{ color:#2E74B5; font-size:13pt; margin:8px 0 24px; }}
table {{ border-collapse:collapse; width:100%; margin:8px 0; font-size:8.8pt; page-break-inside:avoid; }}
th {{ background:#1F3864; color:#fff; padding:5px 7px; text-align:left; }}
td {{ border:1px solid #cdd6e4; padding:4px 7px; vertical-align:top; }}
td.c {{ text-align:center; }}
tr:nth-child(even) td {{ background:#f4f7fb; }}
.sub {{ color:#666; font-size:0.85em; }}
figure {{ margin:10px 0; text-align:center; page-break-inside:avoid; }}
figure img {{ border:1px solid #e2e2e2; }}
figcaption {{ font-size:8.4pt; color:#555; font-style:italic; text-align:justify; margin-top:3px; }}
ul {{ margin:6px 0 6px 0; padding-left:20px; }}
li {{ margin:4px 0; text-align:justify; }}
.box {{ background:#f4f7fb; border-left:4px solid #2E74B5; padding:8px 12px; margin:8px 0; }}
.warn {{ background:#fff8e1; border-left:4px solid #e0a800; padding:8px 12px; margin:8px 0; font-size:9.4pt; }}
.pill-no {{ color:#2e7d32; font-weight:700; }}
.kv td:first-child {{ background:#eaf0f8; font-weight:600; width:26%; }}
.pb {{ page-break-before:always; }}
.foot {{ font-size:8pt; color:#888; }}
</style></head><body>

<div class="title">
<div class="foot">IN SILICO TOXICOLOGICAL QUALIFICATION — REGULATORY-SUPPORT REPORT</div>
<h1>In Silico Toxicological Qualification of Two Vancomycin&nbsp;B Impurities by 3D-QSAR / Structure-Based Read-Across</h1>
<div class="sub">Impurities RRT&nbsp;0.75 [(26R)-DMe-DeLI] and RRT&nbsp;0.87 [(7S,26R)-Vancomycin&nbsp;B]</div>
<table class="kv" style="width:88%;margin:0 auto;text-align:left">
<tr><td>Test items</td><td>Vancomycin B (API) and impurities RRT 0.75 and RRT 0.87</td></tr>
<tr><td>In silico methods</td><td>Toxicophore/pharmacophore profiling; RDKit PAINS/Brenk/NIH alerts; physicochemical &amp; 3D-shape descriptors; 3D conformer generation &amp; MCS/O3A shape read-across; in-house 12-endpoint Tox21 statistical QSAR (RF, 5-fold CV, AD); ICH M7 two-methodology mutagenicity screen (Benigni–Bossa expert rulebase + Ames QSAR); supporting logD/pKa QSAR</td></tr>
<tr><td>Toxicological anchor</td><td>Parent API 90-day i.v. dog study, NOAEL 75/100 mg/kg/day (~2–3× clinical dose)</td></tr>
<tr><td>Date</td><td>2026-07-30</td></tr>
</table>
<p class="warn" style="width:88%;margin:18px auto;text-align:left"><b>IMPORTANT</b> — This report is a computational screening and structure-based read-across assessment intended to support expert toxicological review. It is not, by itself, a regulatory toxicology conclusion or a substitute for expert sign-off or, where required, experimental testing (see §8).</p>
</div>

<h1 class="pb">1.&nbsp; Executive summary</h1>
<p>Two related substances of the glycopeptide antibiotic vancomycin&nbsp;B require in silico toxicological qualification: <b>RRT 0.87 = (7S,26R)-vancomycin B</b>, a heat/low-water-activity degradation product, and <b>RRT 0.75 = (26R)-DMe-DeLI</b>, a fermentation-derived residue-1 variant. Because the differences are largely stereochemical/constitutional and not resolvable by conventional 2D read-across, a three-dimensional, structure-based approach was applied to the parent and both daughters.</p>
<p><b>Principal conclusion. Neither impurity introduces a new toxicophore or a new DNA-reactive (mutagenic) structural alert absent from the parent API.</b> RRT 0.87 is a pure stereoisomer sharing all 101 heavy atoms and the complete functional-group inventory of the parent (differing only in configuration at 2 of 18 stereocentres). RRT 0.75 shares 99 of ~100 heavy atoms; its only constitutional change is a secondary (N-methyl) → primary aliphatic amine at residue 1 — a group class already present in the parent (vancosamine) — with unchanged lipophilicity (logD 0.42) and near-identical basicity (pKa ≈ 8.1 vs 8.2).</p>
<p>A statistical 12-endpoint Tox21 QSAR (mean 5-fold CV-AUC {meanauc:.3f}) reinforces this: the parent is itself a Tox21 training compound (measured profile essentially clean — only a weak NR-AR signal), RRT 0.87 is fingerprint-identical to the parent (Δ ≡ 0 at all endpoints), and RRT 0.75 differs from the parent by ≤ 0.02 at every endpoint. An ICH M7 two-methodology mutagenicity screen (Benigni–Bossa expert rulebase + Ames QSAR, AUC {M7['methodology_2']['cv_auc']:.3f}) returns a concordant <b>negative / Class 5</b> for all three. On a weight-of-evidence basis both impurities qualify for read-across to the parent, whose safety is anchored by the 91-day i.v. dog study (NOAEL 75/100 mg/kg/day) and extensive clinical use, provided each stays within the Ph. Eur. related-substances limits.</p>
<table>
<tr><th>Question</th><th>RRT 0.87</th><th>RRT 0.75</th><th>Basis</th></tr>
<tr><td>New toxicophore vs parent?</td><td class=c class="pill-no" style="color:#2e7d32;font-weight:700">NO</td><td class=c style="color:#2e7d32;font-weight:700">NO</td><td>Identical FG inventory (0.87); iso-functional amine swap (0.75)</td></tr>
<tr><td>New mutagenic (M7) alert?</td><td class=c style="color:#2e7d32;font-weight:700">NO</td><td class=c style="color:#2e7d32;font-weight:700">NO</td><td>Both M7 methodologies negative → Class 5 (§5.7)</td></tr>
<tr><td>Read-across to parent justified?</td><td class=c style="color:#2e7d32;font-weight:700">YES</td><td class=c style="color:#2e7d32;font-weight:700">YES</td><td>≥98% shared skeleton; conserved core toxicophores</td></tr>
</table>

<h1 class="pb">2.&nbsp; Purpose, scope and regulatory framework</h1>
<p>The objective is a toxicological in silico qualification of impurities RRT 0.75 and RRT 0.87, delivering: (i) a 3D-QSAR / structure-based assessment of the parent and each daughter; (ii) identification/assessment of any additional toxophore present in a daughter but absent from the parent; (iii) a documented read-across justification; and (iv) this regulatory-style report.</p>
<p><b>Framework.</b> Non-mutagenic impurity qualification follows ICH Q3A(R2)/Q3B(R2). Mutagenic potential is screened per ICH M7(R2), which anticipates two complementary (Q)SAR methodologies (expert rule-based + statistical) plus expert review. In silico model use is framed by the OECD (Q)SAR validation principles. Because vancomycin B is a large (MW ≈ 1449) flexible glycopeptide outside the applicability domain of many small-molecule tox models, the assessment is weighted toward mechanistic, structure- and shape-based read-across, with statistical models used where their domain is respected (and flagged where not).</p>

<h1 class="pb">3.&nbsp; Test items and structural relationships</h1>
<p>Structures were taken from the sponsor information (IUPAC/SMILES/InChI) and independently validated with RDKit; computed formulae and exact masses reproduce the source document.</p>
<table>
<tr><th>Attribute</th><th>Parent API</th><th>RRT 0.87</th><th>RRT 0.75</th></tr>
<tr><td>Xellia name</td><td>Vancomycin B</td><td>(7S,26R)-Vancomycin B</td><td>Vancomycin B impurity (26R)-DMe-DeLI</td></tr>
<tr><td>Molecular formula</td><td>C66H75Cl2N9O24</td><td>C66H75Cl2N9O24</td><td>C65H73Cl2N9O24</td></tr>
<tr><td>Monoisotopic MW</td><td>1447.43</td><td>1447.43</td><td>1433.41</td></tr>
<tr><td>Defined stereocentres</td><td>18</td><td>18</td><td>19</td></tr>
<tr><td>Origin</td><td>Fermentation (API)</td><td>Degradation (low water activity + heat)</td><td>Fermentation</td></tr>
<tr><td>Relationship to API</td><td>—</td><td>Diastereomer (C7 &amp; C26 epimer)</td><td>Residue-1 variant (N-desMe; Leu→Ile) + C26 epimer</td></tr>
</table>
{fig("fig1_structures.png","Figure 1. 2D structures of the parent API (left) and the two impurities. RRT 0.87 is drawn identically to the parent because it differs only in stereochemistry; RRT 0.75 differs at N-terminal residue 1 (free primary amine + isoleucine side chain instead of N-methyl-leucine).")}
<p><b>RRT 0.87 — (7S,26R)-vancomycin B.</b> Same molecular formula as the API; configuration inverted at C7 (benzylic carbinol of residue 3) and C26 (a peptide α-carbon) — consistent with thermally-driven epimerisation under low water activity. No atoms added or removed.</p>
<p><b>RRT 0.75 — (26R)-DMe-DeLI.</b> Formula C65H73Cl2N9O24 (one CH2 less). Residue 1 changes from N-methyl-D-leucine to des-N-methyl-isoleucine (secondary → primary amine; isobutyl → sec-butyl), with epimerisation at C26. All other residues, both sugars, both aryl chlorides and all phenols are retained.</p>
<div style="display:flex;gap:10px">{fig("fig2_parent_residue1.png","Parent: N-methyl-D-leucine (secondary amine).","48%")}{fig("fig2_rrt075_residue1.png","RRT 0.75: des-methyl isoleucine (primary amine).","48%")}</div>
<figcaption style="text-align:center">Figure 2. The only site of constitutional change between RRT 0.75 and the API.</figcaption>

<h1 class="pb">4.&nbsp; Computational methods</h1>
<h2>4.1 Toxicophore / structural-alert profiling</h2>
<p>Each structure was profiled for the vancomycin pharmacophore (heptapeptide aglycone D-Ala-D-Ala cleft, H-bonding network, vancosamine–glucose disaccharide) and the recognised toxicophores (phenolic/aromatic rings → renal accumulation; basic amine centres → histamine/Red-Man; high MW/hydrophilicity). A quantitative functional-group inventory (SMARTS) was differenced against the parent. Liabilities were screened with RDKit PAINS(A/B/C), Brenk and NIH catalogues and reviewed against ICH M7 DNA-reactive alert classes.</p>
<h2>4.2 Physicochemical &amp; 3D descriptors; 3D shape read-across</h2>
<p>RDKit descriptors (MW, cLogP, TPSA, HBD/HBA, rot. bonds, QED, Lipinski/Veber); Inductive Bio logD/pKa QSAR. ETKDGv3 conformers, MMFF-minimised; 3D shape descriptors (RoG, asphericity, spherocity, eccentricity, NPR1/2). 3D read-across quantified by MCS core-restricted alignment (RMSD) and Crippen Open3DAlign / shape-Tanimoto. A formal CoMFA/CoMSIA 3D-QSAR is not built (see §8).</p>
<h2>4.3 Statistical QSAR — Tox21 and Ames</h2>
<p>A 12-endpoint Tox21 QSAR (class-balanced RandomForest on 2048-bit Morgan fingerprints; public MoleculeNet Tox21; 5-fold OOF AUC) and an Ames QSAR (Hansen benchmark N≈6500) were trained and applied with applicability-domain checks. Parent–impurity differences are the interpretable read-out.</p>

<h1 class="pb">5.&nbsp; Results</h1>
<h2>5.1 Parent (Vancomycin B) — baseline toxicophore profile</h2>
<p>The parent presents the expected glycopeptide profile: 6 phenolic OH, 5 aromatic rings, 2 aryl chlorides, 3 aryl ethers, 2 benzylic alcohols, 1 carboxylic acid, 1 primary amide, 6 secondary amides, 1 primary amine (vancosamine), 1 secondary N-methyl amine. It carries no PAINS/Brenk alert and no DNA-reactive alert; the single NIH flag (≥7 aliphatic OH) is a polarity descriptor. This defines the read-across baseline.</p>
{fig("fig3_toxicophore_inventory.png","Figure 3. Toxicophore / functional-group inventory (counts). Parent and RRT 0.87 are identical across every category; RRT 0.75 differs only in the amine columns (red) — one secondary (N-methyl) amine replaced by one primary amine.")}
<h2>5.2 RRT 0.87 — (7S,26R)-vancomycin B</h2>
<p>Differing from the API only in configuration, its 2D structure, functional-group inventory and alert profile are identical to the parent (ECFP4/MACCS Tanimoto = 1.00; 101/101 shared heavy atoms; 0 alerts). <b>No additional toxicophore.</b> Conventional 2D read-across cannot distinguish it from the API — which is why a 3D/mechanistic assessment is required. 3D shape descriptors are near-congruent (RoG 7.45 vs 7.50 Å; NPR1 0.529 vs 0.528). Stereochemical inversion can modulate target binding (efficacy) but creates no new reactive/DNA-reactive functionality.</p>
<h2>5.3 RRT 0.75 — (26R)-DMe-DeLI</h2>
<p>Retains the entire aglycone, both sugars, both aryl chlorides, all phenols, the carboxyl and the primary amide; 99/~100 heavy atoms shared (ECFP4 0.86 / MACCS 0.94). The single change (secondary N-methyl → primary amine) is <b>not a new functional-group class</b> (parent already carries a primary amine) and is <b>not a DNA-reactive alert</b>. logD unchanged (0.42 vs 0.42), most-basic pKa near-identical (8.08 vs 8.21) — the basic-amine/Red-Man liability is conserved, not amplified. TPSA +14 Å²; 3D shape close to parent (RoG 7.46 Å).</p>
<h2>5.4 Differential toxicophore analysis — is any new toxophore present?</h2>
<table>
<tr><th>Potential toxophore</th><th>RRT 0.87 vs parent</th><th>RRT 0.75 vs parent</th></tr>
<tr><td>Phenolic / aromatic rings</td><td>Identical (6 phenol, 5 aryl)</td><td>Identical (6 phenol, 5 aryl)</td></tr>
<tr><td>Basic amine centres (histamine/Red-Man)</td><td>Identical</td><td style="background:#fff3e0">Conserved: sec→prim amine, iso-basic (pKa 8.1); not amplified</td></tr>
<tr><td>Aryl chlorides</td><td>Identical (2)</td><td>Identical (2)</td></tr>
<tr><td>DNA-reactive alert (arom. amine, nitro, epoxide, Michael acceptor, alkyl halide…)</td><td>None (as parent)</td><td>None (as parent)</td></tr>
<tr><td>PAINS / Brenk liability</td><td>None</td><td>None</td></tr>
<tr><td><b>NET: new toxophore introduced?</b></td><td class=c style="color:#2e7d32;font-weight:700">NO</td><td class=c style="color:#2e7d32;font-weight:700">NO</td></tr>
</table>
<p>The only element warranting comment is the residue-1 amine in RRT 0.75 — a modification of an existing toxicophore class (secondary→primary amine of equal basicity), not a new one. No additional toxophore is identified in either daughter.</p>
<h2>5.5 3D-QSAR / shape read-across metrics</h2>
<div style="display:flex;gap:10px;align-items:flex-start">
<div style="flex:1">{fig("fig5_similarity.png","Figure 4. 2D similarity &amp; shared-skeleton read-across. RRT 0.87 is fingerprint-identical; RRT 0.75 shares ≥98% of the heavy-atom skeleton.")}</div>
</div>
{fig("fig4_descriptors.png","Figure 5. 3D whole-molecule shape descriptors (left) and physicochemical descriptors (right, symlog). The three occupy essentially the same shape/property space.")}

<h2 class="pb">5.6 Statistical (2D) QSAR — 12-endpoint Tox21 profiling</h2>
<p>A 12-endpoint Tox21 QSAR (class-balanced RandomForest on 2048-bit Morgan fingerprints; public MoleculeNet Tox21, n={Q['model'].split('n=')[1].split(')')[0]}) was trained: mean 5-fold OOF ROC-AUC = {meanauc:.3f} (range {min(aucs):.3f}–{max(aucs):.3f}), then each molecule scored with an applicability-domain check.</p>
<div class="box"><b>Applicability domain — a favourable finding.</b> Vancomycin is itself in the Tox21 collection (training compound {', '.join(QA['identical_mol_ids'])}; nearest-neighbour Tanimoto = 1.00 to the parent and to RRT 0.87). The QSAR is therefore in-domain for the parent; RRT 0.75 sits close by (max Tanimoto 0.863).</div>
<p><b>Measured vancomycin profile (ground truth):</b> active only at NR-AR (weak), inactive at every other tested endpoint (incl. genotoxicity SR-ATAD5, SR-p53). This is the toxicological baseline.</p>
{fig("fig6_qsar.png","Figure 6. In-house 12-endpoint Tox21 QSAR. (A) Predicted probabilities — parent and RRT 0.87 rows identical (identical fingerprints). (B) 5-fold OOF AUC per endpoint (mean 0.827). (C) RRT 0.75 − held-out parent Δ — every endpoint within ±0.02; RRT 0.87 Δ ≡ 0 by construction.")}
<p><b>Impurity vs parent.</b> RRT 0.87 returns identical probabilities at all 12 endpoints (Δ ≡ 0): a concrete demonstration that a 2D-QSAR is structurally blind to this stereoisomer (hence the 3D read-across carries the weight). For RRT 0.75, every endpoint lies within ±0.02 of the held-out parent (max |Δ| = 0.02), none crossing the activity threshold — including the genotoxicity-associated SR-ATAD5 and SR-p53.</p>
<table>
<tr><th>Endpoint</th><th>CV-AUC</th><th>Measured (parent)</th><th>Held-out parent P</th><th>RRT 0.75 P (Δ)</th></tr>
{qrows}
</table>
<p class="foot"><i>All predicted probabilities are below the 0.5 threshold (parent NR-AR is a measured weak active the model under-calls — expected at that endpoint's lower AUC; the measured label governs). RRT 0.87 omitted from the table: identical fingerprint ⇒ identical predictions (Δ ≡ 0).</i></p>

<h2 class="pb">5.7 ICH M7 mutagenicity assessment — two complementary (Q)SAR methodologies</h2>
<p>ICH M7(R2) requires potential bacterial (Ames) mutagenicity to be assessed by two complementary methodologies — one <b>expert rule-based</b> and one <b>statistical</b> — with any positive resolved by expert review. Both were applied:</p>
<ul>
<li><b>Methodology 1 (expert rule-based).</b> A Benigni–Bossa structural-alert rulebase for mutagenicity/genotoxic carcinogenicity ({n_alerts} alert classes: aromatic amines, nitro/nitroso, azo/azoxy, N-nitroso, epoxides/aziridines, Michael acceptors, acyl/alkyl halides, hydrazines, aldehydes, quinones, β-lactams, …) — the public scientific rulebase underlying Derek-type systems (Toxtree/VEGA).</li>
<li><b>Methodology 2 (statistical).</b> An Ames QSAR (class-balanced RandomForest, 2048-bit Morgan) trained on the Hansen benchmark (N = {M7['methodology_2']['n_train']}; {round(M7['methodology_2']['pos_frac']*100)}% mutagens); 5-fold OOF AUC = {M7['methodology_2']['cv_auc']:.3f}; decision threshold {M7['methodology_2']['threshold_youden']:.2f}.</li>
</ul>
{fig("fig7_m7.png","Figure 7. ICH M7 two-methodology assessment. (A) Concordance matrix — both methodologies NEGATIVE for all three molecules. (B) Statistical Ames probabilities, all below threshold; amber note flags the low applicability-domain similarity, which is why the domain-independent expert rulebase is the primary M7 call.")}
<table>
<tr><th>Molecule</th><th>Expert rule-based</th><th>Statistical Ames QSAR</th><th>M7</th></tr>
{m7rows}
</table>
<p><b>Outcome.</b> The two methodologies are concordant and negative for the parent and both impurities: no mutagenicity structural alert, and the statistical Ames prediction below threshold in every case. Under ICH M7 this corresponds to <b>Class 5 (no structural alert)</b> — treated as non-mutagenic and controlled as ordinary ICH Q3A related substances, with no dedicated Ames study required. For transparency: (i) RRT 0.87 is fingerprint-identical to the parent, so both methodologies return identical results by construction; (ii) these glycopeptides lie outside the Ames model's applicability domain (nearest-neighbour Tanimoto ≈ 0.22), so the statistical negative is low-confidence and the domain-independent expert rulebase — unambiguously alert-free — is the primary determinant, corroborated by the parent's experimentally inactive Tox21 DNA-damage endpoints (§5.6). A licensed Derek/Sarah Nexus run is expected to reproduce this alert-free, negative call.</p>

<h1 class="pb">6.&nbsp; Read-across justification</h1>
<p><b>Source (analogue):</b> the parent API, vancomycin B (full toxicological package incl. 91-day i.v. dog study, NOAEL 75/100 mg/kg/day, plus decades of clinical exposure). <b>Targets:</b> RRT 0.87 and RRT 0.75.</p>
<ul>
<li><b>Structural similarity.</b> RRT 0.87 shares 100% of the heavy-atom skeleton and functional groups (configurational isomer). RRT 0.75 shares ≥98% (99/100), with a single localised, non-alerting modification at a peptide N-terminus.</li>
<li><b>Common mechanism / conserved toxicophores.</b> All phenolic, aromatic, aryl-chloride, amide, carboxyl and basic-amine features are conserved. No new reactive centre.</li>
<li><b>Metabolic / physicochemical congruence.</b> Vancomycin is renally cleared largely unchanged; near-identical size, polarity, logD (0.42) and pKa (≈8.1) predict the same ADME/disposition.</li>
<li><b>No mutagenic concern.</b> Both ICH M7 methodologies negative (§5.7); QSAR shows no activation at SR-ATAD5/SR-p53 (Δ ≤ 0.02).</li>
<li><b>Acknowledged uncertainty (stereochemistry).</b> Stereoisomers can differ in target binding — bearing on potency/efficacy, not on de novo toxicity; the 3D analysis specifically addresses this dimension.</li>
</ul>

<h1>7.&nbsp; Qualification conclusion</h1>
<p>On a weight-of-evidence basis the in silico assessment supports toxicological qualification of both impurities by read-across to the parent:</p>
<ul>
<li><b>RRT 0.87 [(7S,26R)-vancomycin B]</b> — a pure stereoisomer with no new toxicophore and no new mutagenic alert; read-across fully justified.</li>
<li><b>RRT 0.75 [(26R)-DMe-DeLI]</b> — a single, non-alerting, iso-functional residue-1 modification (plus C26 epimerisation) with no new toxicophore and no new mutagenic alert; read-across justified.</li>
</ul>
<p>Provided both remain within the Ph. Eur. related-substances limits, the parent's NOAEL (75/100 mg/kg/day, ≈2–3× clinical) and clinical record provide adequate toxicological cover. Formal sign-off should be issued by a qualified toxicologist, with the mutagenic dimension documented under the company's ICH M7 two-(Q)SAR-methodology workflow.</p>

<h1 class="pb">8.&nbsp; Applicability domain, honesty statement and limitations</h1>
<ul>
<li><b>Nature of the "3D-QSAR".</b> No congeneric training set with measured toxicity endpoints exists for these compounds, so a formal field-based 3D-QSAR (CoMFA/CoMSIA) predictive model was not — and cannot legitimately be — built. What is reported is 3D-structure/shape read-across plus mechanistic toxicophore profiling.</li>
<li><b>Statistical QSAR — scope of validity.</b> The 12-endpoint Tox21 QSAR is in-domain for the parent (vancomycin is a training compound). Two cautions: (i) the parent's in-set predictions reflect partial memorisation and are reported instead as measured labels + held-out predictions; (ii) the Tox21 panel covers 12 specific mechanisms only — it does not model vancomycin's clinical effects (nephro-/ototoxicity, Red-Man), which are handled mechanistically and by read-across to the parent's in vivo package.</li>
<li><b>ICH M7 statistical model domain.</b> The Ames QSAR is out-of-domain for these glycopeptides (Tanimoto ≈ 0.22); its negative is low-confidence and the expert rulebase is the primary call. Licensed Derek/Sarah Nexus should confirm in the validated workflow.</li>
<li><b>Conformational sampling.</b> 3D descriptors/overlaps derive from single low-energy conformers of a large flexible macrocycle → corroborating, not decisive; the decisive evidence is the shared skeleton and conserved toxicophore inventory.</li>
<li><b>Scope.</b> A hazard-identification / structure-based screening and read-across exercise that supports — does not replace — expert toxicological review and any experimental testing a competent authority may require.</li>
<li><b>Impurity control limits.</b> Ph. Eur. numerical limits for RRT 0.75/0.87 were not provided; the conclusion assumes control within the applicable monograph limits.</li>
</ul>

<h1 class="pb">9.&nbsp; References and tools</h1>
<ul class="foot" style="font-size:9pt">
<li>ICH Q3A(R2) / Q3B(R2) — Impurities in New Drug Substances/Products.</li>
<li>ICH M7(R2) — Assessment and Control of DNA-Reactive (Mutagenic) Impurities.</li>
<li>OECD (2007) — Guidance on the Validation of (Q)SAR Models.</li>
<li>Benigni R., Bossa C. (2008) — Structural alerts of mutagens and carcinogens (Toxtree/VEGA rulebase; Methodology 1).</li>
<li>Hansen K. et al. (2009) J. Chem. Inf. Model. 49:2077 — Benchmark data set for Ames mutagenicity (Methodology 2 training data).</li>
<li>Huang R. et al. — Tox21 / MoleculeNet (Wu Z. et al. 2018, Chem. Sci.) — statistical QSAR training data.</li>
<li>RDKit (2024) — descriptors, PAINS/Brenk/NIH catalogs, ETKDGv3, MMFF94, Open3DAlign, MCS. Inductive Bio logD/pKa QSAR. Ph. Eur. Vancomycin monograph.</li>
<li>Sponsor source: "Vanco_RRT_075_and_087_Info" (Xellia) — identities, SMILES/InChI, origin, 91-day dog NOAEL 75/100 mg/kg/day.</li>
</ul>
<p class="foot" style="margin-top:20px;border-top:1px solid #ccc;padding-top:6px">Computational screening &amp; structure-based read-across — supports, does not replace, expert toxicological review.</p>
</body></html>"""
open(f"{D}/report.html","w").write(HTML)
print("wrote report.html", len(HTML), "chars")
