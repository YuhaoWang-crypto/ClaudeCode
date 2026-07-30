const fs = require("fs");
const {
  Document, Packer, Paragraph, TextRun, HeadingLevel, AlignmentType,
  Table, TableRow, TableCell, WidthType, BorderStyle, ShadingType,
  ImageRun, PageBreak, Header, Footer, PageNumber, LevelFormat
} = require("docx");

const DIR = "/tmp/claude-0/-home-user-ClaudeCode/0e1c2284-3690-5d68-a96a-480ebe908c1c/scratchpad";
const NAVY = "1F3864", BLUE = "2E74B5", GREY = "595959", LIGHT = "DCE6F1", AMBER = "FBE4D5";

function img(file, w, h) {
  return new Paragraph({ alignment: AlignmentType.CENTER, spacing:{before:120,after:60},
    children:[ new ImageRun({ type:"png", data: fs.readFileSync(`${DIR}/${file}`),
      transformation:{ width:w, height:h } }) ]});
}
function cap(t){ return new Paragraph({ alignment:AlignmentType.CENTER, spacing:{after:160},
  children:[ new TextRun({ text:t, italics:true, size:16, color:GREY }) ]}); }
function h1(t){ return new Paragraph({ heading:HeadingLevel.HEADING_1, spacing:{before:280,after:120}, children:[ new TextRun({text:t, bold:true, color:NAVY}) ]}); }
function h2(t){ return new Paragraph({ heading:HeadingLevel.HEADING_2, spacing:{before:200,after:100}, children:[ new TextRun({text:t, bold:true, color:BLUE}) ]}); }
function p(runs, opts={}){ if(typeof runs==="string") runs=[new TextRun({text:runs, size:20})];
  return new Paragraph({ spacing:{after:120, line:276}, alignment: opts.align||AlignmentType.JUSTIFIED, children:runs }); }
function t(text, o={}){ return new TextRun({ text, size:o.size||20, bold:o.bold, italics:o.i, color:o.color, break:o.break }); }
function bullet(text, level=0){ return new Paragraph({ bullet:{level}, spacing:{after:80}, children: Array.isArray(text)?text:[new TextRun({text,size:20})] }); }

// ---- table helper ----
function cell(content, o={}) {
  const kids = Array.isArray(content)? content : [ new Paragraph({ children:[ new TextRun({ text:String(content), size:o.size||18, bold:o.bold, color:o.color }) ], alignment:o.align||AlignmentType.LEFT, spacing:{before:40,after:40} }) ];
  return new TableCell({ width:{size:o.w, type:WidthType.DXA}, shading:o.shade?{type:ShadingType.CLEAR, fill:o.shade, color:"auto"}:undefined,
    margins:{top:40,bottom:40,left:80,right:80}, verticalAlign:"center", children:kids });
}
function hdrCell(text,w){ return cell(text,{w,bold:true,color:"FFFFFF",shade:NAVY,align:AlignmentType.CENTER}); }
function table(colWidths, rows){ return new Table({ columnWidths:colWidths, width:{size:colWidths.reduce((a,b)=>a+b,0),type:WidthType.DXA},
  rows }); }

// identity data
const IDrows = [
  ["Attribute","Parent API","RRT 0.87","RRT 0.75"],
  ["Xellia name","Vancomycin B","(7S,26R)-Vancomycin B","Vancomycin B impurity (26R)-DMe-DeLI"],
  ["Molecular formula","C66H75Cl2N9O24","C66H75Cl2N9O24","C65H73Cl2N9O24"],
  ["Monoisotopic MW","1447.43","1447.43","1433.41"],
  ["[M+H]+ (source doc)","1448.4","1448.3","1435.23"],
  ["Avg MW (calc.)","1449.27","1449.27","1435.24"],
  ["Defined stereocentres","18","18","19"],
  ["Origin","Fermentation (API)","Degradation (low water\nactivity + heat; stable on\nstability)","Fermentation (stable on\nstability)"],
  ["Relationship to API","—","Diastereomer\n(epimerised at C7 & C26)","Residue-1 variant\n(N-desMe; Leu→Ile) +\nC26 epimer"],
];

// physchem table data
const PCrows = [
  ["Descriptor","Parent","RRT 0.87","RRT 0.75","Interpretation"],
  ["cLogP (Crippen)","0.11","0.11","-0.15","Unchanged / lower"],
  ["logD7.4 (Inductive Bio)","0.42","n/c*","0.42","Identical - low lipophilicity"],
  ["Most-basic pKa (pred.)","8.21","~8.2*","8.08","Iso-basic (moderate)"],
  ["TPSA (Å²)","530.5","530.5","544.5","+14 (extra NH)"],
  ["H-bond donors / acceptors","19 / 25","19 / 25","19 / 25","Unchanged"],
  ["Rotatable bonds","13","13","12","Unchanged"],
  ["Aromatic rings","5","5","5","Unchanged"],
  ["PAINS / Brenk alerts","0 / 0","0 / 0","0 / 0","No reactive/interference alert"],
  ["NIH alert","≥7 aliph-OH","≥7 aliph-OH","≥7 aliph-OH","Polarity flag only (non-genotoxic)"],
  ["ECFP4 / MACCS vs parent","1.00 / 1.00","1.00 / 1.00","0.86 / 0.94","Very high 2D similarity"],
  ["Shared heavy atoms /101","101","101","99","≥98% skeleton common"],
];

const doc = new Document({
  creator:"In silico toxicology (Claude Code)",
  title:"3D-QSAR / read-across toxicological qualification of vancomycin impurities RRT 0.75 and RRT 0.87",
  styles:{ default:{ document:{ run:{ font:"Calibri", size:20 } } } },
  numbering:{ config:[] },
  sections:[{
    properties:{ page:{ size:{ width:12240, height:15840 }, margin:{ top:1200, bottom:1200, left:1200, right:1200 } } },
    headers:{ default: new Header({ children:[ new Paragraph({ alignment:AlignmentType.RIGHT, children:[ new TextRun({ text:"In silico Toxicological Qualification — Vancomycin B RRT 0.75 & RRT 0.87", size:14, color:GREY }) ]}) ]}) },
    footers:{ default: new Footer({ children:[ new Paragraph({ alignment:AlignmentType.CENTER, children:[ new TextRun({ text:"Computational screening & structure-based read-across — supports, does not replace, expert toxicological review.  Page ", size:14, color:GREY }), new TextRun({ children:[PageNumber.CURRENT], size:14, color:GREY }) ]}) ]}) },
    children:[
      // ===== TITLE =====
      new Paragraph({ spacing:{before:1200,after:120}, alignment:AlignmentType.CENTER, children:[ new TextRun({ text:"In Silico Toxicological Qualification of Two", bold:true, size:36, color:NAVY }) ]}),
      new Paragraph({ alignment:AlignmentType.CENTER, spacing:{after:120}, children:[ new TextRun({ text:"Vancomycin B Impurities by 3D-QSAR / Structure-Based Read-Across", bold:true, size:36, color:NAVY }) ]}),
      new Paragraph({ alignment:AlignmentType.CENTER, spacing:{after:400}, children:[ new TextRun({ text:"Impurities RRT 0.75  [(26R)-DMe-DeLI]  and  RRT 0.87  [(7S,26R)-Vancomycin B]", size:24, color:BLUE }) ]}),
      table([3000,6600],[
        new TableRow({children:[cell("Document type",{w:3000,bold:true,shade:LIGHT}),cell("In silico toxicological qualification report (research / regulatory-support draft)",{w:6600})]}),
        new TableRow({children:[cell("Test items",{w:3000,bold:true,shade:LIGHT}),cell("Vancomycin B (API) and impurities RRT 0.75 and RRT 0.87",{w:6600})]}),
        new TableRow({children:[cell("In silico methods",{w:3000,bold:true,shade:LIGHT}),cell("Toxicophore/pharmacophore profiling; RDKit structural-alert catalogues (PAINS A/B/C, Brenk, NIH); physicochemical & 3D-shape descriptors; 3D conformer generation & MCS/O3A shape read-across; ICH M7-style expert alert review; supporting (Q)SAR pKa/logD",{w:6600})]}),
        new TableRow({children:[cell("Toxicological anchor",{w:3000,bold:true,shade:LIGHT}),cell("Parent API 90-day i.v. dog study, NOAEL 75/100 mg/kg/day (~2–3× clinical dose)",{w:6600})]}),
        new TableRow({children:[cell("Date",{w:3000,bold:true,shade:LIGHT}),cell("2026-07-30",{w:6600})]}),
      ]),
      new Paragraph({ spacing:{before:300}, children:[ new TextRun({ text:"IMPORTANT — This report is a computational screening and structure-based read-across assessment intended to support expert toxicological review. It is not, by itself, a regulatory toxicology conclusion or a substitute for expert sign-off or, where required, experimental testing. See Section 9 (Limitations).", italics:true, size:18, color:"C00000" }) ]}),
      new Paragraph({children:[new PageBreak()]}),

      // ===== 1. EXECUTIVE SUMMARY =====
      h1("1.  Executive summary"),
      p([t("Two related substances of the glycopeptide antibiotic vancomycin B require in silico toxicological qualification: ",{}),
         t("RRT 0.87 = (7S,26R)-vancomycin B",{bold:true}), t(", a heat/low-water-activity degradation product, and ",{}),
         t("RRT 0.75 = (26R)-DMe-DeLI",{bold:true}), t(", a fermentation-derived residue-1 variant. Because the sponsor correctly notes that the differences are largely stereochemical/constitutional and therefore not resolvable by conventional 2D read-across, a three-dimensional, structure-based approach was applied to each of the parent and the two daughter molecules.",{})]),
      p([t("Principal conclusion. ",{bold:true}),
         t("Neither impurity introduces a new toxicophore or a new DNA-reactive (mutagenic) structural alert that is absent from the parent API. ",{bold:true}),
         t("RRT 0.87 is a pure stereoisomer sharing all 101 heavy atoms and the complete functional-group inventory of the parent; its only difference from vancomycin B is the spatial configuration at 2 of 18 stereocentres. RRT 0.75 shares 99 of ~100 heavy atoms with the parent; its sole structural change is at N-terminal residue 1, where a secondary (N-methyl) amine is replaced by a primary aliphatic amine — a functional-group class already present elsewhere in the parent (vancosamine) — together with epimerisation at C26. Predicted lipophilicity (logD 0.42) and basicity (pKa ≈ 8.1–8.2) are essentially unchanged, so the basic-amine ('Red-Man'/histamine) liability is preserved but not amplified.",{})]),
      p([t("On a weight-of-evidence basis, both impurities qualify for toxicological read-across to the parent vancomycin B, whose safety is anchored by a 91-day intravenous dog study (NOAEL 75/100 mg/kg/day, ≈2–3× the clinical dose) and extensive clinical exposure, provided each impurity is controlled within the Ph. Eur. related-substances limits. The in silico outcome should be confirmed under the company's ICH M7 two-methodology (expert rule-based + statistical) framework and expert review.",{})]),
      table([2600,2400,2400,2240],[
        new TableRow({children:[hdrCell("Question",2600),hdrCell("RRT 0.87",2400),hdrCell("RRT 0.75",2400),hdrCell("Basis",2240)]}),
        new TableRow({children:[cell("New toxicophore vs parent?",{w:2600,bold:true,shade:LIGHT}),cell("NO",{w:2400,align:AlignmentType.CENTER,bold:true,color:"2E7D32"}),cell("NO",{w:2400,align:AlignmentType.CENTER,bold:true,color:"2E7D32"}),cell("Identical FG inventory (0.87); iso-functional amine swap (0.75)",{w:2240})]}),
        new TableRow({children:[cell("New mutagenic (M7) alert?",{w:2600,bold:true,shade:LIGHT}),cell("NO",{w:2400,align:AlignmentType.CENTER,bold:true,color:"2E7D32"}),cell("NO",{w:2400,align:AlignmentType.CENTER,bold:true,color:"2E7D32"}),cell("No DNA-reactive alert in parent or daughters",{w:2240})]}),
        new TableRow({children:[cell("Read-across to parent justified?",{w:2600,bold:true,shade:LIGHT}),cell("YES",{w:2400,align:AlignmentType.CENTER,bold:true,color:"2E7D32"}),cell("YES",{w:2400,align:AlignmentType.CENTER,bold:true,color:"2E7D32"}),cell("≥98% shared skeleton; conserved core toxicophores",{w:2240})]}),
      ]),
      new Paragraph({children:[new PageBreak()]}),

      // ===== 2. PURPOSE & FRAMEWORK =====
      h1("2.  Purpose, scope and regulatory framework"),
      p("The objective is a toxicological in silico qualification of impurities RRT 0.75 and RRT 0.87 in Xellia vancomycin, delivering: (i) a 3D-QSAR / structure-based assessment for the parent molecule and each of the two daughter molecules; (ii) identification and assessment of any additional toxophore present in a daughter but absent from the parent; (iii) a documented read-across justification; and (iv) this regulatory-style report."),
      p([t("Framework. ",{bold:true}), t("Non-mutagenic impurity qualification follows the logic of ICH Q3A(R2)/Q3B(R2) (qualification of impurities in new drug substances/products). Mutagenic-impurity potential is screened using the logic of ICH M7(R2), which anticipates two complementary (Q)SAR methodologies — one expert rule-based, one statistical — plus expert review. In silico model use is framed by the OECD principles for the validation of (Q)SAR models (defined endpoint, unambiguous algorithm, defined applicability domain, appropriate goodness-of-fit/robustness, and mechanistic interpretation where possible). Because vancomycin B is a large (MW ≈ 1449) flexible glycopeptide well outside the applicability domain of conventional small-molecule statistical tox models, the assessment is weighted toward mechanistic, structure- and shape-based read-across, with statistical models used only where their domain is respected (and flagged where not).",{})]),

      // ===== 3. TEST ITEMS =====
      h1("3.  Test items and structural relationships"),
      p("Structures were taken from the sponsor information (IUPAC name, SMILES and InChI) and independently validated with RDKit; computed molecular formulae and exact masses reproduce the values in the source document, confirming faithful structure capture."),
      table([2200,2480,2480,2480], IDrows.map((r,i)=> new TableRow({children:[
        cell(r[0],{w:2200,bold:true,shade: i===0?NAVY:LIGHT, color:i===0?"FFFFFF":undefined}),
        cell(r[1],{w:2480,shade:i===0?NAVY:undefined,color:i===0?"FFFFFF":undefined,bold:i===0,align:i===0?AlignmentType.CENTER:undefined}),
        cell(r[2],{w:2480,shade:i===0?NAVY:undefined,color:i===0?"FFFFFF":undefined,bold:i===0,align:i===0?AlignmentType.CENTER:undefined}),
        cell(r[3],{w:2480,shade:i===0?NAVY:undefined,color:i===0?"FFFFFF":undefined,bold:i===0,align:i===0?AlignmentType.CENTER:undefined}),
      ]}))),
      img("fig1_structures.png", 640, 210),
      cap("Figure 1.  2D structures of the parent API (left) and the two impurities. RRT 0.87 is drawn identically to the parent because it differs only in stereochemistry; RRT 0.75 differs at N-terminal residue 1 (free primary amine + isoleucine side chain instead of N-methyl-leucine)."),
      p([t("RRT 0.87 — (7S,26R)-vancomycin B. ",{bold:true}), t("Same molecular formula as the API (C66H75Cl2N9O24). Configuration is inverted at C7 (the benzylic carbinol carbon of residue 3) and at C26 (a peptide α-carbon). This is consistent with thermally-driven epimerisation under low water activity — a stereochemical, not constitutional, change. No atoms are added or removed.",{})]),
      p([t("RRT 0.75 — (26R)-DMe-DeLI. ",{bold:true}), t("Molecular formula C65H73Cl2N9O24 (one CH2 less than the API). N-terminal residue 1 changes from N-methyl-D-leucine to des-N-methyl-isoleucine: the N-methyl group is lost (secondary → primary amine) and the side chain changes from isobutyl (leucine) to sec-butyl (isoleucine), with additional epimerisation at C26. All other residues, both sugars, both aryl chlorides and all phenolic rings are retained.",{})]),
      img("fig2_parent_residue1.png", 300, 250),
      img("fig2_rrt075_residue1.png", 300, 250),
      cap("Figure 2.  Residue-1 difference highlighted: parent N-methyl-D-leucine (secondary amine) vs RRT 0.75 des-methyl isoleucine (primary amine). This is the only site of constitutional change between RRT 0.75 and the API."),
      new Paragraph({children:[new PageBreak()]}),

      // ===== 4. METHODS =====
      h1("4.  Computational methods"),
      h2("4.1  Toxicophore / pharmacophore and structural-alert profiling"),
      p("Each structure was profiled for (a) the vancomycin pharmacophore elements (heptapeptide aglycone D-Ala-D-Ala binding cleft, H-bonding network, vancosamine–glucose disaccharide) and (b) the recognised vancomycin toxicophores (phenolic/aromatic rings linked to renal-tubular accumulation; basic amine centres linked to non-IgE mast-cell/histamine release / Red-Man syndrome; high MW / hydrophilicity). A quantitative functional-group inventory (SMARTS substructure counting) was computed for all three molecules and differenced against the parent. Structural liabilities were screened with the RDKit FilterCatalog PAINS (A/B/C), Brenk and NIH catalogues, and reviewed against ICH M7 DNA-reactive alert classes (aromatic amine, nitro/nitroso, epoxide, Michael acceptor, aldehyde, alkyl-/benzyl-halide, etc.)."),
      h2("4.2  Physicochemical and 3D descriptors; 3D-QSAR-style shape read-across"),
      p("Physicochemical descriptors (MW, cLogP, TPSA, HBD/HBA, rotatable bonds, aromatic-ring count, QED, Lipinski/Veber) were computed with RDKit; logD7.4 and acidic/basic pKa were predicted with the Inductive Bio (Q)SAR models. Three-dimensional conformers were generated (ETKDG v3) and MMFF-minimised; 3D whole-molecule shape descriptors (radius of gyration, asphericity, spherocity, eccentricity, normalised principal moments NPR1/NPR2) were computed. Read-across in three dimensions was quantified by (i) maximum-common-substructure (MCS) determination and core-restricted 3D alignment (RMSD), and (ii) Crippen Open3DAlign / shape-Tanimoto comparison to the parent. This constitutes a 3D-structure/field-based read-across; it is the applicable in silico 3D approach for a three-compound impurity qualification (a formal CoMFA/CoMSIA 3D-QSAR model is not constructed — see Section 9)."),

      // ===== 5. RESULTS - parent =====
      h1("5.  Results"),
      h2("5.1  Parent (Vancomycin B) — baseline toxicophore profile"),
      p("The parent presents the expected glycopeptide profile: 6 phenolic OH, 5 aromatic rings, 2 aryl chlorides, 3 aryl ethers, 2 benzylic alcohols, 1 carboxylic acid, 1 primary amide (asparagine), 6 secondary amides (peptide backbone), 1 primary amine (vancosamine) and 1 secondary N-methyl amine (residue-1 N-methyl-leucine). It carries no PAINS and no Brenk alert and no DNA-reactive (ICH M7) alert; the single NIH-catalogue flag (≥ 7 aliphatic hydroxyls) is a polarity/solubility descriptor, not a genotoxic alert. Physicochemically it is a large, highly polar, low-lipophilicity molecule (MW 1449; TPSA 531 Å²; logD7.4 ≈ 0.42) — properties that themselves constrain passive membrane permeation and contribute to the known renal-handling / distribution behaviour of vancomycin. These features define the read-across baseline against which the two daughters are compared."),
      img("fig3_toxicophore_inventory.png", 620, 330),
      cap("Figure 3.  Toxicophore / functional-group inventory (counts). Parent and RRT 0.87 are identical across every category; RRT 0.75 differs only in the amine columns (red) — one secondary (N-methyl) amine is replaced by one primary amine."),

      h2("5.2  RRT 0.87 — (7S,26R)-vancomycin B"),
      p([t("Because RRT 0.87 differs from the API only in configuration, its 2D structure, functional-group inventory, structural-alert profile and all constitution-based descriptors are identical to the parent (ECFP4 and MACCS Tanimoto = 1.00; 101/101 shared heavy atoms; 0 PAINS/Brenk/M7 alerts). ",{}),
         t("No additional toxicophore is present. ",{bold:true}),
         t("Conventional 2D QSAR/read-across cannot distinguish this impurity from the API at all — which is precisely why a 3D/mechanistic assessment is required. The 3D shape descriptors of the two are near-congruent (radius of gyration 7.45 vs 7.50 Å; NPR1 0.529 vs 0.528; NPR2 0.615 vs 0.603), i.e. the molecule occupies essentially the same shape space. Stereochemical inversion at C7/C26 can modulate D-Ala-D-Ala target engagement (an efficacy/potency consideration raised by the sponsor) but does not create any new reactive, electrophilic or DNA-reactive functionality, and therefore does not introduce a new toxicological hazard relative to the parent.",{})]),

      h2("5.3  RRT 0.75 — (26R)-DMe-DeLI"),
      p([t("RRT 0.75 retains the entire aglycone, both sugars, both aryl chlorides, all phenols, the carboxylic acid and the primary amide; 99 of ~100 heavy atoms are shared with the API (ECFP4 0.86 / MACCS 0.94). The single constitutional change is at residue 1: the secondary N-methyl amine becomes a primary aliphatic amine. In the functional-group inventory this appears solely as secondary-amine 1→0 and primary-amine 1→2 (Figure 3, red cells). ",{}),
         t("The primary aliphatic amine is not a new functional-group class — the parent already carries one on vancosamine — and a primary aliphatic amine is not a DNA-reactive ICH M7 alert. ",{bold:true}),
         t("No PAINS/Brenk/M7 alert is introduced. Predicted lipophilicity is unchanged (logD7.4 0.42 vs 0.42) and predicted most-basic pKa is essentially unchanged (8.08 vs 8.21, both 'moderate basicity'), so the basic-amine liability associated with histamine release / Red-Man syndrome is conserved in magnitude rather than amplified. TPSA rises only marginally (+14 Å², reflecting the extra N–H). 3D shape remains close to the parent (radius of gyration 7.46 Å); the residue-1 arm is slightly more extended (NPR1 0.42), consistent with loss of the N-methyl and the Leu→Ile change, but the toxicophore-bearing rigid core is unchanged.",{})]),

      h2("5.4  Differential toxicophore analysis — is any new toxophore present?"),
      table([3200,3200,3239],[
        new TableRow({children:[hdrCell("Potential toxophore",3200),hdrCell("RRT 0.87 vs parent",3200),hdrCell("RRT 0.75 vs parent",3239)]}),
        new TableRow({children:[cell("Phenolic / aromatic rings (renal accumulation)",{w:3200}),cell("Identical (6 phenol, 5 aryl)",{w:3200}),cell("Identical (6 phenol, 5 aryl)",{w:3239})]}),
        new TableRow({children:[cell("Basic amine centres (histamine / Red-Man)",{w:3200}),cell("Identical",{w:3200}),cell("Conserved: sec→prim amine, iso-basic (pKa 8.1); not amplified",{w:3239,shade:AMBER})]}),
        new TableRow({children:[cell("Aryl chlorides",{w:3200}),cell("Identical (2)",{w:3200}),cell("Identical (2)",{w:3239})]}),
        new TableRow({children:[cell("DNA-reactive alert (aromatic amine, nitro, epoxide, Michael acceptor, alkyl halide…)",{w:3200}),cell("None (as parent)",{w:3200}),cell("None (as parent)",{w:3239})]}),
        new TableRow({children:[cell("PAINS / Brenk liability",{w:3200}),cell("None",{w:3200}),cell("None",{w:3239})]}),
        new TableRow({children:[cell("NET: new toxophore introduced?",{w:3200,bold:true}),cell("NO",{w:3200,bold:true,color:"2E7D32",align:AlignmentType.CENTER}),cell("NO",{w:3239,bold:true,color:"2E7D32",align:AlignmentType.CENTER}),]}),
      ]),
      p([t("The only element warranting comment is the residue-1 amine in RRT 0.75 (amber). It is a modification of an existing toxicophore class rather than a new one: a secondary amine becomes a primary amine of equal basicity and equal lipophilicity contribution. There is no mechanistic basis for a qualitatively different immunogenic or reactive behaviour, and the change removes (rather than adds) an N-alkyl substituent. Accordingly, no additional toxophore is identified in either daughter molecule.",{})]),

      h2("5.5  3D-QSAR / shape read-across metrics"),
      img("fig5_similarity.png", 470, 300),
      cap("Figure 4.  2D similarity and shared-skeleton read-across metrics vs the parent. RRT 0.87 is fingerprint-identical (stereo-blind methods); RRT 0.75 is highly similar with ≥98% of the heavy-atom skeleton in common."),
      img("fig4_descriptors.png", 640, 250),
      cap("Figure 5.  3D whole-molecule shape descriptors (left) and physicochemical descriptors (right, symlog). The three molecules occupy essentially the same shape and property space; the only visible physchem shift is RRT 0.75's marginally lower MW and cLogP from loss of the N-methyl group."),

      // ===== 6. READ-ACROSS =====
      h1("6.  Read-across justification"),
      p([t("Source (analogue): ",{bold:true}), t("the parent API, vancomycin B, for which a full toxicological package exists (including the 91-day intravenous dog study, NOAEL 75/100 mg/kg/day, and decades of human clinical exposure).",{})]),
      p([t("Target(s): ",{bold:true}), t("RRT 0.87 and RRT 0.75.",{})]),
      bullet([t("Structural similarity. ",{bold:true}), t("RRT 0.87 shares 100% of the heavy-atom skeleton and functional groups (a configurational isomer). RRT 0.75 shares ≥98% (99/100 heavy atoms), with a single localised, non-alerting modification at a peptide N-terminus.",{})]),
      bullet([t("Common mechanism / conserved toxicophores. ",{bold:true}), t("All phenolic, aromatic, aryl-chloride, amide, carboxyl and basic-amine features that drive vancomycin's known effects are conserved. No new reactive/electrophilic centre is created in either daughter.",{})]),
      bullet([t("Metabolic / physicochemical congruence. ",{bold:true}), t("Vancomycin is renally cleared essentially unchanged and is minimally metabolised; the impurities' near-identical size, polarity (TPSA), lipophilicity (logD 0.42) and ionisation (pKa ≈ 8.1) predict the same ADME behaviour, supporting common internal exposure and disposition.",{})]),
      bullet([t("No mutagenic concern. ",{bold:true}), t("Under ICH M7 logic, absence of any DNA-reactive structural alert in the parent or the daughters places them in the non-alerting category; they are controlled as ordinary (non-mutagenic) related substances.",{})]),
      bullet([t("Acknowledged uncertainty (stereochemistry). ",{bold:true}), t("As the sponsor notes, stereoisomers can differ in target binding; here this bears on antibacterial potency/selectivity (efficacy), not on the creation of a new toxicological hazard. The 3D analysis was used specifically to address this dimension that 2D methods cannot.",{})]),

      // ===== 7. CONCLUSION =====
      h1("7.  Qualification conclusion"),
      p([t("On a weight-of-evidence basis the in silico assessment supports toxicological qualification of both impurities by read-across to the parent vancomycin B:",{})]),
      bullet([t("RRT 0.87 [(7S,26R)-vancomycin B] ",{bold:true}), t("— a pure stereoisomer with no new toxicophore and no new mutagenic alert; toxicologically read-across to the API is fully justified.",{})]),
      bullet([t("RRT 0.75 [(26R)-DMe-DeLI] ",{bold:true}), t("— a single, non-alerting, iso-functional residue-1 modification (plus C26 epimerisation) with no new toxicophore and no new mutagenic alert; toxicologically read-across to the API is justified.",{})]),
      p([t("Provided both impurities remain controlled within the Ph. Eur. vancomycin related-substances limits, the parent's established NOAEL (75/100 mg/kg/day, ≈2–3× clinical dose) and clinical safety record provide an adequate toxicological cover for the levels at which these related substances occur. Where an impurity would exceed the applicable ICH Q3A qualification threshold, the same read-across supports qualification without a dedicated animal study; formal sign-off should be issued by a qualified toxicologist and, for the mutagenic dimension, documented under the company's ICH M7 two-(Q)SAR-methodology workflow.",{})]),

      // ===== 8. HONESTY / LIMITATIONS =====
      new Paragraph({children:[new PageBreak()]}),
      h1("8.  Applicability domain, honesty statement and limitations"),
      bullet([t("Nature of the '3D-QSAR'. ",{bold:true}), t("No congeneric training set with measured toxicity endpoints exists for these compounds, so a formal field-based 3D-QSAR (CoMFA/CoMSIA) predictive model was not — and cannot legitimately be — built. What is reported is 3D-structure/shape-based read-across plus mechanistic toxicophore profiling, which is the defensible in silico 3D approach for a three-compound impurity qualification.",{})]),
      bullet([t("Statistical tox models deliberately not over-interpreted. ",{bold:true}), t("Small-molecule statistical panels (e.g. Tox21 nuclear-receptor/stress-response RandomForests) were not used as evidence because a 1435–1449 Da glycopeptide with 100+ heavy atoms lies entirely outside their drug-like applicability domain; their output would be uninterpretable and would convey false precision. The logD/pKa (Q)SAR predictions were flagged out-of-domain by the provider and are used only directionally.",{})]),
      bullet([t("Conformational sampling. ",{bold:true}), t("3D descriptors and whole-molecule shape overlaps derive from single low-energy conformers of a large, highly flexible macrocycle; absolute alignment RMSD/shape-Tanimoto values are approximate and dominated by peripheral (sugar/side-arm) flexibility, so they are used as corroborating, not decisive, evidence. The decisive evidence is the near-complete shared skeleton and the conserved toxicophore inventory.",{})]),
      bullet([t("Scope. ",{bold:true}), t("This is a hazard-identification / structure-based screening and read-across exercise. It does not model quantitative potency, does not address the immunological mechanism of Red-Man syndrome at the molecular level, and does not replace expert toxicological review or any experimental testing that a competent authority may require.",{})]),
      bullet([t("Impurity control limits. ",{bold:true}), t("Ph. Eur. numerical related-substances limits for RRT 0.75 and RRT 0.87 were not provided in the source document; the qualification conclusion assumes control within the applicable monograph limits and should be reconciled with the actual specification.",{})]),

      // ===== 9. METHODS/DATA SUMMARY TABLE =====
      h1("9.  Consolidated data table"),
      table([3040,1650,1650,1650,1649], PCrows.map((r,i)=> new TableRow({children:r.map((c,j)=>
        cell(c,{ w:[3040,1650,1650,1650,1649][j], bold:i===0, color:i===0?"FFFFFF":undefined, shade:i===0?NAVY:(j===0?LIGHT:undefined), align:i===0||j>0?AlignmentType.CENTER:undefined })
      )}))),
      new Paragraph({ spacing:{before:80}, children:[ new TextRun({ text:"* RRT 0.87 is 2D/constitution-identical to the parent, so constitution-based predictions (logD, pKa) are identical to the parent by definition.", italics:true, size:16, color:GREY }) ]}),

      // ===== 10. REFERENCES =====
      h1("10.  References and tools"),
      bullet("ICH Q3A(R2) — Impurities in New Drug Substances; ICH Q3B(R2) — Impurities in New Drug Products (qualification thresholds and read-across logic)."),
      bullet("ICH M7(R2) — Assessment and Control of DNA-Reactive (Mutagenic) Impurities (two complementary (Q)SAR methodologies + expert review)."),
      bullet("OECD (2007) — Guidance on the Validation of (Q)SAR Models (five validation principles)."),
      bullet("RDKit (2024) open-source cheminformatics — structure parsing, descriptors, PAINS/Brenk/NIH FilterCatalogs, ETKDGv3 conformers, MMFF94, Open3DAlign / shape metrics, MCS."),
      bullet("Inductive Bio global logD & pKa (Q)SAR models (v1.4–1.7) — supporting physicochemical predictions."),
      bullet("Ph. Eur. Vancomycin (hydrochloride) monograph — nomenclature and related-substances framework."),
      bullet("Sponsor source document: 'Vanco_RRT_075_and_087_Info' (Xellia) — identities, SMILES/InChI, origin, and 91-day dog NOAEL 75/100 mg/kg/day."),
    ]
  }]
});

Packer.toBuffer(doc).then(buf=>{ fs.writeFileSync(`${DIR}/Vancomycin_RRT075_RRT087_3DQSAR_Qualification_Report.docx`, buf); console.log("WROTE report docx", buf.length, "bytes"); });
