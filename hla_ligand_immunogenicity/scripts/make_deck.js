#!/usr/bin/env node
/**
 * Render report.pptx from results/deck_data.json (written by make_deck.py).
 * Every number on every slide comes out of that file.
 */
const fs = require("fs");
const pptxgen = require("pptxgenjs");

const data = JSON.parse(fs.readFileSync(process.argv[2], "utf8"));

const INK = "12212E", PAPER = "FFFFFF", SOFT = "F2F6F9", LINE = "DCE4EC";
const MUTED = "63748A", ACCENT = "2F6F9F", BAD = "B23B3B", WARN = "C9752B",
      GOOD = "3D8A6B", VIOLET = "7A5EA8";
const HEAD = "Cambria", BODY = "Calibri";

const pres = new pptxgen();
pres.layout = "LAYOUT_WIDE";           // 13.3 x 7.5 in
pres.author = "hla_ligand_immunogenicity pipeline";
pres.title = "AAV Ligand Immunogenicity";

const W = 13.3, H = 7.5, M = 0.62;

const pct = (x, d = 1) => (parseFloat(x) * 100).toFixed(d) + "%";
const num = (x, d = 2) => parseFloat(x).toFixed(d);

function light(title, kicker) {
  const s = pres.addSlide();
  s.background = { color: PAPER };
  if (kicker) s.addText(kicker.toUpperCase(), {
    x: M, y: 0.34, w: W - 2 * M, h: 0.24, fontSize: 10.5, bold: true,
    color: ACCENT, charSpacing: 1.6, fontFace: BODY, margin: 0 });
  s.addText(title, { x: M, y: 0.6, w: W - 2 * M, h: 0.6, fontSize: 27, bold: true,
    color: INK, fontFace: HEAD, margin: 0 });
  return s;
}

function dark(title, sub, top) {
  const s = pres.addSlide();
  s.background = { color: INK };
  const ty = top ? 0.5 : 2.25;
  s.addText(title, { x: M, y: ty, w: W - 2 * M, h: top ? 0.72 : 1.55,
    fontSize: top ? 28 : 40, bold: true, color: PAPER, fontFace: HEAD, margin: 0 });
  if (sub) s.addText(sub, { x: M, y: ty + (top ? 0.78 : 1.6), w: W - 2 * M - 1.5,
    h: 1.35, fontSize: 15, color: "AEC2D4", fontFace: BODY, margin: 0, lineSpacing: 22 });
  return s;
}

function statTile(s, x, y, w, n, label, sub, color) {
  s.addShape(pres.ShapeType.roundRect, { x, y, w, h: 1.42, fill: { color: SOFT },
    line: { color: LINE, width: 0.75 }, rectRadius: 0.06 });
  s.addText(n, { x: x + 0.16, y: y + 0.1, w: w - 0.32, h: 0.5, fontSize: 26, bold: true,
    color: color || INK, fontFace: HEAD, margin: 0 });
  s.addText(label.toUpperCase(), { x: x + 0.16, y: y + 0.6, w: w - 0.32, h: 0.24,
    fontSize: 8.5, bold: true, color: MUTED, charSpacing: 0.9, fontFace: BODY, margin: 0 });
  s.addText(sub, { x: x + 0.16, y: y + 0.85, w: w - 0.32, h: 0.5, fontSize: 9.5,
    color: MUTED, fontFace: BODY, margin: 0, lineSpacing: 12 });
}

function tbl(s, rows, opts) {
  s.addTable(rows, Object.assign({
    x: M, y: 1.35, w: W - 2 * M, border: { type: "solid", color: LINE, pt: 0.5 },
    fontFace: BODY, fontSize: 10, color: INK, valign: "middle",
    autoPage: false,
  }, opts));
}

function head(cells) {
  return cells.map((t) => ({ text: t, options: { bold: true, color: MUTED, fill: { color: SOFT },
    fontSize: 8.5, charSpacing: 0.6 } }));
}

function note(s, text, y) {
  s.addText(text, { x: M, y: y, w: W - 2 * M, h: 0.55, fontSize: 9.5, color: MUTED,
    fontFace: BODY, italic: true, margin: 0, lineSpacing: 12 });
}

function figure(s, path, o) {
  s.addImage(Object.assign({ path, x: M, y: 1.35, w: W - 2 * M }, o));
}

// ------------------------------------------------------------------ 1 title
{
  const s = dark("HLA-DR immunogenicity of an AAV affinity ligand",
    `${data.test_id} · ${data.test_length} aa camelid VHH · ${data.test_source}\n` +
    `${data.panel.panel_size_total}-molecule HLA-DR panel · ` +
    `${pct(data.panel.weighted_coverage)} weighted US/EU coverage · ` +
    `benchmarked against ${data.rank.length - 1} comparators and controls`);
  s.addText("PROCESS-RELATED IMPURITY  ·  RESEARCH USE ONLY", {
    x: M, y: 1.75, w: W - 2 * M, h: 0.3, fontSize: 11, bold: true, color: "7FA8C9",
    charSpacing: 2.2, fontFace: BODY, margin: 0 });
  s.addText("NetMHCIIpan EL + BA  ·  human-self tolerance filter  ·  population-weighted scoring  " +
            "·  BepiPred-2.0  ·  exposure banding", {
    x: M, y: 6.5, w: W - 2 * M, h: 0.4, fontSize: 10.5, color: "6E8AA3",
    fontFace: BODY, margin: 0 });
  s.addNotes("Test article is a public stand-in (PDB 9DC3) so the whole pipeline is reproducible. " +
             "Proprietary ligands drop in as de-identified FASTA and every module runs unchanged.");
}

// -------------------------------------------------------------- 2 exec summ
{
  const s = light("Executive summary", "Verdict");
  const tr = data.test_rank, ar = data.anchor_rank, sm = data.summary;
  const w = (W - 2 * M - 3 * 0.22) / 4;
  const bandColor = { "comparable-to-benchmark": GOOD, "modestly-elevated": WARN,
                      "elevated": BAD, "high": BAD }[tr.risk_band] || ACCENT;
  statTile(s, M, 1.32, w, num(tr.fold_vs_ProteinA_Z, 1) + "×", "vs Protein A Z-domain",
    "the affinity ligand with the longest controlled clinical leachate history", bandColor);
  statTile(s, M + (w + 0.22), 1.32, w, sm.n_foreign_epitopes, "non-self epitopes",
    `of ${sm.n_epitopes} total; ${sm.n_tolerised_epitopes} matched human self and were down-weighted`, BAD);
  statTile(s, M + 2 * (w + 0.22), 1.32, w, Math.round(parseFloat(tr.pop_at_risk_pct)) + "%",
    "US/EU population at risk",
    "carry a DR molecule predicted to present ≥1 non-self epitope", WARN);
  statTile(s, M + 3 * (w + 0.22), 1.32, w, `${sm.max_promiscuity}/${data.panel.panel_size_total}`,
    "peak promiscuity", "DR molecules presenting the dominant epitope", ACCENT);

  const tc = data.top_cluster;
  s.addShape(pres.ShapeType.roundRect, { x: M, y: 3.0, w: W - 2 * M, h: 1.12,
    fill: { color: "FBF1F1" }, line: { color: "E8CDCD", width: 0.75 }, rectRadius: 0.06 });
  s.addText([
    { text: "The risk is concentrated, not diffuse.  ", options: { bold: true, color: BAD } },
    { text: tc
      ? `Residues ${tc.start}–${tc.end}, binding core ${tc.peak_core}, peak 15-mer ` +
        `${tc.peak_peptide} — presented by ${tc.union_sb_alleles} of ` +
        `${data.panel.panel_size_total} DR molecules, reaching ${pct(tc.pop_presenting, 0)} of the ` +
        `weighted US/EU population. That single peptide is what a confirmatory assay should be ` +
        `built around; the rest of the sequence is unremarkable for a non-human protein this size.`
      : "No non-self cluster was called.", options: { color: INK } },
  ], { x: M + 0.2, y: 3.12, w: W - 2 * M - 0.4, h: 0.9, fontSize: 12, fontFace: BODY,
       margin: 0, lineSpacing: 16 });

  const bullets = [
    [`Intrinsic score ${num(tr.pIRS)} population-weighted non-self epitope units / 100 aa, against `
      + `${num(ar.pIRS)} for the Protein A Z-domain — same batch, same panel, same thresholds.`, BAD],
    [`Batch controls ${data.suitability.batch_valid ? "all pass" : "DO NOT all pass"}: the panel `
      + `finds the tetanus universal epitopes, and the tolerance filter suppresses the human `
      + `self controls without flattening the foreign ligands.`, data.suitability.batch_valid ? GOOD : BAD],
    [`Requiring the binding-affinity head to agree with eluted-ligand scoring removed `
      + `${Math.round(100 * (data.n_el_sb - data.n_cons_sb) / data.n_el_sb)}% of strong-binder calls `
      + `across the batch — peptides that look presented but do not measurably bind.`, ACCENT],
    [`Risk is intrinsic potential × dose. At ≤1 ng ligand / mg product and a 10 mg dose the `
      + `exposure is ~0.01 µg — the regime in which qualified Protein A leachate has decades of `
      + `clinical use.`, GOOD],
  ];
  let y = 4.32;
  bullets.forEach(([t, c]) => {
    s.addShape(pres.ShapeType.ellipse, { x: M + 0.04, y: y + 0.09, w: 0.14, h: 0.14,
      fill: { color: c }, line: { color: c } });
    s.addText(t, { x: M + 0.34, y: y - 0.02, w: W - 2 * M - 0.4, h: 0.62, fontSize: 11.5,
      color: INK, fontFace: BODY, margin: 0, lineSpacing: 15 });
    y += 0.66;
  });
}

// ------------------------------------------------- 3 workflow, what changed
{
  const s = light("What a single-tool DR scan cannot tell you", "Workflow review");
  const rows = [
    head(["", "Added step", "The failure mode it fixes"]),
    ...[
      ["M2", "Panel designed against measured coverage",
       `The DR subset of the IEDB class II reference set is widely used as "broadly representative". ` +
       `Measured, it reaches ${pct(data.panel.legacy_weighted_coverage)} weighted US/EU DRB1 coverage — ` +
       `it cannot meet a 95–98% requirement at any size. A greedy build to a stated target reaches ` +
       `${pct(data.panel.weighted_coverage)} at ${data.panel.panel_size_drb1} DRB1 molecules.`],
      ["M3", "Two orthogonal prediction heads",
       `Eluted-ligand scoring alone over-calls: it rewards motif-like peptides with poor measured ` +
       `affinity. Requiring EL and BA agreement dropped ` +
       `${Math.round(100 * (data.n_el_sb - data.n_cons_sb) / data.n_el_sb)}% of strong-binder calls.`],
      ["M4", "Self / pre-existing-tolerance filter",
       `Most DR hits in an antibody-derived ligand are framework cores near-identical to human ` +
       `germline V. Counting them inflates every VHH, scFv and Fab ligand identically and destroys ` +
       `the ranking. Validated here against a shuffled-sequence null.`],
      ["M5–M6", "Population weighting + benchmark calibration with controls",
       `"13 strong binders" is uninterpretable. Weighting each hit by the fraction of the US/EU ` +
       `population carrying the presenting molecule, then expressing it against Protein A Z-domain, ` +
       `makes it a comparison rather than a count.`],
      ["M7–M8", "B-cell / ADA layer and exposure context",
       `The measured endpoint is an anti-drug antibody assay, and impurity risk scales with µg ` +
       `delivered per dose. A T-cell-only, dose-free score cannot reach a risk call.`],
    ].map(([m, a, f]) => ([
      { text: m, options: { bold: true, color: ACCENT, fontSize: 10 } },
      { text: a, options: { bold: true, fontSize: 10 } },
      { text: f, options: { fontSize: 9.5, color: "3B4A5A" } },
    ])),
  ];
  tbl(s, rows, { colW: [0.85, 2.9, W - 2 * M - 3.75], rowH: 0.86, y: 1.3 });
  note(s, "M9 (anchor-position deimmunisation scan) is optional — it only matters if the ligand " +
          "itself can be re-engineered.", 6.72);
}

// ------------------------------------------------------------ 4 panel design
{
  const s = light("The panel had to be designed, not inherited", "Module 2 · panel");
  figure(s, data.figures.panel, { y: 1.26, w: W - 2 * M, h: 3.55 });
  const p = data.panel;
  s.addText([
    { text: `${pct(p.legacy_weighted_coverage)} → ${pct(p.weighted_coverage)}   `,
      options: { bold: true, fontSize: 17, color: ACCENT, fontFace: HEAD } },
    { text: `weighted US/EU DRB1 phenotypic coverage, ` +
            `${p.legacy_drb1_panel.length + 4} → ${p.panel_size_total} DR molecules. ` +
            `The missing coverage in the legacy set sits in alleles it does not contain — ` +
            `DRB1*13:01, *04:04, *11:04, *16:01, *14:01 — so no amount of re-weighting fixes it.`,
      options: { fontSize: 12, color: INK } },
  ], { x: M, y: 5.02, w: W - 2 * M, h: 1.05, fontSize: 12, fontFace: BODY, margin: 0,
       lineSpacing: 16 });
  const drb1 = p.drb1_panel.map((a) => a.replace("HLA-DRB1*", "")).join(" · ");
  const other = p.drb345_panel.map((a) => a.replace("HLA-DR", "DR")).join(" · ");
  s.addText(`Panel — DRB1*  ${drb1}    +    ${other}`,
    { x: M, y: 6.13, w: W - 2 * M, h: 0.38, fontSize: 9.5, color: MUTED, fontFace: BODY,
      margin: 0 });
  note(s, `Coverage is single-locus Hardy–Weinberg over the IEDB tables, reproducing the IEDB CLI ` +
          `to two decimals. DRB3/4/5 carry no frequencies there, so they add presentation breadth ` +
          `without entering this arithmetic. US Asian lands at ` +
          `${pct(p.per_population["United States Asian"], 0)} and US Amerindian at ` +
          `${pct(p.per_population["United States Amerindian"], 0)}: a US/EU-weighted objective does ` +
          `not buy those populations.`, 6.62);
}

// --------------------------------------------------------- 5 binding scan
{
  const s = light("Binding landscape across the full panel", "Module 3 · prediction");
  figure(s, data.figures.landscape, { y: 1.2, w: W - 2 * M, h: 4.35 });
  s.addText([
    { text: `${data.n_el_sb} → ${data.n_cons_sb}   `,
      options: { bold: true, fontSize: 17, color: ACCENT, fontFace: HEAD } },
    { text: `strong-binder calls across the batch once the binding-affinity head has to agree with ` +
            `eluted-ligand scoring (BA %Rank < ${data.ba_confirm_rank}). ` +
            `${Math.round(100 * (data.n_el_sb - data.n_cons_sb) / data.n_el_sb)}% of EL-only calls are ` +
            `peptides that look presented but do not measurably bind; every downstream number uses ` +
            `the consensus call.`, options: { fontSize: 12, color: INK } },
  ], { x: M, y: 5.68, w: W - 2 * M, h: 1.0, fontSize: 12, fontFace: BODY, margin: 0,
       lineSpacing: 16 });
  note(s, "15-mer scan, EL %Rank < 1 = strong, < 5 = weak. Boxes mark consolidated epitope " +
          "clusters; the track below counts DR molecules with a strong call at each frame.", 6.68);
}

// -------------------------------------------------------- 6 tolerance filter
{
  const s = light("Framework hits are not risk", "Module 4 · tolerance filter");
  const v = data.filter_validation;
  const w = (W - 2 * M - 2 * 0.22) / 3;
  statTile(s, M, 1.3, w, pct(v.real_hit_rate_9of9, 1), "cores that are exact human 9-mers",
    "the peptide itself occurs in the human proteome — weight 0", GOOD);
  statTile(s, M + w + 0.22, 1.3, w, pct(v.real_hit_rate_8of9 - v.real_hit_rate_9of9, 1),
    "one substitution from human", "8/9 identity to a human 9-mer — weight 0.35", WARN);
  statTile(s, M + 2 * (w + 0.22), 1.3, w, v.enrichment_8of9_real_over_null + "×",
    "enrichment over null", `real cores ${pct(v.real_hit_rate_8of9, 1)} vs shuffled-sequence null ` +
    `${pct(v.null_hit_rate_8of9, 1)}`, ACCENT);

  s.addText([
    { text: "The threshold is validated, not asserted.\n", options: { bold: true, fontSize: 13 } },
    { text: `Re-running the identical screen on cores drawn from shuffled versions of the same ` +
            `sequences — same composition, no real self-similarity — gives ` +
            `${pct(v.null_hit_rate_8of9, 1)} at the 8/9 cut against ${pct(v.real_hit_rate_8of9, 1)} ` +
            `for real cores. A 5-of-9 "TCR-face" screen, the obvious shortcut, was tried first and ` +
            `discarded: with only five positions specified it matches the human proteome by chance ` +
            `several times per query and flags nearly everything.`, options: { fontSize: 12 } },
  ], { x: M, y: 3.1, w: W - 2 * M, h: 1.5, color: INK, fontFace: BODY, margin: 0,
       lineSpacing: 17 });

  const selfDrop = Math.max(...data.rank.filter((r) => r.role === "negative_control_self")
    .map((r) => parseFloat(r.tolerance_drop_pct)));
  s.addShape(pres.ShapeType.roundRect, { x: M, y: 4.75, w: W - 2 * M, h: 0.95,
    fill: { color: SOFT }, line: { color: LINE, width: 0.75 }, rectRadius: 0.06 });
  s.addText(`The filter discriminates rather than deletes: the human self controls lose up to ` +
            `${selfDrop.toFixed(0)}% of their unfiltered score, while the test article loses ` +
            `${parseFloat(data.test_rank.tolerance_drop_pct).toFixed(0)}%. Every flagged core is ` +
            `written out with the human protein it matched, so each call is checkable.`,
    { x: M + 0.2, y: 4.87, w: W - 2 * M - 0.4, h: 0.75, fontSize: 12, color: INK,
      fontFace: BODY, margin: 0, lineSpacing: 16 });
  note(s, "Reference: UniProt Swiss-Prot Homo sapiens. This is a screen, not JanusMatrix — it " +
          "does not require the human counterpart to bind the same allele, so it errs toward " +
          "calling more peptides tolerised.", 6.0);
}

// --------------------------------------------- 6b boundary of the filter
if (data.promiscuity_calibration && data.promiscuity_calibration.boundary_controls) {
  const bc = data.promiscuity_calibration.boundary_controls;
  const nP = data.promiscuity_calibration.panel_size;
  const s = light("Where the filter is wrong, measured", "Module 4 · boundary controls");
  s.addText("The filter's rule is \u201cself implies tolerated\u201d. That is right often enough to " +
            "be worth applying and wrong often enough to be worth measuring, so the batch carries " +
            "two controls that sit on the boundary.",
    { x: M, y: 1.28, w: W - 2 * M, h: 0.5, fontSize: 12, color: INK, fontFace: BODY,
      margin: 0, lineSpacing: 16 });
  const rows = [
    head(["Control", "Role", "DR at %Rank<1", "at <5", "Cores tolerised", "pIRS unfiltered", "pIRS filtered"]),
    ...Object.entries(bc).map(([k, v]) => ([
      { text: k.replace("_region", ""), options: { fontFace: "Courier New", bold: true, fontSize: 9 } },
      { text: v.role.replace(/_/g, " "), options: { color: WARN, fontSize: 9 } },
      { text: `${v.dr_breadth_sb}/${nP}`, options: { align: "right" } },
      { text: `${v.dr_breadth_wb}/${nP}`, options: { align: "right" } },
      { text: `${v.cores_called_tolerised}/${v.predicted_cores_in_window}`, options: { align: "right" } },
      { text: v.ligand_pIRS_unfiltered.toFixed(2), options: { align: "right" } },
      { text: v.ligand_pIRS.toFixed(2), options: { align: "right", bold: true } },
    ])),
  ];
  tbl(s, rows, { y: 1.95, colW: [2.5, 2.2, 1.5, 1.0, 1.6, 1.6, 1.66], rowH: 0.38 });

  const mbp = bc["MBP_85_99_region"];
  if (mbp) {
    s.addShape(pres.ShapeType.roundRect, { x: M, y: 3.35, w: W - 2 * M, h: 1.75,
      fill: { color: "FBF1F1" }, line: { color: "E8CDCD", width: 0.75 }, rectRadius: 0.06 });
    s.addText([
      { text: "The filter suppresses a real epitope, and this is the clearest example.\n",
        options: { bold: true, color: BAD } },
      { text: `Myelin basic protein 85\u201399 is a human self peptide with positive ` +
              `DR-restricted T-cell assays on 10 distinct DR molecules in IEDB. On this panel it is ` +
              `the broadest binder in the entire batch \u2014 ${mbp.dr_breadth_sb}/${nP} DR molecules ` +
              `at the strong-binder tier, unfiltered pIRS ${mbp.ligand_pIRS_unfiltered.toFixed(2)}, ` +
              `higher than any ligand here. The tolerance filter takes it to ` +
              `${mbp.ligand_pIRS.toFixed(2)}.\n\nThat is the filter working as designed and being ` +
              `wrong. Self-derived epitopes exist \u2014 autoimmunity is what they are \u2014 and a ` +
              `sequence-identity filter cannot see the difference. It is safe for an affinity ligand ` +
              `because a foreign scaffold's framework similarity to human germline is the case it was ` +
              `built for; it is not safe for anything where self-reactivity is the question.`,
        options: { color: INK } },
    ], { x: M + 0.22, y: 3.48, w: W - 2 * M - 0.44, h: 1.5, fontSize: 11.5, fontFace: BODY,
         margin: 0, lineSpacing: 15 });
  }
  const clp = bc["CLIP_87_101_region"];
  if (clp) {
    s.addText(`The opposite control behaves too. CLIP occupies the groove of essentially every DR ` +
              `molecule as the invariant chain's placeholder, yet carries no positive human DR ` +
              `T-cell record in IEDB. Here it reaches ${clp.dr_breadth_wb}/${nP} at the weak tier and ` +
              `${clp.dr_breadth_sb}/${nP} at the strong tier, and ends at pIRS ` +
              `${clp.ligand_pIRS.toFixed(2)} \u2014 binding strength alone is not being read as risk.`,
      { x: M, y: 5.3, w: W - 2 * M, h: 0.7, fontSize: 11.5, color: INK, fontFace: BODY,
        margin: 0, lineSpacing: 15 });
  }
  note(s, "Both controls' roles were checked against IEDB before assignment. They do not pass or " +
          "fail the batch \u2014 they put a number on two things the filter cannot do.", 6.35);
}

// ------------------------------------------------ 6c measured accuracy
if (data.accuracy && data.benchmark && data.figures.calibration) {
  const a = data.accuracy, bm = data.benchmark;
  const s = light("Does the rule actually work?", "Modules 10\u201311 · measured accuracy");
  s.addText("Every choice above this point was made by argument. The EL+BA consensus entered the " +
            "pipeline because eluted-ligand scoring over-calls \u2014 a reasonable argument and no " +
            "evidence. So the rule was scored against ground truth: every HLA-DR-restricted CD4 " +
            "T-cell assay result IEDB holds for these " + a.n_pairs.toLocaleString() +
            " peptide-allele pairs, in human hosts.",
    { x: M, y: 1.26, w: W - 2 * M, h: 0.55, fontSize: 12, color: INK, fontFace: BODY,
      margin: 0, lineSpacing: 16 });
  s.addImage({ path: data.figures.calibration, x: M, y: 1.92, w: 7.5, h: 3.05 });

  const w2 = W - M - (M + 7.8);
  const best = (a.sweep || {}).max_mcc;
  const cur = Object.entries(a.operating_points).find(([k]) => k.startsWith("AND"));
  const geo = (a.comparisons || {}).GEO_vs_EL;
  let yy = 1.95;
  const facts = [
    [`${bm.labelled_pairs.toLocaleString()} labelled outcomes`,
     `${bm.positives.toLocaleString()} positive / ${bm.negatives.toLocaleString()} negative pairs in ` +
     `${bm.clusters.toLocaleString()} 9-mer clusters, so overlapping peptides from one study count once`],
    [`Best rule: ${a.best_continuous_rule}, AUC ${a.rules[a.best_continuous_rule].auc.toFixed(3)}`,
     `against ${a.rules.EL.auc.toFixed(3)} for eluted-ligand scoring alone` +
     (geo ? `; paired cluster bootstrap puts the difference at ${geo.delta >= 0 ? "+" : ""}` +
            `${geo.delta.toFixed(4)}, 95% CI [${geo.ci95[0].toFixed(4)}, ${geo.ci95[1].toFixed(4)}]` : "")],
    cur ? [`Current rule: ${(cur[1].sensitivity * 100).toFixed(0)}% sensitivity, ` +
           `${(cur[1].specificity * 100).toFixed(0)}% specificity`,
           `a flagged peptide is real ${(cur[1].ppv_at_scan_prevalence * 100).toFixed(0)}% of the time ` +
           `at ${(a.scan_prevalence * 100).toFixed(0)}% assumed scan prevalence`] : null,
    best ? [`Calibrated cut: %Rank < ${best.threshold_rank.toFixed(2)}`,
            `${(best.sensitivity * 100).toFixed(0)}% / ${(best.specificity * 100).toFixed(0)}%, ` +
            `MCC ${best.mcc.toFixed(3)}` + (cur ? ` against ${cur[1].mcc.toFixed(3)} for the current rule` : "")] : null,
  ].filter(Boolean);
  facts.forEach(([h, t]) => {
    s.addText([{ text: h + "\n", options: { bold: true, color: ACCENT } },
               { text: t, options: { color: INK } }],
      { x: M + 7.8, y: yy, w: w2, h: 0.95, fontSize: 10.5, fontFace: BODY,
        margin: 0, lineSpacing: 13 });
    yy += 0.82;
  });

  s.addShape(pres.ShapeType.roundRect, { x: M, y: 5.15, w: W - 2 * M, h: 1.15,
    fill: { color: SOFT }, line: { color: LINE, width: 0.75 }, rectRadius: 0.06 });
  s.addText([
    { text: "Read the absolute numbers as an upper bound.  ", options: { bold: true } },
    { text: "NetMHCIIpan is trained on IEDB binding and mass-spec data; the labels here are T-cell " +
            "assay outcomes, a different endpoint, but partial training-set overlap is certain and " +
            "cannot be excluded from outside. The comparison between rules is far more robust \u2014 " +
            "every rule uses the same possibly-leaky predictors on the same peptides, so leakage " +
            "inflates them together and largely cancels in the difference.", options: {} },
  ], { x: M + 0.22, y: 5.28, w: W - 2 * M - 0.44, h: 0.95, fontSize: 10.5, color: INK,
       fontFace: BODY, margin: 0, lineSpacing: 13 });
  note(s, "IEDB's negatives carry their own bias: a peptide tested and found negative in one donor " +
          "set is not a clean non-binder, which inflates apparent specificity.", 6.5);
}

// ----------------------------------------------------------- 7 calibration
{
  const s = light("Calibrated against benchmarks and controls in the same batch",
                  "Modules 5–6 · scoring");
  figure(s, data.figures.ranking, { y: 1.2, w: W - 2 * M, h: 3.9 });
  const checks = data.suitability.checks || [];
  const rows = [
    head(["System suitability check", "Observed", ""]),
    ...checks.map((c) => ([
      { text: c.check, options: { fontSize: 9.5, bold: true } },
      { text: c.detail, options: { fontSize: 9, color: "3B4A5A" } },
      { text: c.pass ? "PASS" : "FAIL",
        options: { fontSize: 9.5, bold: true, color: c.pass ? GOOD : BAD, align: "center" } },
    ])),
  ];
  tbl(s, rows, { y: 5.25, colW: [3.4, W - 2 * M - 4.3, 0.9], rowH: 0.34 });
}

// -------------------------------------------------------------- 8 clusters
{
  const s = light("Where the risk actually sits", "Module 5 · epitope clusters");
  const rows = [
    head(["From", "To", "Core", "Peak 15-mer", "EL %Rank", "DR molecules", "% US/EU presenting", "Class"]),
    ...data.clusters.map((c) => ([
      { text: c.start }, { text: c.end },
      { text: c.peak_core, options: { fontFace: "Courier New", bold: true } },
      { text: c.peak_peptide, options: { fontFace: "Courier New" } },
      { text: c.peak_el_rank, options: { align: "right" } },
      { text: `${c.union_sb_alleles}/${data.panel.panel_size_total}`, options: { align: "right" } },
      { text: pct(c.pop_presenting, 1), options: { align: "right", bold: true } },
      { text: c.tolerance_class.replace(/_/g, " "),
        options: { color: c.tolerance_class === "foreign" ? BAD
                        : c.tolerance_class === "mixed" ? WARN : GOOD, bold: true } },
    ])),
  ];
  tbl(s, rows, { y: 1.3, colW: [0.6, 0.6, 1.25, 1.95, 0.95, 1.25, 1.6, 1.86], rowH: 0.34 });

  const y2 = 1.3 + 0.34 * (data.clusters.length + 1) + 0.35;
  s.addText("Top epitopes by promiscuity", { x: M, y: y2 - 0.32, w: 6, h: 0.28,
    fontSize: 13, bold: true, color: INK, fontFace: HEAD, margin: 0 });
  const erows = [
    head(["Pos", "Core", "15-mer", "EL %Rank", "DR", "% pop", "Class", "Presenting molecules"]),
    ...data.epitopes.map((e) => ([
      { text: e.pos },
      { text: e.core, options: { fontFace: "Courier New", bold: true } },
      { text: e.peptide, options: { fontFace: "Courier New" } },
      { text: e.best_el_rank, options: { align: "right" } },
      { text: e.n_sb_alleles, options: { align: "right" } },
      { text: pct(e.pop_presenting, 0), options: { align: "right" } },
      { text: e.tolerance_class.replace(/_/g, " "),
        options: { color: e.tolerance_class === "foreign" ? BAD : GOOD } },
      { text: e.sb_alleles.replace(/;/g, " "), options: { fontSize: 7.5, color: MUTED } },
    ])),
  ];
  tbl(s, erows, { y: y2, colW: [0.55, 1.25, 1.95, 0.95, 0.55, 0.8, 1.3, W - 2 * M - 7.35],
                  rowH: 0.29, fontSize: 9 });
}

// --------------------------------------------------------------- 9 T/B
{
  const s = light("The measured endpoint is an antibody, not a T cell", "Module 7 · ADA layer");
  figure(s, data.figures.tb, { y: 1.22, w: W - 2 * M, h: 3.1 });
  s.addText("A B cell only class-switches with help from a CD4 T cell recognising a peptide from " +
            "the same protein. The regions worth taking into wet-lab work first are those where a " +
            "non-self DR cluster and a predicted B-cell epitope coincide — shaded above.",
    { x: M, y: 4.45, w: W - 2 * M, h: 0.55, fontSize: 12, color: INK, fontFace: BODY,
      margin: 0, lineSpacing: 16 });
  if (data.coincidence.length) {
    const rows = [
      head(["T-cell cluster", "Core", "% pop presenting", "B-cell region", "Overlap (aa)", "Region sequence"]),
      ...data.coincidence.map((c) => ([
        { text: c.t_cluster },
        { text: c.t_peak_core, options: { fontFace: "Courier New", bold: true } },
        { text: pct(c.t_pop_presenting, 1), options: { align: "right", bold: true } },
        { text: c.b_region },
        { text: c.overlap_aa, options: { align: "right" } },
        { text: c.region_peptide, options: { fontFace: "Courier New", fontSize: 8.5 } },
      ])),
    ];
    tbl(s, rows, { y: 5.1, colW: [1.5, 1.3, 1.7, 1.5, 1.2, W - 2 * M - 7.2], rowH: 0.32 });
  }
  note(s, "BepiPred-2.0 linear B-cell propensity is the weakest model in this pipeline — most " +
          "real ADA epitopes are conformational. This layer prioritises regions; it is never a " +
          "standalone claim.", 6.75);
}

// ---------------------------------------------------------- 10 exposure
{
  const s = light("Intrinsic risk × dose", "Module 8 · exposure context");
  const doses = [...new Set(data.exposure.map((g) => g.dose_mg))];
  const ppms = [...new Set(data.exposure.map((g) => g.leachate_ng_per_mg))];
  const bandColor = { negligible: GOOD, low: GOOD, moderate: WARN, elevated: BAD };
  const rows = [
    head(["Leachate \\ dose", ...doses.map((d) => `${d} mg`)]),
    ...ppms.map((p) => ([
      { text: `${p} ng/mg`, options: { bold: true, fill: { color: SOFT }, color: MUTED } },
      ...doses.map((d) => {
        const g = data.exposure.find((x) => x.dose_mg === d && x.leachate_ng_per_mg === p);
        return { text: `${g.ug_ligand_per_dose} µg\n${g.exposure_band}`,
                 options: { align: "center", color: bandColor[g.exposure_band] || INK,
                            bold: g.exposure_band === "elevated" } };
      }),
    ])),
  ];
  tbl(s, rows, { y: 1.35, colW: [2.2, ...doses.map(() => (W - 2 * M - 2.2) / doses.length)],
                 rowH: 0.62, fontSize: 11 });
  s.addText("A leached ligand is an impurity, and impurity risk is intrinsic potential multiplied " +
            "by how much a patient receives. The intrinsic score places this ligand at " +
            `${num(data.test_rank.fold_vs_ProteinA_Z, 1)}× the Protein A Z-domain; the grid above ` +
            "decides whether that multiple matters for a given product.",
    { x: M, y: 1.4 + 0.62 * (ppms.length + 1) + 0.3, w: W - 2 * M, h: 0.8, fontSize: 12,
      color: INK, fontFace: BODY, margin: 0, lineSpacing: 16 });
  note(s, "Bands are an internal triage convention, not a regulatory threshold. No agency " +
          "publishes a numeric leachate immunogenicity limit; the ICH Q6B / EMA expectation is " +
          "control to a justified, consistently achieved level. Replace this grid with your " +
          "product's measured leachate and dose.", 6.55);
}

// -------------------------------------------------------- 11 deimmunisation
if (data.deimm.length && data.figures.deimm) {
  const s = light("If the ligand can be re-engineered", "Module 9 · deimmunisation");
  s.addImage({ path: data.figures.deimm, x: M, y: 1.3, w: 6.0, h: 3.3 });
  const wt = data.deimm_wt;
  s.addText([
    { text: "Anchor-position saturation scan\n", options: { bold: true, fontSize: 13, fontFace: HEAD } },
    { text: `Every substitution at each MHC-II anchor pocket (P1, P4, P6, P9) of the dominant ` +
            `core was re-scored across the full panel. Wild type is presented by ` +
            `${wt.n_sb_alleles} DR molecules reaching ${pct(wt.pop_presenting, 0)} of the weighted ` +
            `population; the best variants drop that to single digits.\n\n` +
            `BLOSUM62 and the human germline residue at the aligned position are reported because ` +
            `a substitution that abolishes DR binding and breaks the fold is not a design. ` +
            `Germline-matching changes are the ones least likely to create a new epitope or ` +
            `destabilise the framework.`, options: { fontSize: 11.5 } },
  ], { x: M + 6.3, y: 1.3, w: W - M - (M + 6.3), h: 3.3, color: INK, fontFace: BODY,
       margin: 0, lineSpacing: 16 });

  const rows = [
    head(["Variant", "Substitution", "Core", "DR", "% pop", "Δ % pop", "BLOSUM62", "Germline aa"]),
    ...data.deimm.filter((r) => r.variant !== "WT").slice(0, 7).map((r) => ([
      { text: r.variant, options: { bold: true } },
      { text: r.substitution },
      { text: r.core, options: { fontFace: "Courier New" } },
      { text: r.n_sb_alleles, options: { align: "right" } },
      { text: pct(r.pop_presenting, 1), options: { align: "right" } },
      { text: (parseFloat(r.delta_pop_presenting) * 100).toFixed(1),
        options: { align: "right", color: BAD } },
      { text: r.blosum62, options: { align: "right",
        color: parseInt(r.blosum62, 10) >= 0 ? GOOD : WARN, bold: true } },
      { text: r.germline_residue, options: { align: "center" } },
    ])),
  ];
  tbl(s, rows, { y: 4.85, colW: [1.5, 1.4, 1.5, 0.7, 1.0, 1.0, 1.1, W - 2 * M - 8.2],
                 rowH: 0.28, fontSize: 9.5 });
  note(s, "DR presentation only. Effect on ligand–target affinity, resin capacity and alkaline " +
          "stability is not modelled; any candidate goes back through binding and stability " +
          "screens before it means anything.", 6.88);
}

// ------------------------------------------------------------- 12 limits
{
  const s = dark("What this is, and what it is not", null, true);
  s.addText("In-silico DR screening ranks and localises risk. It does not measure it.",
    { x: M, y: 1.3, w: W - 2 * M, h: 0.4, fontSize: 15, color: "AEC2D4", fontFace: BODY,
      margin: 0, italic: true });

  const limits = [
    ["Prediction, not presentation", "NetMHCIIpan scores peptide–MHC binding. It does not model " +
      "antigen uptake, endosomal proteolysis, HLA-DM editing or endosomal stability."],
    ["DR only", "DP and DQ contribute to CD4 responses and are excluded by the panel specification. " +
      "DQ is implicated in several biologic ADA responses."],
    ["pIRS is relative", "Readable only against the benchmarks in this batch. It is not a predicted " +
      "ADA incidence, and no in-silico method available today predicts one."],
    ["No aggregation or adjuvant effect", "Aggregated or particulate impurity is substantially more " +
      "immunogenic than monomer; sequence-based methods cannot see it."],
    ["The benchmark is an argument", "“Comparable predicted epitope content to Protein A” is an " +
      "inference about risk, not evidence of it."],
  ];
  let y = 2.35;
  limits.forEach(([h, t]) => {
    s.addText([
      { text: h + "  —  ", options: { bold: true, color: "E7B76B" } },
      { text: t, options: { color: "C6D5E3" } },
    ], { x: M, y, w: 6.5, h: 0.68, fontSize: 11, fontFace: BODY, margin: 0, lineSpacing: 14 });
    y += 0.72;
  });

  s.addShape(pres.ShapeType.roundRect, { x: 7.55, y: 2.05, w: W - M - 7.55, h: 4.45,
    fill: { color: "1B2E40" }, line: { color: "1B2E40", width: 0 }, rectRadius: 0.05 });
  s.addText("Recommended confirmatory work — all RUO", { x: 7.85, y: 2.25, w: 4.6, h: 0.35,
    fontSize: 14, bold: true, color: PAPER, fontFace: HEAD, margin: 0 });
  const steps = [
    ["1", "HLA-DR competitive binding", "on the flagged peptides against the panel's dominant " +
      "molecules — confirms the prediction directly. Days, lowest cost."],
    ["2", "MAPPs", "on monocyte-derived dendritic cells from HLA-typed donors — shows what is " +
      "actually processed and presented from the intact ligand."],
    ["3", "Ex-vivo PBMC / CD4 proliferation", "across ~50 HLA-typed donors matched to this panel — " +
      "the closest available surrogate for clinical ADA risk."],
  ];
  let sy = 2.85;
  steps.forEach(([n, h, t]) => {
    s.addShape(pres.ShapeType.ellipse, { x: 7.85, y: sy + 0.03, w: 0.3, h: 0.3,
      fill: { color: ACCENT }, line: { color: ACCENT } });
    s.addText(n, { x: 7.85, y: sy + 0.04, w: 0.3, h: 0.28, fontSize: 11, bold: true,
      color: PAPER, align: "center", fontFace: BODY, margin: 0 });
    s.addText([
      { text: h + "  ", options: { bold: true, color: PAPER } },
      { text: t, options: { color: "9FB4C8" } },
    ], { x: 8.3, y: sy, w: W - M - 8.5, h: 0.95, fontSize: 10.5, fontFace: BODY,
         margin: 0, lineSpacing: 13 });
    sy += 1.03;
  });
  s.addText("All three are scoped by the peptide list this pipeline produces — which is the " +
            "practical point of running it.", { x: 7.85, y: 6.05, w: W - M - 8.05, h: 0.45,
    fontSize: 10, color: "7FA8C9", italic: true, fontFace: BODY, margin: 0 });

  s.addText("Research use only — not for regulatory submission without confirmatory wet-lab data.",
    { x: M, y: 6.6, w: W - 2 * M, h: 0.3, fontSize: 10, color: "6E8AA3", fontFace: BODY, margin: 0 });
}

pres.writeFile({ fileName: data.out }).then(() => console.log("wrote " + data.out));
