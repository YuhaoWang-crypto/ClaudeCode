"""Build the campaign PDF report.

Neutral professional styling. A clearly-marked branding slot sits in the header
block so a brand mark can be dropped in later without touching the layout.
"""

from __future__ import annotations

import json
import pathlib
from datetime import date

from reportlab.lib import colors
from reportlab.lib.enums import TA_JUSTIFY, TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.platypus import (BaseDocTemplate, Frame, Image, KeepTogether,
                                NextPageTemplate, PageBreak, PageTemplate,
                                Paragraph, Spacer, Table, TableStyle)

INK = colors.HexColor("#0b0b0b")
INK2 = colors.HexColor("#52514e")
MUTED = colors.HexColor("#898781")
RULE = colors.HexColor("#c3c2b7")
HAIR = colors.HexColor("#e1e0d9")
ACCENT = colors.HexColor("#2a78d6")
WARN = colors.HexColor("#eb6834")
SURFACE = colors.HexColor("#fcfcfb")

PAGE_W, PAGE_H = A4
MARGIN = 18 * mm


def styles():
    ss = getSampleStyleSheet()
    s = {}
    s["title"] = ParagraphStyle("title", parent=ss["Title"], fontName="Helvetica-Bold",
                                fontSize=21, leading=25, textColor=INK,
                                alignment=TA_LEFT, spaceAfter=2)
    s["subtitle"] = ParagraphStyle("subtitle", parent=ss["Normal"],
                                   fontName="Helvetica", fontSize=10.5, leading=14,
                                   textColor=INK2, spaceAfter=2)
    s["slot"] = ParagraphStyle("slot", parent=ss["Normal"], fontName="Helvetica",
                               fontSize=7.5, leading=10, textColor=MUTED)
    s["h1"] = ParagraphStyle("h1", parent=ss["Heading1"], fontName="Helvetica-Bold",
                             fontSize=14, leading=17, textColor=INK,
                             spaceBefore=14, spaceAfter=6)
    s["h2"] = ParagraphStyle("h2", parent=ss["Heading2"], fontName="Helvetica-Bold",
                             fontSize=11, leading=14, textColor=INK,
                             spaceBefore=10, spaceAfter=4)
    s["body"] = ParagraphStyle("body", parent=ss["BodyText"], fontName="Helvetica",
                               fontSize=9.4, leading=13.6, textColor=INK,
                               alignment=TA_JUSTIFY, spaceAfter=6)
    s["small"] = ParagraphStyle("small", parent=ss["BodyText"], fontName="Helvetica",
                                fontSize=8.1, leading=11.4, textColor=INK2,
                                alignment=TA_JUSTIFY, spaceAfter=4)
    s["caption"] = ParagraphStyle("caption", parent=ss["BodyText"],
                                  fontName="Helvetica-Oblique", fontSize=8,
                                  leading=10.8, textColor=INK2, spaceBefore=3,
                                  spaceAfter=10)
    s["ref"] = ParagraphStyle("ref", parent=ss["BodyText"], fontName="Helvetica",
                              fontSize=8.1, leading=11.4, textColor=INK,
                              alignment=TA_LEFT, spaceAfter=4,
                              leftIndent=14, firstLineIndent=-14)
    s["cell"] = ParagraphStyle("cell", parent=ss["BodyText"], fontName="Helvetica",
                               fontSize=7.2, leading=8.8, textColor=INK,
                               alignment=TA_LEFT, spaceAfter=0)
    s["cellh"] = ParagraphStyle("cellh", parent=s["cell"],
                                fontName="Helvetica-Bold", textColor=INK)
    s["callout"] = ParagraphStyle("callout", parent=ss["BodyText"],
                                  fontName="Helvetica", fontSize=8.8, leading=12.4,
                                  textColor=INK, alignment=TA_JUSTIFY,
                                  leftIndent=8, rightIndent=8, spaceBefore=4,
                                  spaceAfter=8, borderPadding=7,
                                  backColor=colors.HexColor("#fdf2ec"),
                                  borderColor=WARN, borderWidth=0)
    return s


def _footer(canvas, doc):
    canvas.saveState()
    canvas.setFillColor(SURFACE)
    canvas.rect(0, 0, PAGE_W, PAGE_H, fill=1, stroke=0)
    canvas.setStrokeColor(HAIR)
    canvas.setLineWidth(0.5)
    canvas.line(MARGIN, 13 * mm, PAGE_W - MARGIN, 13 * mm)
    canvas.setFont("Helvetica", 7.4)
    canvas.setFillColor(MUTED)
    canvas.drawString(MARGIN, 9 * mm,
                      "GSK3B de novo design campaign - computational designs, "
                      "not experimentally validated")
    canvas.drawRightString(PAGE_W - MARGIN, 9 * mm, f"Page {doc.page}")
    canvas.restoreState()


def stat_row(items, width):
    """A row of headline numbers."""
    cells = []
    for value, label in items:
        cells.append([Paragraph(
            f'<font size="15" color="#2a78d6"><b>{value}</b></font><br/>'
            f'<font size="7.2" color="#52514e">{label}</font>',
            ParagraphStyle("st", fontName="Helvetica", leading=13.5))])
    t = Table([[c[0] for c in cells]], colWidths=[width / len(cells)] * len(cells))
    t.setStyle(TableStyle([
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("TOPPADDING", (0, 0), (-1, -1), 8),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 8),
        ("LINEABOVE", (0, 0), (-1, 0), 0.8, RULE),
        ("LINEBELOW", (0, -1), (-1, -1), 0.8, RULE),
        ("LEFTPADDING", (0, 0), (-1, -1), 2),
    ]))
    return t


def fig(path, width, caption, s, max_h=None):
    from PIL import Image as PILImage
    iw, ih = PILImage.open(path).size
    h = width * ih / iw
    if max_h and h > max_h:
        h = max_h
        width = h * iw / ih
    return KeepTogether([Image(str(path), width=width, height=h),
                         Paragraph(caption, s["caption"])])


def build(result_dir: str, out_pdf: str):
    R = pathlib.Path(result_dir)
    figs = R / "figures"
    s = styles()
    D = json.load(open(R / "top10_designs.json"))
    log = json.load(open(R / "cascade_log.json"))
    params = json.load(open(R / "ga_params.json"))
    hist = json.load(open(R / "ga_history.json"))
    retro = {r["design_id"]: r for r in json.load(open(R / "retro_results.json"))}
    meta = json.load(open(R / "run_meta.json"))

    avail = PAGE_W - 2 * MARGIN
    doc = BaseDocTemplate(out_pdf, pagesize=A4,
                          leftMargin=MARGIN, rightMargin=MARGIN,
                          topMargin=MARGIN, bottomMargin=20 * mm,
                          title="De novo small-molecule design campaign: GSK3B",
                          author="Computational design report")
    frame = Frame(MARGIN, 20 * mm, avail, PAGE_H - MARGIN - 20 * mm, id="main")
    doc.addPageTemplates([PageTemplate(id="all", frames=[frame], onPage=_footer)])

    E = []
    P = lambda t, st="body": Paragraph(t, s[st])

    # ---------------- header ----------------
    slot = Table([[Paragraph(
        "[ BRANDING SLOT ]<br/>logo / wordmark", s["slot"])]],
        colWidths=[38 * mm], rowHeights=[13 * mm])
    slot.setStyle(TableStyle([
        ("BOX", (0, 0), (-1, -1), 0.6, HAIR),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("ALIGN", (0, 0), (-1, -1), "CENTER"),
    ]))
    head = Table([[Paragraph(
        "De novo small-molecule design campaign", s["title"]), slot]],
        colWidths=[avail - 40 * mm, 40 * mm])
    head.setStyle(TableStyle([("VALIGN", (0, 0), (-1, -1), "TOP"),
                              ("LEFTPADDING", (0, 0), (-1, -1), 0),
                              ("RIGHTPADDING", (0, 0), (-1, -1), 0)]))
    E += [head]
    E += [P("Target: glycogen synthase kinase-3 beta (GSK3B, UniProt P49841, "
            "ChEMBL CHEMBL262)", "subtitle")]
    E += [P(f"Generated {date.today().isoformat()} &nbsp;·&nbsp; graph genetic "
            f"algorithm &nbsp;·&nbsp; {meta['n_unique']:,} molecules evaluated "
            f"&nbsp;·&nbsp; 10 designs selected", "subtitle")]
    E += [Spacer(1, 7)]

    E += [stat_row([
        (f"{meta['n_unique']:,}", "unique molecules<br/>generated &amp; scored"),
        (f"{log[0]['out']:,}", "valid, un-gated"),
        (f"{log[-1]['out']:,}", "passed full<br/>filter cascade"),
        ("10", "designs selected<br/>(distinct chemotypes)"),
        (f"{meta['n_retro_solved']}/10", "with a proposed<br/>synthetic route"),
    ], avail)]
    E += [Spacer(1, 4)]

    E += [Paragraph(
        "<b>What this report is.</b> Ten computationally designed molecules, "
        "generated de novo and ranked by a surrogate activity model combined with "
        "drug-likeness and synthetic-accessibility heuristics. <b>None of these "
        "compounds has been synthesised or tested.</b> Every activity number here "
        "is a model prediction, not a measurement, and the model is the same kind "
        "of object the generator was optimising against - see Limitations before "
        "acting on any ranking.", s["callout"])]

    # ---------------- introduction ----------------
    E += [P("1. Introduction", "h1")]
    E += [P(
        "GSK3B is a constitutively active serine/threonine kinase that sits at the "
        "intersection of Wnt/beta-catenin signalling, insulin signalling and "
        "microtubule regulation. Its dysregulation is implicated in the principal "
        "hallmarks of Alzheimer's disease - it phosphorylates tau, modulates "
        "amyloid-beta production, and affects neurogenesis and synaptic function - "
        "which has made it a long-standing and heavily prosecuted "
        "central-nervous-system target [6]. That history cuts both ways for a "
        "design exercise: the target is well validated and richly annotated, so a "
        "data-driven generator has plenty to learn from, but the accessible "
        "chemical space around it has also been explored extensively by medicinal "
        "chemists, so genuinely novel chemotypes are correspondingly harder to "
        "reach.")]
    E += [P(
        "This campaign asks a narrow, tractable question: starting from known "
        "potent GSK3B binders, can a graph genetic algorithm produce molecules "
        "that (a) score well on a pretrained GSK3B activity surrogate, (b) are "
        "drug-like and synthetically plausible, and (c) are structurally distinct "
        "from the compounds it started from? The deliverable is a ranked shortlist "
        "of ten designs with properties, novelty measurements, structural alerts "
        "and proposed retrosynthetic routes - a starting point for triage by a "
        "medicinal chemist, not a set of conclusions.")]

    # ---------------- methods ----------------
    E += [P("2. Methods", "h1")]
    E += [P("2.1 Reference and seed data", "h2")]
    E += [P(
        f"Known actives were retrieved from the ChEMBL REST API (release "
        f"<b>{meta['chembl_release']}</b>, retrieved {meta['retrieved']}) for "
        f"target CHEMBL262, restricted to IC50/Ki/Kd measurements with "
        f"pChEMBL &ge; 7 (that is, 100 nM or more potent). This returned "
        f"{meta['n_activity_records']:,} activity records spanning "
        f"<b>{meta['n_reference_actives']:,} unique compounds</b>, each represented "
        f"by its most potent measurement. This set served two distinct roles: the "
        f"source of the GA's starting population, and the reference against which "
        f"novelty was later measured. The starting population was the "
        f"{params['pop_size']} most structurally diverse members, chosen by MaxMin "
        f"picking over Morgan fingerprints so the GA did not begin inside a single "
        f"chemotype.")]

    E += [P("2.2 Generation", "h2")]
    E += [P(
        f"Molecules were generated by a graph genetic algorithm operating directly "
        f"on molecular graphs, following the operator family of Jensen [1] and the "
        f"graph-GA baseline of GuacaMol [2], implemented here on RDKit rather than "
        f"taken from either codebase. Crossover fragments both parents at a random "
        f"acyclic single bond and re-joins one fragment from each. Mutation applies "
        f"one local graph edit: atom append (from a library of common "
        f"medicinal-chemistry rings and substituents), element substitution, "
        f"terminal-atom deletion, bond-order change, or atom insertion into a bond. "
        f"Every child is sanitised and round-tripped through SMILES; anything "
        f"RDKit rejects is discarded. Selection is elitist truncation with "
        f"rank-weighted parent sampling. Parameters: population "
        f"{params['pop_size']}, {params['n_generations']} generations, "
        f"{params['offspring_per_gen']} offspring per generation, mutation rate "
        f"{params['mutation_rate']}, random seed {params['seed']}. Wall-clock "
        f"runtime was {params['runtime_s']:.0f} s on CPU.")]

    E += [P("2.3 Objective function", "h2")]
    E += [P(
        "Each molecule was scored as the geometric mean of three normalised terms - "
        "surrogate GSK3B activity, QED [3], and a synthetic-accessibility term "
        "(10 - SA_Score)/9 [4]. The geometric mean is deliberately conjunctive: a "
        "zero in any term zeroes the score, so the GA cannot trade away "
        "drug-likeness for predicted potency. Activity is the pretrained GSK3B "
        "oracle distributed with the Therapeutics Data Commons [8] - a random "
        "forest over ECFP features returning a value in [0,1]. A hard gate "
        "additionally zeroed any molecule outside MW 150-600, 10-50 heavy atoms, "
        "cLogP &le; 6, SA_Score &le; 6, zero formal charge, single fragment, "
        "no ring larger than 7, and elements restricted to C/N/O/S/F/Cl/Br. "
        "Without this gate a GA optimising a surrogate reliably drifts into "
        "high-scoring but unmakeable space.")]

    E += [P("2.4 Filtering and selection", "h2")]
    E += [P(
        "All molecules seen during the run - not only the final population - formed "
        "the selection pool. The cascade applied, in order: validity, uniqueness "
        "and gate survival; novelty (maximum Tanimoto &lt; 0.40 against every one "
        "of the reference actives, Morgan radius 2, 2048 bits); activity &ge; 0.50 "
        "and QED &ge; 0.60; freedom from PAINS A/B/C substructures [5]; ring "
        "sanity; and SA_Score &le; 4.5. "
        "Selecting the ten highest-scoring survivors directly produced ten "
        "decorations of a single core, so final selection instead clustered the "
        "survivors (Butina, Tanimoto 0.4) and took the best-scoring representative "
        "of each of the ten highest-scoring clusters. Bemis-Murcko scaffold "
        "uniqueness was tried first and rejected as too weak a criterion: swapping "
        "a pendant ring changes the Murcko scaffold while leaving the recognition "
        "motif intact.")]

    E += [P("2.5 Retrosynthesis and structural alerts", "h2")]
    E += [P(
        "Routes were proposed with AiZynthFinder [7] using the public USPTO "
        "template library, expansion and filter policies, and the ZINC stock, with "
        "a 120 s search budget per molecule. Designs were additionally screened "
        "against the Brenk reactive/unsuitable-group alerts and an explicit list of "
        "Michael-acceptor SMARTS; these were recorded as annotations rather than "
        "used as filters, because several validated GSK3B chemotypes are themselves "
        "maleimides.")]


    # ---------------- results ----------------
    E += [P("3. Results", "h1")]
    E += [P("3.1 Convergence", "h2")]
    first, last = hist[0], hist[-1]
    E += [P(
        f"The GA converged smoothly. Population mean score rose from "
        f"{first['mean']:.3f} to {last['mean']:.3f} over "
        f"{params['n_generations']} generations while the best individual improved "
        f"from {first['best']:.3f} to {last['best']:.3f}; the best score among "
        f"molecules that were also novel by the Tanimoto &lt; 0.4 criterion rose "
        f"from {first['best_novel']:.3f} to {last['best_novel']:.3f}. The best-score "
        f"curve is a staircase, which is expected under elitist truncation. The "
        f"mean curve flattening after roughly generation 25 indicates the "
        f"population had largely converged; the run was not extended because "
        f"additional generations were producing decorations of chemotypes already "
        f"present rather than new ones.")]
    E += [fig(figs / "fig1_convergence.png", avail, "<b>Figure 1.</b> Objective "
              "score by generation. The novel-only curve sits persistently below "
              "the unconstrained best, which is the central tension of this "
              "campaign: the surrogate rewards resemblance to known actives, while "
              "the novelty criterion penalises it.", s)]

    E += [P("3.2 Filter attrition", "h2")]
    E += [P(
        f"Of {log[0]['in']:,} unique molecules generated, {log[0]['out']:,} were "
        f"valid, unique and within the property gate. Novelty was by far the "
        f"harshest filter, removing {log[1]['removed']:,} molecules "
        f"({100*log[1]['removed']/log[1]['in']:.0f}% of those entering it) - the "
        f"generator spends most of its effort near its starting material. The "
        f"combined activity and QED gate removed a further "
        f"{log[2]['removed']:,}. PAINS, ring sanity and synthetic accessibility "
        f"together removed only {log[3]['removed']+log[4]['removed']+log[5]['removed']}, "
        f"which says less about the designs' quality than about the property gate "
        f"having already excluded most of what those filters target. "
        f"{log[-1]['out']} molecules survived the full cascade.")]
    E += [fig(figs / "fig3_cascade.png", avail, "<b>Figure 2.</b> Molecules "
              "remaining after each cascade stage, applied in the order shown.", s)]

    E += [P("3.3 Novelty", "h2")]
    nd = meta["novelty_bands"]
    E += [P(
        f"Novelty here is a measured quantity, and the measurement is sobering. "
        f"Across the whole generated set the modal similarity to the nearest known "
        f"active is around 0.45-0.50. Among the {log[-1]['out']} cascade survivors "
        f"the distribution is compressed hard against the 0.40 cutoff: "
        f"{nd['0.35-0.40']} of {log[-1]['out']} sit between 0.35 and 0.40, "
        f"{nd['0.30-0.35']} between 0.30 and 0.35, {nd['0.20-0.30']} between 0.20 "
        f"and 0.30, and <b>none below 0.20</b>. The ten selected designs span "
        f"{meta['sim_min']:.2f}-{meta['sim_max']:.2f}. These molecules are "
        f"therefore novel only in the specific, threshold-relative sense defined in "
        f"Methods. They are not new chemotypes; they are new points near the edge "
        f"of a well-populated region. A different cutoff would have produced a "
        f"different answer, and that sensitivity is itself the finding.")]
    E += [fig(figs / "fig5_novelty.png", avail, "<b>Figure 3.</b> Maximum Tanimoto "
              "similarity of each generated molecule to the nearest of the "
              f"{meta['n_reference_actives']:,} reference actives. The spike at 1.0 "
              "is the seed molecules, which are themselves known actives.", s)]

    E += [P("3.4 Property profile", "h2")]
    E += [P(
        "The surviving designs are systematically smaller and less polar than the "
        "known actives they were derived from - the generated molecular-weight "
        "distribution centres near 250-350 Da against roughly 350-450 Da for the "
        "ChEMBL set, with correspondingly lower TPSA. This is a known bias of "
        "scoring functions containing a QED term, since QED peaks in the middle of "
        "the property ranges rather than at the upper end where many optimised "
        "kinase inhibitors sit. Synthetic accessibility overlaps well with the "
        "reference set. The practical implication is that these designs have "
        "headroom to grow during optimisation, which is a reasonable place for a "
        "starting point to be.")]
    E += [fig(figs / "fig4_properties.png", avail, "<b>Figure 4.</b> Property "
              "distributions for cascade survivors (filled) against a random "
              "400-compound sample of the ChEMBL reference actives (outline).", s)]

    E += [P("3.5 Selected designs", "h2")]
    E += [P(
        f"Ten designs were selected as the best-scoring representatives of the ten "
        f"highest-scoring chemotype clusters (of {meta['n_clusters']} clusters among "
        f"the survivors). Chemotype-aware selection cost roughly 0.03 in mean "
        f"objective score relative to taking the top ten outright "
        f"({meta['mean_score_selected']:.3f} against "
        f"{meta['mean_score_naive']:.3f}) and bought a portfolio spanning "
        f"recognisably different recognition motifs rather than ten variations on "
        f"one.")]
    E += [fig(figs / "fig6_top10.png", avail * 0.82,
              "<b>Figure 5.</b> The ten selected designs, one per chemotype "
              "cluster. Scores, properties and alerts are in Table 1; they are "
              "surrogate-model outputs and heuristics, not measurements.", s)]

    # table
    rows = [[Paragraph(h, s["cellh"]) for h in
             ["ID", "Score", "Act.", "QED", "SA", "Sim.", "MW", "cLogP", "TPSA",
              "Route", "Alerts"]]]
    for d in D:
        r = retro.get(d["design_id"], {})
        n_steps = r.get("n_steps_best")
        route = (f"{n_steps} step{'' if n_steps == 1 else 's'}"
                 if r.get("solved") else "unsolved")
        alerts = ", ".join(d.get("reactive_motifs", [])) or "-"
        alerts = alerts.replace("michael_acceptor", "Michael").replace("_", " ")
        rows.append([Paragraph(x, s["cell"]) for x in [
            d["design_id"].replace("GSK3B-DN-", "DN-"),
            f"{d['score']:.3f}", f"{d['activity']:.2f}", f"{d['qed']:.2f}",
            f"{d['sa']:.2f}", f"{d['max_sim_to_known']:.2f}", f"{d['mw']:.0f}",
            f"{d['logp']:.1f}", f"{d['tpsa']:.0f}", route, alerts]])
    t = Table(rows, colWidths=[14 * mm, 12 * mm, 10 * mm, 10 * mm, 9 * mm,
                               10 * mm, 11 * mm, 12 * mm, 11 * mm, 15 * mm,
                               avail - 114 * mm], repeatRows=1)
    t.setStyle(TableStyle([
        ("LINEBELOW", (0, 0), (-1, 0), 0.8, RULE),
        ("LINEBELOW", (0, 1), (-1, -2), 0.3, HAIR),
        ("LINEBELOW", (0, -1), (-1, -1), 0.8, RULE),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("TOPPADDING", (0, 0), (-1, -1), 3.5),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 3.5),
        ("LEFTPADDING", (0, 0), (-1, -1), 3),
    ]))
    E += [KeepTogether([t, Paragraph(
        "<b>Table 1.</b> Selected designs. Act. = surrogate GSK3B oracle output; "
        "Sim. = maximum Tanimoto to any reference active; Route = shortest "
        "AiZynthFinder route found within the 120 s budget. Full SMILES, scaffolds "
        "and nearest-neighbour references are in the accompanying CSV.",
        s["caption"])])]

    E += [P("3.6 Retrosynthesis", "h2")]
    solved = [d for d in D if retro.get(d["design_id"], {}).get("solved")]
    E += [P(
        f"AiZynthFinder returned at least one fully-solved route (all leaves in the "
        f"ZINC stock) for <b>{len(solved)} of the 10</b> designs, with shortest "
        f"routes of 1-6 steps. Several proposals are chemically conventional and "
        f"plausible on their face - DN-01 resolves to a Suzuki coupling of an "
        f"aryl boronic acid with a bromo-indazole, and DN-08 to an SNAr between a "
        f"methylthio-pyrimidine and an aniline. The three unsolved cases "
        f"(DN-02, DN-06, DN-07) returned partial trees whose leaves were not all "
        f"purchasable within the budget; that is a statement about the template "
        f"library and the search budget, not proof that the molecules are "
        f"inaccessible. All routes are template-based proposals and carry no "
        f"information about yield, selectivity or protecting-group strategy.")]

    # ---------------- conclusions ----------------
    E += [P("4. Conclusions", "h1")]
    for txt in [
        f"<b>The pipeline works end to end and is cheap.</b> {meta['n_unique']:,} "
        f"molecules were generated, scored on three objectives, filtered, "
        f"clustered and routed in under {int(params['runtime_s']//60)+3} minutes of "
        f"CPU time, with a fixed seed and reproducible outputs.",
        "<b>Optimisation succeeded; novelty largely did not.</b> The objective "
        "improved substantially and the designs are drug-like and synthetically "
        "tractable. But no survivor fell below Tanimoto 0.20 of a known active, and "
        "the survivors pile up immediately beneath the 0.40 cutoff. The honest "
        "reading is that the GA rediscovered and recombined established GSK3B "
        "pharmacophore families - maleimides, anilinopyrimidines, azaindoles, "
        "indazoles - rather than inventing new ones.",
        "<b>That rediscovery is itself a weak positive control.</b> A generator "
        "optimising a GSK3B surrogate converging onto chemotypes that real GSK3B "
        "programmes have converged on suggests the surrogate is tracking something "
        "real. It is not evidence that the specific molecules are active.",
        "<b>Chemotype collapse is the default and must be designed against.</b> "
        "The naive top ten were ten decorations of one core. Only explicit "
        "clustering produced a portfolio. Any campaign reporting a top-N list "
        "without a diversity criterion is probably reporting one bet N times.",
    ]:
        E += [P(txt)]

    # ---------------- limitations ----------------
    E += [P("5. Limitations", "h1")]
    E += [Paragraph(
        "These limitations are not boilerplate; several of them materially affect "
        "how the ranking above should be read.", s["callout"])]
    lim = [
        ("The activity score is a surrogate, and the generator optimised against "
         "it directly.", "This is the dominant caveat. The GSK3B oracle is a "
         "random forest over fingerprints trained on public data; a genetic "
         "algorithm given thousands of queries will find inputs that score highly "
         "partly because of model artefacts. High oracle output is evidence of "
         "resemblance to the model's training actives, not of binding. No docking, "
         "free-energy or structural validation was performed."),
        ("Novelty is threshold-relative and marginal.",
         "The 0.40 Tanimoto cutoff is a convention, not a property of chemistry. "
         "Survivors cluster at 0.35-0.40, so modest changes to the cutoff, the "
         "fingerprint or the radius would materially change which molecules qualify."),
        ("Six of ten designs carry structural alerts; five are maleimides.",
         "Maleimides and the thiophene-dione in DN-07 are Michael acceptors with "
         "covalent-reactivity and promiscuity liabilities. PAINS filtering does not "
         "catch them. They were retained because validated GSK3B inhibitors "
         "(SB-216763, bisindolylmaleimides) share the motif - but a maleimide-rich "
         "shortlist is a real liability, not a neutral observation, and DN-01, "
         "DN-06, DN-08 and DN-09 are the alert-free designs."),
        ("QED and SA_Score are heuristics.",
         "QED encodes an aggregate notion of drug-likeness that penalises the "
         "larger, more polar molecules common among optimised kinase inhibitors, "
         "which is visible in Figure 4. SA_Score is a fragment-frequency proxy and "
         "systematically flatters molecules built from common fragments - which is "
         "precisely what this generator builds."),
        ("Retrosynthetic routes are proposals from USPTO templates.",
         "A solved route means the search reached purchasable leaves. It implies "
         "nothing about yield, selectivity, scale or protecting groups, and the "
         "template library reflects historically published reactions."),
        ("Single run, single seed.",
         "One GA run at seed 42. No replicate runs, so the variance of the outcome "
         "is unmeasured and the specific ten molecules should be treated as one "
         "sample from a distribution, not as the optimum."),
        ("Selectivity was not assessed at all.",
         "GSK3A shares close to identical ATP-site architecture with GSK3B, and "
         "nothing in this pipeline addresses kinome-wide selectivity, cellular "
         "activity, permeability, CNS penetration or toxicity."),
    ]
    for head_txt, body in lim:
        E += [Paragraph(f"<b>{head_txt}</b> {body}", s["small"])]

    # ---------------- next steps ----------------
    E += [P("6. Next steps", "h1")]
    nxt = [
        "<b>Break the surrogate's grip.</b> Re-score the 248 survivors with an "
        "orthogonal method the GA never saw - docking into a GSK3B ATP-site "
        "structure (several are available, e.g. 1Q5K, 6B8J) or a co-folding model - "
        "and keep only molecules that both methods like. Disagreement between the "
        "two is the most informative signal available without experiments.",
        "<b>Run replicates and report variance.</b> Five to ten seeds, then report "
        "how often each chemotype recurs. Chemotypes that survive across seeds are "
        "worth more than the single best score from one run.",
        "<b>Put novelty in the objective, not just the filter.</b> Adding an "
        "explicit novelty or diversity reward would let the GA trade predicted "
        "potency for genuine structural departure, instead of being filtered after "
        "the fact - which is what produced the pile-up at the cutoff.",
        "<b>Triage the maleimides deliberately.</b> Decide as a programme whether "
        "covalent/reactive chemotypes are in scope. If not, re-run with a Michael-"
        "acceptor exclusion in the gate; the alert-free designs (DN-01, DN-06, "
        "DN-08, DN-09) indicate the pipeline can produce them.",
        "<b>Add a GSK3A counter-screen and an ADMET pass</b> before any synthesis "
        "decision, and have a medicinal chemist review the proposed routes - "
        "template plausibility is not a synthesis plan.",
    ]
    for n in nxt:
        E += [P(n)]

    # ---------------- references ----------------
    E += [P("7. References", "h1")]
    refs = [
        "Jensen JH. A graph-based genetic algorithm and generative model/Monte "
        "Carlo tree search for the exploration of chemical space. "
        "<i>Chem Sci</i> 2019;10(12):3567-3572. doi:10.1039/c8sc05372c",
        "Brown N, Fiscato M, Segler MHS, Vaucher AC. GuacaMol: Benchmarking Models "
        "for de Novo Molecular Design. <i>J Chem Inf Model</i> "
        "2019;59(3):1096-1108. doi:10.1021/acs.jcim.8b00839",
        "Bickerton GR, Paolini GV, Besnard J, Muresan S, Hopkins AL. Quantifying "
        "the chemical beauty of drugs. <i>Nat Chem</i> 2012;4(2):90-98. "
        "doi:10.1038/nchem.1243",
        "Ertl P, Schuffenhauer A. Estimation of synthetic accessibility score of "
        "drug-like molecules based on molecular complexity and fragment "
        "contributions. <i>J Cheminform</i> 2009;1(1):8. doi:10.1186/1758-2946-1-8",
        "Baell JB, Holloway GA. New substructure filters for removal of pan assay "
        "interference compounds (PAINS) from screening libraries and for their "
        "exclusion in bioassays. <i>J Med Chem</i> 2010;53(7):2719-2740. "
        "doi:10.1021/jm901137j",
        "Lauretti E, Dincer O, Pratic&ograve; D. Glycogen synthase kinase-3 "
        "signaling in Alzheimer's disease. <i>Biochim Biophys Acta Mol Cell Res</i> "
        "2020;1867(5):118664. doi:10.1016/j.bbamcr.2020.118664",
        "Genheden S, Thakkar A, Chadimov&aacute; V, Reymond JL, Engkvist O, "
        "Bjerrum E. AiZynthFinder: a fast, robust and flexible open-source software "
        "for retrosynthetic planning. <i>J Cheminform</i> 2020;12(1):70. "
        "doi:10.1186/s13321-020-00472-1",
        "Huang K, Fu T, Gao W, et al. Therapeutics Data Commons: machine learning "
        "datasets and tasks for drug discovery and development. "
        "tdcommons.ai (software and pretrained GSK3B oracle; accessed "
        f"{meta['retrieved']}).",
        f"ChEMBL database, release {meta['chembl_release']}, EMBL-EBI. "
        f"Target CHEMBL262 bioactivity data retrieved via the ChEMBL REST API on "
        f"{meta['retrieved']}. ebi.ac.uk/chembl",
        "RDKit: Open-source cheminformatics. rdkit.org",
    ]
    for i, r in enumerate(refs, 1):
        E += [Paragraph(f"[{i}]&nbsp;&nbsp;{r}", s["ref"])]

    E += [Spacer(1, 10)]
    E += [Paragraph(
        "Literature references [1]-[7] were retrieved and verified via PubMed. "
        "Bioactivity data were retrieved from ChEMBL. All ten designs are "
        "computational proposals that have not been synthesised or assayed.",
        s["small"])]

    doc.build(E)
    return out_pdf
