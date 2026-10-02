"""
bioif.datasheet -- per-pairing deployment spec: inputs, reliable outputs,
throughput, and the binding limitation.

`map21.py` audits whether a pairing is *buildable*. This module answers the
operational questions that come next: **what must I supply per query, what do
I get back that I can trust, and how many queries per day?**

The organising finding is that **throughput is set by the input you must
supply, not by the model**. A pairing whose input is a string you already
have runs at screening scale; a pairing whose input is "a paired measurement
in the same cell line at the same dose and time" has the throughput of your
wet lab, however fast its model is. That is why the throughput column below
is a property of the *interface*, not of anyone's GPU.

Four throughput classes, with measured anchors where this repo has them:

  S  screening       input is a structure/sequence you already hold.
                     ✅ measured here: 10-13 M SMILES/hour on CPU for
                     featurize + predict (Morgan 2048 + logistic/copula).
  B  batch           needs a reference-dataset lookup, an external API, or
                     GPU inference. 10^3-10^5 queries/day.
  T  targeted        needs per-target setup: a structure, a cell-line panel,
                     a rate-limited API. 10^1-10^2/day.
  X  experiment      needs a NEW paired measurement per query. Throughput is
                     your laboratory's, not your cluster's.

Numbers are attributed: [here] measured in this repo, [rep] from the
independent strategy report's demos, [both] independently reproduced.

Run: python3 -m bioif.datasheet
"""
from __future__ import annotations

from dataclasses import dataclass

from . import map21

S, B, T, X = "S screening", "B batch", "T targeted", "X experiment"

#: Evidence grades, strongest first (the blueprint's ladder).
GRADES = ("measured", "inferred_association", "calibrated_prediction",
          "mechanistic_hypothesis", "none")


@dataclass(frozen=True)
class Deploy:
    key: str
    ref: str                 # the strategy report's pair id
    inputs: str              # what you must supply PER QUERY
    output: str              # what comes back that can be trusted
    grade: str
    throughput: str
    number: str              # the headline measured result, attributed
    limit: str               # the binding limitation


DEPLOY: list[Deploy] = [
    # ---------------- screening-capable, reliable -----------------------
    Deploy("C->P", "P01",
           "SMILES; the target's assay id you want the number expressed in",
           "pIC50 for ONE named assay, with a conformal interval and a "
           "refusal outside the calibrated assay pair",
           "calibrated_prediction", S,
           "[rep] scaffold-split Spearman 0.815 on 4,338 PARP1 compounds, "
           "but AD gating leaves 3.5% of a screening library usable. "
           "[here] assay identity alone moves the same compound's potency a "
           "median 1.00 log (10x); cross-assay conformal coverage collapses "
           "to 0.61 against a nominal 0.80, worst case 0.00",
           "The number is only meaningful per assay. Pooling assays is the "
           "default behaviour of every public source and is the single "
           "largest error here. Applicability-domain gating is what makes it "
           "usable and it discards ~96% of a library"),

    Deploy("C->F", "P05",
           "SMILES only",
           "P(assay endpoint active) for ONE named endpoint and protocol, as "
           "a conformal prediction set {0}/{1}/abstain",
           "calibrated_prediction", S,
           "[here] Tox21 SR-p53 scaffold-split AUROC 0.718 / AP 0.286; "
           "label-conditional conformal abstains on 54% of compounds. "
           "[rep] AUROC 0.765 / AP 0.191 on the same endpoint",
           "Endpoint-specific: a p53 reporter model says nothing about Comet "
           "or micronucleus. Minority-class conformal coverage falls to 0.832 "
           "against a nominal 0.90 BECAUSE a scaffold split deliberately "
           "breaks the exchangeability the guarantee needs"),

    Deploy("C->G", "P04",
           "SMILES only",
           "P(Ames positive)",
           "calibrated_prediction", S,
           "[here] Hansen N=6512 scaffold-split AUROC 0.787 / AP 0.820. "
           "[rep] AUROC 0.817 / AP 0.836",
           "Ames is not adduct formation, not Comet, not in-vivo "
           "genotoxicity. Compounds needing S9 metabolic activation are a "
           "known blind spot and are not flagged by the model"),

    Deploy("G->R", "P11",
           "Variant in HGVS + an explicit genome build; the MANE transcript",
           "Splice-disruption score, and a transcript consequence "
           "(LoF-via-NMD / in-frame skip) under a MANE contract",
           "calibrated_prediction", B,
           "[rep] SpliceAI AUROC 0.796 / Pangolin 0.866 / CADD 0.884 against "
           "3,912 functionally measured MPSA variants; canonical sites ~0.95 "
           "but DEEP EXONIC 0.62-0.78, and deep-exonic P/LP recall only "
           "0.018 at a 0.5 threshold. Hop 2-3 interface loss ~0 under MANE "
           "discipline",
           "The deep-exonic blind spot is structural, not a tuning problem. "
           "Build and schema drift are the live traps: a naive parser of "
           "Pangolin's DS_SG/DS_SL fields returns silent zeros against "
           "SpliceAI's DS_AG/AL/DG/DL. The public lookup API was NOT "
           "reachable from this container"),

    # ---------------- reliable but reference-bound -----------------------
    Deploy("G->P", "P08",
           "Variant + the protein isoform it is to be scored on",
           "Missense pathogenicity / functional-disruption score",
           "calibrated_prediction", B,
           "[rep] AlphaMissense AUROC 0.883, ESM-1v 0.873, PolyPhen-2 0.841, "
           "SIFT 0.774 on 103 classified ClinVar missense. Models disagree on "
           "21% of variants and the disagreements are ENRICHED for P/LP "
           "(18/113) -- usable as a review trigger",
           "Isoform mismatch silently drops variants: all 179 WT1 missense "
           "lack an AlphaMissense score because MANE encodes 497aa and AM's "
           "canonical is 449aa. Splice variants are outside the model's "
           "representation entirely"),

    Deploy("P->F", "P09",
           "A gene identifier and a cell-line panel",
           "Gene-effect / dependency score per cell line",
           "measured", B,
           "[rep] DepMap Chronos is a measured resource, but the median "
           "single-drug viability-to-dependency correlation across 1,201 "
           "compound-target pairs is 0.017 -- the endpoint, not the "
           "interface, is the noise source",
           "Gene-level aggregation destroys domain and isoform resolution. A "
           "dependency score is a statement about cell lines, not patients"),

    Deploy("F->G", "(Demo C)",
           "A measured reporter call, PLUS the cytotoxicity call -- the "
           "second one is not optional",
           "A likelihood-ratio update on mutagenicity risk, published ONLY "
           "in the non-cytotoxic stratum",
           "inferred_association", S,
           "[both] p53+ compounds are 55.2% Ames+ vs 30.5% ([here], "
           "n=1,908, OR 2.81, Fisher p=5e-9); [rep] OR 2.70 on n=1,565. "
           "Sensitivity only 0.13. Stratified: the association VANISHES "
           "among cytotoxic compounds ([here] OR 1.13 p=0.77; [rep] OR 0.98) "
           "and is strongest outside that stratum ([here] OR 2.90)",
           "Ranks candidates, cannot clear them: as a screen it misses 87% "
           "of mutagens. And it is conditional -- a marginal odds ratio "
           "averages two regimes, one of which carries no signal at all"),

    # ---------------- partial: resolution-dependent ---------------------
    Deploy("G->F", "P15",
           "An sgRNA plus its cut-site context; a matched screen in the SAME "
           "editing system",
           "Relative depletion ranking among guides WITHIN one gene",
           "calibrated_prediction (within gene) / mechanistic_hypothesis "
           "(between genes)", B,
           "[rep] within-gene, essential genes: Spearman(frameshift rate, "
           "depletion) = -0.142, p=4e-21; non-essential +0.007 (null). "
           "BETWEEN genes: mean frameshift vs Chronos rho=0.008, p=0.32 -- "
           "no signal. Cell-type-specific repair model gave NO gain "
           "(p=0.767)",
           "The edge connects at guide resolution and BREAKS at gene "
           "resolution, because essentiality varies far more than repair "
           "outcome does. Any chain through it must declare which resolution "
           "it claims. HCT116 could not be paired at all: only a Cas12a "
           "screen exists against a Cas9 repair model"),

    Deploy("C->R", "P02",
           "SMILES + dose + time + cell line; a reference perturbation "
           "signature library",
           "Direction of pathway-level change, not per-gene magnitudes",
           "inferred_association", T,
           "[rep] olaparib up-signature enriches p53 signalling at "
           "adj p=7.8e-12. But chemical and genetic perturbation signatures "
           "agree in the CROSS direction too (OR 6.9, p=1.8e-5) -- "
           "inhibition is not knockout. The data API returned 503; interface "
           "availability is itself the bottleneck",
           "Dose, time and cell line determine the signature and are usually "
           "dropped. 978 landmark genes are extrapolated to the "
           "transcriptome. Trust pathway direction; do not trust gene-level "
           "magnitude"),

    Deploy("R->P", "P06",
           "Matched RNA and protein on the SAME samples",
           "Nothing reliable without that pairing",
           "mechanistic_hypothesis", X,
           "[here] this hop contributes ~50% of endpoint variance in a 4-hop "
           "stub chain -- the largest single contributor -- and no public "
           "paired source was reachable to calibrate it",
           "THE most load-bearing edge in the whole map and the least "
           "verifiable one. Translation, half-life, localisation and PTM all "
           "live in this residual. Buy this measurement before buying any "
           "model"),

    Deploy("R->F", "P12",
           "An expression signature + a same-batch functional readout to "
           "calibrate against",
           "Pathway/mechanism class ranking only",
           "mechanistic_hypothesis", X,
           "[rep] batch effects and cohort shift dominate; signatures "
           "transfer poorly across platforms",
           "Expression is a state proxy; secretion and survival have their "
           "own kinetics. Correlation and causation are routinely conflated "
           "here"),

    # ---------------- blocked: no reachable paired labels ----------------
    Deploy("C->E", "P03", "Dose-time-matched ATAC/CUT&Tag in the right cell type",
           "Nothing", "none", X,
           "[rep] no LINCS-scale epigenomic perturbation library exists",
           "No public paired set. Correctly a B-grade edge requiring new "
           "experiments"),
    Deploy("G->E", "P13", "Same-cell CRISPRi/a plus ATAC", "Nothing", "none", X,
           "[rep] cell-type-specific caQTL data is scarce; LD confounds the "
           "causal variant",
           "cis/trans ambiguity, allele specificity and repair side-effects "
           "all unresolved without matched perturbation"),
    Deploy("E->R", "P10",
           "Genome sequence + matched chromatin tracks in the right cell type",
           "Direction of effect for PROXIMAL elements at best",
           "mechanistic_hypothesis", T,
           "[rep] distal-enhancer prediction is weak (Karollus 2023); "
           "cell-type match is a hard constraint",
           "Co-localisation is not regulation and the distal target gene is "
           "not unique. Needs a GPU and the weights are heavy"),
    Deploy("E->P", "P?", "Paired ATAC/RNA/proteome time series", "Nothing",
           "none", X, "[rep] no paired source",
           "Inherits the whole R->P problem and adds chromatin timing"),
    Deploy("E->F", "P14", "Chromatin-regulator perturbation + phenotype",
           "Nothing", "none", X,
           "[rep] the causal chain enhancer->gene->phenotype is long and the "
           "ABC contact approximation is coarse",
           "State and function share latent causes, so association is "
           "uninformative without perturbation"),
    Deploy("G->D", "P?", "Fine-mapped locus + matched-tissue QTL evidence",
           "A ranked candidate gene list", "mechanistic_hypothesis", T,
           "[rep] locus-to-gene scores exist publicly",
           "Ranking only. A GWAS edge is not a causal edge; LD, "
           "directionality and tissue mismatch all intervene"),

    # ---------------- refused -------------------------------------------
    Deploy("C->D", "P?", "n/a", "REFUSE", "none", X, "",
           "Not a computational edge. PK/PD, immunity, microenvironment and "
           "selection bias; the correct interface behaviour is a refusal"),
    Deploy("R->D", "P?", "n/a", "REFUSE", "none", X, "", "As C->D"),
    Deploy("P->D", "P?", "n/a", "REFUSE", "none", X, "", "As C->D"),
    Deploy("E->D", "P?", "n/a", "REFUSE", "none", X, "", "As C->D"),
    Deploy("F->D", "P?", "n/a", "REFUSE", "none", X, "",
           "The largest context gap in the map. No contract fixes dose "
           "translation, exposure and species; it is a study-design question"),
]

BY_KEY = {d.key: d for d in DEPLOY}


def screening_capable() -> list[Deploy]:
    """Pairings whose per-query input is something you already hold."""
    return [d for d in DEPLOY if d.throughput in (S, B)
            and d.grade not in ("none",)]


def experiment_bound() -> list[Deploy]:
    return [d for d in DEPLOY if d.throughput == X]


def report() -> str:
    L = ["Per-pairing deployment datasheet",
         "  input you must supply -> output you can trust -> throughput",
         "  [here] measured in this repo | [rep] from the independent report "
         "| [both] reproduced twice", ""]

    order = {S: 0, B: 1, T: 2, X: 3}
    for d in sorted(DEPLOY, key=lambda x: (order[x.throughput],
                                           GRADES.index(x.grade.split()[0])
                                           if x.grade.split()[0] in GRADES
                                           else 9, x.key)):
        st = map21.by_status()
        status = next((p.status for p in map21.MAP if p.key == d.key), "?")
        L.append("=" * 76)
        L.append(f"{d.key:<7} [{d.ref}]  {d.throughput:<12} grade: {d.grade}"
                 f"   audit: {status}")
        L.append(f"  INPUT   {d.inputs}")
        L.append(f"  OUTPUT  {d.output}")
        if d.number:
            L.append("  EVIDENCE")
            for line in _wrap(d.number, 70):
                L.append(f"          {line}")
        L.append("  LIMIT")
        for line in _wrap(d.limit, 70):
            L.append(f"          {line}")

    sc = screening_capable()
    xb = experiment_bound()
    L += ["", "=" * 76, "SUMMARY", "=" * 76,
          f"  screening- or batch-capable, with a usable output : "
          f"{len(sc)}/{len(DEPLOY)}",
          f"     {', '.join(d.key for d in sc)}",
          f"  experiment-bound (throughput = your wet lab)       : "
          f"{len(xb)}/{len(DEPLOY)}",
          f"     {', '.join(d.key for d in xb)}",
          "",
          "  The two sets barely overlap, and that is the whole answer to",
          "  'can this be high-throughput?'. Throughput is decided by the",
          "  INPUT the interface demands. Every screening-capable pairing is",
          "  one whose input is a string you already have -- which is also",
          "  why those pairings have large training sets and therefore",
          "  usable models. Reliability and throughput are not a trade-off",
          "  here; they have a common cause, which is cheap input."]
    return "\n".join(L)


CHAIN_EVIDENCE = [
    ("[here] §8.4  direct beats the 2-hop chain",
     "On 1,908 compounds with labels at BOTH ends, direct C->G beats "
     "C->F->G by +0.186 AP [0.166, 0.223] over 3 splits. Decomposed: only "
     "-0.056 is the C->F model's error; +0.279 is information the "
     "intermediate THROWS AWAY. A MEASURED intermediate adds +0.001 "
     "[-0.003, +0.005] on top of structure -- nothing."),
    ("[rep] Demo D  the 4-hop gain is not from depth",
     "Routing AM+ESM reaches AUROC 0.992 against 0.963 for serial-max. But "
     "their own ablation shows routing with AM ALONE also scores 0.963 -- "
     "identical to serial max. So the gain came from adding a SECOND, "
     "complementary model at one node, not from the extra hop. Their own "
     "conclusion: 'the value of routing is mechanistic interpretability, "
     "not accuracy.'"),
    ("[rep] Demo D  what depth DID buy: reach",
     "Deep-exonic P/LP recall 0.018 -> 0.544 (AM route) -> 0.860 (AM+ESM). "
     "WT1 P/LP recall 0.200 -> 0.900. These are variants a splice model "
     "structurally cannot see. That is genuine new capability, and it comes "
     "from covering a blind spot with a different modality."),
    ("[rep] Demo E  a chain breaks at the wrong resolution",
     "The G->F chain connects at guide level (within-gene rho=-0.142, "
     "p=4e-21) and is null at gene level (rho=0.008, p=0.32). Guide-level "
     "accuracy does not extrapolate to gene-level prediction."),
    ("[rep] Demo F  a generative edge loses to a linear one",
     "On 273 wholly unseen perturbations, a conditional flow-matching "
     "(diffusion-equivalent) model scores pearson_delta 0.123 against 0.387 "
     "for a linear co-expression kernel using the SAME conditioning "
     "information -- beaten 3x. Absolute pearson 0.96-0.99 for every method, "
     "confirming that metric is misleading."),
    ("[rep] Demo F / [here] §7  what distribution-valued edges DO buy",
     "Propagating 64 posterior samples instead of a point: 15.4% of "
     "condition x pathway direction calls are FLIPPED by the posterior "
     "majority, and 0% reach 90% confidence. Point-estimate chains "
     "systematically overstate confidence. [here] the same conclusion from "
     "the other side: a properly conditioned joint model does gain from a "
     "measured endpoint (+0.021 [0.015, 0.028], every seed positive) and "
     "still loses to a 66-parameter copula (+0.028)."),
    ("[here] §7  why the big model loses, measured",
     "Under joint-error selection the diffusion checkpoint BEATS the copula "
     "on calibration (0.862 vs 0.968) and loses on test (1.222 vs 1.070), "
     "degrading ~0.36 while the copula barely moves. Since the selection "
     "objective IS joint error, that cannot be a selection artefact: the "
     "correlation structure itself overfits at 8-18 parameters per observed "
     "label. The binding constraint is SAMPLE SIZE, not architecture."),
]


def chain_verdict() -> str:
    """Does chaining buy more biology? Assembled from both efforts."""
    L = ["Can chains cross short causal links and predict more biology?", ""]
    for title, body in CHAIN_EVIDENCE:
        L.append(f"  {title}")
        for line in _wrap(body, 70):
            L.append(f"      {line}")
        L.append("")
    L += [
        "VERDICT, and it is the same from two independent efforts:",
        "",
        "  Chains do NOT buy accuracy. Four demos across two groups say so",
        "  independently -- a direct model beating a 2-hop chain, a 4-hop",
        "  routing gain that ablation shows came from model fusion rather",
        "  than depth, a chain that is null at the resolution people would",
        "  deploy it at, and a generative edge beaten 3x by a linear one.",
        "",
        "  Chains DO buy three things, all measured:",
        "    1. REACH -- deep-exonic recall 0.018 -> 0.860 covers variants a",
        "       splice model cannot represent at all. This is the strongest",
        "       result in either effort and it is a capability, not an",
        "       accuracy gain.",
        "    2. INTERPRETABILITY -- a mechanism attached to each call",
        "       (LoF-via-NMD / in-frame skip / missense disruption).",
        "    3. HONEST UNCERTAINTY -- distribution-valued edges downgrade",
        "       15.4% of direction calls to coin flips that a point chain",
        "       reported as confident.",
        "",
        "  So the architecture that pays is WIDE, NOT DEEP: several",
        "  complementary models at ONE node, fused, with the fusion gated on",
        "  domain and the output distribution-valued. Adding nodes in series",
        "  multiplies unvalidated maps; adding models in parallel at a node",
        "  covers blind spots. Demo D is the existence proof and its own",
        "  ablation is the control.",
        "",
        "  Practical rule from the numbers: build a chain where end-to-end",
        "  labels are MISSING (most of this map), and expect a direct model",
        "  to win wherever they exist. Sell a chain on coverage and",
        "  auditability, never on accuracy.",
    ]
    return "\n".join(L)


def _wrap(text: str, width: int) -> list[str]:
    import textwrap
    return textwrap.wrap(text, width) or [""]


if __name__ == "__main__":
    import sys
    if "--chains" in sys.argv:
        print(chain_verdict())
    else:
        print(report())
        print()
        print(chain_verdict())
