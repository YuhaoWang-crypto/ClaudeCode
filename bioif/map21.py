"""
bioif.map21 -- the interface map's 21 cross-domain pairings, audited.

The blueprint lays out seven object classes and the 21 directed pairings
worth attempting between them, each with an evidence level: **A** verifiable
in quantified form under stated conditions, **B** needs paired experimental
calibration, **C** fit only for mechanistic hypothesis or candidate ranking.

This module is the audit of that map against what can actually be built. For
every pairing it records four things, and keeps them separate on purpose:

  claimed    the blueprint's evidence level
  data       whether a labelled source is reachable from this environment
  built      whether bioif implements the pairing as a real adapter
  verdict    what can honestly be done here, now

The distinction that matters is between **data reachable** and **paired
labels reachable**. Many pairings have abundant data on each side and almost
no compounds/genes/cells measured on BOTH, which is precisely why their
evidence level is B. An audit that only asked "is there a dataset" would
mark them green and be wrong.

Run: python3 -m bioif.map21
"""
from __future__ import annotations

from dataclasses import dataclass, field

# Object classes, in the blueprint's notation.
CLASSES = {
    "C": "compound / exposure",
    "G": "DNA, variant, CRISPR",
    "E": "chromatin / epigenome",
    "R": "transcript / splicing",
    "P": "protein, binding, PTM",
    "F": "cell phenotype / in-vitro assay",
    "D": "disease / donor / clinical",
}

#: Audit status for the "built" column.
BUILT = "built"              # a real adapter on real data in this repo
PARTIAL = "partial"          # infrastructure present, model is a stub
REACHABLE = "reachable"      # labelled data verified reachable, not built
BLOCKED = "blocked"          # no reachable paired labels from here
REFUSED = "refused"          # should not be built as a computational edge


@dataclass
class Pairing:
    key: str
    claimed: str                 # blueprint evidence level: A / A-B / B / B-C / C
    output: str                  # what the edge is supposed to produce
    breakpoint: str              # the dominant failure mode
    status: str
    source: str = ""             # the data source, if any
    paired_labels: str = ""      # what is measured on BOTH sides
    note: str = ""
    module: str = ""             # where it lives in bioif


# --------------------------------------------------------------------------
# The audit. `status` reflects what this container can actually reach, which
# was probed, not assumed.
# --------------------------------------------------------------------------

MAP: list[Pairing] = [
    Pairing("C->P", "A/B", "SMILES -> target potency / occupancy",
            "salt form, stereochemistry, assay type, Ki vs IC50 mixing",
            BUILT, "ChEMBL REST (live + committed snapshot)",
            "5,000 pChEMBL activities, 2,689 compounds, 206 assays on KRAS",
            "Readout-type separation (pIC50/pKi/pKd/pEC50 as distinct "
            "quantities) and assay-transfer conformal calibration both land "
            "here. Measured: assay identity moves potency a median 1.00 log.",
            "bioif/real/sources.py, conformal.py, transfer_adapter.py"),

    Pairing("C->F", "A", "compound -> reporter / viability readout",
            "mixed endpoints, curve QC, general-cytotoxicity false positives",
            BUILT, "Tox21 (MoleculeNet release), committed snapshot",
            "7,453 compounds x 12 endpoints; SR-p53 6,774 labelled (6.2% pos)",
            "Scaffold-split QSAR with label-conditional conformal. Measured: "
            "AUROC 0.718 / AP 0.286 on held-out scaffolds; SR-MMP kept as the "
            "cytotoxicity control the blueprint asks for.",
            "bioif/real/qsar.py"),

    Pairing("C->G", "B", "compound -> DNA adduct / mutation",
            "metabolism, S9, reactive intermediates, different damage endpoints",
            BUILT, "Hansen Ames benchmark N=6512, committed snapshot",
            "6,499 compounds with bacterial reverse-mutation labels",
            "Direct structure->Ames QSAR: AUROC 0.787 / AP 0.820 on held-out "
            "scaffolds. Note this is mutagenicity, NOT the adduct or Comet "
            "endpoint, which remain unbuildable from public data.",
            "bioif/real/qsar.py"),

    Pairing("F->G", "B", "reporter result -> mutagenicity",
            "a reporter is a stress signal, not a mutation count",
            BUILT, "Tox21 x Ames, joined on InChIKey skeleton",
            "2,064 compounds measured in BOTH; 1,908 with an SR-p53 label",
            "The blueprint's own `blocked` edge, now graded. Measured: p53+ "
            "compounds are 55.2% Ames+ vs 30.5% for p53- (RR 1.81, OR 2.81, "
            "Fisher p=5e-9) but sensitivity is only 0.13. Upgraded from "
            "blocked to inferred_association -- usable for ranking, useless "
            "as a screen. And a MEASURED p53 adds -0.000 AP on top of "
            "structure for predicting Ames.",
            "bioif/real/tox.py, chain_vs_direct.py"),

    Pairing("G->F", "A/B", "CRISPR KO -> viability / function",
            "KO and CRISPRi are not equivalent; copy-number artefacts; "
            "genetic background",
            PARTIAL, "DepMap portal API reachable (probed 200)",
            "gene effect x cell line, but the matrices are ~500 MB",
            "The residual_activity->fitness adapter exists with a declared "
            "calibration set, but its transform is a stub. Making it real "
            "needs a DepMap slice, which is a download-size problem rather "
            "than an availability one.",
            "bioif/adapters_demo.py (stub)"),

    Pairing("G->R", "A/B", "variant -> splicing / expression",
            "transcript version, gene dosage, off-target, editing efficiency",
            BLOCKED, "SpliceAI lookup API NOT reachable from here "
            "(tunnel closed); Ensembl REST also blocked",
            "MFASS / Vex-seq style minigene sets exist publicly",
            "The dPSI source in the stub chain stands in for this. UCSC "
            "sequence API IS reachable, so a from-scratch motif model is "
            "possible, but that would be building a model rather than "
            "interfacing one, and would not be competitive.",
            "bioif/adapters_demo.py (stub)"),

    Pairing("G->P", "B", "variant -> structure / stability / protein level",
            "isoform, folding, degradation, compensation",
            REACHABLE, "UniProt REST + EBI Proteins API reachable (200)",
            "domain annotations; no paired variant->protein-level labels",
            "Annotation-rule territory (domain disruption, LoF/NMD "
            "inference), which is what the plan document's Demo A node 2-3 "
            "did. Honest ceiling without paired proteomics: a rule engine, "
            "not an ML functional prediction.", ""),

    Pairing("R->P", "B", "mRNA -> protein abundance / activity",
            "translation, half-life, localisation, phosphorylation",
            PARTIAL, "ProteomicsDB / CPTAC endpoints not usable from here "
            "(probed 400)",
            "needs RNA and protein on the SAME samples",
            "This is the hop the variance attribution repeatedly names as "
            "dominant (~50% of endpoint variance in the stub chain). It is "
            "also the one with no reachable paired labels, so it stays a "
            "flagged bridge. That combination -- largest contributor, least "
            "verifiable -- is the single most important thing the audit "
            "surfaces.",
            "bioif/adapters_demo.py (stub, flagged bridge)"),

    Pairing("P->F", "B/C", "target occupancy / complex -> function",
            "inhibition direction, feedback, ubiquitination, duration",
            PARTIAL, "ChEMBL: a biochemical and a cellular assay on the same "
            "compounds",
            "CHEMBL4357258 (nucleotide exchange) x CHEMBL4357259 (pERK in "
            "MIA PaCa-2): 22 shared compounds",
            "A genuine biochemical->cellular transfer with paired labels, "
            "already inside the conformal assay-transfer machinery (bias "
            "-0.02, sd 0.32, r=0.95 on n=22). Small but real.",
            "bioif/real/conformal.py"),

    Pairing("C->R", "B", "drug -> expression delta",
            "dose/time, batch, cytotoxicity, unmeasured genes",
            BLOCKED, "LINCS / SigCom APIs not probed reachable from here",
            "L1000 signatures keyed by pert_id",
            "The plan document reports a 503 from SigCom in its own run, "
            "which is itself the interface-availability bottleneck. Not "
            "attempted here.", ""),

    Pairing("C->E", "B", "drug -> ATAC / histone marks",
            "time and cell-type dependence between binding and chromatin",
            BLOCKED, "", "needs paired dose-time ATAC/CUT&Tag",
            "No public paired set reachable; correctly a B-level edge "
            "requiring new experiments.", ""),

    Pairing("G->E", "B", "variant / CRISPR -> accessibility",
            "cis vs trans, allele specificity, feedback, repair side-effects",
            BLOCKED, "", "needs same-cell CRISPRi/a + ATAC",
            "Same situation as C->E.", ""),

    Pairing("E->R", "B", "enhancer / methylation -> gene expression",
            "co-localisation is not regulation; the distal target is not unique",
            BLOCKED, "UCSC sequence API reachable, but Enformer/Borzoi-class "
            "weights are a heavy GPU dependency",
            "enhancer-CRISPRi + RNA sets exist publicly",
            "Buildable in principle via a sequence model plus GTEx eQTL "
            "sign-concordance; days of work and a GPU, not hours.", ""),

    Pairing("E->P", "B/C", "chromatin change -> protein level / modification",
            "RNA-to-protein buffering, translation, PTM",
            BLOCKED, "", "needs paired ATAC/RNA/proteome time series", "", ""),

    Pairing("E->F", "B", "chromatin state -> differentiation / resistance",
            "state and function share latent causes",
            BLOCKED, "", "needs chromatin-regulator perturbation + phenotype",
            "", ""),

    Pairing("R->F", "B", "expression programme -> cytokine / viability",
            "expression is a state proxy; secretion and survival have their "
            "own kinetics",
            BLOCKED, "", "needs same-batch expression + functional readout",
            "", ""),

    Pairing("G->D", "B/C", "genetic signal -> disease target",
            "LD, locus-to-gene, directionality, tissue mismatch",
            BLOCKED, "Open Targets GraphQL endpoint reachable but rejected a "
            "bare GET (400); not pursued",
            "L2G scores exist publicly",
            "Ranking-level only. The blueprint's own warning stands: a GWAS "
            "edge is not a causal edge.", ""),

    Pairing("C->D", "C", "drug -> patient benefit / adverse effect",
            "PK/PD, immunity and microenvironment, selection bias",
            REFUSED, "", "",
            "Should not be built as a computational edge. The honest "
            "interface behaviour is REFUSE.", ""),
    Pairing("R->D", "C", "transcriptional signature -> patient outcome",
            "donor and platform drift, confounding, tissue composition",
            REFUSED, "", "", "As C->D.", ""),
    Pairing("P->D", "C", "protein / pathway -> clinical phenotype",
            "tissue distribution, safety window, polypharmacology, feedback",
            REFUSED, "", "", "As C->D.", ""),
    Pairing("E->D", "C", "epigenetic mark -> disease risk",
            "cell composition, reverse causation along disease course",
            REFUSED, "", "", "As C->D.", ""),
    Pairing("F->D", "C", "in-vitro assay -> patient safety / efficacy",
            "dose translation, exposure, species and microenvironment",
            REFUSED, "", "",
            "The largest context gap in the whole map. No contract fixes it; "
            "it is a study-design question. bioif's answer is a REFUSE, and "
            "that is the correct answer rather than a gap.", ""),
]

ORDER = {BUILT: 0, PARTIAL: 1, REACHABLE: 2, BLOCKED: 3, REFUSED: 4}


def by_status() -> dict[str, list[Pairing]]:
    out: dict[str, list[Pairing]] = {}
    for p in sorted(MAP, key=lambda q: (ORDER[q.status], q.key)):
        out.setdefault(p.status, []).append(p)
    return out


def report() -> str:
    g = by_status()
    L = ["The interface map's 21 pairings, audited against what is buildable",
         f"  {len(MAP)} rows: the map's 21 directed pairings, plus F->G, which",
         "  the blueprint raises in its own Demo B and marks `blocked`.",
         "  Object classes: " + ", ".join(f"{k}={v}" for k, v in CLASSES.items()),
         ""]
    L.append(f"  {'status':<11}{'n':>3}   pairings")
    for st in sorted(g, key=lambda s: ORDER[s]):
        L.append(f"  {st:<11}{len(g[st]):>3}   "
                 + " ".join(p.key for p in g[st]))
    L.append("")
    for st in sorted(g, key=lambda s: ORDER[s]):
        L.append("=" * 74)
        L.append(f"{st.upper()}  ({len(g[st])})")
        L.append("=" * 74)
        for p in g[st]:
            L.append(f"\n{p.key}   [blueprint level {p.claimed}]   {p.output}")
            L.append(f"   breakpoint     : {p.breakpoint}")
            if p.source:
                L.append(f"   data source    : {p.source}")
            if p.paired_labels:
                L.append(f"   paired labels  : {p.paired_labels}")
            if p.note:
                L.append("   finding        : " + p.note.replace(
                    "\n", "\n                    "))
            if p.module:
                L.append(f"   implemented in : {p.module}")
    L += ["", "=" * 74,
          "Two observations the audit produces that the map alone does not:",
          "",
          "1. Every pairing that is BUILT here is one where some public source",
          f"   happens to hold labels on BOTH sides. That is 4 of 21. The",
          "   bottleneck across the map is not models and not formats -- it is",
          "   paired measurement, which is exactly why the blueprint grades so",
          "   much of the map B.",
          "",
          "2. R->P is simultaneously the hop that contributes the most endpoint",
          "   variance in our chains (~50%) and the hop with no reachable",
          "   paired labels. The most load-bearing edge is the least",
          "   verifiable one. Any honest long-chain programme should buy that",
          "   measurement first; no modelling choice substitutes for it."]
    return "\n".join(L)


if __name__ == "__main__":
    print(report())
