"""Figures and the per-species markdown report.

Every claim in the generated report carries one of three labels:

  RIGOROUS      derived from database sequences by set arithmetic or alignment;
                as reliable for a cat as for a human.
  MEASURED      a property of the workflow that was actually measured here
                (panel size, applicability-domain distance, control outcomes).
  ILLUSTRATIVE  produced by the untrained surrogate scorer; plumbing only.

Nothing in the report is allowed to state an immunogenicity risk from an
ILLUSTRATIVE number.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Dict, List, Sequence

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .groove import DomainCall
from .insulin import Insulin, diff

LABELS = {
    "rigorous": "**RIGOROUS**",
    "measured": "**MEASURED**",
    "illustrative": "**⚠️ ILLUSTRATIVE**",
}


def write_csv(path: Path, rows: Sequence[Dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


# --------------------------------------------------------------------------
# figures
# --------------------------------------------------------------------------

def figure_difference_map(self_insulin: Insulin, products: Sequence[Insulin],
                          path: Path) -> Path:
    """Which residues of each product the recipient has never seen."""
    max_a = max(len(p.A) for p in products)
    max_b = max(len(p.B) for p in products)
    width = max_a + max_b
    n = len(products)

    def column(chain: str, pos: int) -> int:
        return (pos - 1) if chain == "A" else max_a + pos - 1

    fig, ax = plt.subplots(figsize=(0.30 * width + 4.5, 0.42 * n + 1.9))
    for row, prod in enumerate(products):
        y = n - 1 - row                                  # first product on top
        foreign = {(d.chain, d.position): d for d in diff(self_insulin, prod)}
        for chain, length in (("A", max_a), ("B", max_b)):
            seq = prod.chain(chain)
            for i in range(length):
                x = column(chain, i + 1)
                present = i < len(seq)
                d = foreign.get((chain, i + 1))
                if not present:
                    face, edge = "#f0f0f0", "#f0f0f0"    # residue absent (des-B30)
                elif d is None:
                    face, edge = "#ffffff", "#dddddd"    # identical to recipient
                else:
                    face, edge = "#c0392b", "#7b241c"    # foreign residue
                ax.add_patch(plt.Rectangle((x, y), 1, 1, facecolor=face,
                                           edgecolor=edge, linewidth=0.4))
                if present and d is not None:
                    ax.text(x + 0.5, y + 0.5, seq[i], ha="center", va="center",
                            fontsize=7, color="white", fontweight="bold")

    ax.set_xlim(0, width)
    ax.set_ylim(0, n)
    ax.set_yticks([n - 0.5 - i for i in range(n)])
    ax.set_yticklabels([p.name for p in products], fontsize=8)
    ticks = list(range(0, max_a, 2)) + list(range(max_a, width, 2))
    ax.set_xticks([t + 0.5 for t in ticks])
    ax.set_xticklabels(
        [f"A{t + 1}" for t in range(0, max_a, 2)] +
        [f"B{t - max_a + 1}" for t in range(max_a, width, 2)],
        fontsize=6.5, rotation=90)
    ax.axvline(max_a, color="k", lw=1.4)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(length=0)
    species = self_insulin.name.split()[0]
    ax.set_title(f"Residues foreign to a {species}\n"
                 "red = not present in this species' own insulin; "
                 "grey = residue absent from the product", fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)
    return path


def figure_applicability(loo: Sequence[float], calls: Sequence[DomainCall],
                         species: str, path: Path) -> Path:
    """Where the species panel sits relative to the model's own training space."""
    ident = [c.identity for c in calls]
    bins = np.linspace(min(min(loo), min(ident)) - 0.02, 1.005, 45)

    fig, ax = plt.subplots(figsize=(7.6, 4.4))
    ax.hist(loo, bins=bins, color="#9ecae1", edgecolor="white",
            label=f"NetMHCIIpan training molecules (n={len(loo)}),\n"
                  f"leave-one-out nearest neighbour")
    ax.set_ylabel("training molecules", color="#3182bd")
    ax.tick_params(axis="y", labelcolor="#3182bd")

    # Panel counts live on their own axis: a 20-molecule panel would be
    # invisible against a 1700-molecule training histogram.
    ax2 = ax.twinx()
    ax2.hist(ident, bins=bins, color="#d62728", alpha=0.75, edgecolor="white",
             label=f"{species} panel (n={len(ident)}),\ndistance to nearest training molecule")
    ax2.set_ylabel(f"{species} panel molecules", color="#d62728")
    ax2.tick_params(axis="y", labelcolor="#d62728")

    p5 = float(np.percentile(loo, 5))
    ax.axvline(p5, color="k", ls="--", lw=1.2)
    ax.annotate(f"training 5th percentile = {p5:.2f}", xy=(p5, ax.get_ylim()[1] * 0.92),
                xytext=(-6, 0), textcoords="offset points", ha="right", fontsize=8)
    ax.set_xlabel("groove-contact identity to nearest training molecule")
    ax.set_title(f"Applicability domain: is a {species} allele something the model has seen?",
                 fontsize=11)
    handles = ax.get_legend_handles_labels()[0] + ax2.get_legend_handles_labels()[0]
    labels = ax.get_legend_handles_labels()[1] + ax2.get_legend_handles_labels()[1]
    ax.legend(handles, labels, fontsize=7.5, loc="upper left")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)
    return path


def figure_core_landscape(drug: Insulin, coverage: Dict[str, List[int]], path: Path) -> Path:
    """Per-residue count of non-self 9-mer cores."""
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.2),
                             gridspec_kw={"width_ratios": [len(drug.A), len(drug.B)]})
    for ax, ch in zip(axes, ("A", "B")):
        vals = coverage[ch]
        ax.bar(range(1, len(vals) + 1), vals, color="#d62728", width=0.85)
        ax.set_title(f"{ch} chain", fontsize=10)
        ax.set_xlabel("residue")
        ax.set_xticks(range(1, len(vals) + 1, 2))
        ax.tick_params(labelsize=7)
    axes[0].set_ylabel("non-self cores spanning residue")
    fig.suptitle(f"{drug.name}: where the non-self binding cores are", fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)
    return path


# --------------------------------------------------------------------------
# markdown report
# --------------------------------------------------------------------------

def _table(headers: Sequence[str], rows: Sequence[Sequence[object]]) -> str:
    out = ["| " + " | ".join(headers) + " |",
           "|" + "|".join(["---"] * len(headers)) + "|"]
    for row in rows:
        out.append("| " + " | ".join(str(c) for c in row) + " |")
    return "\n".join(out)


def build_report(ctx: dict) -> str:
    species = ctx["species"]
    common = ctx["common_name"]
    self_ins = ctx["self_insulin"]
    parts: List[str] = []

    parts.append(f"# MHC-II immunogenicity workflow — {common}\n")
    parts.append(f"*Generated by `vetimmuno`. Backend: `{ctx['backend']}`. "
                 f"Every number below is labelled {LABELS['rigorous']}, "
                 f"{LABELS['measured']} or {LABELS['illustrative']}.*\n")

    # ---- 1. sequence layer
    parts.append("## 1. What is actually foreign to this animal\n")
    parts.append(f"{LABELS['rigorous']} — mature A/B chains taken from the UniProt "
                 f"records' own chain features; comparison is positional (insulin "
                 f"chains are colinear across mammals).\n")
    parts.append(f"Recipient self-antigen: **{self_ins.name}** ({self_ins.source})\n")
    rows = []
    for burden in ctx["burdens"]:
        rows.append([burden["product"], burden["n_differences"],
                     ", ".join(burden["differences"]) or "—",
                     burden["n_neo_cores"], f"{burden['frac_neo_cores']:.2f}"])
    parts.append(_table(
        ["product", "foreign residues", "positions", "non-self 9-mer cores",
         "fraction of all cores"], rows))
    parts.append("")
    parts.append("A *non-self core* is a 9-mer register in the administered insulin that "
                 "does not occur anywhere in this animal's own insulin. Cores the animal "
                 "already carries are removed by central tolerance and are dropped here — "
                 "that filter is set arithmetic, and it is exactly as reliable for a cat "
                 "as for a human.\n")

    # ---- 2. panel
    parts.append("## 2. The MHC-II panel that exists for this species\n")
    parts.append(f"{LABELS['measured']}\n")
    parts.append(_table(
        ["locus", "records considered", "rejected", "molecules kept",
         "unique grooves", "curated nomenclature", "source"],
        [[s["locus"], s["records_considered"], s["records_rejected"], s["molecules"],
          s["unique_grooves"], "yes" if s["curated_nomenclature"] else "**no**", s["source"]]
         for s in ctx["panel_summaries"]]))
    parts.append("")
    parts.append("> There is no allele-**frequency** database for either DLA or FLA, and "
                 "IEDB's population-coverage tool is human-only. A panel here is a "
                 "diversity sample of known sequences — it cannot be turned into a "
                 "\"% of the population covered\" number the way an HLA panel can.\n")

    # ---- 3. applicability domain
    ad = ctx["domain_summary"]
    parts.append("## 3. Is a pan-specific predictor even allowed to score these molecules?\n")
    parts.append(f"{LABELS['measured']}\n")
    parts.append(
        f"NetMHCIIpan-4.3 is trained on human HLA-DR/DQ/DP, mouse H-2 and bovine "
        f"BoLA-DRB3. We measured how far each {species} molecule sits from that "
        f"training space, in the same representation the model uses (groove-contact "
        f"residues), and calibrated the scale against the training set's own "
        f"leave-one-out nearest-neighbour distribution.\n")
    parts.append(_table(
        ["locus", "training molecules (unique grooves)", "training NN identity: median",
         "training NN identity: 5th pct", f"{species} panel: median",
         f"{species} panel: range", "out-of-domain", "marginal"],
        [[r["locus"], r["training_molecules"], f"{r['train_median']:.3f}",
          f"{r['train_p5']:.3f}", f"{r['panel_median']:.3f}",
          f"{r['panel_min']:.3f} – {r['panel_max']:.3f}",
          f"{r['n_out']}/{r['n_total']}", f"{r['n_marginal']}/{r['n_total']}"]
         for r in ad["rows"]]))
    parts.append("")
    parts.append("Each locus is compared against its own training pool, so the rows are "
                 "not on a common scale and are not pooled. \"Training NN identity\" is "
                 "the leave-one-out nearest-neighbour identity *within* the training set "
                 "— the calibration curve everything else is read against.\n")
    parts.append(ad["verdict_text"] + "\n")

    # ---- 4. binding predictions
    parts.append("## 4. Binding predictions\n")
    if ctx["backend_trained"]:
        parts.append(f"{LABELS['measured']} — produced by `{ctx['backend']}` in "
                     "custom-molecule mode. Read together with section 3.\n")
    else:
        parts.append(
            f"{LABELS['illustrative']} — produced by the untrained "
            "pocket-complementarity surrogate, which exists so this pipeline can run "
            "without a licensed NetMHCIIpan install. **These numbers rank peptides "
            "plausibly and predict nothing.** No risk statement in this report is "
            "derived from them. Install NetMHCIIpan-4.3 and re-run with "
            "`--backend netmhciipan` to make section 4 evidential.\n\n"
            "It also reads only the four anchor positions, so cores differing "
            "solely at a non-anchor residue come out with the same score. "
            "Identical scores in the table below are that blindness, not a "
            "finding.\n")
    if ctx["top_hits"]:
        parts.append(_table(
            ["products carrying this core", "chain", "core", "position",
             "foreign residues in core", "strongest molecule", "%rank (background)",
             "class"],
            ctx["top_hits"]))
        parts.append("")
        parts.append("%rank is computed against a per-molecule background distribution of "
                     "20 000 random peptides, because NetMHCIIpan emits no %Rank for "
                     "custom molecules. Lower is stronger; the MHC-II convention is "
                     "strong ≤ 2, weak ≤ 10.\n")

    # ---- 5. validation
    parts.append("## 5. Validation harness\n")
    parts.append(f"{LABELS['measured']} — {ctx['validation'].passed} of "
                 f"{len(ctx['validation'].checks)} checks passed.\n")
    parts.append(_table(
        ["id", "check", "status", "observed", "expected"],
        [[c.id, c.name, c.status, c.observed, c.expected or "—"]
         for c in ctx["validation"].checks]))
    parts.append("")
    parts.append("No check here tests *immunogenicity prediction accuracy* — there is no "
                 "dog or cat benchmark to test it against, which is the honest headline "
                 "of this whole exercise. What is tested is that the sequence layer "
                 "reproduces independently known facts, that the tolerance filter fires "
                 "exactly where it should, that positive and negative controls behave, "
                 "that the rank calibration is uniform, and that the out-of-domain "
                 "guard-rail fires instead of quietly passing an extrapolation through.\n")

    # ---- 6. what to do
    parts.append("## 6. Reading this for a real programme\n")
    parts.append(ctx["interpretation"])
    parts.append("")
    parts.append("### Figures\n")
    for caption, rel in ctx["figures"]:
        parts.append(f"![{caption}]({rel})\n\n*{caption}*\n")
    return "\n".join(parts)
