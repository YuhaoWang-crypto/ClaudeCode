"""Tables, figures and a labelled markdown report.

Every claim in the generated report carries one of three labels, following the
convention already used elsewhere in this repository:

* OK   - a direct measurement of the input structure
* HEUR - an uncalibrated heuristic; ranks candidates, predicts nothing
* TODO - requires a tool or experiment this package does not run
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .constants import AA20
from .liability import summarise_gates
from .offtarget import ProteomeHit
from .scan import WEIGHTS, Mutation
from .structure import Complex

LABELS = {
    "measured": "`[measured]`",
    "heuristic": "`[heuristic]`",
    "external": "`[needs external tool]`",
}


def write_tables(
    complex_obj: Complex,
    mutations: list[Mutation],
    hits: list[ProteomeHit],
    directory: Path,
) -> dict[str, Path]:
    directory.mkdir(parents=True, exist_ok=True)
    paths: dict[str, Path] = {}

    interface = pd.DataFrame(
        [
            {
                "residue": r.label,
                "chain": r.chain,
                "number": r.number,
                "aa": r.aa,
                "loop": r.loop,
                "contacts_peptide": r.contacts_peptide,
                "contacts_mhc": r.contacts_mhc,
                "mhc_contact_fraction": round(r.mhc_contact_fraction, 3),
                "buried_fraction": round(r.buried_fraction, 3),
                "delta_sasa_A2": round(r.sasa_free - r.sasa_complex, 1),
                "void_slots": r.void_slots,
                "growth_reaches_partner": r.growth_reaches_partner,
                "partner_charge": round(r.partner_charge, 3),
                "partner_polar_atoms": r.partner_polar_atoms,
            }
            for r in complex_obj.paratope()
        ]
    )
    paths["interface"] = directory / "interface_residues.csv"
    interface.to_csv(paths["interface"], index=False)

    scan = pd.DataFrame([m.as_row() for m in mutations])
    if not scan.empty:
        scan = scan.sort_values(
            ["verdict", "peptide_z", "z_score"], ascending=[True, False, False]
        )
    paths["scan"] = directory / "mutation_scan.csv"
    scan.to_csv(paths["scan"], index=False)

    offtarget = pd.DataFrame([h.as_row() for h in hits])
    paths["offtarget"] = directory / "offtarget_candidates.csv"
    offtarget.to_csv(paths["offtarget"], index=False)
    return paths


def plot_interface(complex_obj: Complex, contact_profile: list[int], path: Path) -> Path:
    paratope = complex_obj.paratope()
    fig, axes = plt.subplots(2, 1, figsize=(11, 8), height_ratios=[1, 1.6])

    peptide = complex_obj.peptide
    positions = [f"P{i + 1}\n{a}" for i, a in enumerate(peptide)]
    axes[0].bar(positions, contact_profile, color="#2b6cb0")
    axes[0].set_ylabel("TCR heavy-atom contacts")
    axes[0].set_title(
        f"{complex_obj.pdb_id}: which peptide positions the TCR reads "
        f"(peptide {peptide})"
    )
    axes[0].grid(axis="y", alpha=0.3)

    labels = [r.label for r in paratope]
    pep = np.array([r.contacts_peptide for r in paratope])
    mhc = np.array([r.contacts_mhc for r in paratope])
    x = np.arange(len(labels))
    axes[1].bar(x, pep, label="contacts to peptide", color="#2f855a")
    axes[1].bar(x, mhc, bottom=pep, label="contacts to MHC", color="#c05621")
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(labels, rotation=90, fontsize=7)
    axes[1].set_ylabel("heavy-atom contacts")
    axes[1].set_title(
        "TCR paratope: peptide contacts are the affinity opportunity, "
        "MHC contacts are the cross-reactivity risk"
    )
    axes[1].legend()
    axes[1].grid(axis="y", alpha=0.3)

    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def plot_mutation_heatmap(mutations: list[Mutation], path: Path) -> Path:
    if not mutations:
        return path
    positions: list[str] = []
    for mutation in mutations:
        key = f"{mutation.chain}:{mutation.wild_type}{mutation.number}"
        if key not in positions:
            positions.append(key)

    grid = np.full((len(positions), len(AA20)), np.nan)
    index = {key: i for i, key in enumerate(positions)}
    for mutation in mutations:
        key = f"{mutation.chain}:{mutation.wild_type}{mutation.number}"
        grid[index[key], AA20.index(mutation.mutant)] = mutation.triage_score

    limit = np.nanmax(np.abs(grid))
    fig, ax = plt.subplots(figsize=(10, 0.42 * len(positions) + 2.6))
    image = ax.imshow(grid, cmap="RdBu_r", vmin=-limit, vmax=limit, aspect="auto")
    ax.set_xticks(range(len(AA20)))
    ax.set_xticklabels(list(AA20))
    ax.set_yticks(range(len(positions)))
    ax.set_yticklabels(positions, fontsize=8)
    ax.set_xlabel("substituted to")
    ax.set_title(
        "Triage score per substitution (arbitrary units, uncalibrated)\n"
        "red = predicted to help, blue = predicted to hurt",
        fontsize=11,
    )
    fig.colorbar(image, ax=ax, label="triage score (a.u.)", shrink=0.7)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def plot_gate_funnel(mutations: list[Mutation], path: Path) -> Path:
    counts = summarise_gates(mutations)
    if not counts:
        return path
    order = [
        "advance",
        "watch-specificity",
        "below-threshold",
        "reject-specificity",
        "reject-liability",
    ]
    keys = [k for k in order if k in counts] + [k for k in counts if k not in order]
    values = [counts[k] for k in keys]
    colours = {
        "advance": "#2f855a",
        "watch-specificity": "#b7791f",
        "below-threshold": "#a0aec0",
        "reject-specificity": "#c05621",
        "reject-liability": "#9b2c2c",
    }
    fig, ax = plt.subplots(figsize=(8, 3.6))
    bars = ax.barh(keys, values, color=[colours.get(k, "#718096") for k in keys])
    ax.bar_label(bars, padding=3)
    ax.invert_yaxis()
    ax.set_xlabel("single mutants")
    ax.set_title(f"Gate outcome for {len(mutations)} scanned single mutants")
    ax.grid(axis="x", alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def plot_additivity(results: dict, path: Path) -> Path:
    """Sum-of-parts against measured, and how the error grows with stack size."""
    pooled = results["pooled"]["decompositions"]
    if not pooled:
        return path

    predicted = np.array([r["sum_of_parts"] for r in pooled], dtype=float)
    measured = np.array([r["measured"] for r in pooled], dtype=float)
    sizes = np.array([r["n_mutations"] for r in pooled], dtype=float)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.4))

    scatter = axes[0].scatter(
        predicted, measured, c=sizes, cmap="viridis", s=55, alpha=0.85,
        edgecolors="white", linewidths=0.5,
    )
    limits = [
        min(predicted.min(), measured.min()) - 0.5,
        max(predicted.max(), measured.max()) + 0.5,
    ]
    axes[0].plot(limits, limits, "--", color="#c05621", lw=1.2, label="perfect additivity")
    axes[0].set_xlim(limits)
    axes[0].set_ylim(limits)
    axes[0].set_xlabel("sum of the measured parts (kcal/mol)")
    axes[0].set_ylabel("measured whole (kcal/mol)")
    axes[0].set_title(
        "Does stacking mutations add up?\n"
        "points above the line underdelivered against the sum"
    )
    axes[0].legend(loc="upper left")
    axes[0].grid(alpha=0.3)
    fig.colorbar(scatter, ax=axes[0], label="mutations stacked")

    within = results["within_study"]["summary"].get("by_stack_size", {})
    pooled_sizes = results["pooled"]["summary"].get("by_stack_size", {})
    keys = sorted({int(k) for k in pooled_sizes} | {int(k) for k in within})
    x = np.arange(len(keys))
    width = 0.38
    pooled_values = [pooled_sizes.get(k, {}).get("mean_error", np.nan) for k in keys]
    within_values = [within.get(k, {}).get("mean_error", np.nan) for k in keys]

    axes[1].bar(x - width / 2, pooled_values, width, label="pooled across studies",
                color="#a0aec0")
    axes[1].bar(x + width / 2, within_values, width, label="within one study",
                color="#2b6cb0")
    axes[1].axhline(0, color="black", lw=0.9)
    axes[1].axhspan(-0.5, 0.5, color="#2f855a", alpha=0.12,
                    label="within 0.5 kcal/mol")
    axes[1].set_xticks(x)
    axes[1].set_xticklabels([str(k) for k in keys])
    axes[1].set_xlabel("mutations stacked")
    axes[1].set_ylabel("mean error, measured minus predicted (kcal/mol)")
    axes[1].set_title(
        "Pairs add up; deep stacks fall short\n"
        "positive means the stack underdelivered"
    )
    axes[1].legend(fontsize=8)
    axes[1].grid(axis="y", alpha=0.3)

    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def plot_calibration(entries, predictions, results: dict, path: Path) -> Path:
    """Three panels: the unfitted score, the held-out fit, and per-structure spread."""
    measured = np.array([e.ddg for e in entries], dtype=float)
    triage = np.array([e.triage_score for e in entries], dtype=float)
    predictions = np.asarray(predictions, dtype=float)
    usable = ~np.isnan(predictions)

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    unfitted = results["unfitted"]
    axes[0].scatter(triage, measured, s=18, alpha=0.6, color="#718096")
    axes[0].axhline(0, color="black", lw=0.8)
    axes[0].set_xlabel("triage score, as shipped (a.u.)")
    axes[0].set_ylabel("measured ddG (kcal/mol)")
    axes[0].set_title(
        f"Hand-weighted score\nSpearman {unfitted['spearman']:+.3f} "
        "(negative would be correct)"
    )
    axes[0].grid(alpha=0.3)

    held = results["fitted_extended"]["held_out"]
    axes[1].scatter(
        predictions[usable], measured[usable], s=18, alpha=0.6, color="#2b6cb0"
    )
    limits = [
        min(predictions[usable].min(), measured[usable].min()) - 0.5,
        max(predictions[usable].max(), measured[usable].max()) + 0.5,
    ]
    axes[1].plot(limits, limits, color="#c05621", lw=1, ls="--", label="y = x")
    axes[1].set_xlim(limits)
    axes[1].set_ylim(limits)
    axes[1].set_xlabel("predicted ddG, held out by structure (kcal/mol)")
    axes[1].set_ylabel("measured ddG (kcal/mol)")
    spearman_text = (
        "n/a" if held["spearman"] is None else f"{held['spearman']:+.3f}"
    )
    axes[1].set_title(
        f"Fitted on structural measurements\n"
        f"held-out Spearman {spearman_text}, RMSE {held['rmse']:.2f}"
    )
    axes[1].legend()
    axes[1].grid(alpha=0.3)

    per_structure = unfitted["per_structure"]
    names = list(per_structure)
    values = [
        per_structure[n]["spearman"] if per_structure[n]["spearman"] is not None else 0.0
        for n in names
    ]
    colours = ["#2f855a" if v < 0 else "#9b2c2c" for v in values]
    axes[2].barh(names, values, color=colours)
    axes[2].axvline(0, color="black", lw=0.8)
    axes[2].set_xlabel("Spearman of shipped score vs ddG")
    axes[2].set_title(
        "Per structure (green = correct sign)\nthe score does not transfer"
    )
    axes[2].invert_yaxis()
    axes[2].grid(axis="x", alpha=0.3)

    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def _mutation_table(rows: list[Mutation], limit: int = 12) -> str:
    if not rows:
        return "_none_\n"
    header = (
        "| mutant | loop | predicted ddG | score (a.u.) | z | pep contacts | "
        "MHC contacts | MHC frac | buried | why |\n"
        "|---|---|---|---|---|---|---|---|---|---|\n"
    )
    lines = []
    for mutation in rows[:limit]:
        z = mutation.peptide_z if mutation.peptide_z is not None else mutation.z_score
        why = mutation.notes[0] if mutation.notes else "-"
        predicted = (
            "-" if mutation.predicted_ddg is None else f"{mutation.predicted_ddg:+.2f}"
        )
        lines.append(
            f"| `{mutation.name}` | {mutation.loop or '-'} | {predicted} | "
            f"{mutation.triage_score:+.2f} | {z:+.2f} | "
            f"{mutation.contacts_peptide} | {mutation.contacts_mhc} | "
            f"{mutation.mhc_contact_fraction:.2f} | "
            f"{mutation.buried_fraction:.2f} | {why} |"
        )
    return header + "\n".join(lines) + "\n"


def write_report(
    complex_obj: Complex,
    mutations: list[Mutation],
    hits: list[ProteomeHit],
    motif: dict,
    contact_profile: list[int],
    figures: dict[str, Path],
    tables: dict[str, Path],
    handoff_dir: Path,
    path: Path,
) -> Path:
    """Write the labelled markdown report."""
    counts = summarise_gates(mutations)
    advance = sorted(
        [m for m in mutations if m.verdict == "advance"],
        key=lambda m: -(m.peptide_z if m.peptide_z is not None else m.z_score),
    )
    watch = sorted(
        [m for m in mutations if m.verdict == "watch-specificity"],
        key=lambda m: -(m.peptide_z if m.peptide_z is not None else m.z_score),
    )
    rejected = sorted(
        [m for m in mutations if m.verdict == "reject-specificity"],
        key=lambda m: -m.z_score,
    )
    paratope = complex_obj.paratope()
    buried = sum(r.sasa_free - r.sasa_complex for r in paratope)

    peptide = complex_obj.peptide
    read_positions = [
        f"P{i + 1}{peptide[i]} ({contact_profile[i]})"
        for i in np.argsort(contact_profile)[::-1][:4]
        if contact_profile[i] > 0
    ]

    lines = [
        f"# TCR interface triage: {complex_obj.pdb_id}",
        "",
        "Generated by `tcr_design`. Every claim below is labelled:",
        "",
        f"* {LABELS['measured']} direct measurement of the deposited structure",
        f"* {LABELS['heuristic']} uncalibrated heuristic: it ranks, it does not predict",
        f"* {LABELS['external']} needs a tool or experiment this package does not run",
        "",
        "## 1. System",
        "",
        f"| role | chain | length |",
        "|---|---|---|",
    ]
    for role in ("mhc", "b2m", "peptide", "tcr_alpha", "tcr_beta"):
        chain_id = complex_obj.chains.get(role)
        if chain_id:
            lines.append(
                f"| {role} | {chain_id} | {len(complex_obj.sequences[chain_id])} |"
            )
    lines += [
        "",
        f"Peptide: `{peptide}`",
        "",
        "Chain typing decisions:",
        "",
    ]
    lines += [f"* {note}" for note in complex_obj.typing_notes]
    lines += [
        "",
        "CDR loops located from the conserved V-domain anchors "
        "(IMGT junction definition for CDR3):",
        "",
        "| loop | chain | sequence |",
        "|---|---|---|",
    ]
    for name, (chain_id, start, end) in sorted(complex_obj.loops.items()):
        lines.append(f"| {name} | {chain_id} | `{complex_obj.sequences[chain_id][start:end]}` |")

    lines += [
        "",
        "## 2. Interface geometry " + LABELS["measured"],
        "",
        f"* Paratope: **{len(paratope)} TCR residues** contact the pMHC "
        f"(heavy-atom cutoff 5 A).",
        f"* Buried TCR surface on binding: **{buried:.0f} A^2**. "
        "Typical TCR-pMHC interfaces bury roughly 600-1000 A^2 per side, so "
        "this is a normal interface.",
        f"* Peptide positions the receptor reads hardest: "
        f"{', '.join(read_positions)}.",
        "",
        "The peptide/MHC split of each paratope residue's contacts is the key "
        "number for everything that follows. Affinity added at peptide-contacting "
        "positions sharpens recognition; affinity added at MHC-contacting positions "
        "broadens it.",
        "",
        f"![interface]({figures['interface'].name})",
        "",
        f"Full per-residue measurements: `{tables['interface'].name}`",
        "",
        "## 3. Mutation scan " + LABELS["heuristic"],
        "",
        f"Scanned **{len(mutations)} single mutants** across the contacting CDR "
        "positions.",
        "",
        "The score is a weighted sum of six terms. The weights were set by hand "
        "from interface-energetics principles and **were not fitted to any "
        "affinity data**:",
        "",
        "| term | weight | what it measures |",
        "|---|---|---|",
        f"| packing | {WEIGHTS['packing']} | side-chain volume change against measured local void |",
        f"| apolar_burial | {WEIGHTS['apolar_burial']} | hydropathy change scaled by measured burial |",
        f"| electrostatic | {WEIGHTS['electrostatic']} | charge complementarity to nearby partner side chains |",
        f"| hbond | {WEIGHTS['hbond']} | donor/acceptor change against partner polar atoms |",
        f"| conformational | {WEIGHTS['conformational']} | proline/glycine backbone penalties |",
        f"| substitution_prior | {WEIGHTS['substitution_prior']} | BLOSUM62 plausibility, as a fold/expression proxy |",
        "",
        "Scores are reported in arbitrary units plus a z-score. Positions that "
        "touch the peptide are z-scored against each other, because across a whole "
        "TCR-pMHC interface the MHC-contacting positions would otherwise dominate "
        "the ranking.",
        "",
        f"![heatmap]({figures['heatmap'].name})",
        "",
        "### Gate outcomes",
        "",
        f"![funnel]({figures['funnel'].name})",
        "",
        "| verdict | count |",
        "|---|---|",
    ]
    for key, value in counts.items():
        lines.append(f"| {key} | {value} |")

    lines += [
        "",
        "### Advance",
        "",
        "Clean, peptide-focused, and predicted to tighten binding. These are the "
        "candidates worth spending real compute on. When this list is empty, the "
        "best surviving candidates are still handed to tier 2: the calibrated "
        "model is known to under-call real improvements, so an empty list is a "
        "reason to spend the compute rather than to stop.",
        "",
        _mutation_table(advance),
        "### Watch for specificity",
        "",
        "Predicted gain, but the position leans on MHC contacts.",
        "",
        _mutation_table(watch),
        "### Rejected on specificity",
        "",
        "These score well and are still rejected: their contacts are almost "
        "entirely to MHC rather than the peptide, so any affinity they buy is "
        "peptide-independent. This is the dominant outcome of the scan and it is "
        "the main result worth taking away.",
        "",
        _mutation_table(rejected, limit=8),
        "## 4. Off-target screen " + LABELS["heuristic"],
        "",
        "Peptide positions are weighted by measured TCR contact, filtered to those "
        "the HLA allele could present, and matched against the reviewed human "
        "proteome. Identity is required at the most heavily contacted position(s).",
        "",
        "| parameter | value |",
        "|---|---|",
        f"| target peptide | `{motif['target_peptide']}` |",
        f"| weighted positions | {', '.join(motif['weighted_positions'])} |",
        f"| identity required at | {', '.join(motif['identity_required_at'])} |",
        f"| HLA anchor filter | {motif['anchor_filter']} |",
        f"| proteome windows scanned | {motif['windows_scanned']:,} |",
        f"| passing all filters | {motif['windows_passing_all_filters']:,} |",
        f"| above similarity {motif['min_similarity']} | {motif['hits_above_threshold']:,} |",
        "",
        "Top candidate human mimics:",
        "",
        "| peptide | similarity | identical at contacts | pattern | gene | protein |",
        "|---|---|---|---|---|---|",
    ]
    for hit in hits[:15]:
        lines.append(
            f"| `{hit.peptide}` | {hit.similarity:.2f} | {hit.identical_contacts} | "
            f"`{hit.contact_positions_matched}` | {hit.gene} | {hit.protein[:48]} |"
        )

    lines += [
        "",
        f"Full list: `{tables['offtarget'].name}`",
        "",
        "**What this list is.** A prioritised shortlist of human peptides that "
        "resemble the target at the positions this receptor actually reads. The "
        "pattern column shows per-position agreement: uppercase identical, "
        "lowercase conservative, `-` non-conservative, `.` unweighted.",
        "",
        "**What this list is not.** A cross-reactivity prediction. The motif comes "
        "from one crystal structure rather than from an experimental scan of the "
        "receptor; presentation is filtered by a coarse anchor motif rather than "
        "scored with NetMHCpan; and the screen covers only the canonical reviewed "
        "proteome, so isoforms, modified peptides and non-canonical products are "
        "invisible to it. " + LABELS["external"],
        "",
        "This step is included because it is the one that cannot be skipped. An "
        "affinity-enhanced MAGE-A3 TCR cross-reacted with a titin peptide and "
        "killed two patients; the peptide was found only on retrospective review, "
        "and it was invisible both in cardiac cell lines and in mouse models. "
        "A published tool, CrossDome, implements the same contact-weighting idea "
        "more thoroughly and moved that titin peptide from rank 27 to rank 6 of "
        "36,000 candidates.",
        "",
        "## 5. What to run next " + LABELS["external"],
        "",
        f"Tier-2 inputs for every advancing candidate are in `{handoff_dir.name}/`, "
        f"with exact commands in `{handoff_dir.name}/RUN_TIER2.md`:",
        "",
        "* AlphaFold3 and Boltz co-folding inputs for each mutant complex",
        "* Rosetta resfiles for flex ddG",
        "* a FoldX `individual_list.txt`",
        "* ProteinMPNN chain and fixed-position files for CDR loop redesign",
        "",
        "Order of trust for this problem, based on published TCR-specific "
        "validation: MM/PBSA against ATLAS mutants (r = 0.74-0.9, small n) > "
        "Rosetta flex ddG (well validated on general interfaces, no TCR benchmark) "
        "> FoldX > co-folding confidence scores (which are not affinity). "
        "None of these replace a binding measurement.",
        "",
        "## 6. Honest limits of this pipeline",
        "",
        "**Both the heuristic and the calibrated model fail the retrospective "
        "benchmark.** The one test buildable from public structures is the DMF5 "
        "TCR: `3QDG` is wild type and `4L3E` is a designed variant reported at "
        "roughly 400-fold improved affinity with retained specificity. The two "
        "substitutions that achieved it (alpha D26Y in CDR1, beta L98W in CDR3) "
        "rank in the 47th and 44th percentile of a blind scan with the raw "
        "heuristic, and the 25th and 26th with the calibrated model refitted "
        "without those structures. Calibration roughly halves the percentile, and "
        "L98W becomes the top-ranked substitution of 19 at its own position, but "
        "both are still predicted to *weaken* binding when D26Y in fact tightens "
        "it by about 1.8 kcal/mol. Neither would have been shortlisted. Run "
        "`python -m tcr_design validate 3QDG 4L3E --calibrated` to reproduce it.",
        "",
        "* Against pooled ATLAS and SKEMPI measurements the calibrated model "
        "reaches a held-out Spearman of about +0.43 and an RMSE near 0.88 "
        "kcal/mol, versus 0.98 for predicting the mean. That is real but weak, and "
        "its predictions are compressed towards the middle of the range, so it "
        "rarely calls anything strongly stabilising.",
        "* The benchmark behind it is dominated by alanine scanning, so the model "
        "mostly learned that deleting a buried side chain weakens binding. Adding "
        "volume to gain affinity is the rarer direction and the one design needs.",
        "* The raw triage score is unfitted. Against the same measurements its "
        "correlation with ddG is about -0.08, which is no signal at all; refitting "
        "its six weights does not beat predicting the mean. Treat it as a way to "
        "allocate compute, never as a result.",
        "* The specificity gate uses contact partitioning in a single wild-type "
        "structure as a proxy for cross-reactivity risk. The DMF5 case above shows "
        "that proxy is crude: an MHC-contacting position can be enhanced without "
        "losing specificity. The gate raises a flag for triage, it does not "
        "establish risk.",
        "* Single point mutations only. Real affinity maturation explores "
        "combinations, and the best published TCR design result came from "
        "structure-based physics plus wet-lab iteration, not from a one-shot scan.",
        "* One static structure. Side-chain and loop rearrangement on binding is "
        "not modelled, and TCR CDR3 loops are flexible.",
        "* No expression, pairing, stability or immunogenicity assessment.",
        "* Higher affinity is not automatically better: very high affinity can "
        "reduce T cell activity, and affinity enhancement is what produced the "
        "best-known off-target fatality in this field.",
        "",
    ]
    path.write_text("\n".join(lines) + "\n")
    return path


def write_summary_json(
    complex_obj: Complex,
    mutations: list[Mutation],
    hits: list[ProteomeHit],
    motif: dict,
    path: Path,
) -> Path:
    payload = {
        "pdb_id": complex_obj.pdb_id,
        "chains": complex_obj.chains,
        "peptide": complex_obj.peptide,
        "loops": {k: complex_obj.sequences[v[0]][v[1] : v[2]] for k, v in complex_obj.loops.items()},
        "paratope_size": len(complex_obj.paratope()),
        "buried_surface_A2": round(
            sum(r.sasa_free - r.sasa_complex for r in complex_obj.paratope()), 1
        ),
        "mutations_scanned": len(mutations),
        "gate_counts": summarise_gates(mutations),
        "advance": [m.name for m in mutations if m.verdict == "advance"],
        "offtarget_motif": motif,
        "offtarget_top": [h.as_row() for h in hits[:10]],
        "disclaimer": (
            "Triage scores are uncalibrated arbitrary units, not affinity or ddG "
            "predictions. The off-target list is a similarity shortlist, not a "
            "cross-reactivity prediction."
        ),
    }
    path.write_text(json.dumps(payload, indent=2) + "\n")
    return path
