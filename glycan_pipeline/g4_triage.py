"""
Module G4 - Design triage: sequons, confidence gates, and what the filter
cascade actually buys you.

Three things a glycoprotein design campaign has to get right before any
protein is expressed, all reproduced here from the paper's own numbers.

1. SEQUON PLACEMENT. An N-glycan needs N-X-S/T with X != Pro. The paper
   biases +2 positions to threonine because NXT sequons are occupied more
   often than NXS. Occupancy is not a detail: a beautifully shielded site
   that is only 5% occupied tells you nothing, which is why the paper gates
   on occupancy > 5% before analysing composition at all, and > 50% before
   selecting a design.

2. THE CONFIDENCE GATE DEPENDS ON THE PREDICTOR. The paper's first-generation
   filter used AF2 without a glycan; the second used glycan-aware RF3/AF3.
   Extended Data Table 2 gives both for O1-O4, and it is the single most
   useful table in the paper for anyone adapting the method, because it shows
   a structure-prediction filter rejecting designs that were in fact correct.
   AF2, blind to the glycan, reads a deep glycan-housing pocket as a large
   empty cavity and reports a broken model.

3. THE DELTA-RASA FILTER IS WHAT MOVED THE NUMBERS. The paper reports 20/43
   designs over 25% oligomannose in the DeltaRASA-filtered set against 1/97
   in the unfiltered first generation. That contingency table is reproduced
   here with an exact Fisher test computed from scratch.

Rigor
-----
(rigorous) ED Table 2 values are transcribed from the paper. The gate
evaluations, the Fisher test and the sequon scan are deterministic.
(caution) The AF2 gate thresholds are taken from the paper's description of
its first-generation filter (pLDDT > 80, pAE < 10, Ca-RMSD < 3 A). The paper
states the conclusion - only O2 would have passed - but does not print the
gate alongside the table, so the mapping of threshold to column is inferred.
The reproduction agreeing with the paper's stated conclusion is the check.
"""

from __future__ import annotations

import math
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .g1_delta_rasa import FIG_DIR

# ---------------------------------------------------------------------------
# Extended Data Table 2, transcribed. "Model" rows are model indices
# (AF2 indexed 1-5, AF3 from 0), not scores. RF3 reports no pTM.
# ---------------------------------------------------------------------------
ED_TABLE_2: dict[str, dict] = {
    "O1": {"af2_plddt": 42.7, "af2_ptm": 0.3, "af2_pae": 21.5, "af2_rmsd": 21.8,
           "rf3_plddt": 87.2, "rf3_pae": 6.5, "rf3_rmsd": 1.50,
           "af3_plddt": 84.0, "af3_pae": 7.8, "af3_ptm": 0.85, "af3_rmsd": 1.2},
    "O2": {"af2_plddt": 85.6, "af2_ptm": 0.8, "af2_pae": 7.5, "af2_rmsd": 0.7,
           "rf3_plddt": 82.2, "rf3_pae": 8.8, "rf3_rmsd": 2.82,
           "af3_plddt": 81.9, "af3_pae": 9.3, "af3_ptm": 0.83, "af3_rmsd": 3.1},
    "O3": {"af2_plddt": 67.7, "af2_ptm": 0.7, "af2_pae": 12.6, "af2_rmsd": 16.8,
           "rf3_plddt": 81.2, "rf3_pae": 8.8, "rf3_rmsd": 2.26,
           "af3_plddt": 82.8, "af3_pae": 8.6, "af3_ptm": 0.83, "af3_rmsd": 2.1},
    "O4": {"af2_plddt": 77.5, "af2_ptm": 0.8, "af2_pae": 9.2, "af2_rmsd": 6.1,
           "rf3_plddt": 75.4, "rf3_pae": 11.4, "rf3_rmsd": 7.60,
           "af3_plddt": 75.4, "af3_pae": 10.3, "af3_ptm": 0.78, "af3_rmsd": 5.7},
}

#: Confidence gates as described in the paper.
GATES = {
    # First generation: AF2 in single-sequence mode, no glycan modelled.
    "af2_gen1": {"plddt_min": 80.0, "pae_max": 10.0, "rmsd_max": 3.0},
    # Second generation: glycan-aware RF3 or AF3 server.
    "rf3_gen2": {"plddt_min": 75.0, "pae_max": 15.0, "rmsd_max": 8.0},
    "af3_gen2": {"plddt_min": 75.0, "pae_max": 15.0, "rmsd_max": 8.0},
    # THRgly redesign campaign used a tighter AF2 gate.
    "af2_thr": {"plddt_min": 90.0, "pae_max": 5.0, "rmsd_max": 2.0},
}

#: Outcome of the DeltaRASA filter, as reported in the paper's Results.
DELTA_RASA_CONTINGENCY = {
    "filtered_oligomannose": 20, "filtered_total": 43,
    "unfiltered_oligomannose": 1, "unfiltered_total": 97,
    "criterion": ">25% oligomannose among designs with >5% glycan occupancy",
}


def evaluate_gate(metrics: dict, gate: str, prefix: str) -> dict:
    """Apply one confidence gate to one design, returning per-term verdicts."""
    g = GATES[gate]
    plddt = metrics[f"{prefix}_plddt"]
    pae = metrics[f"{prefix}_pae"]
    rmsd = metrics[f"{prefix}_rmsd"]
    terms = {
        "pLDDT": (plddt, f">{g['plddt_min']:g}", plddt > g["plddt_min"]),
        "pAE": (pae, f"<{g['pae_max']:g}", pae < g["pae_max"]),
        "RMSD": (rmsd, f"<{g['rmsd_max']:g}", rmsd < g["rmsd_max"]),
    }
    return {"passes": all(t[2] for t in terms.values()), "terms": terms}


SEQUON_RE = re.compile(r"(?=(N[^P][ST]))")


def scan_sequons(sequence: str) -> list[dict]:
    """Find every N-X-S/T sequon (X != Pro) in a protein sequence.

    Returns 1-based Asn positions, the motif, and whether it is the
    higher-occupancy NXT form.
    """
    seq = sequence.strip().upper().replace(" ", "")
    out = []
    for m in SEQUON_RE.finditer(seq):
        motif = m.group(1)
        out.append({
            "asn_position": m.start() + 1,
            "motif": motif,
            "x_residue": motif[1],
            "plus2": motif[2],
            "is_nxt": motif[2] == "T",
        })
    return out


def _hypergeom_pmf(k, K, n, N):
    """P(X = k) for drawing n from N with K successes."""
    return (math.comb(K, k) * math.comb(N - K, n - k)) / math.comb(N, n)


def fisher_exact_two_sided(a, b, c, d) -> dict:
    """Two-sided Fisher exact test on [[a, b], [c, d]] (no scipy).

    Sums the probability of every table at least as extreme as observed,
    the standard "sum of small p-values" convention.
    """
    row1, row2 = a + b, c + d
    col1 = a + c
    total = row1 + row2
    p_obs = _hypergeom_pmf(a, row1, col1, total)
    lo = max(0, col1 - row2)
    hi = min(row1, col1)
    p = 0.0
    for k in range(lo, hi + 1):
        pk = _hypergeom_pmf(k, row1, col1, total)
        if pk <= p_obs * (1 + 1e-9):
            p += pk
    odds = float("inf") if (b == 0 or c == 0) else (a * d) / (b * c)
    return {"p_value": p, "odds_ratio": odds,
            "table": [[a, b], [c, d]], "p_observed_table": p_obs}


def analyse():
    designs = {}
    for name, m in ED_TABLE_2.items():
        designs[name] = {
            "metrics": m,
            "af2_gen1": evaluate_gate(m, "af2_gen1", "af2"),
            "rf3_gen2": evaluate_gate(m, "rf3_gen2", "rf3"),
            "af3_gen2": evaluate_gate(m, "af3_gen2", "af3"),
        }

    af2_pass = [n for n, d in designs.items() if d["af2_gen1"]["passes"]]
    af3_pass = [n for n, d in designs.items() if d["af3_gen2"]["passes"]]
    rf3_pass = [n for n, d in designs.items() if d["rf3_gen2"]["passes"]]

    c = DELTA_RASA_CONTINGENCY
    fisher = fisher_exact_two_sided(
        c["filtered_oligomannose"],
        c["filtered_total"] - c["filtered_oligomannose"],
        c["unfiltered_oligomannose"],
        c["unfiltered_total"] - c["unfiltered_oligomannose"],
    )

    return {
        "designs": designs,
        "af2_gen1_pass": af2_pass,
        "rf3_gen2_pass": rf3_pass,
        "af3_gen2_pass": af3_pass,
        "paper_claim_only_O2_passes_af2": af2_pass == ["O2"],
        "fisher": fisher,
        "contingency": c,
        "sequon_audit": sequon_audit(),
    }


#: Structures used for the sequon audit. 5FYL needs its biological assembly
#: so all three protomers' sequons are present.
AUDIT_STRUCTURES = [("5FYL", 1), ("7PFP", None), ("6VXX", None)]


def sequon_audit(structures=None) -> dict:
    """Cross-check sequence-level sequon rules against observed glycans.

    Two independent code paths meet here: ``sequons_in_structure`` reads
    N-X-S/T off the amino-acid sequence, and ``find_glycan_trees`` finds
    sugars covalently bonded to an Asn by geometry alone. Every glycan ought
    to land on a sequon. Anything else means a bug in one of them, or a
    genuinely non-canonical site.
    """
    from .glycans import find_glycan_trees, load_model, sequons_in_structure
    from .pdbcache import fetch

    structures = structures or AUDIT_STRUCTURES
    per_structure, occupied, empty = [], [], []
    for pdb_id, assembly in structures:
        model = load_model(fetch(pdb_id, assembly=assembly))
        sq = sequons_in_structure(model)
        glycosylated = {(t.asn_chain, t.asn_resseq)
                        for t in find_glycan_trees(model)
                        if t.asn_resseq is not None}
        occ = [v for k, v in sq.items() if k in glycosylated]
        emp = [v for k, v in sq.items() if k not in glycosylated]
        orphans = sorted(glycosylated - set(sq))
        occupied += occ
        empty += emp
        per_structure.append({
            "pdb_id": pdb_id, "assembly": assembly,
            "n_sequons": len(sq),
            "n_with_glycan": len(occ),
            "n_without_glycan": len(emp),
            "n_glycans_outside_a_sequon": len(orphans),
            "orphans": orphans[:5],
        })

    def n_nxt(seqs):
        return sum(1 for s in seqs if s.endswith("T"))

    a, b = n_nxt(occupied), len(occupied) - n_nxt(occupied)
    c, d = n_nxt(empty), len(empty) - n_nxt(empty)
    test = fisher_exact_two_sided(a, b, c, d) if (c + d) else None

    return {
        "per_structure": per_structure,
        "total_glycans_outside_a_sequon":
            sum(s["n_glycans_outside_a_sequon"] for s in per_structure),
        "nxt_occupied": a, "nxs_occupied": b,
        "nxt_empty": c, "nxs_empty": d,
        "nxt_test": test,
        "powered": bool((c + d) >= 30),
    }


def make_figure(r, path=None):
    path = path or os.path.join(FIG_DIR, "g4_triage.png")
    os.makedirs(os.path.dirname(path), exist_ok=True)

    names = list(r["designs"])
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))

    # pLDDT / RMSD by predictor
    ax = axes[0]
    x = np.arange(len(names))
    for i, (pref, col, lab) in enumerate([
            ("af2", "#C44E52", "AF2 (no glycan)"),
            ("rf3", "#4C72B0", "RF3 (glycan-aware)"),
            ("af3", "#55A868", "AF3 (glycan-aware)")]):
        vals = [r["designs"][n]["metrics"][f"{pref}_plddt"] for n in names]
        ax.bar(x + (i - 1) * 0.27, vals, width=0.27, color=col, label=lab)
    ax.axhline(80, ls="--", c="gray", lw=1)
    ax.text(len(names) - 0.4, 81, "AF2 gate", fontsize=7, color="gray", ha="right")
    ax.set_xticks(x)
    ax.set_xticklabels(names)
    ax.set_ylabel("pLDDT")
    ax.set_title("Glycan-blind AF2 under-scores\nglycan-housing pockets", fontsize=10)
    ax.legend(fontsize=7, frameon=False)
    ax.grid(axis="y", alpha=0.3)

    ax = axes[1]
    for i, (pref, col, lab) in enumerate([
            ("af2", "#C44E52", "AF2"),
            ("rf3", "#4C72B0", "RF3"),
            ("af3", "#55A868", "AF3")]):
        vals = [r["designs"][n]["metrics"][f"{pref}_rmsd"] for n in names]
        ax.bar(x + (i - 1) * 0.27, vals, width=0.27, color=col, label=lab)
    ax.axhline(3.0, ls="--", c="gray", lw=1)
    ax.set_xticks(x)
    ax.set_xticklabels(names)
    ax.set_ylabel(r"C$\alpha$-RMSD to design model ($\AA$)")
    ax.set_title("AF2 calls O1/O3 broken;\nglycan-aware models do not", fontsize=10)
    ax.legend(fontsize=7, frameon=False)
    ax.grid(axis="y", alpha=0.3)

    ax = axes[2]
    c = r["contingency"]
    fracs = [c["unfiltered_oligomannose"] / c["unfiltered_total"],
             c["filtered_oligomannose"] / c["filtered_total"]]
    bars = ax.bar(["first generation\n(no $\\Delta$RASA)",
                   "$\\Delta$RASA-filtered"], fracs,
                  color=["#C44E52", "#55A868"])
    for bar, num, den in zip(bars,
                             [c["unfiltered_oligomannose"], c["filtered_oligomannose"]],
                             [c["unfiltered_total"], c["filtered_total"]]):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.01,
                f"{num}/{den}", ha="center", fontsize=9)
    ax.set_ylabel("fraction >25% oligomannose")
    ax.set_title(r"$\Delta$RASA filter enrichment"
                 f"\nFisher exact p = {r['fisher']['p_value']:.2e}", fontsize=10)
    ax.grid(axis="y", alpha=0.3)

    fig.suptitle("Design triage - sequons, confidence gates, and filter payoff",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def report():
    r = analyse()
    print("=" * 68)
    print("MODULE G4 - Design triage (Extended Data Table 2 + Results counts)")
    print("=" * 68)

    print("confidence gates applied to O1-O4:")
    print(f"    {'design':>7} {'AF2 gen-1':>11} {'RF3 gen-2':>11} {'AF3 gen-2':>11}")
    for n, d in r["designs"].items():
        print(f"    {n:>7} "
              f"{'PASS' if d['af2_gen1']['passes'] else 'FAIL':>11} "
              f"{'PASS' if d['rf3_gen2']['passes'] else 'FAIL':>11} "
              f"{'PASS' if d['af3_gen2']['passes'] else 'FAIL':>11}")

    print(f"\n    AF2 gen-1 survivors : {r['af2_gen1_pass']}")
    print(f"    RF3 gen-2 survivors : {r['rf3_gen2_pass']}")
    print(f"    AF3 gen-2 survivors : {r['af3_gen2_pass']}")
    print(f"    paper states only O2 would have passed AF2 -> reproduced: "
          f"{r['paper_claim_only_O2_passes_af2']}")

    print("\n    why each AF2 rejection happened:")
    for n, d in r["designs"].items():
        if d["af2_gen1"]["passes"]:
            continue
        bad = [f"{k}={v[0]:g} (need {v[1]})"
               for k, v in d["af2_gen1"]["terms"].items() if not v[2]]
        print(f"        {n}: " + "; ".join(bad))

    c, f = r["contingency"], r["fisher"]
    print(f"\nDeltaRASA filter payoff ({c['criterion']}):")
    print(f"    first generation : {c['unfiltered_oligomannose']}/"
          f"{c['unfiltered_total']} "
          f"({c['unfiltered_oligomannose'] / c['unfiltered_total']:.1%})")
    print(f"    DeltaRASA-filtered: {c['filtered_oligomannose']}/"
          f"{c['filtered_total']} "
          f"({c['filtered_oligomannose'] / c['filtered_total']:.1%})")
    print(f"    Fisher exact (two-sided) p = {f['p_value']:.3e}, "
          f"odds ratio = {f['odds_ratio']:.1f}")
    print(f"    paper reports p < 1e-10 -> reproduced: "
          f"{f['p_value'] < 1e-10}")

    a = r["sequon_audit"]
    print("\nsequon audit (sequence rule vs observed glycans):")
    print(f"    {'structure':>16} {'sequons':>8} {'glycosylated':>13} "
          f"{'empty':>6} {'glycan off-sequon':>18}")
    for s in a["per_structure"]:
        tag = f"{s['pdb_id']}" + (f" asm{s['assembly']}" if s["assembly"] else "")
        print(f"    {tag:>16} {s['n_sequons']:>8} {s['n_with_glycan']:>13} "
              f"{s['n_without_glycan']:>6} {s['n_glycans_outside_a_sequon']:>18}")
    print(f"    -> {a['total_glycans_outside_a_sequon']} glycans sat outside a "
          f"canonical sequon, across "
          f"{sum(s['n_with_glycan'] for s in a['per_structure'])} observed glycans.")

    print(f"\n    NXT vs NXS   glycosylated {a['nxt_occupied']}/{a['nxs_occupied']}"
          f"   empty {a['nxt_empty']}/{a['nxs_empty']}")
    if not a["powered"]:
        print("    (caution) too few unoccupied sequons to test the NXT occupancy")
        print("              preference here"
              + (f" (p = {a['nxt_test']['p_value']:.2f})" if a["nxt_test"] else "")
              + "; and an unmodelled glycan is")
        print("              not an unoccupied one. That preference has to come")
        print("              from LC-MS/MS occupancy, which is what the paper used.")

    p = make_figure(r)
    print(f"\nfigure written to: {p}")
    return r


if __name__ == "__main__":
    report()
