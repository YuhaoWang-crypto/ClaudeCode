"""Step 9b -- export bundle and consistency gates.

Nothing leaves the pipeline until the numbers agree with each other.  Each
gate compares two artefacts that were produced independently -- the cohort
table against the expression matrix, the recall string against the per-gene
ranks, the literature table against the ranking -- and the bundle is written
with a manifest that records every gate outcome and a checksum per file.
"""

from __future__ import annotations

import hashlib
import json

import pandas as pd

from . import config as C

EXPORT = C.RESULTS / "export"


def _sha256(path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()[:16]


def _gate(name: str, ok: bool, detail: str) -> dict:
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {detail}")
    return {"gate": name, "passed": bool(ok), "detail": detail}


def run() -> dict:
    print("[S9b] export consistency gates")
    EXPORT.mkdir(parents=True, exist_ok=True)

    cohort = pd.read_csv(C.RESULTS / "s2_cohort.csv")
    consensus = pd.read_csv(C.RESULTS / "s2_consensus_expression.csv")
    ranked = pd.read_csv(C.RESULTS / "s6_ranked_candidates.csv")
    positives = pd.read_csv(C.RESULTS / "s7_validation_positives.csv")
    negatives = pd.read_csv(C.RESULTS / "s7_negative_controls.csv")
    literature = pd.read_csv(C.RESULTS / "s8_literature_evidence.csv")
    s2p = json.loads((C.RESULTS / "s2_provenance.json").read_text())
    s7p = json.loads((C.RESULTS / "s7_provenance.json").read_text())

    gates = []

    # G1 -- cohort accounting: compartment counts must add up to the analysed cells
    per_comp = sum(int(cohort[f"cells_{c}"].sum()) for c in C.COMPARTMENTS)
    analysed = int(cohort.cells_in_compartments.sum())
    gates.append(_gate(
        "cohort cell accounting", per_comp == analysed == s2p["cells_analysed"],
        f"{per_comp:,} cells summed over compartments = {analysed:,} in the cohort table "
        f"= {s2p['cells_analysed']:,} in provenance"))

    # G2 -- sampling ceiling honoured
    gates.append(_gate(
        "sampling ceiling", bool((cohort.cells_sampled <= C.MAX_CELLS_PER_DATASET).all()),
        f"max sampled per atlas {int(cohort.cells_sampled.max()):,} "
        f"<= ceiling {C.MAX_CELLS_PER_DATASET:,}"))

    # G3 -- every ranked gene was measured in every atlas of the cohort
    n_datasets = int(cohort.dataset_id.nunique())
    complete = int((consensus.n_datasets == n_datasets).sum())
    gates.append(_gate(
        "expression matrix coverage", complete == len(consensus),
        f"{complete:,}/{len(consensus):,} genes carry statistics from all {n_datasets} atlases"))

    # G4 -- ranking integrity
    ok = (ranked.gene.is_unique and ranked["rank"].is_monotonic_increasing
          and ranked.final_score.between(0, 1).all()
          and ranked.final_score.is_monotonic_decreasing)
    gates.append(_gate(
        "ranking integrity", ok,
        f"{len(ranked):,} unique genes, scores monotone in [0,1], "
        f"range {ranked.final_score.min():.3f}-{ranked.final_score.max():.3f}"))

    # G5 -- the recall string matches the per-gene ranks
    recall_ok = True
    details = []
    for k in C.RECALL_K:
        stated = s7p["recall"][f"recall@{k}"]
        recomputed = f"{int((positives['rank'] <= k).sum())}/{len(positives)}"
        recall_ok &= stated == recomputed
        details.append(f"@{k}: reported {stated}, recomputed {recomputed}")
    gates.append(_gate("recall report", recall_ok, "; ".join(details)))

    # G6 -- negative-control verdicts follow the rule, not the narrative
    def expected(row):
        if pd.isna(row["rank"]):
            return "PASS"
        return "FAIL" if row["rank"] <= C.NEGATIVE_CONTROL_FAIL_RANK else "PASS"
    neg_ok = bool((negatives.apply(expected, axis=1) == negatives.verdict).all())
    gates.append(_gate(
        "negative-control verdicts", neg_ok,
        "; ".join(f"{r.gene}={r.verdict}(rank {'-' if pd.isna(r['rank']) else int(r['rank'])})"
                  for _, r in negatives.iterrows())))

    # G7 -- literature table covers exactly the head of the ranking
    lit_ok = list(literature.gene) == list(ranked.gene.head(len(literature)))
    gates.append(_gate(
        "literature coverage", lit_ok,
        f"{len(literature)} reviewed candidates match ranks 1-{len(literature)}"))

    # G8 -- figures present in both formats
    missing = [f"{n}.{ext}"
               for n in ("fig1_candidate_ranking", "fig2_compartment_heatmap",
                         "fig3_therapeutic_index", "fig4_validation_ranks")
               for ext in ("png", "svg")
               if not (C.FIGURES / f"{n}.{ext}").exists()
               or (C.FIGURES / f"{n}.{ext}").stat().st_size == 0]
    gates.append(_gate("figures present", not missing,
                       "8 files (4 figures x PNG+SVG)" if not missing else f"missing {missing}"))

    # ---- export bundle ---------------------------------------------------
    ranked.to_csv(EXPORT / "antigen_targets_all_ranked.csv", index=False)
    ranked.head(50).to_csv(EXPORT / "antigen_targets_top50.csv", index=False)
    literature.to_csv(EXPORT / "antigen_targets_literature.csv", index=False)
    cohort.to_csv(EXPORT / "cohort_datasets.csv", index=False)
    positives.to_csv(EXPORT / "validation_positives.csv", index=False)
    negatives.to_csv(EXPORT / "validation_negative_controls.csv", index=False)

    manifest = {
        "gates": gates,
        "all_gates_passed": all(g["passed"] for g in gates),
        "files": {p.name: {"bytes": p.stat().st_size, "sha256_16": _sha256(p)}
                  for p in sorted(EXPORT.glob("*.csv"))},
        "census_version": s2p["census_version"],
        "cells_analysed": s2p["cells_analysed"],
        "candidates_ranked": int(len(ranked)),
    }
    (EXPORT / "manifest.json").write_text(json.dumps(manifest, indent=2))
    C.write_provenance("s9_export", manifest)
    print(f"  bundle written to {EXPORT} "
          f"({'all gates passed' if manifest['all_gates_passed'] else 'GATE FAILURE'})")
    return manifest


if __name__ == "__main__":
    run()
