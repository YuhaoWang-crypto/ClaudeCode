"""Stage 8 - Boltz-2 as a screening filter, on the same compounds as docking.

The question people actually want answered about co-folding models is whether
their scores can stand in for an enzymatic assay as a first-pass filter. That
is a ranking question, so it is asked the same way the docking baseline is
asked: rank the compounds by the model's score and correlate with measured
pIC50, on exactly the compounds the docking ensemble also scored.

Two honest caveats are built into this stage rather than left out of it.

**What the library-screen endpoint returns is confidence, not affinity.** The
per-molecule output carries `iptm`, `binding_confidence`, `structure_confidence`
and `optimization_score` - interface and structure confidences - plus ADME
predictions. There is no calibrated affinity or predicted pIC50 in it. So every
number below is "does a Boltz-2 confidence score rank CDK2 potency", not "how
accurate is Boltz-2's predicted affinity". Reporting it the second way would be
an overclaim about what was measured.

**The default SMARTS alert filter silently dropped 35 of 100 compounds.** That
is worth knowing on its own if you plan to run a real library through this
endpoint - a third of a literature-derived CDK2 set never gets scored unless
the filter is turned off. Whether it also biases the comparison is testable, so
this stage tests it (Mann-Whitney on the potency of kept vs dropped) instead of
assuming either way.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, pearsonr, spearmanr

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
RESULTS = ROOT / "results"

ACTIVE_THRESHOLD = 7.0

BOLTZ_METRICS = [
    "boltz_binding_confidence",
    "boltz_iptm",
    "boltz_optimization_score",
    "boltz_structure_confidence",
    "boltz_complex_plddt",
]


def rank_report(name: str, y: np.ndarray, score: np.ndarray) -> dict:
    """Spearman / Pearson / EF5% for one ranking, higher score = better binder."""
    ok = np.isfinite(score) & np.isfinite(y)
    y, score = y[ok], score[ok]
    if len(y) < 10:
        return {"name": name, "n": int(len(y)), "note": "too few points"}
    actives = y >= ACTIVE_THRESHOLD
    n_top = max(1, int(round(0.05 * len(y))))
    ef5 = (
        float(actives[np.argsort(-score)[:n_top]].mean() / actives.mean())
        if actives.any()
        else None
    )
    rho, p_rho = spearmanr(y, score)
    return {
        "name": name,
        "n": int(len(y)),
        "spearman_rho": round(float(rho), 3),
        "spearman_p": float(f"{p_rho:.3g}"),
        "pearson_r": round(float(pearsonr(y, score)[0]), 3),
        "ef5_percent": round(ef5, 2) if ef5 is not None else None,
        "n_actives": int(actives.sum()),
    }


def main() -> None:
    boltz = pd.read_csv(RESULTS / "boltz" / "affinity_screen_metrics.csv")
    lig = pd.read_csv(DATA / "screening_set.csv")
    requested = lig.head(100)  # the set submitted to Boltz

    # --- what the filter removed, and whether it mattered ------------------
    kept = requested[requested["ligand_id"].isin(boltz["ligand_id"])]
    dropped = requested[~requested["ligand_id"].isin(boltz["ligand_id"])]
    u, p_filter = mannwhitneyu(kept["pchembl"], dropped["pchembl"])
    print(
        f"Boltz default SMARTS filter: {len(kept)} scored, {len(dropped)} dropped "
        f"of {len(requested)} submitted"
    )
    print(
        f"  kept    pIC50 mean {kept['pchembl'].mean():.2f} sd {kept['pchembl'].std():.2f}\n"
        f"  dropped pIC50 mean {dropped['pchembl'].mean():.2f} sd {dropped['pchembl'].std():.2f}\n"
        f"  Mann-Whitney p = {p_filter:.3f} -> "
        + (
            "the dropped set is NOT a potency-biased subset"
            if p_filter >= 0.05
            else "WARNING: the filter is potency-biased, the comparison is confounded"
        )
    )

    df = boltz.merge(lig[["ligand_id", "pchembl", "scaffold"]], on="ligand_id")
    y = df["pchembl"].to_numpy(dtype=float)

    # --- Boltz-2 confidence scores as rankers ------------------------------
    print("\nBoltz-2 confidence scores as rankers of measured pIC50")
    rows = [rank_report(m, y, df[m].to_numpy(dtype=float)) for m in BOLTZ_METRICS]
    for r in rows:
        if "spearman_rho" in r:
            print(
                f"  {r['name']:<30s} rho {r['spearman_rho']:+.3f} "
                f"(p={r['spearman_p']:<8}) EF5% {r['ef5_percent']}  n={r['n']}"
            )

    # --- the docking ensemble on exactly the same compounds ----------------
    dock_rows = []
    scores_path = RESULTS / "docking_scores.csv"
    overlap = pd.DataFrame()
    if scores_path.exists():
        scores = pd.read_csv(scores_path, index_col=0)
        overlap = df[df["ligand_id"].isin(scores.index)].copy()
        print(
            f"\ncompounds with both a Boltz score and a full docking ensemble: "
            f"{len(overlap)}"
        )
        if len(overlap) >= 10:
            yo = overlap["pchembl"].to_numpy(dtype=float)
            S = scores.loc[overlap["ligand_id"]].to_numpy(dtype=float)
            # smina scores are negative-is-better, so negate to rank.
            dock_rows.append(rank_report("dock_best_of_ensemble", yo, -np.nanmin(S, axis=1)))
            dock_rows.append(rank_report("dock_mean_of_ensemble", yo, -np.nanmean(S, axis=1)))
            for j, col in enumerate(scores.columns):
                dock_rows.append(rank_report(f"dock_{col}", yo, -S[:, j]))
            for m in BOLTZ_METRICS[:2]:
                dock_rows.append(
                    rank_report(f"{m}__same_subset", yo, overlap[m].to_numpy(dtype=float))
                )
            print("  paired comparison on that subset:")
            for r in dock_rows:
                if "spearman_rho" in r:
                    print(
                        f"    {r['name']:<34s} rho {r['spearman_rho']:+.3f} "
                        f"(p={r['spearman_p']:<8}) EF5% {r['ef5_percent']}"
                    )
    else:
        print("\nno docking_scores.csv yet - run s05 first for the paired comparison")

    out = {
        "what_boltz_returned": (
            "iptm / binding_confidence / structure_confidence / optimization_score "
            "plus ADME. No calibrated affinity or predicted pIC50 is available from "
            "the library-screen endpoint, so these are confidence rankers."
        ),
        "smarts_filter": {
            "n_submitted": int(len(requested)),
            "n_scored": int(len(kept)),
            "n_dropped": int(len(dropped)),
            "kept_pchembl_mean": round(float(kept["pchembl"].mean()), 3),
            "dropped_pchembl_mean": round(float(dropped["pchembl"].mean()), 3),
            "mannwhitney_p": float(f"{p_filter:.3g}"),
            "potency_biased": bool(p_filter < 0.05),
        },
        "boltz_rankers": rows,
        "paired_with_docking": dock_rows,
        "n_paired": int(len(overlap)),
        "active_threshold_pchembl": ACTIVE_THRESHOLD,
    }
    (RESULTS / "boltz_comparison.json").write_text(json.dumps(out, indent=2))
    print(f"\nwrote {RESULTS / 'boltz_comparison.json'}")


if __name__ == "__main__":
    main()
