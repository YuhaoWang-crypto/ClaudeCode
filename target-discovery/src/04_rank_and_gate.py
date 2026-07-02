"""
Step 4 — Rank perturbation hits and apply the novelty gate.

Combines per-state goal-shift tables from Step 3 into a single ranked list, then
DOWN-WEIGHTS genes whose mechanism is already crowded in the IPF clinic
(data/crowded_targets.csv, derived from ClinicalTrials.gov). The point of the
whole exercise is first-in-class novelty, so being crowded is a penalty.

No GPU. Output: results/ranked_targets.csv
"""
from __future__ import annotations
import argparse
import os

import pandas as pd
import yaml


def load_config(path: str) -> dict:
    with open(path) as fh:
        return yaml.safe_load(fh)


def load_shift_tables(cfg: dict) -> pd.DataFrame:
    """Read each state's goal_state_shift stats and stack them.

    Verified against geneformer's InSilicoPerturberStats(mode="goal_state_shift")
    output columns: Gene, Gene_name, Ensembl_ID, Shift_to_goal_end,
    Goal_end_vs_random_pval, Goal_end_FDR, N_Detections, Sig.
    """
    frames = []
    for state in cfg["target_cell_states"]:
        name = state["name"]
        path = os.path.join(
            cfg["perturbation"]["out_dir"], name, f"stats_{name}.csv"
        )
        if not os.path.exists(path):
            print(f"[warn] missing {path}; skipping {name}")
            continue
        df = pd.read_csv(path)
        df = _normalize_columns(df)
        df = _apply_significance(df, cfg)
        df["cell_state"] = name
        frames.append(df)
    if not frames:
        raise SystemExit("No shift tables found — run Step 3 first.")
    return pd.concat(frames, ignore_index=True)


# Confirmed geneformer goal_state_shift column names -> our canonical names.
_COL_MAP = {
    "Gene_name": "gene_symbol",
    "Ensembl_ID": "ensembl_id",
    "Shift_to_goal_end": "goal_shift",
    "Goal_end_FDR": "fdr",
    "Sig": "sig",
}


def _normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    # Case-insensitive rename so minor version casing differences still map.
    lower = {c.lower(): c for c in df.columns}
    ren = {lower[k.lower()]: v for k, v in _COL_MAP.items() if k.lower() in lower}
    df = df.rename(columns=ren)
    if "gene_symbol" not in df or "goal_shift" not in df:
        raise ValueError(
            f"Unexpected stats columns: {list(df.columns)} "
            f"(expected geneformer goal_state_shift output)"
        )
    keep = [c for c in ("gene_symbol", "ensembl_id", "goal_shift", "fdr", "sig")
            if c in df.columns]
    return df[keep]


def _apply_significance(df: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    """Keep only significant knockouts before ranking, if the flags are present."""
    fdr_max = cfg["ranking"].get("fdr_max", 0.05)
    if "sig" in df.columns:
        df = df[df["sig"] == 1]
    elif "fdr" in df.columns:
        df = df[df["fdr"] <= fdr_max]
    return df


def apply_novelty_gate(df: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    crowded = pd.read_csv(
        cfg["ranking"]["crowded_targets_csv"], comment="#"
    )
    crowd_map = dict(zip(crowded["gene_symbol"], crowded["crowding"]))
    penalty = cfg["ranking"]["novelty_penalty"]
    weights = {"high": penalty, "medium": (1 + penalty) / 2, "low": 0.9}

    def novelty_weight(g: str) -> float:
        return weights.get(crowd_map.get(g), 1.0)  # unknown => fully novel (1.0)

    df["novelty_weight"] = df["gene_symbol"].map(novelty_weight)
    df["crowding"] = df["gene_symbol"].map(lambda g: crowd_map.get(g, "novel"))
    df["gated_score"] = df["goal_shift"] * df["novelty_weight"]
    return df


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config.yaml")
    args = ap.parse_args()
    cfg = load_config(args.config)

    df = load_shift_tables(cfg)
    # Take the strongest shift across states per gene, then gate for novelty.
    df = (
        df.sort_values("goal_shift", ascending=False)
        .drop_duplicates("gene_symbol", keep="first")
        .reset_index(drop=True)
    )
    df = apply_novelty_gate(df, cfg)
    df = df.sort_values("gated_score", ascending=False).head(cfg["ranking"]["top_n"])

    out = cfg["ranking"]["out_csv"]
    os.makedirs(os.path.dirname(out), exist_ok=True)
    df.to_csv(out, index=False)
    print(f"[done] top {len(df)} gated targets -> {out}")
    print(df.head(15).to_string(index=False))


if __name__ == "__main__":
    main()
