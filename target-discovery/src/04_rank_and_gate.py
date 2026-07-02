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

    Expects columns including a gene identifier and a shift score. Column names
    from geneformer stats vary by version; we normalize the common ones.
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
        df["cell_state"] = name
        frames.append(df)
    if not frames:
        raise SystemExit("No shift tables found — run Step 3 first.")
    return pd.concat(frames, ignore_index=True)


def _normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    ren = {}
    for c in df.columns:
        cl = c.lower()
        if cl in ("gene_name", "gene", "gene_symbol"):
            ren[c] = "gene_symbol"
        elif "shift" in cl and "goal" in cl:
            ren[c] = "goal_shift"
        elif cl in ("shift_to_goal_end", "test_stat"):
            ren.setdefault(c, "goal_shift")
    df = df.rename(columns=ren)
    if "gene_symbol" not in df or "goal_shift" not in df:
        raise ValueError(f"Unexpected stats columns: {list(df.columns)}")
    return df[["gene_symbol", "goal_shift"]]


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
