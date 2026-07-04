"""Expanded combinatorial library (v2): adds Su-2026-style degradable tails
(internal esters, disulfides) and more amine heads, for a deeper organ-focused
screen. Reuses the v1 blocks + the aza-Michael enumerator.

18 heads x {ester, amide} x 20 tails. Degradable tails matter because internal
esters and disulfides give reduction-/hydrolysis-triggered cargo release and
better tolerability (as in SM-102 / ALC-0315).
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from lipidlib.lion_library import enumerate_michael_lipids, build_library, BASELINES
from analysis.enumerate_library import HEADS as HEADS_V1, TAILS as TAILS_V1

# --- new amine heads --------------------------------------------------------
NEW_HEADS = {
    "H15_bHEAEA": "OCCN(CCO)CCN",   # bis-hydroxyethyl aminoethylamine (polar -OH head)
    "H16_DMABA":  "CN(C)CCCCN",     # dimethylaminobutyl-amine (longer spacer)
    "H17_CHDA":   "NC1CCCCC1N",     # cyclohexane-1,2-diamine (rigid, 4-tail)
    "H18_TETA":   "NCCNCCNCCN",     # triethylenetetramine (large polyamine)
}
# --- new degradable / reduction-cleavable tails -----------------------------
NEW_TAILS = {
    "T17_SS_C6":  "CCCCCCSSCCCCCC",         # disulfide (glutathione-cleavable)
    "T18_SS_C8":  "CCCCCCCCSSCCCCCCCC",     # longer disulfide
    "T19_estC8":  "CCCCCCCC(=O)OCCCCCC",    # internal ester (degradable)
    "T20_estBr":  "CCCCCC(=O)OCCCCCCCC",    # branched internal ester
}

HEADS = {**HEADS_V1, **NEW_HEADS}
TAILS = {**TAILS_V1, **NEW_TAILS}


def main() -> int:
    df = enumerate_michael_lipids(HEADS, TAILS, linkers=("ester", "amide"),
                                  mw_range=(400.0, 1500.0))
    # tag degradable chemistry
    deg_tails = set(NEW_TAILS) | {"T14_thioC8", "T15_OC8"}
    df["degradable"] = df["tail"].isin(deg_tails)
    df["disulfide"] = df["tail"].isin({"T17_SS_C6", "T18_SS_C8"})
    print(f"enumerated {len(df)} valid lipids from {len(HEADS)}x2x{len(TAILS)} "
          f"= {len(HEADS)*2*len(TAILS)} combos")
    print(f"  degradable-tail lipids: {int(df['degradable'].sum())} "
          f"(disulfide: {int(df['disulfide'].sum())})")
    print(df["n_tails"].value_counts().rename("count").to_string())

    out = ROOT / "data" / "libraries" / "combo_v2"
    out.mkdir(parents=True, exist_ok=True)
    df.to_csv(out / "combo_v2_enumerated.csv", index=False)
    build_library(df["smiles"].tolist(), out.parent, name="combo_v2",
                  formulation=BASELINES["spleen_IV"])  # spleen-focused context
    print(f"wrote {len(df)} lipids + spleen screen input -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
