"""Enumerate a Su-2026-style combinatorial ionizable-lipid library and write the
LiON screen input, so the lite model can score it on Modal.

Design (aza-Michael chemistry, the real ionizable-lipid route):
  14 amine heads  x  {ester, amide} linkers  x  16 tails
Every head N-H is alkylated with an acrylate/acrylamide tail, so the tail count =
the head's reactive amine count (2-4 here), mirroring Su 2026's 2/3/4-tail lipids.
These are *inspired by* (not identical to) Su's H1-H14 / T1-T16 — we don't have
their exact block SMILES, so we use a chemically valid, representative set.

Outputs:
  data/libraries/combo_v1/combo_v1.csv                 SMILES (LiON screen input)
  data/libraries/combo_v1/combo_v1_extra_x.csv         formulation features
  data/libraries/combo_v1/combo_v1_metadata.csv
  data/libraries/combo_v1/combo_v1_enumerated.csv      full block annotation
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from lipidlib.lion_library import enumerate_michael_lipids, build_library, Formulation

# --- 14 amine heads (id -> SMILES). Reactive amine N-H count sets the tail number.
HEADS = {
    "H1_DMAEA":   "CN(C)CCN",            # dimethylaminoethyl-amine   -> 2 tails
    "H2_DEAEA":   "CCN(CC)CCN",          # diethylaminoethyl-amine    -> 2
    "H3_MPZ":     "CN1CCN(CCN)CC1",      # methylpiperazinyl-ethylamine -> 2
    "H4_AEPZ":    "C1CNCCN1CCN",         # aminoethyl-piperazine      -> 3
    "H5_bAE":     "NCCOCCN",             # bis-aminoethyl ether       -> 4
    "H6_DAP":     "NCCCN",               # 1,3-diaminopropane         -> 4
    "H7_MePr":    "CNCCCN",              # N-methyl-1,3-diaminopropane -> 3
    "H8_AEEA":    "NCCNCCO",             # aminoethyl-ethanolamine    -> 3
    "H9_PZE":     "OCCN1CCN(CCN)CC1",    # hydroxyethyl-piperazinyl-ethylamine -> 2
    "H10_DMPDA":  "CN(C)CCCN",           # dimethylaminopropyl-amine  -> 2
    "H11_bDMA":   "CN(C)CCNCCN(C)C",     # bis(dimethylaminoethyl)amine -> 1(skip if <2)
    "H12_TRIS":   "NCC(CO)(CO)N",        # diamino-diol core          -> 4
    "H13_MXA":    "NCc1ccc(CN)cc1",      # xylylenediamine            -> 4
    "H14_AEM":    "CNCCN(C)C",           # N-methyl-dimethylaminoethylamine -> 1(skip)
}

# --- 16 tails (id -> alkyl/hetero R attached via the ester -O-R / amide -NH-R)
TAILS = {
    "T1_C6":   "CCCCCC",
    "T2_C8":   "CCCCCCCC",
    "T3_C10":  "CCCCCCCCCC",
    "T4_C12":  "CCCCCCCCCCCC",
    "T5_C14":  "CCCCCCCCCCCCCC",
    "T6_C16":  "CCCCCCCCCCCCCCCC",
    "T7_C18":  "CCCCCCCCCCCCCCCCCC",
    "T8_C8br": "CCCCCCC(C)",             # branched
    "T9_C10br":"CCCCCCCCC(C)",
    "T10_C12i":"CCCCCCCCCCC(C)C",        # iso
    "T11_oleyl":"CCCCCCCC/C=C\\CCCCCCCC", # unsaturated (oleyl)
    "T12_linoleyl":"CCCCC/C=C\\C/C=C\\CCCCCCCC",  # di-unsaturated
    "T13_C9":  "CCCCCCCCC",
    "T14_thioC8":"CCCCCCCCSCC",          # thioether (S)
    "T15_OC8": "CCCCCCCCOCC",            # ether (O)
    "T16_C11": "CCCCCCCCCCC",
}


def main() -> int:
    df = enumerate_michael_lipids(HEADS, TAILS, linkers=("ester", "amide"))
    print(f"enumerated {len(df)} valid lipids "
          f"(from {len(HEADS)}x2x{len(TAILS)} = {len(HEADS)*2*len(TAILS)} combos)")
    print(df["n_tails"].value_counts().rename("count").to_string())
    print(f"MW range: {df['mw'].min():.0f}-{df['mw'].max():.0f} (median {df['mw'].median():.0f})")

    out = ROOT / "data" / "libraries" / "combo_v1"
    out.mkdir(parents=True, exist_ok=True)
    df.to_csv(out / "combo_v1_enumerated.csv", index=False)

    # LiON screen input at the KK/HeLa baseline (stable predictions per the authors)
    build_library(df["smiles"].tolist(), out.parent, name="combo_v1",
                  formulation=Formulation())
    print(f"wrote screen input + annotation -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
