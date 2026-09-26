"""Stage 3 - the screening library: ChEMBL compounds with clinical history.

One library, screened against every target in the panel. Using the same
compounds everywhere is what makes the cross-target selectivity view in p05
meaningful - a per-target library would confound "selective" with "only tested
here".

max_phase >= 1 keeps anything that has reached human trials, which is the
population a repurposing screen is actually about. Phase is retained per
compound so marketed (4) can be separated from early clinical (1) later rather
than being collapsed now.

Standardisation is byte-identical to p01's, because the applicability domain in
p04 measures similarity between library and training compounds. If the two were
standardised differently, every Tanimoto would be measuring a formatting
difference on top of a chemical one.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from p01_curate_panel import api, standardize  # noqa: E402  (shared, deliberately)

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
RESULTS = ROOT / "results"

MIN_PHASE = 1
MAX_MW = 650.0
PAGE = 1000


def fetch_clinical_molecules() -> pd.DataFrame:
    rows, offset = [], 0
    while True:
        res = api(
            "molecule",
            max_phase__gte=MIN_PHASE,
            limit=PAGE,
            offset=offset,
            only="molecule_chembl_id,pref_name,max_phase,molecule_type,"
            "molecule_structures,first_approval,withdrawn_flag",
        )
        batch = res.get("molecules", [])
        for m in batch:
            struct = m.get("molecule_structures") or {}
            rows.append(
                {
                    "chembl_id": m.get("molecule_chembl_id"),
                    "pref_name": m.get("pref_name"),
                    "max_phase": m.get("max_phase"),
                    "molecule_type": m.get("molecule_type"),
                    "first_approval": m.get("first_approval"),
                    "withdrawn_flag": m.get("withdrawn_flag"),
                    "smiles": struct.get("canonical_smiles"),
                }
            )
        meta = res.get("page_meta", {})
        offset += PAGE
        print(f"  fetched {len(rows)}/{meta.get('total_count', '?')}", flush=True)
        if not batch or offset >= meta.get("total_count", 0):
            break
    return pd.DataFrame(rows)


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    print(f"fetching ChEMBL molecules with max_phase >= {MIN_PHASE} ...")
    df = fetch_clinical_molecules()
    n_raw = len(df)

    audit = {"n_raw": n_raw}
    df = df[df["molecule_type"].eq("Small molecule")]
    audit["after_small_molecule_only"] = int(len(df))
    df = df[df["smiles"].notna()]
    audit["after_has_structure"] = int(len(df))

    std = df["smiles"].map(standardize)
    df["std_smiles"] = [s[0] for s in std]
    df["scaffold"] = [s[1] for s in std]
    df["mw"] = [s[2] for s in std]
    df = df[df["std_smiles"].notna()]
    audit["after_parseable"] = int(len(df))
    df = df[df["mw"] <= MAX_MW]
    audit["after_mw_filter"] = int(len(df))

    # One row per distinct standardised structure; keep the most advanced phase
    # and a readable name. Salts and formulations collapse onto their parent.
    df["max_phase"] = pd.to_numeric(df["max_phase"], errors="coerce")
    df = (
        df.sort_values(["max_phase", "pref_name"], ascending=[False, True])
        .drop_duplicates("std_smiles", keep="first")
        .reset_index(drop=True)
    )
    audit["after_dedup_by_structure"] = int(len(df))

    df.to_csv(DATA / "library_clinical.csv", index=False)
    phase = df["max_phase"].value_counts().sort_index(ascending=False).to_dict()
    print(
        f"\nlibrary: {len(df)} unique small molecules, "
        f"{df['scaffold'].nunique()} scaffolds"
    )
    print(f"  by max_phase: {phase}")
    print(f"  withdrawn: {int(df['withdrawn_flag'].fillna(False).astype(bool).sum())}")
    (RESULTS / "library_provenance.json").write_text(
        json.dumps(
            {
                "min_max_phase": MIN_PHASE,
                "max_mw": MAX_MW,
                "n_compounds": int(len(df)),
                "n_scaffolds": int(df["scaffold"].nunique()),
                "by_max_phase": {str(k): int(v) for k, v in phase.items()},
                "audit": audit,
                "licence": "ChEMBL data is CC BY-SA 3.0 - attribute and share alike",
            },
            indent=2,
        )
    )
    print(f"wrote {DATA / 'library_clinical.csv'}")


if __name__ == "__main__":
    main()
