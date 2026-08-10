#!/usr/bin/env python
"""Build the I-SPY2 transcriptomic immune labels and verify the identifier link.

This is the label half of the sample-size step: 717 patients across ~20 trial
sites, against 91 in the TCGA-BRCA imaging cohort.

The identifier link needs checking rather than assuming.  The gene-level
expression matrix from GSE194040 has columns like ``756412``, which look like a
separate research ID and not the ``100899``-style patient IDs that TCIA uses.
They are in fact the same identifiers in a different order -- 986 of 988 column
headers appear in the GEO patient-id list -- and this script re-verifies that
before using it, because a silently wrong join would attach every label to the
wrong patient while looking perfectly healthy.

Usage::

    python scripts/build_ispy2_labels.py
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
warnings.filterwarnings("ignore")

import numpy as np  # noqa: E402

from mri_cd8_til.labels.transcriptome import SIGNATURES  # noqa: E402

EXPRESSION_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE194nnn/GSE194040/suppl/"
    "GSE194040_ISPY2ResID_AgilentGeneExp_990_FrshFrzn_meanCol_geneLevel_n988.txt.gz"
)
CONTROL_GENES = ("EPCAM", "PTPRC")


def read_expression(path: Path, wanted: set[str]) -> tuple[list[str], list[str], np.ndarray]:
    genes, rows = [], []
    with gzip.open(path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if fields[0] in wanted:
                genes.append(fields[0])
                rows.append(
                    [float(x) if x not in ("", "NA") else np.nan for x in fields[1:]]
                )
    matrix = np.array(rows)
    patients = header[1:] if len(header) == matrix.shape[1] + 1 else header
    return patients, genes, matrix


def zscore_signature(matrix: np.ndarray, genes: list[str], names) -> np.ndarray:
    rows = [genes.index(g) for g in names if g in genes]
    if not rows:
        raise KeyError(f"none of {names} present")
    sub = matrix[rows]
    sd = np.nanstd(sub, axis=1, ddof=1, keepdims=True)
    centred = (sub - np.nanmean(sub, axis=1, keepdims=True)) / np.where(sd > 0, sd, 1.0)
    return np.nanmean(centred, axis=0)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--expression", default="results/ispy2_expr.txt.gz")
    parser.add_argument("--out", default="results/ispy2_immune_scores.json")
    args = parser.parse_args()

    path = Path(args.expression)
    if not path.exists():
        print(f"expression matrix not found at {path}")
        print(f"download it with:\n  curl -sSL '{EXPRESSION_URL}' -o {path}")
        return 1

    wanted = {g for genes in SIGNATURES.values() for g in genes} | set(CONTROL_GENES)
    patients, genes, matrix = read_expression(path, wanted)
    print(f"expression: {matrix.shape[0]} signature genes x {matrix.shape[1]} samples")
    print(f"  genes found: {len(genes)}/{len(wanted)}")

    from idc_index import IDCClient

    index = IDCClient().index
    mri_ids = {
        str(p).split("-")[-1]
        for p in index[index.collection_id.astype(str) == "ispy2"].PatientID.unique()
    }
    linked = [i for i, p in enumerate(patients) if p in mri_ids]
    print("\nidentifier check:")
    print(f"  TCIA ispy2 patients with imaging: {len(mri_ids)}")
    print(f"  expression columns matching a TCIA patient ID: {len(linked)}")
    if len(linked) < 0.9 * len(mri_ids):
        print("  REFUSING: fewer than 90% link. The join is probably on the wrong key.")
        return 1

    scores = {name: zscore_signature(matrix, genes, gs) for name, gs in SIGNATURES.items()}
    for control in CONTROL_GENES:
        scores[f"{control.lower()}_control"] = zscore_signature(matrix, genes, [control])

    cd8 = scores["cd8_core"][linked]
    tis = scores["tis"][linked]
    epcam = scores["epcam_control"][linked]
    ptprc = scores["ptprc_control"][linked]

    print("\nsanity checks on the scores:")
    checks = [
        ("corr(CD8 core, TIS)", float(np.corrcoef(cd8, tis)[0, 1]), "high", lambda v: v > 0.5),
        (
            "corr(CD8 core, EPCAM epithelial control)",
            float(np.corrcoef(cd8, epcam)[0, 1]),
            "near zero or negative",
            lambda v: v < 0.2,
        ),
        (
            "corr(CD8 core, PTPRC pan-leukocyte)",
            float(np.corrcoef(cd8, ptprc)[0, 1]),
            "positive",
            lambda v: v > 0.2,
        ),
    ]
    all_pass = True
    for label, value, expectation, test in checks:
        ok = test(value)
        all_pass &= ok
        print(f"  {label:<44} {value:+.3f}   (expect {expectation}) {'ok' if ok else 'FAIL'}")

    payload = {
        "n_samples": int(matrix.shape[1]),
        "n_linked_to_mri": len(linked),
        "patients": [patients[i] for i in linked],
        "scores": {k: [float(v[i]) for i in linked] for k, v in scores.items()},
        "sanity_checks": {label: value for label, value, _, _ in checks},
        "all_checks_pass": bool(all_pass),
        "caveat": (
            "Transcriptomic labels support the inflamed_vs_rest contrast only. "
            "I-SPY2 has no public H&E whole-slide images, so no TIL maps, no geometry "
            "descriptors, and no desert-vs-excluded distinction exist for this cohort."
        ),
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2))
    print(f"\nwrote {out}  ({len(linked)} patients)")
    print(
        "\nNote: the imaging half of this cohort is blocked. I-SPY2's public SEG objects "
        "are\nbreast analysis masks, not tumour segmentations -- see docs/ISPY2.md."
    )
    return 0 if all_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
