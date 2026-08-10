#!/usr/bin/env python
"""Is the H&E TIL label any good? Check it against CD8 mRNA in the same patients.

This is the diagnostic that decides how to read the null result in
``run_real_experiment.py``.  Two very different conclusions fit an AUC of 0.54:

* the **imaging** carries no signal for lymphocyte burden, or
* the **label** is too noisy to detect signal with.

They are separable, because the same 91 patients have RNA-seq. If the H&E-derived
TIL fraction tracks CD8A/CD8B and a cytolytic signature, the label is measuring
what it claims to and the null belongs to the imaging. If it does not, the null
says nothing about MRI at all.

Uses the GDC STAR-Counts tables (TPM), CC-BY open access.
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
from scipy.stats import spearmanr  # noqa: E402

from mri_cd8_til.labels.transcriptome import SIGNATURES  # noqa: E402

#: Genes read out of every STAR-Counts table.
WANTED = sorted({g for genes in SIGNATURES.values() for g in genes} | {"PTPRC", "MS4A1", "EPCAM"})


def read_star_counts(path: Path, genes: set[str]) -> dict[str, float]:
    """Pull TPM for the requested genes out of one GDC STAR-Counts table.

    The format carries four ``N_*`` summary rows before the header line, so the
    parser keys on the column names rather than on a fixed row offset.
    """
    values: dict[str, float] = {}
    with gzip.open(path, "rt") as handle:
        header_columns: list[str] | None = None
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if header_columns is None:
                if fields[0] == "gene_id":
                    header_columns = fields
                continue
            if fields[0].startswith("N_"):
                continue
            name = fields[header_columns.index("gene_name")]
            if name in genes:
                tpm = fields[header_columns.index("tpm_unstranded")]
                try:
                    values[name] = float(tpm)
                except ValueError:
                    continue
    return values


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rna-dir", required=True, help="directory of <PATIENT>.tsv.gz files")
    parser.add_argument("--labels", default="results/til_labels.json")
    parser.add_argument("--out", default="results/label_validation.json")
    args = parser.parse_args()

    labels = {
        r["patient_id"]: r for r in json.loads(Path(args.labels).read_text())["records"]
    }
    wanted = set(WANTED)

    expression: dict[str, dict[str, float]] = {}
    for path in sorted(Path(args.rna_dir).glob("*.tsv.gz")):
        patient = path.name.replace(".tsv.gz", "")
        if patient not in labels:
            continue
        values = read_star_counts(path, wanted)
        if values:
            expression[patient] = values

    patients = sorted(expression)
    print(f"patients with both a TIL label and RNA-seq: {len(patients)}")
    if len(patients) < 20:
        print("too few to correlate")
        return 1

    til = np.array([labels[p]["til_fraction_tissue"] for p in patients])
    moran = np.array([labels[p]["til_moran_i"] for p in patients])

    def tpm(gene: str) -> np.ndarray:
        return np.array([expression[p].get(gene, np.nan) for p in patients])

    def signature(genes) -> np.ndarray:
        matrix = np.vstack([np.log2(tpm(g) + 1) for g in genes if np.isfinite(tpm(g)).any()])
        centred = (matrix - matrix.mean(axis=1, keepdims=True)) / np.where(
            matrix.std(axis=1, ddof=1, keepdims=True) > 0,
            matrix.std(axis=1, ddof=1, keepdims=True),
            1.0,
        )
        return centred.mean(axis=0)

    print("\nSpearman correlation of the H&E TIL label with transcriptome immune measures")
    print("(the label is derived from pathology; the transcriptome is an independent assay)\n")

    checks = {
        "CD8A": np.log2(tpm("CD8A") + 1),
        "CD8B": np.log2(tpm("CD8B") + 1),
        "PTPRC (CD45, pan-leukocyte)": np.log2(tpm("PTPRC") + 1),
        "CD8 core signature": signature(SIGNATURES["cd8_core"]),
        "cytolytic (GZMA, PRF1)": signature(SIGNATURES["cytolytic"]),
        "TIS (18-gene inflamed)": signature(SIGNATURES["tis"]),
        "EPCAM (epithelial, negative control)": np.log2(tpm("EPCAM") + 1),
    }

    results = {}
    for name, values in checks.items():
        ok = np.isfinite(values) & np.isfinite(til)
        rho_amount = spearmanr(til[ok], values[ok])
        rho_arrangement = spearmanr(moran[ok], values[ok])
        results[name] = {
            "n": int(ok.sum()),
            "rho_til_fraction": float(rho_amount.statistic),
            "p_til_fraction": float(rho_amount.pvalue),
            "rho_moran": float(rho_arrangement.statistic),
            "p_moran": float(rho_arrangement.pvalue),
        }
        star = "***" if rho_amount.pvalue < 0.001 else "**" if rho_amount.pvalue < 0.01 else "*" if rho_amount.pvalue < 0.05 else "  "
        print(
            f"  {name:<38} TIL amount rho={rho_amount.statistic:+.3f} "
            f"p={rho_amount.pvalue:.2e} {star}   "
            f"| arrangement rho={rho_arrangement.statistic:+.3f} p={rho_arrangement.pvalue:.3f}"
        )

    cd8 = results["CD8 core signature"]["rho_til_fraction"]
    epcam = results["EPCAM (epithelial, negative control)"]["rho_til_fraction"]
    print("\nReading:")
    if cd8 > 0.3 and results["CD8 core signature"]["p_til_fraction"] < 0.01:
        print(
            f"  The H&E TIL label tracks the CD8 signature at rho={cd8:+.3f}, while the "
            f"epithelial\n  negative control sits at {epcam:+.3f}. The label measures lymphocyte "
            "burden as claimed,\n  so the null in the imaging experiment belongs to the imaging, "
            "not to the label."
        )
    else:
        print(
            f"  The H&E TIL label tracks the CD8 signature at only rho={cd8:+.3f}. The label is "
            "too\n  noisy to support a strong conclusion about the imaging."
        )

    payload = {"n_patients": len(patients), "correlations": results}
    path = Path(args.out)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2))
    print(f"\nwrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
