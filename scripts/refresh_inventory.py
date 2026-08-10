#!/usr/bin/env python
"""Re-measure the public-data inventory against the live TCIA/GDC/GEO APIs.

Everything in ``mri_cd8_til.cohorts.REGISTRY`` is a *claim*.  This script turns
each claim back into a measurement and reports drift, so the design document
never silently rots.

Usage::

    python scripts/refresh_inventory.py                 # full refresh
    python scripts/refresh_inventory.py --skip-linkage  # TCIA counts only (fast)
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from mri_cd8_til.clients.gdc import GDCClient  # noqa: E402
from mri_cd8_til.clients.geo import GEOClient  # noqa: E402
from mri_cd8_til.clients.http import PublicDataSession  # noqa: E402
from mri_cd8_til.clients.tcia import TCIAClient  # noqa: E402
from mri_cd8_til.cohorts import REGISTRY  # noqa: E402

TCIA_COLLECTIONS = {
    "tcga_brca": "TCGA-BRCA",
    "ispy1": "ISPY1",
    "ispy2": "ISPY2",
    "acrin_6698": "ACRIN-6698",
    "duke": "Duke-Breast-Cancer-MRI",
    "advanced_mri_breast": "Advanced-MRI-Breast-Lesions",
    "nact_pilot": "Breast-MRI-NACT-Pilot",
    "qin_breast": "QIN-BREAST",
}


def _drift(label: str, claimed: int, measured: int) -> str:
    flag = "ok" if claimed == measured else f"DRIFT {measured - claimed:+d}"
    return f"  {label:<28} claimed={claimed:<6} measured={measured:<6} {flag}"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default="results/cohort_inventory.json")
    parser.add_argument("--skip-linkage", action="store_true")
    args = parser.parse_args()

    session = PublicDataSession()
    tcia = TCIAClient(session)
    report: dict[str, object] = {}

    print("=== TCIA collection counts ===")
    tcia_report = {}
    for key, collection in TCIA_COLLECTIONS.items():
        try:
            summary = tcia.summarize(collection)
        except Exception as exc:  # noqa: BLE001
            print(f"  {collection:<28} ERROR {exc}")
            continue
        cohort = REGISTRY[key]
        print(_drift(collection, cohort.n_with_mri, summary.n_patients_with_mr))
        tcia_report[key] = {
            "collection": collection,
            "n_patients": summary.n_patients,
            "n_patients_with_mr": summary.n_patients_with_mr,
            "n_mr_series": summary.n_mr_series,
            "modalities": summary.modalities,
            "n_distinct_scanners": summary.n_distinct_scanners,
            "scanner_counts": summary.scanner_counts,
        }
    report["tcia"] = tcia_report

    if args.skip_linkage:
        _write(report, args.out)
        return 0

    print("\n=== TCGA-BRCA pairing (tier A) ===")
    gdc = GDCClient(session)
    availability = gdc.availability("TCGA-BRCA")
    mri_ids = tcia.patients_with_modality("TCGA-BRCA", "MR")
    paired = mri_ids & availability.rna_seq & availability.diagnostic_slides
    print(f"  GDC availability: {availability.summary()}")
    print(_drift("MRI & RNA-seq & diag WSI", REGISTRY["tcga_brca"].n_paired, len(paired)))
    report["tcga_brca_linkage"] = {
        "n_mri": len(mri_ids),
        "gdc": availability.summary(),
        "n_mri_and_rna": len(mri_ids & availability.rna_seq),
        "n_mri_and_diag_wsi": len(mri_ids & availability.diagnostic_slides),
        "n_fully_paired": len(paired),
        "paired_ids": sorted(paired),
    }

    print("\n=== I-SPY2 pairing (tier B) ===")
    geo = GEOClient(session)
    series = geo.series("GSE194040")
    geo_ids = series.patient_ids("patient id")
    mri_numeric = {pid.split("-")[-1] for pid in tcia.patients("ISPY2")}
    linked = mri_numeric & geo_ids
    print(f"  GSE194040 samples: {series.n_samples}, distinct patient ids: {len(geo_ids)}")
    print(_drift("ISPY2 MRI & expression", REGISTRY["ispy2"].n_paired, len(linked)))
    report["ispy2_linkage"] = {
        "n_mri": len(mri_numeric),
        "n_geo_samples": series.n_samples,
        "n_geo_patients": len(geo_ids),
        "n_linked": len(linked),
        "linked_ids": sorted(linked),
    }

    _write(report, args.out)
    return 0


def _write(report: dict[str, object], out: str) -> None:
    path = Path(out)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2, sort_keys=True))
    print(f"\nwrote {path}")


if __name__ == "__main__":
    raise SystemExit(main())
