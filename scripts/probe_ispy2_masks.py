#!/usr/bin/env python
"""Measure what I-SPY2's published SEG objects actually contain, and save it.

The claim in ``docs/ISPY2.md`` is that these are breast analysis masks rather
than tumour segmentations, and that two ways of deriving a tumour ROI from them
both fail.  This script produces the measurements that claim rests on, so the
figures and the prose read from one file instead of from a transcript.

For each patient it records, per bit plane: volume, connected-component count,
and the share of the plane held by its largest component.  A tumour-shaped plane
would be a small fraction of the breast, few components, one of them dominant.
Then it records the two candidate derivations so their failure is quantified
rather than asserted.

Usage::

    python scripts/probe_ispy2_masks.py --n 20
"""

from __future__ import annotations

import argparse
import json
import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
warnings.filterwarnings("ignore")

import numpy as np  # noqa: E402

from mri_cd8_til.clients.idc import IDCClient  # noqa: E402
from mri_cd8_til.experiments import ispy2  # noqa: E402

BIT_PLANES = (0, 1, 4, 5)
#: Published I-SPY2 pre-treatment functional tumour volume, for comparison.
PUBLISHED_FTV_MEDIAN_CM3 = 27.5


def profile(mask: np.ndarray, voxel_cm3: float) -> dict:
    from scipy.ndimage import label

    if not mask.any():
        return {"volume_cm3": 0.0, "n_components": 0, "largest_share": 0.0}
    labelled, n = label(mask)
    sizes = np.bincount(labelled.ravel())[1:]
    return {
        "volume_cm3": float(mask.sum() * voxel_cm3),
        "n_components": int(n),
        "largest_share": float(sizes.max() / mask.sum()) if len(sizes) else 0.0,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n", type=int, default=20)
    parser.add_argument("--out", default="results/ispy2_mask_probe.json")
    args = parser.parse_args()

    from idc_index import IDCClient as Index

    index = Index().index
    cohort = index[index.collection_id.astype(str) == "ispy2"]
    client = IDCClient()

    records = []
    attempted = 0
    for patient in sorted(cohort.PatientID.unique()):
        if len(records) >= args.n:
            break
        attempted += 1
        rows = ispy2.pretreatment_series(cohort, patient)
        if rows is None:
            continue
        try:
            case = ispy2.load_case(client, rows)
        except Exception as exc:  # noqa: BLE001
            print(f"  {patient}: load failed ({type(exc).__name__})")
            continue
        if case is None:
            continue

        voxel_cm3 = float(np.prod(case.spacing) / 1000.0)
        # Re-read the raw mask bits: load_case keeps only the VOI plane.
        mask_row = rows[rows.SeriesDescription == ispy2.ANALYSIS_MASK].iloc[0]
        raw = ispy2._read_series(client, mask_row.crdc_series_uuid)
        import SimpleITK as sitk

        resampled = sitk.Resample(
            raw,
            sitk.GetImageFromArray(case.ser),
            sitk.Transform(),
            sitk.sitkNearestNeighbor,
            0,
        )
        del resampled  # geometry already validated in load_case; use its VOI below
        bits = {}
        raw_array = sitk.GetArrayFromImage(raw).astype(int)
        raw_voxel_cm3 = float(np.prod(raw.GetSpacing()) / 1000.0)
        for bit in BIT_PLANES:
            bits[f"bit{bit}"] = profile(((raw_array >> bit) & 1).astype(bool), raw_voxel_cm3)

        enhancing = case.voi & (case.pe2 >= ispy2.PE_THRESHOLD)
        ftv = case.functional_tumour_volume()
        record = {
            "patient_id": case.patient_id,
            "study_date": case.study_date,
            "spacing": list(case.spacing),
            "bits": bits,
            "voi_volume_cm3": case.volume_cm3(case.voi),
            "enhancing_volume_cm3": case.volume_cm3(enhancing),
            "derived_ftv_cm3": case.volume_cm3(ftv),
            "derived_ftv_share_of_enhancing": float(
                ftv.sum() / enhancing.sum() if enhancing.any() else 0.0
            ),
        }
        records.append(record)
        print(
            f"  [{len(records)}/{args.n}] {case.patient_id}  "
            f"VOI {record['voi_volume_cm3']:6.0f} cm3  "
            f"derived FTV {record['derived_ftv_cm3']:7.1f} cm3"
        )

    if not records:
        print("no cases loaded")
        return 1

    ftv = np.array([r["derived_ftv_cm3"] for r in records])
    print(f"\nprobed {len(records)} patients (from {attempted} examined)")
    print(
        f"  derived FTV: median {np.median(ftv):.1f} cm3  "
        f"IQR [{np.quantile(ftv, .25):.1f}, {np.quantile(ftv, .75):.1f}]  max {ftv.max():.1f}"
    )
    print(f"  published I-SPY2 pre-treatment FTV median ~{PUBLISHED_FTV_MEDIAN_CM3} cm3")
    implausible = int((ftv > 100).sum())
    print(
        f"  {implausible}/{len(ftv)} patients exceed 100 cm3 -- lesion merged with "
        "background parenchymal enhancement"
    )

    payload = {
        "n_patients": len(records),
        "published_ftv_median_cm3": PUBLISHED_FTV_MEDIAN_CM3,
        "derived_ftv_median_cm3": float(np.median(ftv)),
        "derived_ftv_iqr": [float(np.quantile(ftv, 0.25)), float(np.quantile(ftv, 0.75))],
        "derived_ftv_max_cm3": float(ftv.max()),
        "n_implausible_over_100cm3": implausible,
        "records": records,
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2))
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
