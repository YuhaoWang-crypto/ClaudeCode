#!/usr/bin/env python
"""Run the whole imaging path on one real TCGA-BRCA series, end to end.

Everything else in this project is either an API query or a simulation.  This
script closes the loop on real data: it downloads an open-access DCE
subtraction series from TCIA, preprocesses it, builds all four ROI variants and
extracts the full feature block, so the pipeline is demonstrated against real
DICOM geometry rather than against synthetic arrays.

What it does **not** do is segment the tumour.  The ROI here is a
top-percentile threshold on the subtraction image -- a surrogate whose only job
is to exercise the ROI and feature code on a real voxel grid.  A study would use
an expert contour, or MAMA-MIA's pretrained nnU-Net followed by human review.

Usage::

    python scripts/run_real_dicom_demo.py
    python scripts/run_real_dicom_demo.py --patient TCGA-AR-A1AQ --keep
"""

from __future__ import annotations

import argparse
import sys
import tempfile
import time
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
warnings.filterwarnings("ignore")

import numpy as np  # noqa: E402
import SimpleITK as sitk  # noqa: E402

from mri_cd8_til.clients.tcia import TCIAClient  # noqa: E402
from mri_cd8_til.features.extractor import ExtractionConfig, RadiomicFeatureExtractor  # noqa: E402
from mri_cd8_til.imaging.preprocess import PreprocessConfig, preprocess_volume  # noqa: E402
from mri_cd8_til.imaging.roi import (  # noqa: E402
    DEFAULT_REGIONS,
    body_mask_from_image,
    build_region,
    region_is_viable,
)

DEFAULT_PATIENT = "TCGA-AR-A1AQ"
#: TCIA encodes GE subtraction series as "(study/series/N)-(study/series/1)".
SUBTRACTION_MARKER = ")-("


def pick_series(client: TCIAClient, patient: str, max_images: int):
    series = [s for s in client.series("TCGA-BRCA", modality="MR") if s.patient_id == patient]
    if not series:
        raise SystemExit(f"no MR series found for {patient}")
    subtraction = [
        s for s in series if SUBTRACTION_MARKER in s.description and s.image_count <= max_images
    ]
    chosen = min(subtraction or series, key=lambda s: s.image_count)
    kind = "subtraction" if SUBTRACTION_MARKER in chosen.description else "fallback"
    print(f"   chose {kind} series '{chosen.description}' ({chosen.image_count} images)")
    print(f"   scanner: {chosen.scanner_key}")
    return chosen


def surrogate_lesion_mask(image: sitk.Image, body: sitk.Image, percentile: float) -> sitk.Image:
    """Largest connected component of the top-percentile enhancing voxels."""
    array = sitk.GetArrayFromImage(image)
    threshold = float(np.percentile(array[sitk.GetArrayFromImage(body) > 0], percentile))
    seed = sitk.And(sitk.Cast(image > threshold, sitk.sitkUInt8), body)
    seed = sitk.BinaryMorphologicalClosing(seed, [2, 2, 1])
    labelled = sitk.RelabelComponent(sitk.ConnectedComponent(seed), sortByObjectSize=True)
    return sitk.Cast(labelled == 1, sitk.sitkUInt8)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--patient", default=DEFAULT_PATIENT)
    parser.add_argument("--max-images", type=int, default=140)
    parser.add_argument("--percentile", type=float, default=99.5)
    parser.add_argument("--keep", action="store_true", help="keep the downloaded DICOM")
    parser.add_argument("--out", default="results/real_dicom_endtoend.txt")
    args = parser.parse_args()

    lines: list[str] = []

    def emit(text: str = "") -> None:
        print(text)
        lines.append(text)

    client = TCIAClient()
    emit("=" * 74)
    emit(f"END-TO-END ON REAL TCIA DATA -- {args.patient} (TCGA-BRCA, open access)")
    emit("=" * 74)
    chosen = pick_series(client, args.patient, args.max_images)

    workdir = Path(tempfile.mkdtemp(prefix="tcia_"))
    try:
        emit(f"\n   downloading into {workdir} ...")
        dicom_dir = client.download_series(chosen.series_uid, workdir / "dcm")

        reader = sitk.ImageSeriesReader()
        reader.SetFileNames(reader.GetGDCMSeriesFileNames(str(dicom_dir)))
        image = reader.Execute()
        emit(
            f"   native : size={image.GetSize()} "
            f"spacing={tuple(round(s, 3) for s in image.GetSpacing())}"
        )

        config = PreprocessConfig(
            target_spacing_mm=(0.6, 0.6, 2.0), apply_n4=True, normalization="percentile"
        )
        start = time.time()
        processed = preprocess_volume(image, config)
        emit(
            f"   after  : size={processed.GetSize()} "
            f"spacing={tuple(round(s, 2) for s in processed.GetSpacing())} "
            f"({time.time() - start:.1f}s)"
        )

        body = sitk.Cast(body_mask_from_image(processed), sitk.sitkUInt8)
        mask = surrogate_lesion_mask(processed, body, args.percentile)
        n_voxels = int(sitk.GetArrayFromImage(mask).sum())
        emit(f"\n   surrogate enhancing-lesion ROI: {n_voxels} voxels")
        emit("   (a threshold surrogate, NOT a clinical segmentation)")

        for spec in DEFAULT_REGIONS:
            region = build_region(mask, spec, body_mask=body)
            count = int(sitk.GetArrayFromImage(region).sum())
            emit(f"     {spec.name:<18} voxels={count:>7}  viable={region_is_viable(region)}")

        extractor = RadiomicFeatureExtractor(ExtractionConfig(bin_width=0.05))
        start = time.time()
        features = extractor.extract(processed, mask, args.patient, "whole_tumour", 1)
        finite = sum(1 for v in features.values.values() if np.isfinite(v))
        emit(
            f"\n   extracted {features.n_features} features in {time.time() - start:.0f}s "
            f"({finite} finite, expected {extractor.expected_n_features()})"
        )
        for key in (
            "shape_VoxelVolume",
            "original_glcm_Contrast",
            "original_ngtdm_Busyness",
            "wavelet-sym4-HHH_glrlm_RunEntropy",
        ):
            if key in features.values:
                emit(f"     {key:<40} = {features.values[key]:.5g}")

        emit(
            "\n   Note on runtime: extraction cost scales with ROI size, and this "
            "surrogate\n   ROI is far larger than a real lesion. Expect a genuine "
            "tumour VOI to be\n   several times faster."
        )
    finally:
        if not args.keep:
            import shutil

            shutil.rmtree(workdir, ignore_errors=True)
        else:
            emit(f"\n   DICOM kept at {workdir}")

    path = Path(args.out)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")
    print(f"\nwrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
