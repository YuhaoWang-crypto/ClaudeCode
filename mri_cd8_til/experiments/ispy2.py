"""I-SPY2 cohort assembly: 717 patients with MRI and a transcriptomic immune label.

This is the sample-size step. The TCGA-BRCA cohort tops out at 91 patients with
imaging features; I-SPY2 has 719 with MRI on TCIA, 717 of which link to the
GSE194040 expression resource by patient ID, across roughly 20 trial sites.

**Two things change, and both are losses that have to be stated.**

*The label stops being spatial.* I-SPY2 has no H&E whole-slide images in a public
repository, so the Saltz TIL maps -- and therefore the geometry descriptors and
the desert/excluded/inflamed trichotomy -- do not exist for it. What exists is a
pre-treatment biopsy transcriptome, which supports the ``inflamed_vs_rest``
contrast and nothing finer. This is the constraint ``labels.schema`` has encoded
from the start, arriving in practice.

*The tumour ROI stops being an expert contour.* The trial publishes a
``VOLSER Analysis Mask`` per study, and it is tempting to read that as a tumour
segmentation. It is not: its DICOM segment is labelled "VOLSER Analysis Mask"
with segmented property *Breast*, it covers ~4,400 cm3, and none of its bit
planes is tumour-shaped -- the sparsest one scatters 2.7 cm3 across 2,164
disconnected components. It marks the breast region and threshold flags that
feed the trial's functional tumour volume computation.

Two ways of deriving a ROI from what is published were tried and **both were
rejected**, so this module deliberately stops short of producing one:

* *FTV by the trial's own rule* (VOI and PE above threshold, largest component)
  gives a median of 30.3 cm3, which matches published I-SPY2 functional tumour
  volume closely -- but an IQR of [5.3, 225.8] and a maximum of 493 cm3. Where
  background parenchymal enhancement is strong the lesion merges with the
  parenchyma into a single component, and the extracted features become features
  of the whole breast. Nothing errors; the corruption is silent and affects only
  some patients, which is the worst case for a downstream model;
* *bit 1 as the trial's own FTV flag* tracks that derived volume across patients
  (r = +0.855) but is a scattered voxel flag -- thousands of components, largest
  holding under 1% -- and a texture matrix over disconnected single voxels is
  not a measurement of anything.

:func:`functional_tumour_volume` implements the first route because it is worth
having reproducible, and because someone should be able to re-derive the numbers
above rather than take them on trust. **It should not be used as a ROI for
feature extraction.** Doing that properly needs the MAMA-MIA expert
segmentations (Synapse ``syn60868042``); see ``docs/ISPY2.md``.
"""

from __future__ import annotations

import shutil
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np

#: Series produced by the trial's VOLSER pipeline, cropped and co-registered.
ANALYSIS_MASK = "ISPY2: VOLSER: uni-lateral cropped: Analysis Mask"
SER_SERIES = "ISPY2: VOLSER: uni-lateral cropped: SER"
PE2_SERIES = "ISPY2: VOLSER: uni-lateral cropped: PE2"

#: Bit plane of the analysis mask that carries the breast analysis VOI. Measured
#: rather than assumed: it is the only plane forming a single connected
#: component of breast-scale volume across the patients checked.
VOI_BIT = 5

#: The trial's own percent-enhancement threshold, taken from the mask's own
#: DerivationDescription ("VOLSER,BACKGROUND=475,SER PE THRESHOLD=70").
PE_THRESHOLD = 70


@dataclass
class ISPY2Case:
    """One patient's pre-treatment imaging, already reduced to arrays."""

    patient_id: str
    study_date: str
    ser: np.ndarray
    pe2: np.ndarray
    voi: np.ndarray
    spacing: tuple[float, float, float]

    @property
    def shape(self) -> tuple[int, ...]:
        return self.ser.shape

    def functional_tumour_volume(
        self, pe_threshold: int = PE_THRESHOLD, min_component_fraction: float = 0.05
    ) -> np.ndarray:
        """Derived tumour ROI: enhancing voxels inside the VOI, largest component.

        Background parenchymal enhancement scatters hundreds of small components
        through the breast; keeping only the largest is what separates the lesion
        from it.  ``min_component_fraction`` guards the degenerate case where the
        largest component is itself tiny, which means no lesion was found rather
        than a small lesion.
        """
        from scipy.ndimage import label

        enhancing = self.voi & (self.pe2 >= pe_threshold)
        if not enhancing.any():
            return np.zeros_like(enhancing)
        labelled, n = label(enhancing)
        if n == 0:
            return np.zeros_like(enhancing)
        sizes = np.bincount(labelled.ravel())[1:]
        largest = int(sizes.argmax()) + 1
        if sizes.max() / enhancing.sum() < min_component_fraction:
            return np.zeros_like(enhancing)
        return labelled == largest

    def volume_cm3(self, mask: np.ndarray) -> float:
        return float(mask.sum() * np.prod(self.spacing) / 1000.0)


def pretreatment_series(index, patient_id: str):
    """The earliest study for a patient, restricted to the three VOLSER series.

    I-SPY2 images each patient at up to four timepoints; the GSE194040
    transcriptome comes from a **pre-treatment** biopsy, so anything later would
    pair imaging after therapy with a label from before it.
    """
    rows = index[index.PatientID == patient_id].copy()
    if rows.empty:
        return None
    rows["StudyDate"] = rows["StudyDate"].astype(str)
    first = sorted(rows.StudyDate.unique())[0]
    study = rows[rows.StudyDate == first]
    wanted = study[study.SeriesDescription.isin([ANALYSIS_MASK, SER_SERIES, PE2_SERIES])]
    if wanted.SeriesDescription.nunique() < 3:
        return None
    return wanted.drop_duplicates(subset="SeriesDescription")


def load_case(client, series_rows) -> ISPY2Case | None:
    """Download and decode one patient's three series onto a common grid.

    Geometry matters here.  The analysis mask is stored with more slices than
    the SER and PE volumes, and centre-cropping it to match -- the obvious
    shortcut -- is a guess that silently misregisters every feature downstream.
    Instead each series is read with its DICOM geometry and the mask is resampled
    onto the SER grid through the shared frame of reference.
    """
    import SimpleITK as sitk

    volumes: dict[str, "sitk.Image"] = {}
    for _, row in series_rows.iterrows():
        image = _read_series(client, row.crdc_series_uuid)
        if image is None:
            return None
        volumes[row.SeriesDescription] = image

    if len(volumes) < 3:
        return None

    reference = volumes[SER_SERIES]
    mask_image = sitk.Resample(
        volumes[ANALYSIS_MASK], reference, sitk.Transform(), sitk.sitkNearestNeighbor, 0
    )
    pe_image = sitk.Resample(
        volumes[PE2_SERIES], reference, sitk.Transform(), sitk.sitkLinear, 0
    )

    mask = sitk.GetArrayFromImage(mask_image).astype(np.int32)
    voi = ((mask >> VOI_BIT) & 1).astype(bool)
    patient_id = str(series_rows.iloc[0].PatientID)
    return ISPY2Case(
        patient_id=patient_id,
        study_date=str(series_rows.iloc[0].StudyDate),
        ser=sitk.GetArrayFromImage(reference).astype(np.float32),
        pe2=sitk.GetArrayFromImage(pe_image).astype(np.float32),
        voi=voi,
        spacing=tuple(float(s) for s in reference.GetSpacing()[:3]),
    )


def _read_series(client, crdc_series_uuid: str):
    """Fetch one IDC series and read it as a SimpleITK image, geometry intact."""
    import pydicom
    import SimpleITK as sitk

    keys = client._object_keys(crdc_series_uuid)
    if not keys:
        return None
    directory = Path(tempfile.mkdtemp())
    try:
        for key in keys:
            payload = client.session.get(
                f"https://idc-open-data.s3.amazonaws.com/{key}"
            ).content
            (directory / key.split("/")[-1]).write_bytes(payload)

        files = sorted(directory.iterdir())
        if len(files) == 1:
            # Multi-frame object (the SEG). ImageSeriesReader mishandles these,
            # so build the image from the pixel array plus its stated geometry.
            dataset = pydicom.dcmread(str(files[0]))
            return _image_from_multiframe(dataset)

        reader = sitk.ImageSeriesReader()
        reader.SetFileNames(reader.GetGDCMSeriesFileNames(str(directory)))
        return reader.Execute()
    finally:
        shutil.rmtree(directory, ignore_errors=True)


def _image_from_multiframe(dataset):
    """Build a positioned SimpleITK image from an enhanced multi-frame object.

    Three details decide whether this lands on the right voxels, and getting any
    of them wrong misregisters silently rather than failing:

    * **orientation** lives in the *shared* functional groups for these masks,
      not the per-frame ones. Reading only per-frame leaves SimpleITK's identity
      default, which flips x and y relative to the SER volume's (-1, -1, 1) and
      moves the mask almost entirely off the target grid;
    * **frame order** is whatever the writer chose, so frames are sorted by their
      position along the slice normal rather than trusted as stored;
    * **origin** is then the position of the first frame *after* sorting.
    """
    import SimpleITK as sitk

    array = np.asarray(dataset.pixel_array)
    shared = dataset.get("SharedFunctionalGroupsSequence")
    per_frame = dataset.get("PerFrameFunctionalGroupsSequence")

    cosines = _first_orientation(shared, per_frame)
    row, col = np.array(cosines[:3]), np.array(cosines[3:])
    normal = np.cross(row, col)

    positions = _frame_positions(per_frame, array.shape[0])
    if positions is not None:
        order = np.argsort(positions @ normal)
        array = array[order]
        origin = tuple(float(x) for x in positions[order][0])
    else:
        origin = (0.0, 0.0, 0.0)

    image = sitk.GetImageFromArray(array)
    image.SetSpacing(_multiframe_spacing(shared, positions, normal))
    image.SetOrigin(origin)
    image.SetDirection(tuple(np.column_stack([row, col, normal]).ravel()))
    return image


def _first_orientation(shared, per_frame) -> list[float]:
    for source in (shared, per_frame):
        if not source:
            continue
        orientation = source[0].get("PlaneOrientationSequence")
        if orientation:
            return [float(x) for x in orientation[0].ImageOrientationPatient]
    return [1.0, 0.0, 0.0, 0.0, 1.0, 0.0]


def _frame_positions(per_frame, n_frames: int) -> np.ndarray | None:
    if not per_frame or len(per_frame) != n_frames:
        return None
    positions = []
    for frame in per_frame:
        plane = frame.get("PlanePositionSequence")
        if not plane:
            return None
        positions.append([float(x) for x in plane[0].ImagePositionPatient])
    return np.array(positions)


def _multiframe_spacing(shared, positions, normal) -> tuple[float, float, float]:
    in_plane = (1.0, 1.0)
    thickness = None
    if shared:
        measures = shared[0].get("PixelMeasuresSequence")
        if measures:
            pixel_spacing = measures[0].get("PixelSpacing", [1.0, 1.0])
            in_plane = (float(pixel_spacing[1]), float(pixel_spacing[0]))
            thickness = measures[0].get("SpacingBetweenSlices") or measures[0].get(
                "SliceThickness"
            )
    if thickness is None and positions is not None and len(positions) > 1:
        # Derive slice spacing from the frame positions themselves; the stated
        # SliceThickness is the slab thickness and need not equal the step.
        along = np.sort(positions @ normal)
        thickness = float(np.median(np.diff(along)))
    return (in_plane[0], in_plane[1], float(thickness or 1.0))


def load_immune_scores(path: Path | str) -> dict[str, dict[str, float]]:
    """Per-patient immune-signature scores keyed by TCIA patient ID."""
    import json

    payload = json.loads(Path(path).read_text())
    return {
        f"ISPY2-{pid}": {name: values[i] for name, values in payload["scores"].items()}
        for i, pid in enumerate(payload["patients"])
    }
