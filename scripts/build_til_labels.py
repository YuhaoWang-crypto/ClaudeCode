#!/usr/bin/env python
"""Build the real TIL label set for the 136 paired TCGA-BRCA patients.

This is the label half of the external test cohort, computed from real data
rather than proposed.  For each patient it downloads the Stony Brook fractional
TIL map for their diagnostic slide and reduces it to per-patient measures:

* ``til_fraction_tissue`` -- infiltrated fraction of *tissue* patches;
* ``til_moran_i``         -- spatial autocorrelation of the TIL field, i.e. how
  clustered the infiltrate is rather than how much of it there is;
* ``til_cluster_share``   -- share of infiltrated patches sitting in the single
  largest connected cluster.

The last two matter because the whole question is spatial: two slides with
identical TIL fractions can be diffusely infiltrated or have one dense focus,
and only the second pair of measures can tell them apart.

What this still cannot produce is the desert / excluded / inflamed trichotomy.
That needs a tumour region to separate core from invasive margin, and the
published maps carry no tumour boundary. See ``docs/TEST_SET.md``.

Usage::

    python scripts/build_til_labels.py                  # all 136
    python scripts/build_til_labels.py --limit 10       # quick check
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

from mri_cd8_til.clients.idc import FRACTIONAL_TIL, IDCClient, save_map  # noqa: E402

INVENTORY = "results/cohort_inventory.json"


#: Fractional maps encode the TIL probability as uint8 0-255, with 0 doubling as
#: "no tissue".  128 is the 0.5 probability cut that reproduces the companion
#: binary map -- verified per patient in ``--verify`` mode.
TIL_PROB_CUT = 128


def morans_i(field: np.ndarray, mask: np.ndarray) -> float:
    """Moran's I of the TIL field over *tissue patches only*.

    Restricting to ``mask`` is not a refinement, it is the difference between
    measuring two different things.  Run over the full grid, the statistic is
    dominated by the background zeros and mostly reports that tissue is one
    contiguous blob -- a number near 0.6 for every slide, regardless of how the
    lymphocytes are arranged.

    Pairs where either patch is outside the mask are dropped, so the result is
    the spatial autocorrelation of infiltration within the tissue.
    """
    x = np.where(mask, field.astype(float), np.nan)
    mean = np.nanmean(x)
    deviation = np.where(mask, x - mean, 0.0)
    denominator = (deviation[mask] ** 2).sum()
    if denominator <= 0 or mask.sum() < 4:
        return 0.0

    numerator = 0.0
    n_pairs = 0
    for shift, axis in ((1, 0), (1, 1)):
        head = slice(1, None)
        tail = slice(None, -1)
        if axis == 0:
            a, b = deviation[head, :], deviation[tail, :]
            valid = mask[head, :] & mask[tail, :]
        else:
            a, b = deviation[:, head], deviation[:, tail]
            valid = mask[:, head] & mask[:, tail]
        numerator += (a[valid] * b[valid]).sum()
        n_pairs += int(valid.sum())

    if n_pairs == 0:
        return 0.0
    # Each unordered pair is counted once above, so the weight total is 2*n_pairs.
    return float((int(mask.sum()) / (2 * n_pairs)) * (2 * numerator / denominator))


def largest_cluster_share(binary: np.ndarray) -> float:
    """Share of positive patches in the largest 8-connected component."""
    positives = int(binary.sum())
    if positives == 0:
        return 0.0
    try:
        from scipy.ndimage import label

        labelled, n = label(binary, structure=np.ones((3, 3)))
        if n == 0:
            return 0.0
        sizes = np.bincount(labelled.ravel())[1:]
        return float(sizes.max() / positives)
    except ImportError:  # pragma: no cover
        return float("nan")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--cache", default="results/til_maps")
    parser.add_argument("--out", default="results/til_labels.json")
    parser.add_argument("--inventory", default=INVENTORY)
    args = parser.parse_args()

    inventory = json.loads(Path(args.inventory).read_text())
    paired = set(inventory["tcga_brca_linkage"]["paired_ids"])
    print(f"paired TCGA-BRCA patients (MRI + RNA-seq + diagnostic WSI): {len(paired)}")

    client = IDCClient()
    series = client.til_series(paired, FRACTIONAL_TIL)
    print(f"fractional TIL-map series found: {len(series)} for {series.PatientID.nunique()} patients")

    rows = series.drop_duplicates(subset="PatientID")
    if args.limit:
        rows = rows.head(args.limit)

    records = []
    for i, (_, row) in enumerate(rows.iterrows(), start=1):
        try:
            til_map = client.fetch_map(row.crdc_series_uuid, row.PatientID, FRACTIONAL_TIL)
        except Exception as exc:  # noqa: BLE001
            print(f"  [{i}/{len(rows)}] {row.PatientID}: FAILED {type(exc).__name__} {exc}")
            continue
        save_map(til_map, args.cache)

        grid = til_map.grid
        # Value 0 doubles as "no tissue" in these maps, so the tissue mask and
        # the TIL call come from different thresholds on the same array.
        tissue = grid > 0
        positive = tissue & (grid >= TIL_PROB_CUT)
        n_tissue = int(tissue.sum())
        if n_tissue < 1000:
            print(f"  [{i}/{len(rows)}] {til_map.patient_id}: SKIPPED, only {n_tissue} tissue patches")
            continue

        record = {
            "patient_id": til_map.patient_id,
            "grid_shape": list(grid.shape),
            "n_patches": til_map.n_patches,
            "n_tissue_patches": n_tissue,
            "tissue_fraction": float(tissue.mean()),
            "til_fraction_tissue": float(positive.sum() / n_tissue),
            "til_mean_prob": float(grid[tissue].mean() / 255.0),
            "til_moran_i": morans_i(grid, tissue),
            "til_cluster_share": largest_cluster_share(positive),
        }
        records.append(record)
        if i % 20 == 0 or i == len(rows):
            print(
                f"  [{i}/{len(rows)}] {til_map.patient_id}  "
                f"tissue={n_tissue:6d}  TIL={record['til_fraction_tissue']:.4f}  "
                f"Moran I={record['til_moran_i']:+.3f}"
            )

    if not records:
        print("no TIL maps retrieved")
        return 1

    fractions = np.array([r["til_fraction_tissue"] for r in records])
    moran = np.array([r["til_moran_i"] for r in records])
    cluster = np.array([r["til_cluster_share"] for r in records])

    print(f"\nreal TIL labels for {len(records)} patients")
    print(f"  TIL fraction over tissue  median {np.median(fractions):.4f}  "
          f"IQR [{np.quantile(fractions, .25):.4f}, {np.quantile(fractions, .75):.4f}]  "
          f"range [{fractions.min():.4f}, {fractions.max():.4f}]")
    print(f"  Moran's I      median {np.median(moran):.3f}  "
          f"IQR [{np.quantile(moran, .25):.3f}, {np.quantile(moran, .75):.3f}]")
    print(f"  cluster share  median {np.median(cluster):.3f}")
    print(
        f"\n  {int((fractions >= np.quantile(fractions, 2 / 3)).sum())} patients above the "
        f"upper-tertile cut ({np.quantile(fractions, 2 / 3):.4f}) -- the 'inflamed' group "
        "under a tertile split"
    )

    payload = {
        "n_patients": len(records),
        "map_variant": FRACTIONAL_TIL,
        "records": records,
        "summary": {
            "til_fraction_tissue_median": float(np.median(fractions)),
            "til_fraction_tissue_iqr": [float(np.quantile(fractions, 0.25)), float(np.quantile(fractions, 0.75))],
            "moran_i_median": float(np.median(moran)),
            "cluster_share_median": float(np.median(cluster)),
        },
        "missing_for_supervised_test": [
            "tumour region maps: the published TIL maps carry no tumour boundary, so "
            "core vs invasive margin -- and therefore desert vs excluded -- cannot be derived",
            "expert tumour VOIs on the TCGA-BRCA MRI: no public segmentations exist for "
            "this collection (MAMA-MIA covers ISPY1/ISPY2/NACT/DUKE, not TCGA-BRCA)",
        ],
    }
    path = Path(args.out)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2))
    print(f"\nwrote {path} and {len(records)} cached maps under {args.cache}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
