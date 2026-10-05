"""M2 -- the DC-bin artifact: why bin 0 of the protein spectrum must be dropped.

The finding
-----------
The FFT encoding mean-centres the residue-value vector before transforming, which
sets the DC term (bin 0) to **exactly zero**. The published pipeline nonetheless
keeps bin 0 as a feature, and feeds the matrix to `PLSRegression(scale=True)`,
which divides every column by its standard deviation. For bin 0 that standard
deviation is the floating-point round-off itself (~1e-16), so standardization
rescales pure numerical noise into an O(1) predictor and PLS fits through it.

Consequences, all demonstrated below on the paper's own round-1 training set:

* The cross-validated error depends on the **storage precision** of the features.
  The reference implementation casts to float32 (`np.array(..., dtype="float32")`
  in `FFT.get_aai_encoding`); at 10 components that lands on cvMSE 0.359, while
  the identical matrix in float64 gives 0.900.
* Because descriptors are ranked *by* that error, which of the 566 AAindex entries
  "wins" is decided by round-off.
* More components means more noise to fit, so the artifact pushes selection toward
  the maximum allowed `n_components` -- exactly what the published round-1 table
  shows (k=10 from 12 training points).
* Dropping bin 0 makes float32 and float64 agree to machine precision.

This is why `fftplsr.encode.spectrum` defaults to ``drop_dc=True``.

Run: ``python3 -m fftplsr.m2_dc_artifact``
"""

from __future__ import annotations

import numpy as np
from sklearn.model_selection import LeaveOneOut

from . import datasets
from .encode import MutationEncoder
from .model import cv_score, pls_path_predict, screen_indices
from .variants import mutated_positions

COMPONENTS = (2, 5, 8, 10)


def _round1_matrix(drop_dc: bool):
    parent = datasets.ifrs()
    table, _ = datasets.trainset(1)
    labels = table["Variants"].tolist()
    y = table["Fitness"].to_numpy(dtype=float)
    encoder = MutationEncoder(
        parent, positions=mutated_positions(labels + datasets.ROUND12_SITES)
    )
    return encoder, labels, y, encoder.encode(labels, ["OOBM850103"], drop_dc=drop_dc)


def _loocv_curve(X, y, max_k=10):
    oof = np.empty((max_k, len(y)))
    for train, test in LeaveOneOut().split(X):
        oof[:, test] = pls_path_predict(X[train], y[train], X[test], max_k)
    return ((y[None, :] - oof) ** 2).mean(axis=1)


def main() -> None:
    print("=" * 78)
    print("M2  The DC-bin artifact in the published FFT-PLSR encoding")
    print("=" * 78)

    _, _, y, X_dc = _round1_matrix(drop_dc=False)
    _, _, _, X_no = _round1_matrix(drop_dc=True)
    X_dc32 = X_dc.astype(np.float32).astype(np.float64)

    col0 = X_dc[:, 0]
    print(f"\nBin 0 over the 13 round-1 variants ({X_dc.shape[1]} features total):")
    print(f"  min={col0.min():.3e}  max={col0.max():.3e}  std={col0.std(ddof=1):.3e}")
    print("  -> identically zero up to round-off, as mean-centering guarantees.")
    print(
        f"  after standardization its values span "
        f"{np.ptp((col0 - col0.mean()) / col0.std(ddof=1)):.2f} standard deviations "
        "of pure noise."
    )

    print("\nLeave-one-out MSE on the round-1 training set (descriptor OOBM850103):")
    header = "  {:<30}" + "".join(f"  k={k:<8}" for k in COMPONENTS)
    print(header.format("configuration"))
    for tag, X in [
        ("bin 0 kept, numpy rfft f64", X_dc),
        ("bin 0 kept, numpy rfft f32", X_dc32),
        ("bin 0 dropped, f64", X_no),
        ("bin 0 dropped, f32", X_no.astype(np.float32).astype(np.float64)),
    ]:
        curve = _loocv_curve(X, y)
        print("  {:<30}".format(tag) + "".join(f"  {curve[k - 1]:<10.4f}" for k in COMPONENTS))

    # Any FFT routine is entitled to a different round-off in the DC bin: summation
    # order, precision and the choice of full-complex vs real transform all change
    # it. The published code used scipy.fftpack on a float32 array, whose bin 0
    # spans 0 .. 5.9e-16 rather than numpy rfft's 1.68e-16 .. 1.77e-16. So the
    # honest way to read bin 0 is as a lottery: redraw it at the round-off scale
    # and see how far the "cross-validated error" moves.
    print(
        "\nRe-drawing bin 0 at the round-off scale (uniform 0..6e-16, 200 draws)\n"
        "-- every draw is an equally legitimate FFT implementation:"
    )
    rng = np.random.default_rng(0)
    rows = []
    for _ in range(200):
        X = X_dc.copy()
        X[:, 0] = rng.uniform(0.0, 6e-16, size=len(y))
        curve = _loocv_curve(X, y)
        rows.append((curve.min(), int(curve.argmin()) + 1, curve[9]))
    best = np.array([r[0] for r in rows])
    ks = [r[1] for r in rows]
    k10 = np.array([r[2] for r in rows])
    print(f"  best cvMSE over k:  min={best.min():.4f}  median={np.median(best):.4f}  max={best.max():.4f}")
    print(f"  cvMSE at k=10:      min={k10.min():.4f}  median={np.median(k10):.4f}  max={k10.max():.4f}")
    print(f"  selected k:         {sorted(set(ks))}  (mode {max(set(ks), key=ks.count)})")
    print(
        f"  the published round-1 table reports cvMSE 0.3594 at k=10 -- inside this\n"
        f"  round-off lottery ({(k10 <= 0.3594).mean() * 100:.0f}% of draws reach it), "
        "not a property of the data."
    )

    rows_no = []
    for _ in range(50):
        X = X_no.copy()
        X[:, 0] += rng.uniform(-6e-16, 6e-16, size=len(y))
        rows_no.append(_loocv_curve(X, y))
    spread = np.ptp(np.array(rows_no), axis=0).max()
    print(f"\n  Same perturbation with bin 0 dropped: cvMSE moves by at most {spread:.2e}.")

    print("\nEffect on descriptor selection (full 566-entry screen, round-1 data):")
    encoder, labels, y, _ = _round1_matrix(drop_dc=False)
    for drop in (False, True):
        sel = screen_indices(
            lambda codes, d=drop: encoder.encode(labels, codes, drop_dc=d),
            y,
            n_rounds=1,
            cv=None,
            verbose=False,
        )
        chose = cv_score(encoder.encode(labels, sel.indices, drop_dc=drop), y, cv=None)
        print(
            f"  drop_dc={str(drop):<5} -> winner {sel.label:<12} "
            f"cvMSE={sel.cv_mse:.4f}  k={chose.n_components}"
        )

    print(
        "\nConclusion: bin 0 carries no information by construction, but standardizing it\n"
        "turns round-off into a fitted predictor. `fftplsr` drops it by default."
    )


if __name__ == "__main__":
    main()
