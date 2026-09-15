"""
Adapter for the real inputs, for use once access is granted.

This is the switch between the demonstration and an actual reproduction.
Nothing else in the pipeline changes: M3 through M7 consume the same
`Preprocessed` object either way.

UNTESTED. The real files are behind access control and have never been run
through this code, so the column names below are inferred from the accession
metadata and the Olink Explore export format, not observed. Expect to adjust
`load_omix008341` on first contact with the real spreadsheet. It is written
to fail loudly with a readable message rather than silently mis-parse.

WHAT YOU NEED
-------------
1. TRIAL PROTEOME - CNCB OMIX008341, file OMIX008341-01, XLSX, 3.9 MB.
   "OLINK data for IPF patient plasma", Olink Explore 3072, from the
   INS018_055 (rentosertib) phase 2a trial.
   CONTROLLED ACCESS: submit a request at
   https://ngdc.cncb.ac.cn/omix/release/OMIX008341
   Note the accession is registered against the 2025 Nature Medicine trial
   paper, so confirm it carries all four visits and not just baseline.

2. REFERENCE COHORT - UK Biobank Olink (the paper uses n=55,319).
   Requires an approved UK Biobank application. Any large Olink Explore
   cohort with age and sex can substitute, at the cost of comparability.

3. CLOCK WEIGHTS - drop CSVs into a directory and pass it as `weight_dir`
   to m3_clocks.fit_all. Format: columns `protein`, `coefficient`, with an
   optional `(Intercept)` row. Name each file after the clock, e.g.
   `OrganAge_chrono.csv`. Availability per clock is recorded in m3_clocks.

IDENTIFIER NAMESPACES
---------------------
The single most likely failure is a namespace mismatch: Olink exports use
assay names (e.g. "IL6"), the organAging tables use UniProt accessions, and
UK Biobank uses its own assay identifiers. m3_clocks.load_linear_weights
raises if nothing matches, but a PARTIAL match is worse than none because it
silently scores a clock on a fraction of its proteins. Check the reported
match count against the clock's expected `n_features` every time.
"""
import os

import numpy as np
import pandas as pd

from protclock_pipeline.m2_preprocess import Preprocessed

# Columns the trial export is expected to carry, in some spelling.
SAMPLE_KEYS = ("SampleID", "sample_id", "Sample")
PATIENT_KEYS = ("SubjectID", "patient_id", "Subject", "USUBJID")
ASSAY_KEYS = ("Assay", "protein", "OlinkID", "UniProt")
NPX_KEYS = ("NPX", "npx", "value")
VISIT_KEYS = ("Visit", "week", "Timepoint", "AVISIT")
ARM_KEYS = ("Arm", "arm", "Treatment", "TRT01P")


def _pick(df, candidates, what):
    for c in candidates:
        if c in df.columns:
            return c
    raise KeyError(
        f"could not find the {what} column; tried {candidates}. "
        f"Columns present: {list(df.columns)[:25]}")


def load_omix008341(path, sheet=0):
    """
    Read the trial Olink export into (npx wide matrix, sample metadata).

    Accepts either long format (one row per sample x assay, with an NPX
    column) or wide format (one row per sample, one column per assay).
    Returns a wide matrix and a metadata frame with the columns the rest of
    the pipeline needs: sample_id, patient_id, arm, week, age, sex, plate.
    """
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"{path} not found. This file is CONTROLLED ACCESS; request it at "
            f"https://ngdc.cncb.ac.cn/omix/release/OMIX008341")

    raw = pd.read_excel(path, sheet_name=sheet)
    if isinstance(raw, dict):
        raw = next(iter(raw.values()))

    npx_col = next((c for c in NPX_KEYS if c in raw.columns), None)
    if npx_col is not None:
        # Long format: pivot to samples x assays.
        s = _pick(raw, SAMPLE_KEYS, "sample id")
        a = _pick(raw, ASSAY_KEYS, "assay")
        wide = raw.pivot_table(index=s, columns=a, values=npx_col,
                               aggfunc="mean")
        meta_cols = [c for c in raw.columns if c not in (a, npx_col)]
        meta = (raw[meta_cols].drop_duplicates(subset=[s])
                .set_index(s).loc[wide.index].reset_index())
    else:
        s = _pick(raw, SAMPLE_KEYS, "sample id")
        meta_cols = [c for c in raw.columns
                     if c in SAMPLE_KEYS + PATIENT_KEYS + VISIT_KEYS
                     + ARM_KEYS + ("age", "Age", "sex", "Sex", "plate",
                                   "PlateID")]
        wide = raw.set_index(s).drop(columns=[c for c in meta_cols
                                              if c != s], errors="ignore")
        meta = raw[meta_cols].copy()

    meta = _normalise_meta(meta)
    if len(meta) != len(wide):
        raise ValueError(f"metadata rows ({len(meta)}) do not match NPX rows "
                         f"({len(wide)})")
    return wide, meta


def _normalise_meta(meta):
    """Rename whatever the export calls things to the pipeline's names."""
    out = pd.DataFrame(index=range(len(meta)))
    ren = {"sample_id": SAMPLE_KEYS, "patient_id": PATIENT_KEYS,
           "arm": ARM_KEYS, "week": VISIT_KEYS,
           "age": ("age", "Age", "AGE"), "sex": ("sex", "Sex", "SEX"),
           "plate": ("plate", "PlateID", "Plate")}
    for target, cands in ren.items():
        col = next((c for c in cands if c in meta.columns), None)
        if col is not None:
            out[target] = meta[col].to_numpy()

    for required in ("sample_id", "patient_id", "arm", "week"):
        if required not in out.columns:
            raise KeyError(
                f"the export has no column this adapter recognises as "
                f"'{required}'. Add its spelling to the *_KEYS tuples in "
                f"real_data.py.")

    # Visits may be labelled "Week 12" or "V4" rather than a number.
    out["week"] = pd.to_numeric(
        out["week"].astype(str).str.extract(r"(\d+)")[0], errors="coerce")
    if out["week"].isna().any():
        raise ValueError("could not parse a week number from every visit "
                         "label; map them manually")

    if "sex" in out.columns and out["sex"].dtype == object:
        out["sex"] = (out["sex"].astype(str).str.upper()
                      .map({"M": 1.0, "MALE": 1.0, "F": 0.0, "FEMALE": 0.0}))
    if "plate" not in out.columns:
        out["plate"] = 1
    if "age" not in out.columns:
        raise KeyError("age is required to compute age acceleration and is "
                       "absent; it may be in a separate clinical file")
    return out


def build_preprocessed(trial_xlsx, reference_npx, reference_meta,
                       min_detection=None):
    """
    Assemble a `Preprocessed` from real trial and reference data.

    `reference_npx` is a samples x assays DataFrame whose columns are assay
    identifiers in the SAME namespace as the trial export, and
    `reference_meta` must carry age, sex and died_10y (the last only if you
    intend to fit surrogate mortality clocks).

    The ground-truth properties on `Preprocessed` do not exist for real data.
    Anything in the pipeline that compares against the planted truth (M1's
    truth table, the recovery numbers in M6, M5's power sweep) is meaningless
    here and must be skipped; the measurements themselves are not.
    """
    trial_wide, trial_meta = load_omix008341(trial_xlsx)

    shared = [c for c in trial_wide.columns if c in set(reference_npx.columns)]
    if len(shared) < 100:
        raise ValueError(
            f"only {len(shared)} assays are shared between the trial export "
            f"and the reference cohort. This is almost certainly an "
            f"identifier namespace mismatch rather than a real panel "
            f"difference; see the module docstring.")

    tr = trial_wide[shared].to_numpy(dtype=float)
    rf = reference_npx[shared].to_numpy(dtype=float)

    med = np.nanmedian(rf, axis=0)
    for mat in (rf, tr):
        miss = np.isnan(mat)
        if miss.any():
            mat[miss] = np.take(med, np.where(miss)[1])

    pre = Preprocessed(
        protein_ids=np.array(shared),
        trial_npx=tr, trial_meta=trial_meta.reset_index(drop=True),
        ref_npx=rf, ref_meta=reference_meta.reset_index(drop=True),
        panel=None, keep_idx=np.arange(len(shared)),
        qc={"n_shared": len(shared), "source": "REAL", "n_input": len(shared),
            "n_low_detection": 0, "n_panel_only": 0,
            "trial_missing_imputed": 0, "median_detection_rate": float("nan")})
    return pre


def report():
    print("REAL-DATA ADAPTER")
    print("-" * 68)
    print("  trial proteome  : CNCB OMIX008341 (CONTROLLED ACCESS)")
    print("  reference       : UK Biobank Olink (application required)")
    print("  clock weights   : see m3_clocks.CLOCKS[*].where")
    print()
    print("  This adapter has never been run against the real files, because")
    print("  they are not obtainable without approval. Column names are")
    print("  inferred. Verify the assay-identifier namespace before trusting")
    print("  any number it produces.")
    return None


if __name__ == "__main__":
    report()
