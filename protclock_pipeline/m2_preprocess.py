"""
Module 2 - Olink NPX preprocessing and panel harmonisation.

The steps a real run of this pipeline has to perform before any clock can be
scored, applied identically to the synthetic matrices from M1 and to the real
OMIX008341 export once access is granted:

  1. LOD censoring         - Olink reports values below the limit of detection;
                             they are flagged, not silently kept.
  2. Assay QC              - drop proteins whose detectability is too low to
                             carry signal.
  3. Missing-value imputation - per-protein median over the reference cohort.
  4. Plate/batch centring  - remove the plate offsets that Olink bridging
                             normalisation would otherwise leave behind.
  5. Panel harmonisation   - intersect trial and reference assay lists. The
                             paper scores clocks on the 2,832 proteins present
                             in both the trial panel and UK Biobank.

Every clock in M3 consumes the output of this module, so the trial and the
reference cohort are guaranteed to be on the same protein axis and the same
NPX scale.
"""
import numpy as np
import pandas as pd

from protclock_pipeline import m1_cohort

# Fraction of assays dropped to emulate the trial panel / UKB panel mismatch.
N_SHARED = m1_cohort.N_SHARED

# An assay must be above LOD in at least this fraction of samples to be kept.
MIN_DETECTION_RATE = 0.50


class Preprocessed:
    """Trial and reference NPX on a common, QC'd protein axis."""

    def __init__(self, protein_ids, trial_npx, trial_meta, ref_npx, ref_meta,
                 panel, keep_idx, qc):
        self.protein_ids = protein_ids
        self.trial_npx = trial_npx
        self.trial_meta = trial_meta
        self.ref_npx = ref_npx
        self.ref_meta = ref_meta
        self.panel = panel
        self.keep_idx = keep_idx          # index back into the original panel
        self.qc = qc

    # Ground-truth masks, re-indexed onto the surviving proteins. Used only to
    # score the pipeline against the planted truth, never as a model input.
    @property
    def is_aging_true(self):
        return self.panel.is_aging[self.keep_idx]

    @property
    def is_responsive_true(self):
        return self.panel.is_responsive[self.keep_idx]

    def gene_sets_reindexed(self):
        """Map the panel's pathway sets onto surviving-protein positions."""
        pos = {orig: i for i, orig in enumerate(self.keep_idx)}
        out = {}
        for name, idx in self.panel.gene_sets.items():
            out[name] = np.array(sorted({pos[i] for i in np.unique(idx)
                                         if i in pos}))
        return out


def _apply_lod(npx, lod, rng):
    """Censor values below LOD to NaN and report the detection rate."""
    below = npx < lod[None, :]
    out = npx.copy()
    out[below] = np.nan
    detection_rate = 1.0 - below.mean(axis=0)
    return out, detection_rate


def run(panel=None, n_reference=m1_cohort.N_REFERENCE_DEFAULT,
        effect_scale=1.0, seed=3, verbose=False):
    """Simulate, censor, QC, impute, batch-centre and harmonise."""
    rng = np.random.default_rng(seed)
    if panel is None:
        panel = m1_cohort.build_panel()

    # The seed must reach the cohort simulators, not just the QC rng, or every
    # "replicate" redraws the same patients and the same noise. That makes
    # repeated runs look far more stable than they are and puts a fixed
    # patient-level fluke into every replicate of a power sweep.
    ref_npx, ref_meta = m1_cohort.simulate_reference(panel, n=n_reference,
                                                     seed=seed * 7919 + 11)
    tr_npx, tr_meta = m1_cohort.simulate_trial(panel, seed=seed * 104729 + 7,
                                               effect_scale=effect_scale)

    # --- 1. LOD censoring -------------------------------------------------
    # Most assays sit well above their LOD; a tail of ~200 poorer assays is
    # detected in only 55-99% of samples, which is what makes step 2 and the
    # imputation in step 3 do real work.
    lod = panel.mu - rng.uniform(1.0, 3.0, len(panel))
    poor = rng.choice(len(panel), size=200, replace=False)
    lod[poor] = panel.mu[poor] - rng.uniform(0.15, 0.90, len(poor))

    ref_npx, ref_det = _apply_lod(ref_npx, lod, rng)
    tr_npx, tr_det = _apply_lod(tr_npx, lod, rng)

    # --- 2. assay QC ------------------------------------------------------
    detectable = (ref_det >= MIN_DETECTION_RATE) & (tr_det >= MIN_DETECTION_RATE)

    # --- 5a. panel harmonisation -----------------------------------------
    # Assays present in the trial panel but absent from the reference. The
    # paper scores 2,832 of its 2,841 assays, so 9 are dropped here; when
    # QC has already removed assays, the shortfall is taken from the least
    # detectable survivors so the analysis panel is always the best N_SHARED.
    n_detectable = int(detectable.sum())
    n_drop = max(n_detectable - N_SHARED, 0)
    if n_drop:
        order = np.argsort(tr_det[detectable])          # worst detection first
        panel_only = np.where(detectable)[0][order[:n_drop]]
        shared = detectable.copy()
        shared[panel_only] = False
    else:
        shared = detectable.copy()

    keep_idx = np.where(shared)[0]
    ref_npx = ref_npx[:, keep_idx]
    tr_npx = tr_npx[:, keep_idx]
    protein_ids = panel.ids[keep_idx]

    # --- 3. imputation ----------------------------------------------------
    # Medians come from the reference cohort only, so the trial never informs
    # its own imputation. This is the same discipline a real run needs.
    med = np.nanmedian(ref_npx, axis=0)
    n_missing_trial = int(np.isnan(tr_npx).sum())
    for mat in (ref_npx, tr_npx):
        miss = np.isnan(mat)
        mat[miss] = np.take(med, np.where(miss)[1])

    # --- 4. plate centring ------------------------------------------------
    # Centre each plate on the reference-cohort mean profile.
    ref_mean = ref_npx.mean(axis=0)
    for plate, grp in tr_meta.groupby("plate"):
        rows = grp.index.to_numpy()
        offset = np.median(tr_npx[rows, :] - ref_mean[None, :])
        tr_npx[rows, :] -= offset

    qc = {
        "n_input": len(panel),
        "n_low_detection": int((~detectable).sum()),
        "n_panel_only": int(n_drop),
        "n_shared": len(keep_idx),
        "trial_missing_imputed": n_missing_trial,
        "median_detection_rate": float(np.median(tr_det[keep_idx])),
    }

    if verbose:
        print(qc)

    return Preprocessed(protein_ids, tr_npx, tr_meta, ref_npx, ref_meta,
                        panel, keep_idx, qc)


def report():
    pre = run(n_reference=1500)
    q = pre.qc
    print("M2  NPX PREPROCESSING + PANEL HARMONISATION")
    print("-" * 68)
    print(f"  assays in                   : {q['n_input']}")
    print(f"  dropped, detection < {MIN_DETECTION_RATE:.0%}    : {q['n_low_detection']}")
    print(f"  dropped, trial-panel only   : {q['n_panel_only']}")
    print(f"  shared analysis panel       : {q['n_shared']}  "
          f"(paper: {N_SHARED})")
    print(f"  median detection rate       : {q['median_detection_rate']:.3f}")
    print(f"  trial NPX imputed (<LOD)    : {q['trial_missing_imputed']}")
    print(f"  trial matrix                : {pre.trial_npx.shape}")
    print(f"  reference matrix            : {pre.ref_npx.shape}")
    assert not np.isnan(pre.trial_npx).any(), "trial NPX still has NaN"
    assert not np.isnan(pre.ref_npx).any(), "reference NPX still has NaN"
    print("  no NaN remains in either matrix.")
    return pre


if __name__ == "__main__":
    report()
