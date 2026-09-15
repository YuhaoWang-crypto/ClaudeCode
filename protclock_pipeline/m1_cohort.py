"""
Module 1 - Synthetic Olink-shaped cohorts with a PLANTED ground truth.

Reproduction target
-------------------
Zhavoronkov et al., "Integration of proteomic aging clocks in a phase 2a
clinical trial supports simultaneous geroprotective assessment",
Nat Biotechnol (2026), doi:10.1038/s41587-026-03286-y.

The real inputs are NOT openly available:
  * trial serum proteome  -> CNCB OMIX008341, CONTROLLED access
  * reference cohort      -> UK Biobank (n=55,319), application required

So this module manufactures data with the SAME SHAPE as the real inputs
(Olink Explore 3072 NPX, 4 arms x 4 visits, 42 patients) and with a known
effect written into it on purpose. Everything downstream is then judged by
whether it RECOVERS the planted truth, which tests the code rather than the
biology.

WARNING - nothing in this module is a biological finding. The aging
signal, the drug effect and the dose ordering are all injected here by
construction. Real numbers require the controlled data; see data/README.md.

Generative model (NPX, per protein p, per sample i)
--------------------------------------------------
  NPX[i,p] = mu_p
           + beta_p * (age_i - AGE_CENTER)      # chronological aging
           + gamma_p * sex_i                    # sex effect
           + delta_p * ipf_i                    # disease (fibrosis) shift
           + drug_shift(p, arm_i, week_i)       # treatment, see below
           + u_i                                # per-sample/plate offset
           + eps                                # assay noise

  drug_shift = -beta_p * EFFECT_YEARS * dose_w[arm] * ramp[week]
               for p in the drug-responsive set

The drug term is expressed in "years of proteomic age equivalent": shifting a
protein by -beta_p * Y moves it exactly as far as being Y years younger would.
Because only a minority of the responsive set is aging-associated, the age read
out by a clock is diluted well below EFFECT_YEARS, as in the real trial.
"""
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

# ----------------------------------------------------------------- constants
# Olink Explore 3072 assay count reported in the paper.
N_PROTEINS = 2841
# Proteins shared between the trial panel and the UK Biobank reference.
N_SHARED = 2832

AGE_CENTER = 60.0

# Subject-level biological-age offset, in years. This is what separates
# biological from chronological age: a subject whose offset is +8 carries the
# proteome of someone eight years older. Without it every clock recovers
# chronological age almost exactly, age acceleration has no variance, and the
# trial statistics become meaningless. Mortality risk keys off age + offset,
# which is why the mortality-trained clocks behave differently from the
# chronological ones.
BIO_SD = 6.0
# IPF patients start biologically older than the population reference.
TRIAL_BIO_MEAN = 3.0

# Trial design: 42 IPF patients, 4 arms, 4 visits (baseline, wk2, wk4, wk12).
ARMS = ("placebo", "rento_30mg_QD", "rento_30mg_BID", "rento_60mg_QD")
ARM_N = {"placebo": 11, "rento_30mg_QD": 10, "rento_30mg_BID": 11,
         "rento_60mg_QD": 10}
WEEKS = (0, 2, 4, 12)

# Planted dose-response. The paper reports 30 mg BID as the arm with the most
# consistent signal across clocks, so BID is given the largest weight here.
DOSE_WEIGHT = {"placebo": 0.0, "rento_30mg_QD": 0.60,
               "rento_30mg_BID": 1.00, "rento_60mg_QD": 0.80}
# Exposure ramp; the trial's primary proteomic contrast is baseline -> wk12.
TIME_RAMP = {0: 0.0, 2: 0.25, 4: 0.50, 12: 1.0}

# Planted set sizes. 326 = number of circulating proteins the paper reports as
# altered by rentosertib. The aging fraction and the overlap are chosen so the
# TRUE odds ratio of aging-enrichment lands near the paper's 1.74.
N_RESPONSIVE = 326
AGING_FRACTION = 0.15
N_RESPONSIVE_AGING = 72

# Size of the drug effect in proteomic-age-equivalent years, before dilution.
EFFECT_YEARS = 30.0

# The paper's reference cohort is n=55,319. That is 1.2 GB of dense NPX, so the
# default here is smaller; raise it if you have the memory and the time.
N_REFERENCE_DEFAULT = 6000


# ------------------------------------------------------------------- panel
@dataclass
class ProteinPanel:
    """The assay panel plus the ground truth written into it."""
    ids: np.ndarray                 # protein identifiers, e.g. "OLK0001_PROT1"
    mu: np.ndarray                  # baseline NPX mean
    beta: np.ndarray                # NPX change per year of age
    gamma: np.ndarray               # sex effect (male = 1)
    delta: np.ndarray               # IPF disease shift
    sigma: np.ndarray               # assay noise SD
    is_aging: np.ndarray            # bool, truly age-associated
    is_responsive: np.ndarray       # bool, truly drug-responsive
    drug_delta: np.ndarray = None   # NPX shift at full dose and full exposure
    gene_sets: dict = field(default_factory=dict)

    def __len__(self):
        return len(self.ids)

    def truth_table(self):
        """2x2 contingency of the PLANTED aging x responsive overlap."""
        a = int(np.sum(self.is_aging & self.is_responsive))
        b = int(np.sum(~self.is_aging & self.is_responsive))
        c = int(np.sum(self.is_aging & ~self.is_responsive))
        d = int(np.sum(~self.is_aging & ~self.is_responsive))
        odds_ratio = (a / b) / (c / d)
        return {"responsive_aging": a, "responsive_other": b,
                "other_aging": c, "other_other": d,
                "odds_ratio": odds_ratio}


def real_panel_symbols():
    """
    The real Olink gene symbols, when the paper's library is installed.

    Returns None if proteoclock is absent, in which case the panel falls back
    to invented identifiers and only surrogate clocks can score it. A
    published clock is keyed on gene symbols, so it cannot read a panel whose
    proteins are called "OLK0001_PROT1".
    """
    try:
        from protclock_pipeline import proteoclock_backend
        if proteoclock_backend.available():
            return proteoclock_backend.panel_symbols()
    except Exception:
        pass
    return None


def _published_aging_coefficients(ids):
    """
    Age-slope direction per protein, from a published clock when available.

    Returns an array aligned with `ids`, zero where the clock has no
    coefficient, or None when the paper's library is not installed.
    """
    try:
        from protclock_pipeline import proteoclock_backend
        coefs = proteoclock_backend.reference_aging_coefficients()
    except Exception:
        coefs = None
    if not coefs:
        return None
    out = np.array([float(coefs.get(str(i), 0.0)) for i in ids])
    return out if np.any(out != 0.0) else None


def build_panel(seed=20260907, symbols=None):
    """Construct the protein panel and decide which proteins carry signal."""
    rng = np.random.default_rng(seed)

    if symbols is None:
        symbols = real_panel_symbols()
    if symbols is not None:
        ids = np.asarray(symbols)
        n = len(ids)
    else:
        n = N_PROTEINS
        ids = np.array([f"OLK{i:04d}_PROT{i}" for i in range(n)])
    mu = rng.normal(5.0, 1.5, n)
    sigma = rng.uniform(0.25, 0.60, n)
    gamma = rng.normal(0.0, 0.15, n)

    # Aging-associated proteins get a real age slope; the rest get a slope
    # small enough to be indistinguishable from noise in a finite cohort.
    n_aging = int(round(AGING_FRACTION * n))
    beta = rng.normal(0.0, 0.0008, n)                      # background drift

    coefs = _published_aging_coefficients(ids)
    if coefs is not None:
        # Ground the aging axis in a published clock's coefficients, so the
        # REAL clocks can read this cohort. See the backend's
        # reference_aging_coefficients docstring on why this is circular for
        # those clocks and what it is and is not evidence of.
        strength = np.abs(coefs)
        aging_idx = np.argsort(strength)[::-1][:n_aging]
        scale = np.median(np.abs(coefs[aging_idx]))
        beta[aging_idx] = (coefs[aging_idx] / max(scale, 1e-12)) * 0.015
    else:
        aging_idx = rng.choice(n, size=n_aging, replace=False)
        signs = rng.choice([-1.0, 1.0], size=n_aging)
        beta[aging_idx] = signs * rng.uniform(0.008, 0.030, n_aging)

    is_aging = np.zeros(n, dtype=bool)
    is_aging[aging_idx] = True

    # Drug-responsive set: N_RESPONSIVE_AGING drawn from the aging proteins,
    # the remainder from the rest. This fixes the true enrichment odds ratio.
    resp_aging = rng.choice(aging_idx, size=N_RESPONSIVE_AGING, replace=False)
    other_idx = np.setdiff1d(np.arange(n), aging_idx)
    resp_other = rng.choice(other_idx, size=N_RESPONSIVE - N_RESPONSIVE_AGING,
                            replace=False)
    is_responsive = np.zeros(n, dtype=bool)
    is_responsive[resp_aging] = True
    is_responsive[resp_other] = True

    # Drug effect per protein, in NPX at full dose and full exposure.
    #
    # Aging-associated responders move exactly against their own age slope, so
    # shifting them by -beta_p * Y is worth Y years to any clock reading them.
    # Responders that are NOT aging-associated move by a comparable magnitude
    # in a random direction: they are pharmacologically real but carry no age
    # information. Making their shift proportional to beta_p instead would
    # leave them near zero, so differential-abundance testing would recover
    # only the aging responders and the enrichment odds ratio would be
    # unbounded rather than the modest number the paper reports.
    drug_delta = np.zeros(n)
    resp_aging_mask = is_aging & is_responsive
    drug_delta[resp_aging_mask] = -beta[resp_aging_mask] * EFFECT_YEARS
    typical = float(np.mean(np.abs(drug_delta[resp_aging_mask])))
    resp_other_mask = (~is_aging) & is_responsive
    drug_delta[resp_other_mask] = rng.normal(0.0, typical,
                                             int(resp_other_mask.sum()))

    # IPF shifts a fibrosis programme at baseline, overlapping the drug target.
    n_fibrosis = 240
    fibrosis_idx = np.concatenate([
        rng.choice(resp_aging, size=60, replace=False),
        rng.choice(resp_other, size=90, replace=False),
        rng.choice(other_idx, size=n_fibrosis - 150, replace=False)])
    delta = np.zeros(n)
    delta[fibrosis_idx] = rng.normal(0.45, 0.20, len(fibrosis_idx))

    # Pathway sets for the enrichment module. Senescence deliberately overlaps
    # the aging set, metabolic partially, immune barely.
    gene_sets = {
        "senescence_SASP": rng.choice(aging_idx, size=180, replace=False),
        "metabolic_process": np.concatenate([
            rng.choice(aging_idx, size=70, replace=False),
            rng.choice(other_idx, size=150, replace=False)]),
        "fibrosis_ECM": fibrosis_idx,
        "immune_activation": rng.choice(other_idx, size=200, replace=False),
    }

    return ProteinPanel(ids=ids, mu=mu, beta=beta, gamma=gamma, delta=delta,
                        sigma=sigma, is_aging=is_aging,
                        is_responsive=is_responsive, drug_delta=drug_delta,
                        gene_sets=gene_sets)


# -------------------------------------------------------- reference cohort
def simulate_reference(panel, n=N_REFERENCE_DEFAULT, seed=11):
    """
    Population reference standing in for UK Biobank Olink.

    Returns (npx, meta). npx is n x N_PROTEINS. meta carries age, sex and a
    10-year mortality outcome, which is what the mortality-trained clocks
    (OrganAge_mortality, PAC) need as a fitting target.
    """
    rng = np.random.default_rng(seed)

    age = np.clip(rng.normal(57.0, 8.0, n), 39, 71)
    sex = rng.integers(0, 2, n).astype(float)
    bio = rng.normal(0.0, BIO_SD, n)          # biological-age offset, years
    eff_age = age + bio                       # what the proteome reflects

    npx = (panel.mu[None, :]
           + np.outer(eff_age - AGE_CENTER, panel.beta)
           + np.outer(sex, panel.gamma))
    npx += rng.normal(0.0, 1.0, (n, 1)) * 0.25            # per-subject offset
    npx += rng.normal(0.0, 1.0, (n, len(panel))) * panel.sigma[None, :]

    # 10-year mortality is driven by BIOLOGICAL age, so a mortality-trained
    # clock has a different target from a chronological-age clock even though
    # both read the same proteins.
    # Intercept set for roughly 15-20% 10-year mortality, which is what the
    # mortality-trained clocks need in order to have enough events to fit.
    logit = -1.4 + 0.085 * (eff_age - AGE_CENTER)
    died = rng.random(n) < 1.0 / (1.0 + np.exp(-logit))

    meta = pd.DataFrame({"subject_id": [f"REF{i:06d}" for i in range(n)],
                         "age": age, "sex": sex,
                         "bio_offset_true": bio,
                         "died_10y": died.astype(int)})
    return npx, meta


# ------------------------------------------------------------------- trial
def simulate_trial(panel, seed=7, effect_scale=1.0):
    """
    Phase 2a trial: 42 IPF patients x 4 visits, Olink NPX.

    effect_scale=1.0 plants the full drug effect. effect_scale=0.0 produces a
    NULL trial with no treatment effect at all, which m5/m6 use as a negative
    control: a pipeline that reports significance there is broken.
    """
    rng = np.random.default_rng(seed)

    rows, arms_col, pids = [], [], []
    for arm in ARMS:
        for k in range(ARM_N[arm]):
            pids.append(f"{arm}_{k:02d}")
            arms_col.append(arm)
    n_pat = len(pids)

    pat_age = np.clip(rng.normal(66.0, 7.0, n_pat), 45, 82)
    pat_sex = (rng.random(n_pat) < 0.72).astype(float)      # IPF is male-skewed
    pat_offset = rng.normal(0.0, 0.30, n_pat)               # patient random eff.
    # Biological-age offset is a property of the patient, so it is drawn once
    # and held fixed across all four visits. Only the drug term moves with time.
    pat_bio = rng.normal(TRIAL_BIO_MEAN, BIO_SD, n_pat)

    meta_rows = []
    for i, pid in enumerate(pids):
        for wk in WEEKS:
            base = (panel.mu
                    + panel.beta * (pat_age[i] + pat_bio[i] - AGE_CENTER)
                    + panel.gamma * pat_sex[i]
                    + panel.delta)                          # every patient IPF+

            w = DOSE_WEIGHT[arms_col[i]] * TIME_RAMP[wk] * effect_scale
            shift = panel.drug_delta * w

            vals = (base + shift + pat_offset[i]
                    + rng.normal(0.0, 1.0, len(panel)) * panel.sigma)
            rows.append(vals)
            meta_rows.append({"sample_id": f"{pid}_wk{wk}", "patient_id": pid,
                              "arm": arms_col[i], "week": wk,
                              "age": pat_age[i], "sex": pat_sex[i],
                              "bio_offset_true": pat_bio[i],
                              "planted_effect_years": EFFECT_YEARS * w,
                              "plate": (len(meta_rows) // 24) + 1})

    npx = np.vstack(rows)
    meta = pd.DataFrame(meta_rows)

    # Plate/batch offsets, so that m2's normalisation has something to remove.
    for plate, grp in meta.groupby("plate"):
        npx[grp.index.to_numpy(), :] += rng.normal(0.0, 0.12)

    return npx, meta


def report():
    panel = build_panel()
    truth = panel.truth_table()
    ref_npx, ref_meta = simulate_reference(panel, n=1200)
    tr_npx, tr_meta = simulate_trial(panel)

    print("M1  SYNTHETIC OLINK-SHAPED COHORTS  (planted ground truth)")
    print("-" * 68)
    print(f"  panel proteins              : {len(panel)}  (paper: {N_PROTEINS})")
    print(f"  truly aging-associated      : {int(panel.is_aging.sum())}")
    print(f"  truly drug-responsive       : {int(panel.is_responsive.sum())}"
          f"  (paper: {N_RESPONSIVE})")
    print(f"  overlap (responsive ^ aging): {truth['responsive_aging']}")
    print(f"  TRUE enrichment odds ratio  : {truth['odds_ratio']:.3f}"
          f"  (paper: 1.74)")
    print(f"  reference cohort            : {ref_npx.shape}  "
          f"deaths={int(ref_meta.died_10y.sum())}")
    print(f"  trial matrix                : {tr_npx.shape}  "
          f"patients={tr_meta.patient_id.nunique()} "
          f"visits={sorted(tr_meta.week.unique())}")
    print("  WARNING: every number above is injected by construction.")
    return {"panel": panel, "truth": truth}


if __name__ == "__main__":
    report()
