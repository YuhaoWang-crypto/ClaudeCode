"""
Module 3 - The six proteomic aging clocks.

The paper scores six independently published clocks on the same samples:

  ProtAge             Argentieri 2024   deep/stacked ML, chronological age
  OrganAge_chrono     Goeminne  2025    linear, chronological age
  OrganAge_mortality  Goeminne  2025    linear, log mortality hazard
  PAC                 Kuo       2024    elastic net, mortality-informed age
  PAOPAC              Han       2026    organ-specific, chronological age
  ipfP3GPT            Galkin    2025    transformer (Precious3GPT), age

WEIGHT AVAILABILITY - from the paper's code-availability statement
------------------------------------------------------------------
The paper releases `proteoclock` (https://github.com/Insilico-org/proteoclock),
which ships the REAL published weights for PAC and both OrganAge variants.
When that package is installed those three stop being surrogates. Each
ClockSpec below records its status:

  PROTEOCLOCK      real weights shipped in the paper's own library
  PUBLIC_REPO      a usable repository exists but needs separate work
  ON_REQUEST       authors must be emailed for the trained model
  RESTRICTED       weights exist but may not leave a controlled platform

Two clocks cannot be obtained at all:

  ipfP3GPT  proteoclock carries its feature order and nothing else. The
            weights are usable only inside the UK Biobank Research Analysis
            Platform, so they cannot be downloaded here at any access level.
  ProtAge   not in its repository; the authors must be emailed.

PAOPAC is a third case: the weights exist and are downloadable, but the
release is a Windows-only compiled extension for Python 3.9, so it does not
run on this Linux container without a reimplementation.

When real weights are absent this module fits a SURROGATE of the same model
class on the M2 reference cohort. A surrogate exercises every downstream step
(age acceleration, mixed models, enrichment) but it is NOT the published
clock, and no number it produces is comparable to the paper's.
"""
from dataclasses import dataclass

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.linear_model import ElasticNet, LogisticRegression, Ridge

# Weight-availability states, see module docstring.
PROTEOCLOCK = "PROTEOCLOCK"
PUBLIC_REPO = "PUBLIC_REPO"
ON_REQUEST = "ON_REQUEST"
RESTRICTED = "RESTRICTED"


@dataclass
class ClockSpec:
    name: str
    source: str
    target: str            # "age" or "mortality"
    model: str             # surrogate model class
    n_features: int        # proteins the published clock uses
    availability: str
    where: str             # how to obtain the real weights


# The published clocks, with the real provenance of their weights.
CLOCKS = [
    ClockSpec("ProtAge", "Argentieri 2024 (Nat Med)", "age", "boosting", 204,
              ON_REQUEST,
              "github.com/miargentieri/proteomic-age-ukb - model not in repo; "
              "email aargentieri@mgh.harvard.edu"),
    ClockSpec("OrganAge_chrono", "Goeminne 2025 (Cell Metab)", "age",
              "ridge", 500, PROTEOCLOCK,
              "proteoclock: goeminne_2025_full_chrono. Original code at "
              "github.com/ludgergoeminne/organAging"),
    ClockSpec("OrganAge_mortality", "Goeminne 2025 (Cell Metab)", "mortality",
              "logistic", 500, PROTEOCLOCK,
              "proteoclock: goeminne_2025_full_mortality (CPH model)"),
    ClockSpec("PAC", "Kuo 2024", "mortality", "elasticnet", 204,
              PROTEOCLOCK,
              "proteoclock: kuo_2024 (Gompertz, needs age). Original code at "
              "github.com/kuo-lab-uchc/PAC"),
    ClockSpec("PAOPAC", "Han 2026 (Nat Aging)", "age", "ridge", 300,
              PUBLIC_REPO,
              "github.com/JackieHanLab/PAOPAC - model.bin on Releases, but "
              "the extension is a Windows-only .pyd for Python 3.9"),
    ClockSpec("ipfP3GPT", "Galkin 2025 (Precious3GPT)", "age", "boosting", 800,
              RESTRICTED,
              "osf.io/457w8 for code; proteoclock ships feature order only. "
              "Weights usable only inside the UK Biobank RAP"),
]


class FittedClock:
    """A clock that can score an NPX matrix, real weights or surrogate."""

    def __init__(self, spec, feature_idx, model, is_surrogate, scaler=None):
        self.spec = spec
        self.feature_idx = feature_idx
        self.model = model
        self.is_surrogate = is_surrogate
        self.scaler = scaler          # (mean, sd) used to z-score features

    def predict(self, npx, meta=None):
        # `meta` is accepted and ignored: surrogates read only the proteome,
        # but PAC (a Gompertz model) needs chronological age, so every clock
        # has to take the same arguments.
        x = npx[:, self.feature_idx]
        if self.scaler is not None:
            mean, sd = self.scaler
            x = (x - mean) / sd
        if self.spec.target == "mortality" and hasattr(self.model,
                                                       "decision_function"):
            # Mortality clocks output a relative log hazard, not an age.
            return self.model.decision_function(x)
        return self.model.predict(x)


# --------------------------------------------------------- real weight path
def load_linear_weights(spec, path, protein_ids):
    """
    Build a FittedClock from published linear coefficients.

    Expects a CSV with columns `protein` and `coefficient`, plus an optional
    row named `(Intercept)`. This is the format of the organAging supplementary
    tables and of the coefficient block inside pac_proteomic_age.R, so once
    those files are on disk the published clock is used instead of a surrogate.
    Proteins missing from the local panel are dropped and reported.
    """
    tab = pd.read_csv(path)
    tab["protein"] = tab["protein"].astype(str)
    intercept = 0.0
    mask_int = tab["protein"].str.lower().isin(["(intercept)", "intercept"])
    if mask_int.any():
        intercept = float(tab.loc[mask_int, "coefficient"].iloc[0])
        tab = tab[~mask_int]

    pos = {p: i for i, p in enumerate(protein_ids)}
    hit = tab[tab["protein"].isin(pos)]
    missing = len(tab) - len(hit)
    if len(hit) == 0:
        raise ValueError(f"{spec.name}: no coefficient protein matched the "
                         f"panel; check identifier namespace")

    idx = np.array([pos[p] for p in hit["protein"]])
    coef = hit["coefficient"].to_numpy(dtype=float)

    class _Linear:
        def predict(self, x, meta=None):
            return x @ coef + intercept
        def decision_function(self, x):
            return x @ coef

    print(f"    {spec.name}: loaded {len(idx)} published coefficients "
          f"({missing} not on panel)")
    return FittedClock(spec, idx, _Linear(), is_surrogate=False)


# ----------------------------------------------------------- surrogate path
def _select_features(ref_npx, y, k):
    """Top-k proteins by absolute correlation with the training target."""
    x = ref_npx - ref_npx.mean(0)
    yc = y - y.mean()
    denom = np.sqrt((x ** 2).sum(0) * (yc ** 2).sum())
    r = np.abs((x * yc[:, None]).sum(0) / np.maximum(denom, 1e-12))
    k = min(k, ref_npx.shape[1])
    return np.sort(np.argsort(r)[::-1][:k])


def fit_surrogate(spec, ref_npx, ref_meta, seed=0):
    """Fit a stand-in clock of the published model class on the reference."""
    y = (ref_meta["age"].to_numpy() if spec.target == "age"
         else ref_meta["died_10y"].to_numpy())

    idx = _select_features(ref_npx, y.astype(float), spec.n_features)
    x = ref_npx[:, idx]
    mean, sd = x.mean(0), np.maximum(x.std(0), 1e-8)
    xz = (x - mean) / sd

    if spec.model == "ridge":
        model = Ridge(alpha=50.0, random_state=seed).fit(xz, y)
    elif spec.model == "elasticnet":
        model = ElasticNet(alpha=0.005, l1_ratio=0.5, max_iter=5000,
                           random_state=seed).fit(xz, y)
        if np.allclose(model.coef_, 0.0):
            raise RuntimeError(
                f"{spec.name}: elastic net shrank every coefficient to zero; "
                f"lower alpha or check the training target")
    elif spec.model == "logistic":
        model = LogisticRegression(C=0.05, max_iter=2000,
                                   random_state=seed).fit(xz, y)
    elif spec.model == "boosting":
        model = HistGradientBoostingRegressor(
            max_iter=200, learning_rate=0.08, max_depth=4,
            random_state=seed).fit(xz, y)
    else:
        raise ValueError(spec.model)

    return FittedClock(spec, idx, model, is_surrogate=True, scaler=(mean, sd))


def fit_all(pre, weight_dir=None, seed=0, verbose=True, use_proteoclock=True):
    """
    Assemble all six clocks, preferring real weights in this order:

      1. a CSV in `weight_dir` named `<ClockName>.csv`, which overrides
         everything so a specific coefficient set can be forced;
      2. the paper's own `proteoclock` library, if importable;
      3. a surrogate fitted on the M2 reference cohort.

    Steps 1 and 2 are the switch between a demonstration and a real
    reproduction. Whichever applies, the result reports how many clocks ended
    up real, because that number governs how any downstream figure should be
    read.
    """
    import os

    real = {}
    if use_proteoclock:
        from protclock_pipeline import proteoclock_backend
        real = proteoclock_backend.load_real_clocks(
            CLOCKS, pre.protein_ids, verbose=verbose)

    fitted = {}
    for spec in CLOCKS:
        wpath = (os.path.join(weight_dir, f"{spec.name}.csv")
                 if weight_dir else None)
        if wpath and os.path.exists(wpath):
            fitted[spec.name] = load_linear_weights(spec, wpath,
                                                    pre.protein_ids)
        elif spec.name in real:
            fitted[spec.name] = real[spec.name]
        else:
            fitted[spec.name] = fit_surrogate(spec, pre.ref_npx, pre.ref_meta,
                                              seed=seed)
    if verbose:
        n_sur = sum(f.is_surrogate for f in fitted.values())
        print(f"    {len(fitted)} clocks: {len(fitted) - n_sur} with real "
              f"published weights, {n_sur} surrogate")
    return fitted


def evaluate_holdout(pre, frac=0.25, seed=0):
    """
    Held-out accuracy of each surrogate on the reference cohort.

    Feature selection AND fitting happen on the training split only, so this
    is free of the selection leakage that makes an in-sample correlation
    meaningless. Published clocks report r about 0.9 against chronological
    age; a surrogate far above that means the simulation is too easy.
    """
    rng = np.random.default_rng(seed)
    n = pre.ref_npx.shape[0]
    perm = rng.permutation(n)
    n_test = int(round(frac * n))
    test, train = perm[:n_test], perm[n_test:]

    tr_npx = pre.ref_npx[train]
    tr_meta = pre.ref_meta.iloc[train].reset_index(drop=True)
    te_age = pre.ref_meta["age"].to_numpy()[test]

    te_meta = pre.ref_meta.iloc[test].reset_index(drop=True)
    out = {}
    for spec in CLOCKS:
        clk = fit_surrogate(spec, tr_npx, tr_meta, seed=seed)
        pred = clk.predict(pre.ref_npx[test], te_meta)
        r = float(np.corrcoef(pred, te_age)[0, 1])
        mae = (float(np.mean(np.abs(pred - te_age)))
               if spec.target == "age" else float("nan"))
        out[spec.name] = {"r": r, "mae_years": mae, "n_test": len(test)}
    return out


def report(pre=None):
    from protclock_pipeline import m2_preprocess
    if pre is None:
        pre = m2_preprocess.run(n_reference=2500)

    print("M3  SIX PROTEOMIC AGING CLOCKS")
    print("-" * 68)
    print(f"  {'clock':<19}{'target':<11}{'model':<12}{'n_prot':>7}"
          f"  weights")
    for s in CLOCKS:
        print(f"  {s.name:<19}{s.target:<11}{s.model:<12}{s.n_features:>7}"
              f"  {s.availability}")

    fitted = fit_all(pre)

    # Sanity: an age clock must track age in reference data it never saw.
    # Selection and fitting both happen inside the training split.
    out = evaluate_holdout(pre)
    print(f"\n  HELD-OUT reference accuracy (25% split, no leakage):")
    for name, d in out.items():
        mae = ("      -" if np.isnan(d["mae_years"])
               else f"{d['mae_years']:6.2f}")
        print(f"    {name:<19} r = {d['r']:+.3f}   MAE = {mae} yr")

    n_real = sum(1 for c in fitted.values() if not c.is_surrogate)
    print(f"\n  clocks running on REAL published weights: {n_real}/6")
    for s in CLOCKS:
        if fitted[s.name].is_surrogate:
            print(f"    {s.name:<19} surrogate - {s.availability}")
    return {"fitted": fitted, "holdout": out, "n_real": n_real}


if __name__ == "__main__":
    report()
