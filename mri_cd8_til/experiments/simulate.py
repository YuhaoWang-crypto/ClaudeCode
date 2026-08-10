"""A radiomics simulator with the three properties that make radiomics hard.

Simulation is used here for one specific purpose: to run the reference study's
analysis protocol on data whose true effect size is *known*, so the gap between
what the protocol reports and what is really there can be measured rather than
argued about.  That requires the simulated features to be wrong in the same ways
real ones are:

**1. They are massively correlated.**  A radiomic feature set is not 3,332
independent measurements.  GLCM contrast and GLCM dissimilarity are near
duplicates; a wavelet sub-band of a feature correlates strongly with the
unfiltered version.  Correlated features make selection *more* dangerous, not
less: a spuriously predictive feature drags its whole correlated block along.

**2. Signal lives in few features.**  Whatever real biology couples CD8
infiltration to MR texture, it does not do so through all 3,332 features
equally.  Sparse true signal is the regime where selection-based protocols
break down.

**3. Centres shift them.**  Each centre applies its own additive and
multiplicative offset per feature -- the ComBat model, used generatively -- and
centres differ in class prevalence, so "centre" is a confounder rather than
merely a nuisance.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class SimulationConfig:
    """Parameters of a synthetic multi-centre radiomics cohort.

    Defaults reproduce the reference study's dimensions: 182 patients and
    3,332 candidate features (833 per phase x 4 DCE phases) for one region.
    """

    n_patients: int = 182
    n_features: int = 3332
    #: Features carrying real signal.
    n_informative: int = 12
    #: Separation of the informative features between classes, in SDs. The
    #: default corresponds to a per-feature AUC near 0.62 -- a genuine but
    #: unremarkable biological effect.
    effect_size: float = 0.45
    #: Size of the correlated feature blocks.
    block_size: int = 40
    within_block_correlation: float = 0.75
    positive_fraction: float = 0.35
    n_centers: int = 1
    #: Per-centre additive shift SD (in feature SDs) -- the ComBat gamma.
    center_shift_sd: float = 0.0
    #: Per-centre multiplicative scale SD -- the ComBat delta.
    center_scale_sd: float = 0.0
    #: Spread of class prevalence across centres; >0 makes centre a confounder.
    center_prevalence_spread: float = 0.0
    random_state: int = 0


@dataclass
class SimulatedCohort:
    X: np.ndarray = field(repr=False)
    y: np.ndarray = field(repr=False)
    center: np.ndarray = field(repr=False)
    informative_idx: np.ndarray = field(repr=False)
    config: SimulationConfig

    @property
    def n_patients(self) -> int:
        return self.X.shape[0]

    @property
    def n_features(self) -> int:
        return self.X.shape[1]

    def train_validation_split(
        self, train_fraction: float = 137 / 182, random_state: int = 0
    ) -> tuple[np.ndarray, np.ndarray]:
        """Stratified split reproducing the reference study's 137/45 ratio."""
        rng = np.random.default_rng(random_state)
        train, valid = [], []
        for label in np.unique(self.y):
            idx = np.where(self.y == label)[0]
            rng.shuffle(idx)
            cut = int(round(train_fraction * idx.size))
            train.extend(idx[:cut])
            valid.extend(idx[cut:])
        return np.array(sorted(train)), np.array(sorted(valid))


def simulate_cohort(config: SimulationConfig | None = None) -> SimulatedCohort:
    """Generate one synthetic multi-centre radiomics cohort."""
    config = config or SimulationConfig()
    rng = np.random.default_rng(config.random_state)

    center = _assign_centers(config, rng)
    y = _assign_labels(config, center, rng)
    X = _correlated_noise(config, rng)

    informative_idx = _spread_informative_features(config, rng)
    if config.effect_size != 0.0 and informative_idx.size:
        X[:, informative_idx] += config.effect_size * y[:, None]

    X = _apply_center_effects(X, center, config, rng)
    return SimulatedCohort(X=X, y=y, center=center, informative_idx=informative_idx, config=config)


def _assign_centers(config: SimulationConfig, rng: np.random.Generator) -> np.ndarray:
    if config.n_centers <= 1:
        return np.zeros(config.n_patients, dtype=int)
    # Unequal centre sizes: real consortia never contribute equally, and equal
    # sizes make leave-one-centre-out look more stable than it is.
    weights = rng.dirichlet(np.full(config.n_centers, 2.0))
    return rng.choice(config.n_centers, size=config.n_patients, p=weights)


def _assign_labels(
    config: SimulationConfig, center: np.ndarray, rng: np.random.Generator
) -> np.ndarray:
    if config.center_prevalence_spread <= 0 or config.n_centers <= 1:
        y = (rng.random(config.n_patients) < config.positive_fraction).astype(int)
    else:
        prevalence = np.clip(
            rng.normal(config.positive_fraction, config.center_prevalence_spread, config.n_centers),
            0.05,
            0.95,
        )
        y = (rng.random(config.n_patients) < prevalence[center]).astype(int)

    # Guarantee both classes are present; a degenerate draw would silently
    # invalidate every AUC downstream.
    if len(np.unique(y)) < 2:
        y[rng.integers(config.n_patients)] ^= 1
    return y


def _correlated_noise(config: SimulationConfig, rng: np.random.Generator) -> np.ndarray:
    """Block-equicorrelated Gaussian features, generated in O(n*p).

    Each block uses the one-factor representation
    ``x = sqrt(rho)*z_block + sqrt(1-rho)*e``, which reproduces equicorrelation
    exactly without ever forming a p x p covariance matrix.
    """
    n, p = config.n_patients, config.n_features
    rho = config.within_block_correlation
    n_blocks = int(np.ceil(p / config.block_size))

    block_factor = rng.standard_normal((n, n_blocks))
    block_of = np.repeat(np.arange(n_blocks), config.block_size)[:p]
    idiosyncratic = rng.standard_normal((n, p))
    return np.sqrt(rho) * block_factor[:, block_of] + np.sqrt(1.0 - rho) * idiosyncratic


def _spread_informative_features(
    config: SimulationConfig, rng: np.random.Generator
) -> np.ndarray:
    """Place informative features in distinct blocks.

    Concentrating them in one block would make the problem artificially easy:
    a single lucky selection would capture all of the signal at once.
    """
    if config.n_informative <= 0:
        return np.array([], dtype=int)
    n_blocks = int(np.ceil(config.n_features / config.block_size))
    blocks = rng.choice(n_blocks, size=min(config.n_informative, n_blocks), replace=False)
    idx = [
        min(b * config.block_size + rng.integers(config.block_size), config.n_features - 1)
        for b in blocks
    ]
    return np.unique(np.array(idx, dtype=int))


def _apply_center_effects(
    X: np.ndarray, center: np.ndarray, config: SimulationConfig, rng: np.random.Generator
) -> np.ndarray:
    if config.n_centers <= 1 or (config.center_shift_sd == 0 and config.center_scale_sd == 0):
        return X
    out = X.copy()
    for c in np.unique(center):
        rows = center == c
        gamma = rng.normal(0.0, config.center_shift_sd, X.shape[1])
        delta = np.exp(rng.normal(0.0, config.center_scale_sd, X.shape[1]))
        out[rows] = out[rows] * delta + gamma
    return out


def theoretical_feature_auc(effect_size: float) -> float:
    """AUC of a single informative feature, for a given standardised effect.

    For two unit-variance Gaussians separated by ``d``, ``AUC = Phi(d/sqrt(2))``.
    Lets a simulation be specified by the effect size a reader can interpret.
    """
    from math import erf

    return 0.5 * (1.0 + erf(effect_size / 2.0))


def effect_size_for_auc(target_auc: float) -> float:
    """Inverse of :func:`theoretical_feature_auc`."""
    from mri_cd8_til.modeling.metrics import _erfinv

    return 2.0 * _erfinv(2.0 * target_auc - 1.0)
