"""How much of a 0.99 AUC does the analysis protocol produce by itself?

The reference study reports validation AUCs of 0.985 and 0.984 from 182
single-centre patients and ~3,300 candidate features per region.  Rather than
argue about whether that is plausible, this module measures it: run the
published protocol on synthetic data whose true effect size is known, including
the case where the true effect size is exactly **zero**, and read off what it
reports.

Three protocols are compared on identical data:

``published``
    Rank all features by univariate AUC over the **whole cohort**, keep the top
    50, fit LASSO on the 137 training patients, report AUC on the 45 validation
    patients.  This is the protocol as a reader would implement it from the
    paper's description.

``split_first``
    Identical, except the ranking sees only the 137 training patients.  The
    single defensible change.

``repeated_cv``
    The whole pipeline -- scaling, selection, LASSO -- re-fitted inside every
    fold of repeated stratified cross-validation.

The difference ``published - split_first`` is optimism attributable to one step
in the analysis.  The difference ``published - repeated_cv`` is the total gap
between the headline number and a defensible internal estimate.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..modeling.pipeline import PipelineConfig, build_pipeline, run_reference_protocol
from ..modeling.validation import repeated_cv
from .simulate import SimulationConfig, effect_size_for_auc, simulate_cohort


@dataclass
class ProtocolScores:
    """Replicate-level validation AUCs for one protocol at one effect size."""

    protocol: str
    true_feature_auc: float
    scores: np.ndarray = field(repr=False)

    @property
    def mean(self) -> float:
        return float(np.nanmean(self.scores))

    @property
    def median(self) -> float:
        return float(np.nanmedian(self.scores))

    def quantile(self, q: float) -> float:
        return float(np.nanquantile(self.scores, q))

    def fraction_above(self, threshold: float) -> float:
        return float(np.nanmean(self.scores >= threshold))

    def describe(self) -> str:
        return (
            f"{self.protocol:<12} true_feature_AUC={self.true_feature_auc:.2f}  "
            f"reported={self.median:.3f} "
            f"[{self.quantile(0.05):.3f}, {self.quantile(0.95):.3f}]  "
            f"P(>=0.95)={self.fraction_above(0.95):.0%}"
        )


@dataclass
class OptimismResult:
    """Everything the experiment produced, keyed by (protocol, effect size)."""

    scores: dict[tuple[str, float], ProtocolScores]
    true_feature_aucs: tuple[float, ...]
    protocols: tuple[str, ...]
    config: SimulationConfig

    def get(self, protocol: str, true_auc: float) -> ProtocolScores:
        return self.scores[(protocol, true_auc)]

    def optimism(self, true_auc: float, baseline: str = "split_first") -> float:
        """AUC inflation of ``published`` over ``baseline`` at one effect size."""
        return self.get("published", true_auc).median - self.get(baseline, true_auc).median

    def as_table(self) -> list[dict[str, float | str]]:
        return [
            {
                "protocol": p,
                "true_feature_auc": a,
                "median_reported_auc": self.get(p, a).median,
                "p05": self.get(p, a).quantile(0.05),
                "p95": self.get(p, a).quantile(0.95),
                "frac_ge_0.95": self.get(p, a).fraction_above(0.95),
            }
            for a in self.true_feature_aucs
            for p in self.protocols
        ]


def run_optimism_experiment(
    true_feature_aucs: tuple[float, ...] = (0.50, 0.55, 0.60, 0.65, 0.70),
    n_replicates: int = 100,
    base_config: SimulationConfig | None = None,
    pipeline_config: PipelineConfig | None = None,
    cv_repeats: int = 4,
    include_repeated_cv: bool = True,
    seed: int = 0,
    progress: bool = False,
) -> OptimismResult:
    """Run the three protocols across a grid of true effect sizes."""
    base_config = base_config or SimulationConfig()
    pipeline_config = pipeline_config or PipelineConfig()
    protocols = ("published", "split_first") + (("repeated_cv",) if include_repeated_cv else ())

    collected: dict[tuple[str, float], list[float]] = {
        (p, a): [] for p in protocols for a in true_feature_aucs
    }

    for auc_target in true_feature_aucs:
        effect = effect_size_for_auc(auc_target)
        for replicate in range(n_replicates):
            config = SimulationConfig(
                **{**base_config.__dict__, "effect_size": effect, "random_state": seed + replicate}
            )
            cohort = simulate_cohort(config)
            train_idx, valid_idx = cohort.train_validation_split(random_state=seed + replicate)

            for select_on_all, protocol in ((True, "published"), (False, "split_first")):
                result = run_reference_protocol(
                    cohort.X,
                    cohort.y,
                    train_idx,
                    valid_idx,
                    config=pipeline_config,
                    select_on_all_data=select_on_all,
                )
                collected[(protocol, auc_target)].append(result.validation_auc)

            if include_repeated_cv:
                cv = repeated_cv(
                    build_pipeline(pipeline_config),
                    cohort.X,
                    cohort.y,
                    n_splits=5,
                    n_repeats=cv_repeats,
                    random_state=seed + replicate,
                )
                collected[("repeated_cv", auc_target)].append(cv.mean)

        if progress:
            print(f"  true_feature_AUC={auc_target:.2f} done ({n_replicates} replicates)")

    scores = {
        key: ProtocolScores(protocol=key[0], true_feature_auc=key[1], scores=np.array(values))
        for key, values in collected.items()
    }
    return OptimismResult(
        scores=scores,
        true_feature_aucs=tuple(true_feature_aucs),
        protocols=protocols,
        config=base_config,
    )


@dataclass
class DimensionSweep:
    """How the leak scales with the number of candidate features."""

    n_features: tuple[int, ...]
    median_auc: np.ndarray = field(repr=False)
    p95_auc: np.ndarray = field(repr=False)
    frac_above_95: np.ndarray = field(repr=False)

    def describe(self) -> str:
        lines = ["  n_features   median reported AUC   95th pct   P(>=0.95)"]
        for p, m, u, f in zip(self.n_features, self.median_auc, self.p95_auc, self.frac_above_95):
            lines.append(f"  {p:>10d}   {m:>18.3f}   {u:>8.3f}   {f:>8.0%}")
        return "\n".join(lines)


def sweep_feature_dimension(
    n_features_grid: tuple[int, ...] = (50, 200, 800, 1666, 3332, 6664),
    n_replicates: int = 100,
    base_config: SimulationConfig | None = None,
    pipeline_config: PipelineConfig | None = None,
    seed: int = 0,
) -> DimensionSweep:
    """Under **zero** true signal, how does reported AUC grow with p?

    This isolates the mechanism.  With no signal anywhere in the data, any AUC
    above 0.5 is manufactured, and if it grows with the number of candidate
    features then the cause is the number of chances the selection step gets --
    not the biology, the model, or the imaging.
    """
    base_config = base_config or SimulationConfig()
    pipeline_config = pipeline_config or PipelineConfig()

    median, p95, frac = [], [], []
    for n_features in n_features_grid:
        scores = []
        for replicate in range(n_replicates):
            config = SimulationConfig(
                **{
                    **base_config.__dict__,
                    "n_features": n_features,
                    "effect_size": 0.0,
                    "n_informative": 0,
                    "random_state": seed + replicate,
                }
            )
            cohort = simulate_cohort(config)
            train_idx, valid_idx = cohort.train_validation_split(random_state=seed + replicate)
            scores.append(
                run_reference_protocol(
                    cohort.X,
                    cohort.y,
                    train_idx,
                    valid_idx,
                    config=pipeline_config,
                    select_on_all_data=True,
                ).validation_auc
            )
        arr = np.array(scores)
        median.append(float(np.nanmedian(arr)))
        p95.append(float(np.nanquantile(arr, 0.95)))
        frac.append(float(np.nanmean(arr >= 0.95)))

    return DimensionSweep(
        n_features=tuple(n_features_grid),
        median_auc=np.array(median),
        p95_auc=np.array(p95),
        frac_above_95=np.array(frac),
    )
