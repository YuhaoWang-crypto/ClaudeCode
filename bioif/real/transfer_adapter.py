"""
bioif.real.transfer_adapter -- the first adapters in this package whose
uncertainty is earned rather than asserted, and the first hop with more than
one model competing for it.

`ConformalAssayTransfer` answers "may I reuse this potency, measured in assay
A, where the downstream model wants assay B?" It is fitted on compounds
measured in both, emits a conformally calibrated interval whose coverage has
been checked empirically, and refuses any assay pair it was not calibrated on
-- because leave-one-pair-out coverage for an uncalibrated pair is 0.61
against a nominal 0.80, worst case 0.00.

Several of these are registered at once. They differ in two ways, and the
registry uses both:

  * different SOURCE assays -- only one is in domain for a given claim, so
    the domain test alone picks the right family;
  * different PREDICTORS for the same source assay (identity / shift /
    linear) -- all conformally valid, so the tie is broken on calibrated
    interval width.

Contract details doing real work:

  * consumes pIC50, produces pIC50-in-the-reference-assay. Two different
    quantities, so "which assay is this number from" is a type error rather
    than a comment, and a raw measurement cannot reach the downstream model
    without an explicit, calibrated transfer.
  * kind is EMPIRICAL, so a measured input leaves as `calibrated_prediction`.
    A measurement transferred to another assay is not a measurement of that
    assay, and the demotion is now backed by a coverage number.
"""
from __future__ import annotations

import math

from ..adapter import Adapter, EMPIRICAL
from ..core import Context, EXTRAPOLATE, IN_DOMAIN, REFUSE, Verdict
from ..quantities import PIC50, PIC50_REF, REFERENCE_ASSAY
from . import conformal as C
from .transfer_models import ConformalModel, PREDICTORS, fit_conformal


class ConformalAssayTransfer(Adapter):
    """Express a pIC50 from one assay in the reference assay, with an interval."""

    kind = EMPIRICAL
    consumes = PIC50
    produces = PIC50_REF
    requires_context = ("assay",)

    def __init__(self, pair: C.PairInfo, model: ConformalModel,
                 alpha: float = 0.2):
        self.pair = pair
        self.model = model
        self.alpha = alpha
        self.name = (f"transfer {pair.source[-6:]}->{pair.target[-6:]} "
                     f"[{model.predictor.name}]")
        self.version = f"conformal/a={alpha}"
        self.calibrated_systems = ()

    # -- what the registry ranks on ---------------------------------------
    def calibrated_width(self, alpha: float = 0.2) -> float:
        return self.model.width(alpha)

    def conformal_interval(self, x: float, alpha: float | None = None):
        return self.model.interval(x, self.alpha if alpha is None else alpha)

    def domain(self, claim) -> Verdict:
        base = super().domain(claim)
        if base.status == REFUSE:
            return base
        if claim.context.assay != self.pair.source:
            return Verdict(REFUSE,
                           f"{self.name} is calibrated for source assay "
                           f"{self.pair.source}; this claim is from "
                           f"{claim.context.assay}. Cross-pair coverage is "
                           f"0.61 against a nominal 0.80 (worst 0.00), so "
                           f"there is no interval to offer")
        if not self.model.attainable(self.alpha):
            return Verdict(REFUSE,
                           f"{self.model.n_cal} calibration points cannot "
                           f"support alpha={self.alpha}")
        lo, hi = self.model.x_range
        x = claim.estimate.mean
        if not (lo <= x <= hi):
            return Verdict(EXTRAPOLATE,
                           f"source pIC50 {x:.2f} is outside the calibrated "
                           f"range [{lo:.2f}, {hi:.2f}]", inflate=1.5)
        return Verdict(IN_DOMAIN)

    def map_context(self, c: Context, e) -> Context:
        return Context(system=c.system, dose_uM=c.dose_uM, time_h=c.time_h,
                       assay=self.pair.target,
                       covariates=c.covariates +
                       (("transferred_from", self.pair.source),
                        ("transfer_model", self.model.predictor.name)))

    def _forward(self, xs, claim, rng, noise):
        """
        Resample the calibration residuals around the model's prediction.

        The central (1-alpha) interval of the result agrees with
        `conformal_interval` up to the finite-sample (n+1)/n correction that
        conformal spends to buy the guarantee. With noise off (variance
        attribution) only the point prediction applies.
        """
        if not noise or not self.model.residuals:
            return [self.model.predict(x) for x in xs]
        return [self.model.predict(x) + rng.choice(self.model.residuals)
                for x in xs]


def build_transfer(pair_key: str, predictor: str = "shift", alpha: float = 0.2,
                   seed: int = 11) -> ConformalAssayTransfer:
    """Fit one deployable transfer adapter from the committed snapshot."""
    pairs = C.discover_pairs()
    info = pairs[pair_key]
    if info.status != "continuous":
        raise ValueError(
            f"{pair_key} is {info.status}: {info.reason} -- a conformal "
            f"regression interval on binned data would be meaningless")
    pts = C.build_points(pairs)
    model = fit_conformal(pts, pair_key, predictor, (alpha, 0.1), seed)
    return ConformalAssayTransfer(info, model, alpha)


def build_all_transfers(target: str = REFERENCE_ASSAY, alpha: float = 0.2,
                        seed: int = 11) -> list[ConformalAssayTransfer]:
    """
    Every continuous transfer into the reference assay, times every predictor.

    This is the registry's competition: several source assays (only one of
    which is in domain for any given claim) and several models per source
    (all valid, differing in width).
    """
    pairs = C.discover_pairs()
    pts = C.build_points(pairs)
    out = []
    for key, info in pairs.items():
        if info.target != target or info.status != "continuous":
            continue
        for pname in PREDICTORS:
            m = fit_conformal(pts, key, pname, (alpha, 0.1), seed)
            if m.attainable(alpha):
                out.append(ConformalAssayTransfer(info, m, alpha))
    return out
