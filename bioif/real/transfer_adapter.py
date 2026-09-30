"""
bioif.real.transfer_adapter -- the first adapter in this package whose
uncertainty is earned rather than asserted.

`ConformalAssayTransfer` answers "may I reuse this potency, measured in assay
A, where assay B is needed?" It is fitted on compounds measured in both, it
emits a conformally calibrated interval whose coverage has been checked
empirically (see `conformal.py` and `demo_conformal.py`), and it refuses
outright for any assay pair it was not calibrated on -- because the
leave-one-pair-out experiment shows the guarantee does not survive that move
(mean coverage 0.61 against a nominal 0.80, worst case 0.00).

Three contract details are doing real work here:

  * The adapter consumes and produces the SAME quantity (pIC50) and changes
    only `context.assay`. Routing therefore never inserts it on its own: an
    assay transfer is a deliberate act, and has to be asked for.
  * Its kind is EMPIRICAL, so the evidence level of a measured input drops to
    `calibrated_prediction` on the way through. A measurement transferred to
    another assay is no longer a measurement of that assay. That demotion is
    now backed by a coverage number instead of a convention.
  * Calibration is per assay pair -- Mondrian conformal -- because the
    experiment shows pooled calibration meets its marginal guarantee (0.799
    vs 0.80) while under-covering 11 of 18 individual pairs, worst 0.009.
    The class that restores conditional coverage is the assay: the context
    field the contract already carries.
"""
from __future__ import annotations

import math
import statistics

from ..adapter import Adapter, EMPIRICAL
from ..core import Context, EXTRAPOLATE, IN_DOMAIN, REFUSE, Verdict
from ..quantities import PIC50
from . import conformal as C


class ConformalAssayTransfer(Adapter):
    """Reuse a pIC50 from one assay in another, with a checked interval."""

    kind = EMPIRICAL
    consumes = PIC50
    produces = PIC50
    requires_context = ("assay",)

    def __init__(self, pair: C.PairInfo, cal: C.Calibration,
                 residuals: list[float], alpha: float = 0.2):
        self.pair = pair
        self.cal = cal
        self.residuals = residuals
        self.alpha = alpha
        self.name = f"assay-transfer {pair.source[-6:]}->{pair.target[-6:]}"
        self.version = f"conformal/a={alpha}"
        self.calibrated_systems = ()

    # -- the guarantee, stated exactly ------------------------------------
    def conformal_interval(self, x: float, alpha: float | None = None):
        a = self.alpha if alpha is None else alpha
        return self.cal.interval(x, a)

    def domain(self, claim) -> Verdict:
        base = super().domain(claim)
        if base.status == REFUSE:
            return base
        if claim.context.assay != self.pair.source:
            return Verdict(REFUSE,
                           f"{self.name} is calibrated for source assay "
                           f"{self.pair.source}; this claim is from "
                           f"{claim.context.assay}. Leave-one-pair-out "
                           f"coverage for an uncalibrated pair is 0.61 "
                           f"against a nominal 0.80 (worst 0.00), so there "
                           f"is no interval to offer")
        if not self.cal.attainable(self.alpha):
            return Verdict(REFUSE,
                           f"{self.cal.n_cal} calibration points cannot "
                           f"support alpha={self.alpha}")
        lo, hi = self.cal.x_range
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
                       (("transferred_from", self.pair.source),))

    def _forward(self, xs, claim, rng, noise):
        """
        Resample the calibration residuals around the shifted prediction.

        This gives the full empirical predictive distribution; its central
        (1-alpha) interval agrees with `conformal_interval` up to the
        finite-sample (n+1)/n correction, which conformal spends to buy the
        guarantee. With noise off (variance attribution) only the shift
        applies.
        """
        if not noise or not self.residuals:
            return [x + self.cal.shift for x in xs]
        return [x + self.cal.shift + rng.choice(self.residuals) for x in xs]


def build_transfer(pair_key: str, alpha: float = 0.2,
                   seed: int = 11) -> ConformalAssayTransfer:
    """Fit a deployable transfer adapter for one assay pair from the snapshot."""
    pairs = C.discover_pairs()
    info = pairs[pair_key]
    if info.status != "continuous":
        raise ValueError(
            f"{pair_key} is {info.status}: {info.reason} -- a conformal "
            f"regression interval on binned data would be meaningless")
    pts = C.build_points(pairs)
    cal = C.fit_deployable(pts, pair_key, alphas=(alpha, 0.1), seed=seed)
    import random
    tr, ca, _ = C.split_by_compound([p for p in pts if p.pair == pair_key],
                                    random.Random(seed), fracs=(0.5, 0.5, 0.0))
    residuals = [p.y - C.predict(p.x, cal.shift) for p in ca]
    return ConformalAssayTransfer(info, cal, residuals, alpha)
