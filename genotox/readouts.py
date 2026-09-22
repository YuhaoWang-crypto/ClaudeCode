"""
Readout layer — turns core trajectories into what an instrument reports.

Kept separate from the cores so that one simulated biology can be measured
several ways, which is the situation the real kits are in: the same DNA
damage is read as beta-gal absorbance (umu), tail moment (comet), or
micronucleus frequency, and the differences between those endpoints are
partly differences of *measurement*, not of biology.

A readout consumes the observables dict a core emits and returns a flat dict
of scalars.  It must declare which observables it needs, so a mismatch fails
loudly at wiring time instead of silently reading a zero.

Contract with the decision layer
--------------------------------
Exactly one readout per assay emits ``signal`` and ``biomass``; everything
downstream forms induction ratios and growth factors from those two names and
nothing else.  The instrument-specific names (``A410``, ``fluorescence``, ...)
travel alongside for reporting.  Without this, the decision layer would have
to know that umu means absorbance and a reporter line means fluorescence, and
the layer split would be fictional.
"""
from __future__ import annotations

from dataclasses import dataclass


class Readout:
    name = "abstract"
    requires: tuple = ()

    def check(self, obs: dict) -> None:
        missing = [k for k in self.requires if k not in obs]
        if missing:
            raise KeyError(
                f"readout {self.name!r} needs observables {missing} which "
                f"core does not expose (has: {sorted(obs)})")

    def __call__(self, obs: dict) -> dict:
        raise NotImplementedError


@dataclass
class BetaGalONPG(Readout):
    """The umu test's actual measurement.

    beta-galactosidase hydrolyses ONPG to o-nitrophenol, quantified at
    410 nm.  Absorbance is proportional to *total* enzyme in the well, i.e.
    to (enzyme per cell) x (cell density), which is why the assay must be
    normalised by growth (A600) before any induction ratio is formed.  That
    normalisation is the step where a purely cytotoxic compound would
    otherwise masquerade as an inducer.

    Substrate depletion is not modelled (assumes ONPG saturating over the
    incubation); at very high enzyme this would overestimate A410.
    """

    name = "umu_betagal_ONPG"
    requires = ("reporter_enzyme", "density")
    onpg_min: float = 30.0       # ONPG incubation
    gain: float = 1.0            # lumps kcat, path length, extinction coeff.

    def __call__(self, obs: dict) -> dict:
        self.check(obs)
        B_end = float(obs["reporter_enzyme"][-1])
        N_end = float(obs["density"][-1])
        a410 = self.gain * B_end * N_end * self.onpg_min
        a600 = N_end
        # Miller-style specific activity: absorbance per unit biomass per min
        units = a410 / (a600 * self.onpg_min) if a600 > 0 else 0.0
        return {"A410": a410, "A600": a600, "units": units,
                "enzyme_per_cell": B_end,
                "signal": units, "biomass": a600}


@dataclass
class GrowthReadout(Readout):
    """Culture density, used as the assay's own validity gate.

    The umu protocol discards wells whose growth factor falls below ~0.5,
    because at that point the induction ratio is no longer interpretable.
    Reporting the gate alongside the signal is not optional bookkeeping — an
    unreported failed gate is how a cytotoxic compound gets called genotoxic.
    """

    name = "growth"
    requires = ("density", "viability")

    def __call__(self, obs: dict) -> dict:
        self.check(obs)
        return {"density_end": float(obs["density"][-1]),
                "viability_end": float(obs["viability"][-1])}


@dataclass
class DamageProbe(Readout):
    """Internal state, exposed as a pseudo-measurement.

    Not something the umu kit reports.  It exists because lesion burden is
    the quantity shared with comet (break density) and with gamma-H2AX foci,
    so it is the natural place to cross-calibrate endpoint #3 against this
    core rather than re-deriving damage from scratch.
    """

    name = "damage_probe"
    requires = ("lesions",)
    #: reported when the core happens to expose them; cores differ here and
    #: that is fine, because nothing downstream depends on these.
    optional = ("recA_active", "promoter", "ddr_active", "mitotic_stress")

    def __call__(self, obs: dict) -> dict:
        self.check(obs)
        out = {"lesions_end": float(obs["lesions"][-1]),
               "lesions_max": float(obs["lesions"].max())}
        for k in self.optional:
            if k in obs:
                out[f"{k}_end"] = float(obs[k][-1])
                out[f"{k}_max"] = float(obs[k].max())
        return out


@dataclass
class ReporterFluorescence(Readout):
    """Mammalian reporter-line readout (GADD45a-GFP and relatives).

    Two facts about GFP that the enzyme readout does not share, and that
    change what the assay can resolve:

    * it needs a maturation step, so fluorescence lags synthesis;
    * it is stable, so it *integrates* — but it integrates promoter activity,
      not p53.  Since the GADD45a promoter is a steep switch, p53 excursions
      below its threshold add area to the p53 integral and almost nothing to
      the fluorescence.  The assay is therefore blind to low-level sustained
      p53, and two exposures with equal p53 AUC can read very differently.
      :class:`PulseProbe` exposes all three quantities so the distinction is
      testable rather than asserted.
    """

    name = "reporter_gfp"
    requires = ("mature_gfp", "density")

    def __call__(self, obs: dict) -> dict:
        self.check(obs)
        f_cell = float(obs["mature_gfp"][-1])
        n_end = float(obs["density"][-1])
        return {"fluorescence_per_cell": f_cell,
                "fluorescence_total": f_cell * n_end,
                "cell_density": n_end,
                "signal": f_cell, "biomass": n_end}


@dataclass
class PulseProbe(Readout):
    """Dynamics of the upstream signal, which no plate reader sees.

    Counts p53 peaks and reports peak height versus time-integral.  Used to
    check the claim that a stable reporter reads the integral: two exposures
    with similar AUC but different pulse structure should give similar
    fluorescence, and this is the only place that is visible.
    """

    name = "pulse_probe"
    requires = ("t", "p53")
    #: a trace whose peak-to-trough span is below this fraction of its mean
    #: is treated as flat and scores zero pulses.  Without an absolute floor,
    #: a relative-prominence rule happily counts integrator ripple on a
    #: constant trace, which is how an undamaged control ends up "pulsing".
    min_rel_span: float = 0.10
    rel_prominence: float = 0.20

    def __call__(self, obs: dict) -> dict:
        import numpy as np
        from scipy.signal import find_peaks
        self.check(obs)
        t, p = np.asarray(obs["t"], float), np.asarray(obs["p53"], float)
        span, mean = p.max() - p.min(), p.mean()
        n = 0
        if mean > 0 and span / mean >= self.min_rel_span:
            peaks, _ = find_peaks(p, prominence=self.rel_prominence * span)
            n = int(len(peaks))
        out = {"p53_peak": float(p.max()),
               "p53_mean": float(p.mean()),
               "p53_auc": float(np.trapezoid(p, t)),
               "p53_pulses": n}
        if "gadd_promoter" in obs:
            pr = np.asarray(obs["gadd_promoter"], float)
            out["promoter_auc"] = float(np.trapezoid(pr, t))
            out["promoter_peak"] = float(pr.max())
        return out


#: Readouts available by name.  Endpoint #3 registers here.
REGISTRY = {r.name: r for r in (BetaGalONPG(), GrowthReadout(), DamageProbe(),
                                ReporterFluorescence(), PulseProbe())}


def apply_readouts(obs: dict, names=None) -> dict:
    names = list(REGISTRY) if names is None else list(names)
    out = {}
    for n in names:
        if n not in REGISTRY:
            raise KeyError(f"unknown readout {n!r}; have {sorted(REGISTRY)}")
        out.update(REGISTRY[n](obs))
    return out
