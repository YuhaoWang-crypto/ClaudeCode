"""
bioif.real.tox_adapters -- three pairings of the interface map as typed edges.

    C->F   CompoundToP53        SMILES  -> P(SR-p53 reporter active)
    C->G   CompoundToAmes       SMILES  -> P(Ames positive)
    F->G   P53ToMutagenicity    P(p53)  -> P(Ames positive)

The third is the interesting one. The blueprint's own Demo B marks the
p53 -> mutagenicity edge `blocked`, correctly, because there was no bridge to
justify it. There is now a measured one: 2,064 compounds are assayed in both
Tox21 and the Hansen Ames benchmark, which gives the edge a 2x2 and therefore
an effect size. So the edge can be published -- but only at the evidence
level the 2x2 supports.

That is what the new `ASSOCIATION` edge kind is for. The transform is not a
fitted model and not an algebraic identity: it is a measured co-occurrence
between two readouts, empirically grounded and not causal. It caps evidence
at `inferred_association`, which is stronger than a model prediction and
weaker than a measurement -- exactly the rung the blueprint's ladder has for
it, and a rung bioif did not previously have.

The edge also refuses to pretend it is a screen. `screening_verdict()`
reports sensitivity alongside enrichment, because an edge with RR 1.81 and
sensitivity 0.13 is useful for ordering a list and useless for clearing one.
"""
from __future__ import annotations

import math

from ..adapter import Adapter, ASSOCIATION, EMPIRICAL
from ..core import (Claim, Context, Entity, Estimate, IN_DOMAIN, Provenance,
                    REFUSE, Verdict)
from ..quantities import AMES_POSITIVE, P53_ACTIVE
from . import qsar, tox

#: A compound claim needs a structure to featurize, which lives in context
#: rather than in the identifier: an InChIKey skeleton is not invertible to a
#: structure, so the SMILES has to travel with the claim.
SMILES_COV = "smiles"


class _QsarAdapter(Adapter):
    """Shared machinery: featurize the claim's SMILES, score, conformalise."""

    kind = EMPIRICAL
    requires_covariates = (SMILES_COV,)

    def __init__(self, dataset: str, produces, alpha: float = 0.1,
                 seed: int = 0, name: str | None = None):
        built = qsar.build(dataset, alpha=alpha, seed=seed)
        self.cc = built["cc"]
        self.ev = built["eval"]
        self.dataset = dataset
        self.alpha = alpha
        self.produces = produces
        self.name = name or f"{dataset}-qsar"
        self.version = f"logreg-morgan2048/scaffold-split/a={alpha}"

    def calibrated_width(self, alpha: float = 0.2) -> float:
        """
        A classifier has no interval width, so report the abstention rate
        instead: the fraction of predictions that come back as {0,1}. It is
        the same quantity in spirit -- how much the model declines to say --
        and it lets the registry rank two classifiers for one hop.
        """
        return float(self.ev.set_sizes.get(2, 1.0))

    def domain(self, claim) -> Verdict:
        base = super().domain(claim)
        if base.status == REFUSE:
            return base
        smi = dict(claim.context.covariates).get(SMILES_COV, "")
        if qsar.featurize([smi])[1][0] == False:            # unparseable
            return Verdict(REFUSE,
                           f"{self.name}: RDKit cannot parse the SMILES "
                           f"{smi[:40]!r}; refusing rather than scoring a "
                           f"zero fingerprint")
        return Verdict(IN_DOMAIN)

    def prediction_set(self, smiles: str) -> frozenset:
        X, _ = qsar.featurize([smiles])
        return self.cc.predict_set(X)[0]

    def _forward(self, xs, claim, rng, noise):
        """
        The incoming value is ignored: this is a source-like edge whose real
        input is the structure in context. Samples are drawn as Bernoulli
        around the model's probability so the claim carries a distribution
        rather than a bare score, and so an abstaining prediction set widens
        it to the full [0,1] range.
        """
        smi = dict(claim.context.covariates)[SMILES_COV]
        X, _ = qsar.featurize([smi])
        p = float(self.cc.proba(X)[0])
        pset = self.cc.predict_set(X)[0]
        if not noise:
            return [p] * len(xs)
        if len(pset) != 1:
            # the conformal set abstains (or rejects both labels): the honest
            # output is the whole admissible range, not a confident p
            return [rng.random() for _ in xs]
        return [min(max(p + rng.gauss(0.0, 0.05), 0.0), 1.0) for _ in xs]


class CompoundToP53(_QsarAdapter):
    """C->F : structure -> P(Tox21 SR-p53 reporter active)."""

    def __init__(self, alpha: float = 0.1, seed: int = 0):
        super().__init__("p53", P53_ACTIVE, alpha, seed, name="C->F p53-qsar")
        self.consumes = P53_ACTIVE          # source-like; see _forward


class CompoundToAmes(_QsarAdapter):
    """C->G : structure -> P(Ames positive), one hop, no intermediate."""

    def __init__(self, alpha: float = 0.1, seed: int = 0):
        super().__init__("ames", AMES_POSITIVE, alpha, seed,
                         name="C->G ames-qsar")
        self.consumes = AMES_POSITIVE


#: Optional covariate: the measured general-cytotoxicity call (SR-MMP), as
#: "1" / "0". Supplying it changes the edge materially -- see below.
COV_CYTOTOX = "cytotox_mmp"


class P53ToMutagenicity(Adapter):
    """
    F->G : the blueprint's `blocked` edge, graded by a measured 2x2 -- and
    **conditional on the cytotoxicity readout**, which is the part that
    matters.

    The first version of this adapter carried one global 2x2 (OR 2.81,
    p=5e-9) and was wrong in a way a marginal analysis cannot see. An
    independent report on the same datasets did the sharper thing and
    stratified on a general-cytotoxicity readout; stratifying here, on 32%
    more paired compounds, replicates it:

        marginal          OR 2.81   p=5e-09    informative
        SR-MMP+ stratum   OR 1.13   p=0.77     NOT informative
        SR-MMP- stratum   OR 2.90   p=1.3e-03  informative

    So a p53 reporter positive says something about mutagenicity **only in
    compounds that are not already flagged cytotoxic**. Among cytotoxic
    compounds it says nothing, and the marginal number quietly averages the
    two regimes together.

    The adapter therefore behaves differently depending on what it is told:

      cytotox known, positive -> REFUSES to update. It returns the stratum's
          base rate and flags that the reporter adds nothing here. A 1.13
          odds ratio at p=0.77 is not an edge.
      cytotox known, negative -> uses the MMP- conditionals, the strongest
          and cleanest form of the edge.
      cytotox unknown -> uses the marginal conditionals but flags the claim
          as confounded, and `value_of_knowing_cytotox` reports how far the
          answer would move if the readout were supplied.

    That last method is the useful one for experiment planning: it prices one
    extra measurement in units of the thing being predicted.
    """

    name = "F->G p53->mutagenicity"
    version = "tox21xames-2x2/stratified"
    kind = ASSOCIATION
    consumes = P53_ACTIVE
    produces = AMES_POSITIVE
    assumption = ("a reporter positive is treated as evidence about "
                  "mutagenicity via measured co-occurrence only; no "
                  "mechanism is claimed, and Ames is not the Comet or "
                  "adduct endpoint")

    def __init__(self):
        from .cross_check import stratified
        self.t = tox.two_by_two(tox.P53)
        self.strata = stratified(tox.CYTOTOX_CONTROL)
        self._p_pos = self.t.ppv
        self._p_neg = self.t.base

    # -- which conditionals apply to this claim ---------------------------
    def _regime(self, claim) -> tuple[str, object]:
        v = dict(claim.context.covariates).get(COV_CYTOTOX)
        if v is None:
            return "marginal", self.strata["marginal"]
        return ("cytotoxic", self.strata["+"]) if str(v) in ("1", "True", "true") \
            else ("non-cytotoxic", self.strata["-"])

    def calibrated_width(self, alpha: float = 0.2) -> float:
        """A wider gap between the two conditional risks is a more
        informative edge, so the registry prefers it."""
        return 1.0 - abs(self._p_pos - self._p_neg)

    def screening_verdict(self) -> str:
        t = self.t
        return (f"RR {t.risk_ratio:.2f}, OR {t.odds_ratio:.2f} "
                f"(Fisher p={t.fisher_p():.1e}) on n={t.n} -- but sensitivity "
                f"{t.sensitivity:.2f}: this edge ranks candidates, it does "
                f"not clear them")

    def stratum_table(self) -> str:
        L = [f"  {'stratum':<16}{'n':>6}{'OR':>7}{'p':>10}  informative"]
        for key in ("marginal", "+", "-"):
            s = self.strata[key]
            L.append(f"  {s.name:<16}{s.n:>6}{s.odds_ratio:>7.2f}"
                     f"{s.p:>10.1e}  {'yes' if s.informative else 'NO'}")
        return "\n".join(L)

    def value_of_knowing_cytotox(self, p_p53: float) -> float:
        """
        How far the predicted Ames probability moves between the two strata,
        at a given p53 probability. This is the edge's own value-of-
        information statement: it prices the cytotoxicity measurement in
        units of the endpoint.
        """
        hi = self._apply(p_p53, self.strata["-"])
        lo = self._apply(p_p53, self.strata["+"], refuse_update=True)
        return abs(hi - lo)

    @staticmethod
    def _apply(p, stratum, refuse_update: bool = False) -> float:
        if refuse_update:
            # No usable signal in this stratum, so the reporter result is
            # discarded and the stratum's own Ames prevalence is returned:
            # (a + c) / n, recovered from the 2x2's margins.
            n_pos, n = stratum.n_pos, max(stratum.n, 1)
            return (stratum.ppv * n_pos + stratum.base * (n - n_pos)) / n
        return p * stratum.ppv + (1.0 - p) * stratum.base

    def map_context(self, c: Context, e: Entity) -> Context:
        return c.with_cov("bridge", "tox21xames 2x2 (stratified)")

    def apply(self, claim, rng, noise: bool = True, inflate: float = 1.0):
        out = super().apply(claim, rng, noise=noise, inflate=inflate)
        regime, s = self._regime(claim)
        if regime == "cytotoxic":
            out = out.add_flag(
                f"[no-update] this compound is flagged cytotoxic "
                f"({tox.CYTOTOX_CONTROL}+), where the p53->mutagenicity "
                f"association is absent (OR {s.odds_ratio:.2f}, p={s.p:.2f}, "
                f"n={s.n}); the reporter result was NOT used")
        elif regime == "marginal":
            voi = self.value_of_knowing_cytotox(float(claim.estimate.mean))
            out = out.add_flag(
                f"[confounded] {tox.CYTOTOX_CONTROL} is unknown, so the "
                f"marginal 2x2 was used; measuring it would move this "
                f"probability by up to {voi:.3f}")
        return out

    def _forward(self, xs, claim, rng, noise):
        regime, s = self._regime(claim)
        refuse = (regime == "cytotoxic")
        return [self._apply(p, s, refuse_update=refuse) for p in xs]


def compound_claim(smiles: str, skeleton: str, quantity=P53_ACTIVE) -> Claim:
    """
    A compound as a claim. The identifier is the InChIKey skeleton actually
    used to join the two datasets, and the structure rides in context
    because the identifier cannot be inverted back to it.
    """
    return Claim(
        entity=Entity("inchikey", skeleton + "-UHFFFAOYSA-N"),
        quantity=quantity,
        context=Context(assay="Tox21 qHTS", time_h=None)
        .with_cov(SMILES_COV, smiles),
        estimate=Estimate.point(0.0, 2000),
        provenance=(Provenance("structure", "input", "", smiles[:60]),),
    )
