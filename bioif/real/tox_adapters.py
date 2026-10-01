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


class P53ToMutagenicity(Adapter):
    """
    F->G : the blueprint's `blocked` edge, graded by a measured 2x2.

    P(Ames+ | p53 state) is read straight off compounds assayed on both
    sides. No model sits in between, which is why the kind is ASSOCIATION
    and the evidence ceiling is `inferred_association`.
    """

    name = "F->G p53->mutagenicity"
    version = "tox21xames-2x2"
    kind = ASSOCIATION
    consumes = P53_ACTIVE
    produces = AMES_POSITIVE
    assumption = ("a reporter positive is treated as evidence about "
                  "mutagenicity via measured co-occurrence only; no "
                  "mechanism is claimed, and Ames is not the Comet or "
                  "adduct endpoint")

    def __init__(self):
        self.t = tox.two_by_two(tox.P53)
        self._p_pos = self.t.ppv
        self._p_neg = self.t.base

    def calibrated_width(self, alpha: float = 0.2) -> float:
        """How far apart the two conditional risks are: a wider gap is a
        more informative edge, so the registry prefers it."""
        return 1.0 - abs(self._p_pos - self._p_neg)

    def screening_verdict(self) -> str:
        t = self.t
        return (f"RR {t.risk_ratio:.2f}, OR {t.odds_ratio:.2f} "
                f"(Fisher p={t.fisher_p():.1e}) on n={t.n} -- but sensitivity "
                f"{t.sensitivity:.2f}: this edge ranks candidates, it does "
                f"not clear them")

    def map_context(self, c: Context, e: Entity) -> Context:
        return c.with_cov("bridge", "tox21xames 2x2")

    def _forward(self, xs, claim, rng, noise):
        # total probability over the (uncertain) reporter state
        return [p * self._p_pos + (1.0 - p) * self._p_neg for p in xs]


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
