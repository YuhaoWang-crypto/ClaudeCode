"""
A scoring target for the upstream layer that does not exist yet.

``QSARSource`` is a declared seam with a deliberate ``NotImplementedError``.
A seam with no target is an intention; a seam with a benchmark is a task. This
module supplies the benchmark: reference compounds with structures, the
endpoint outcomes they are known for, and a harness that scores any
:class:`~genotox.damage.DamageSource` against them.

What the labels are
-------------------
The ``expected`` field on each compound is the **published reference
classification** for that chemical — the consensus that makes it a standard
positive or negative control in guideline testing — not a measurement made
here and not an output of this package. Potencies are deliberately absent:
this package's rate constants are illustrative, so scoring it on EC values
would be scoring noise. What is scored is the *call*: positive, negative, and
for micronuclei whether the mechanism is clastogenic or aneugenic.

What passing would and would not mean
-------------------------------------
A source that reproduces these calls has shown it can route known chemistry
into the right channels. It has **not** shown that its channel magnitudes are
real lesion fluxes — these compounds are famous precisely because their
outcomes are widely known, so any model trained on public genotoxicity data
has almost certainly seen them. Treat a good score here as a floor the
upstream must clear, never as validation. The first genuinely informative
test is a compound set the source has not seen, scored before the answers
are looked up.

Why these compounds
-------------------
Each row is here to exercise something the architecture claims:

* direct-acting versus S9-dependent, which the activation layer must split;
* alkylating versus bulky versus strand-breaking, which route to different
  channels and give different comet behaviour;
* aneugens, which must reach micronuclei without touching any break pool;
* the "misleading positives" — chemicals positive in mammalian cell assays
  and not regarded as genotoxic carcinogens — which is what the validity
  gates and the cytotoxicity machinery exist to catch;
* clear negatives, without which a source that says "positive" to everything
  would score well.
"""
from __future__ import annotations

from dataclasses import dataclass, field

CHANNEL_HINT_NOTE = (
    "channels a correct upstream should populate — a hypothesis about "
    "mechanism, not a measured flux"
)


@dataclass(frozen=True)
class Reference:
    """One benchmark compound.

    ``expected`` keys are endpoint names; values are ``"+"``, ``"-"`` or
    ``"?"`` where the published position is genuinely mixed. ``mn_mechanism``
    is ``"clastogenic"``, ``"aneugenic"`` or ``None``.
    """

    name: str
    smiles: str
    formula: str
    role: str
    expected: dict
    channels: tuple            # channel names a correct source should fire
    needs_s9: bool = False
    mn_mechanism: str | None = None
    #: the in-vitro positive for this compound is attributed to cytotoxicity
    #: rather than to a genotoxic mechanism.  Scored separately: a source
    #: that predicts damage channels SHOULD decline to fire here, so counting
    #: the in-vitro "+" as a miss penalises it for being right.
    artifact: bool = False
    note: str = ""


#: Reference set. Structures verified against molecular formula with RDKit;
#: classifications are the published consensus for these control chemicals.
REFERENCES = (
    # ---- direct-acting bacterial mutagens -------------------------------
    Reference(
        "4-nitroquinoline 1-oxide", "O=[N+]([O-])c1cc[n+]([O-])c2ccccc21",
        "C9H6N2O3", "direct-acting positive control",
        {"ames": "+", "umu": "+", "comet": "+", "mn": "+"},
        ("bulky_adduct", "oxidative"), mn_mechanism="clastogenic",
        note="the standard -S9 positive control for the umu test",
    ),
    Reference(
        "methyl methanesulfonate", "COS(C)(=O)=O", "C2H6O3S",
        "direct alkylating agent",
        {"ames": "+", "umu": "+", "comet": "+", "mn": "+"},
        ("alkylation", "ssb"), mn_mechanism="clastogenic",
        note="adducts are alkali-labile, so comet-positive without waiting "
             "for excision",
    ),
    Reference(
        "ethyl methanesulfonate", "CCOS(C)(=O)=O", "C3H8O3S",
        "direct alkylating agent",
        {"ames": "+", "umu": "+", "comet": "+", "mn": "+"},
        ("alkylation", "ssb"), mn_mechanism="clastogenic",
    ),
    Reference(
        "sodium azide", "[N-]=[N+]=[N-].[Na+]", "N3Na",
        "strain-specific positive control",
        {"ames": "+", "umu": "?", "comet": "-", "mn": "-"},
        ("alkylation",),
        note="the classic TA1535 Ames positive — same strain background this "
             "package's SOS core uses — but it acts through a metabolite and "
             "is not a general mammalian clastogen. A source that fires every "
             "channel for it is wrong.",
    ),

    # ---- promutagens: the activation layer decides ----------------------
    Reference(
        "2-aminoanthracene", "Nc1ccc2cc3ccccc3cc2c1", "C14H11N",
        "promutagen", {"ames": "+", "umu": "+", "comet": "+", "mn": "+"},
        ("bulky_adduct",), needs_s9=True, mn_mechanism="clastogenic",
        note="the standard +S9 positive control for the umu test; near-silent "
             "without activation",
    ),
    Reference(
        "benzo[a]pyrene", "c1ccc2c(c1)cc1ccc3cccc4ccc2c1c34", "C20H12",
        "promutagen", {"ames": "+", "umu": "+", "comet": "+", "mn": "+"},
        ("bulky_adduct",), needs_s9=True, mn_mechanism="clastogenic",
    ),
    Reference(
        "cyclophosphamide", "ClCCN(CCCl)P1(=O)NCCCO1", "C7H15Cl2N2O2P",
        "promutagen, clinical alkylating prodrug",
        {"ames": "+", "umu": "+", "comet": "+", "mn": "+"},
        ("alkylation", "icl"), needs_s9=True, mn_mechanism="clastogenic",
        note="the standard +S9 clastogen control in the micronucleus test",
    ),

    # ---- direct clastogens ----------------------------------------------
    Reference(
        "mitomycin C",
        "CC1=C(N)C(=O)C2=C(C1=O)N1CC3NC3C1(COC(N)=O)C2OC", "C15H18N4O5",
        "crosslinking clastogen",
        {"ames": "+", "umu": "+", "comet": "+", "mn": "+"},
        ("icl", "dsb"), mn_mechanism="clastogenic",
        note="the standard -S9 clastogen control in the micronucleus test",
    ),
    Reference(
        "etoposide",
        "COc1cc(cc(OC)c1O)C1c2cc3OCOc3cc2C(OC2OC3COC(C)OC3C(O)C2O)"
        "C2COC(=O)C12", "C29H32O13",
        "topoisomerase II poison, clinical",
        {"ames": "-", "umu": "?", "comet": "+", "mn": "+"},
        ("topo_mammalian", "dsb"), mn_mechanism="clastogenic",
        note="a strong mammalian clastogen that is largely Ames-negative — "
             "the case that a bacterial endpoint alone would miss",
    ),
    Reference(
        "doxorubicin",
        "COc1cccc2C(=O)c3c(O)c4CC(O)(CC(OC5CC(N)C(O)C(C)O5)c4c(O)c3C(=O)c12)"
        "C(=O)CO", "C27H29NO11",
        "anthracycline, topoisomerase II poison, clinical",
        {"ames": "+", "umu": "+", "comet": "+", "mn": "+"},
        ("topo_mammalian", "dsb", "oxidative"), mn_mechanism="clastogenic",
    ),
    Reference(
        "ciprofloxacin", "OC(=O)C1=CN(C2CC2)c2cc(N3CCNCC3)c(F)cc2C1=O",
        "C17H18FN3O3", "bacterial gyrase poison, clinical",
        {"ames": "-", "umu": "+", "comet": "-", "mn": "-"},
        ("topo_bacterial",),
        note="poisons the bacterial enzyme, so it induces SOS while being "
             "unremarkable in mammalian cells — the sharpest test that a "
             "source must not treat 'genotoxic' as species-independent",
    ),

    # ---- aneugens: micronuclei without breaks ----------------------------
    Reference(
        "colchicine",
        "COc1cc2CCC(NC(C)=O)c3cc(=O)c(OC)ccc3-c2c(OC)c1OC", "C22H25NO6",
        "aneugen, clinical",
        {"ames": "-", "umu": "-", "comet": "-", "mn": "+"},
        ("aneugenic",), mn_mechanism="aneugenic",
        note="tubulin binder; micronuclei are centromere-positive",
    ),
    Reference(
        "carbendazim", "COC(=O)Nc1nc2ccccc2[nH]1", "C9H9N3O2",
        "aneugen", {"ames": "-", "umu": "-", "comet": "-", "mn": "+"},
        ("aneugenic",), mn_mechanism="aneugenic",
        note="benzimidazole carbamate; the structural class is a clean "
             "SMARTS target for an aneugenic channel",
    ),
    Reference(
        "nocodazole", "COC(=O)Nc1nc2cc(ccc2[nH]1)C(=O)c1cccs1", "C14H11N3O3S",
        "aneugen", {"ames": "-", "umu": "-", "comet": "-", "mn": "+"},
        ("aneugenic",), mn_mechanism="aneugenic",
    ),

    # ---- misleading positives: what the gates exist for ------------------
    Reference(
        "eugenol", "C=CCc1ccc(O)c(OC)c1", "C10H12O2",
        "misleading positive",
        {"ames": "-", "umu": "-", "comet": "?", "mn": "+"},
        (), mn_mechanism=None, artifact=True,
        note="positive in mammalian cell assays at cytotoxic concentrations; "
             "not regarded as a genotoxic carcinogen. A source should give it "
             "little or nothing and let the validity gate explain the assay "
             "result.",
    ),
    Reference(
        "2-biphenylol", "Oc1ccccc1-c1ccccc1", "C12H10O",
        "misleading positive",
        {"ames": "-", "umu": "-", "comet": "?", "mn": "+"},
        (), artifact=True,
        note="same class: positive in vitro, negative in vivo",
    ),
    Reference(
        "saccharin", "O=C1NS(=O)(=O)c2ccccc21", "C7H5NO3S",
        "misleading positive",
        {"ames": "-", "umu": "-", "comet": "-", "mn": "?"},
        (), artifact=True,
        note="high-dose in-vitro effects with no genotoxic mechanism",
    ),

    # ---- negatives -------------------------------------------------------
    Reference(
        "D-mannitol", "OC[C@@H](O)[C@@H](O)[C@H](O)[C@H](O)CO", "C6H14O6",
        "negative control",
        {"ames": "-", "umu": "-", "comet": "-", "mn": "-"}, (),
    ),
    Reference(
        "caffeine", "Cn1cnc2c1c(=O)n(C)c(=O)n2C", "C8H10N4O2",
        "negative at relevant doses",
        {"ames": "-", "umu": "-", "comet": "-", "mn": "-"}, (),
        note="modifies the response to other agents at high concentrations, "
             "but is not itself a genotoxicant at testable doses",
    ),
    Reference(
        "methotrexate",
        "CN(Cc1cnc2nc(N)nc(N)c2n1)c1ccc(cc1)C(=O)NC(CCC(O)=O)C(O)=O",
        "C20H22N8O5", "cytotoxic, non-DNA-reactive, clinical",
        {"ames": "-", "umu": "-", "comet": "-", "mn": "?"}, (), artifact=True,
        note="an antifolate: strongly antiproliferative with no direct DNA "
             "reactivity, so it probes whether a source confuses cytotoxicity "
             "with genotoxicity",
    ),
)

BY_NAME = {r.name: r for r in REFERENCES}


def validate() -> dict:
    """Parse every structure and check it against its stated formula.

    A benchmark whose structures are wrong is worse than no benchmark, so
    this runs before anything is scored against it.
    """
    from rdkit import Chem, RDLogger
    from rdkit.Chem import rdMolDescriptors as rmd
    RDLogger.DisableLog("rdApp.*")

    rows, bad = [], []
    for r in REFERENCES:
        mol = Chem.MolFromSmiles(r.smiles)
        if mol is None:
            bad.append((r.name, "does not parse"))
            continue
        got = rmd.CalcMolFormula(mol).replace("+", "").replace("-", "")
        if got != r.formula:
            bad.append((r.name, f"formula {got} != {r.formula}"))
        rows.append({"name": r.name, "formula": got,
                     "canonical": Chem.MolToSmiles(mol)})
    return {"rows": rows, "problems": bad, "n": len(REFERENCES)}


# --------------------------------------------------------------------------
@dataclass
class Score:
    endpoint: str
    tp: int = 0
    fp: int = 0
    tn: int = 0
    fn: int = 0
    skipped: int = 0
    #: in-vitro positives attributed to cytotoxicity, scored separately
    artifact_total: int = 0
    artifact_declined: int = 0
    errors: list = field(default_factory=list)

    @property
    def n(self) -> int:
        return self.tp + self.fp + self.tn + self.fn

    @property
    def accuracy(self) -> float:
        return (self.tp + self.tn) / self.n if self.n else float("nan")

    @property
    def sensitivity(self) -> float:
        d = self.tp + self.fn
        return self.tp / d if d else float("nan")

    @property
    def specificity(self) -> float:
        d = self.tn + self.fp
        return self.tn / d if d else float("nan")


def score_source(source, build_assay, protocol, endpoint: str,
                 doses, s9: bool | None = None) -> Score:
    """Run one endpoint over the reference set and compare calls to labels.

    ``build_assay(source)`` must return a configured
    :class:`~genotox.assay.VirtualAssay`. Compounds whose published position
    is genuinely mixed (``"?"``) are skipped rather than counted either way —
    scoring against an ambiguous label manufactures accuracy.

    INCONCLUSIVE is not folded into either class: it is counted as a miss
    against a positive label and as a correct abstention against a negative
    one, because an assay that declines to decide has not produced a false
    positive.
    """
    from .damage import Compound, DamageFlux
    from .doseresponse import call_result, dose_series

    sc = Score(endpoint)
    for ref in REFERENCES:
        want = ref.expected.get(endpoint, "?")
        if want == "?":
            sc.skipped += 1
            continue
        if ref.artifact and want == "+":
            # An in-vitro positive attributed to cytotoxicity.  A
            # damage-channel source that declines to fire is behaving
            # correctly, so this is counted in `artifact_declined`, not as a
            # false negative.  Reproducing the artifact would need a
            # cytotoxicity prediction, which no source here makes.
            sc.artifact_total += 1
            try:
                fired = source.flux_from_smiles(ref.smiles).nonzero()
            except NotImplementedError:
                fired = {}
            sc.artifact_declined += not fired
            continue
        try:
            source.flux_from_smiles(ref.smiles)   # fail fast on an empty seam
        except NotImplementedError as e:
            sc.errors.append((ref.name, str(e).split(".")[0]))
            sc.skipped += 1
            continue
        comp = Compound(name=ref.name, per_uM=DamageFlux(),
                        smiles=ref.smiles,
                        direct_fraction=0.02 if ref.needs_s9 else 1.0,
                        s9_fraction=0.95 if ref.needs_s9 else 0.0)
        use_s9 = ref.needs_s9 if s9 is None else s9
        res = call_result(dose_series(build_assay(source), comp, doses,
                                      s9=use_s9, protocol=protocol))
        got = res["verdict"]
        if want == "+":
            sc.tp += got == "POSITIVE"
            sc.fn += got != "POSITIVE"
        else:
            sc.tn += got != "POSITIVE"
            sc.fp += got == "POSITIVE"
    return sc
