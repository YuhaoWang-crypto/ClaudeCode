"""
A first upstream: structural alerts to damage channels.

The seam this fills has been empty because the obvious way to fill it is
circular — train on genotoxicity endpoints, predict genotoxicity endpoints,
and call the intermediate a lesion flux. These rules avoid that by encoding
*chemistry and pharmacology* rather than assay outcomes: a sulfonate ester
alkylates because of what it is, not because a table says it is an Ames
positive. Nothing here was fitted to the benchmark labels.

That honesty has a price, and the benchmark charges it. Alerts describe
**reactivity**, so they cover the electrophiles well. They describe
**target binding** badly, because binding a pocket is not a substructure
property: the benzimidazole carbamate aneugens are a clean chemotype and are
caught, while colchicine — a tubulin binder with no reactive group — is
invisible to every rule here and always will be. That class needs a
target-binding model, not a SMARTS list. The same applies to the
podophyllotoxin topo-II poisons.

Rate magnitudes are illustrative (H), on the same arbitrary scale as
:data:`genotox.damage.DEMO_COMPOUNDS`. What this source is for is testing
whether structures route into the right *channels*; its numbers are not
potencies.

Metabolic activation is **not** predicted. The benchmark supplies it from
reference metadata, so a score here measures channel routing only.
"""
from __future__ import annotations

from dataclasses import dataclass, field

from .damage import CHANNELS, DamageFlux, DamageSource


@dataclass(frozen=True)
class Alert:
    """One rule: a substructure, the channels it implies, and why."""

    name: str
    smarts: str
    channels: dict
    rationale: str


#: Reactivity alerts. Each cites the chemistry, not an endpoint.
ALERTS = (
    Alert("alkyl sulfonate ester", "[CX4][OX2]S(=O)(=O)[#6]",
          {"alkylation": 0.40, "ssb": 0.25},
          "direct S_N2 alkylating agent; the adducts are alkali-labile, and "
          "base excision leaves transient gaps"),
    Alert("nitrogen mustard", "[NX3]([CH2][CH2][Cl,Br])[CH2][CH2][Cl,Br]",
          {"alkylation": 0.35, "icl": 0.20},
          "bifunctional alkylator via aziridinium; the second arm crosslinks"),
    Alert("epoxide", "[CX4]1[OX2][CX4]1",
          {"alkylation": 0.30},
          "strained electrophile, opens on nucleophilic attack by DNA bases"),
    Alert("aziridine", "[CX4]1[NX3][CX4]1",
          {"alkylation": 0.25, "icl": 0.30},
          "strained electrophile; bifunctional agents of this class "
          "crosslink"),
    Alert("aromatic nitro", "[a][NX3+](=O)[O-]",
          {"bulky_adduct": 0.45, "oxidative": 0.15},
          "nitroreduction gives nitrenium/hydroxylamine species that form "
          "bulky adducts, with redox cycling alongside"),
    Alert("aromatic primary amine", "[NX3;H2][c]",
          {"bulky_adduct": 0.30},
          "N-hydroxylation then esterification gives a nitrenium ion — the "
          "classic promutagen route, which is also this rule's weakness: the "
          "group is common in non-mutagenic drugs"),
    Alert("quinone", "O=C1[#6]=,:[#6]C(=O)[#6]=,:[#6]1",
          {"oxidative": 0.35},
          "redox cycling and Michael acceptor chemistry"),
    Alert("benzimidazole carbamate", "[#6]OC(=O)[NX3]c1nc2ccccc2[nH]1",
          {"aneugenic": 0.90},
          "the tubulin-binding chemotype; a target-class rule rather than a "
          "reactivity rule, and legitimate only because this class is "
          "structurally tight"),
    Alert("quinolone-3-carboxylate", "O=c1c(-C(=O)[OX2H1,OX1-])cn([#6])c2ccccc12",
          {"topo_bacterial": 0.50},
          "traps the BACTERIAL gyrase/topo-IV cleavage complex; deliberately "
          "not the mammalian channel"),
    Alert("anthracycline", "O=C1c2ccccc2C(=O)c2c1cccc2",
          {"topo_mammalian": 0.40, "oxidative": 0.25, "dsb": 0.015},
          "intercalating topo-II poison with a redox-active quinone"),
)

#: A fused polycyclic aromatic system needs a ring count, not a substructure
#: match, so it sits outside the SMARTS list.
PAH_MIN_FUSED_RINGS = 4
PAH_CHANNELS = {"bulky_adduct": 0.40}


@dataclass
class StructuralAlertSource(DamageSource):
    """Reads a structure, fires reactivity alerts, emits a channel vector.

    Channels no rule fired are marked **unknown**, not zero. The distinction
    matters here more than anywhere else in the package: a rule set that
    finds no alert has established that *its rules* found nothing, which is
    not the same as establishing the compound is clean. Downstream, a
    NEGATIVE verdict built on this source therefore states its own coverage.
    """

    provenance = ("structural alerts (SMARTS + fused-ring count), "
                  "illustrative magnitudes (H), not fitted to any labels")
    alerts: tuple = ALERTS
    _patterns: dict = field(default_factory=dict, repr=False)

    def __post_init__(self):
        from rdkit import Chem, RDLogger
        RDLogger.DisableLog("rdApp.*")
        for a in self.alerts:
            pat = Chem.MolFromSmarts(a.smarts)
            if pat is None:
                raise ValueError(f"alert {a.name!r} has invalid SMARTS")
            self._patterns[a.name] = pat

    # -- chemistry ---------------------------------------------------------
    def _fused_aromatic_carbocycles(self, mol) -> int:
        rings = [r for r in mol.GetRingInfo().AtomRings()
                 if all(mol.GetAtomWithIdx(i).GetIsAromatic()
                        and mol.GetAtomWithIdx(i).GetSymbol() == "C"
                        for i in r)]
        return sum(1 for r in rings
                   if any(len(set(r) & set(o)) == 2 for o in rings if o != r))

    def explain(self, smiles: str) -> dict:
        """Which rules fired, and what each contributed."""
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"cannot parse SMILES: {smiles!r}")
        fired = []
        for a in self.alerts:
            if mol.HasSubstructMatch(self._patterns[a.name]):
                fired.append((a.name, dict(a.channels)))
        if self._fused_aromatic_carbocycles(mol) >= PAH_MIN_FUSED_RINGS:
            fired.append(("fused polycyclic aromatic", dict(PAH_CHANNELS)))
        return {"smiles": smiles, "fired": fired}

    def flux(self, compound, dose_uM, s9=False) -> DamageFlux:
        """Predict from the structure when the compound carries one.

        Falling back to the tabulated value keeps mixed panels runnable, but
        the fallback is a table lookup and must not be mistaken for a
        prediction -- which is why it is only reached when there is no
        structure to read.
        """
        if not compound.smiles:
            return compound.flux(dose_uM, s9=s9)
        eff = compound.direct_fraction + (compound.s9_fraction if s9 else 0.0)
        return self.flux_from_smiles(compound.smiles).scaled(dose_uM * eff)

    def flux_from_smiles(self, smiles: str) -> DamageFlux:
        totals: dict = {}
        for _, chans in self.explain(smiles)["fired"]:
            for c, v in chans.items():
                totals[c] = totals.get(c, 0.0) + v
        unknown = frozenset(c for c in CHANNELS if c not in totals)
        return DamageFlux(unknown=unknown, **totals)
