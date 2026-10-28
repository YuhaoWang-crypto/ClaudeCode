"""Saturation mutagenesis triage over the TCR paratope.

READ THIS BEFORE USING THE NUMBERS
----------------------------------
The score computed here is a *triage heuristic in arbitrary units*. It is not
a free-energy prediction. The weights below were set by hand from
well-established interface-energetics principles; they were not fitted to any
affinity dataset, and this package ships no calibration set.

What the score is for: shrinking a space of ~20 x 19 = 380+ single mutants
down to a few dozen that are worth spending real compute on. What it is not
for: deciding that one mutant binds more tightly than another. That decision
needs Rosetta flex ddG, FoldX, an MD/FEP protocol, or a co-folding model such
as Boltz-2 run on each candidate, and ultimately a binding measurement.

The structural quantities the score consumes (contact counts, buried surface
area, side-chain void, partner charge and polar-atom counts) are direct
measurements of the input structure and are reliable to the extent the
structure is.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict

import numpy as np
from Bio.Align import substitution_matrices

from .constants import AA20, AROMATIC, CHARGE, HBOND_CAPACITY, HYDROPATHY, VOLUME
from .structure import Complex, Residue

BLOSUM62 = substitution_matrices.load("BLOSUM62")

# Hand-set term weights. Unfitted: see the module docstring.
WEIGHTS = {
    "packing": 1.0,
    "apolar_burial": 1.5,
    "electrostatic": 1.0,
    "hbond": 0.8,
    "conformational": 1.0,
    "substitution_prior": 0.3,
}

# A position whose partner contacts are this heavily skewed towards MHC is
# treated as a cross-reactivity liability rather than an affinity opportunity.
MHC_FRACTION_WATCH = 0.5
MHC_FRACTION_REJECT = 0.8


@dataclass
class Mutation:
    """One scored single-point substitution."""

    chain: str
    number: int
    wild_type: str
    mutant: str
    loop: str | None

    packing: float = 0.0
    apolar_burial: float = 0.0
    electrostatic: float = 0.0
    hbond: float = 0.0
    conformational: float = 0.0
    substitution_prior: float = 0.0
    triage_score: float = 0.0
    z_score: float = 0.0
    peptide_z: float | None = None
    # Filled in when the calibrated model is available: predicted change in
    # binding free energy, kcal/mol, negative meaning tighter binding.
    predicted_ddg: float | None = None

    contacts_peptide: int = 0
    contacts_mhc: int = 0
    mhc_contact_fraction: float = 0.0
    buried_fraction: float = 0.0
    void_slots: int = 0

    specificity_flag: str = "ok"
    liabilities: tuple[str, ...] = ()
    verdict: str = "pending"
    notes: tuple[str, ...] = ()

    @property
    def name(self) -> str:
        return f"{self.chain}:{self.wild_type}{self.number}{self.mutant}"

    def as_row(self) -> dict:
        row = asdict(self)
        row["name"] = self.name
        row["liabilities"] = "; ".join(self.liabilities)
        row["notes"] = "; ".join(self.notes)
        return row


def ranking_key(mutation: Mutation):
    """Sort key putting the most promising candidate first.

    Predicted ddG wins when the calibrated model is available, because it is
    the only quantity here with a measured relationship to binding. Otherwise
    fall back to the heuristic's z-score.
    """
    if mutation.predicted_ddg is not None:
        return mutation.predicted_ddg  # lower is tighter binding
    z = mutation.peptide_z if mutation.peptide_z is not None else mutation.z_score
    return -z


def _packing_term(residue: Residue, wt: str, mut: str) -> tuple[float, list[str]]:
    """Reward side-chain growth only where there is room and a partner to touch.

    ``void_slots`` counts how many of three probe points along the side-chain
    direction are free of neighbouring atoms, so it measures real space in the
    deposited structure. ``growth_reaches_partner`` distinguishes growth into
    the interface from growth into bulk solvent.
    """
    notes: list[str] = []
    delta_volume = (VOLUME[mut] - VOLUME[wt]) / 100.0
    room = residue.void_slots / 3.0
    reach = 1.0 if residue.growth_reaches_partner else 0.25

    if delta_volume > 0:
        # 2*room-1 flips sign when the position is crowded, turning growth
        # into a steric penalty.
        value = delta_volume * (2.0 * room - 1.0) * reach
        if room < 0.34:
            notes.append("no room to enlarge side chain (steric clash likely)")
        elif not residue.growth_reaches_partner:
            notes.append("added volume points at solvent, not the interface")
    else:
        # Shrinking a residue at a packed interface costs contacts.
        value = delta_volume * 0.5 * (1.0 if residue.contacts_partner else 0.2)
        if residue.contacts_partner > 15:
            notes.append("removes volume from a well-packed contact")

    if mut in AROMATIC and wt not in AROMATIC and room >= 0.67 and residue.growth_reaches_partner:
        value += 0.3
        notes.append("aromatic insertion into available interface volume")
    return value, notes


def _apolar_term(residue: Residue, wt: str, mut: str) -> tuple[float, list[str]]:
    """Apolar burial, the dominant driver of real interface affinity gains."""
    notes: list[str] = []
    delta = (HYDROPATHY[mut] - HYDROPATHY[wt]) / 9.0
    coupling = residue.buried_fraction
    value = delta * coupling
    if delta > 0 and coupling > 0.6:
        notes.append("adds apolar surface at a well-buried position")
    if delta > 0 and coupling < 0.25:
        notes.append("exposed position: new apolar surface stays solvated")
        value *= 0.3
    return value, notes


def _electrostatic_term(residue: Residue, wt: str, mut: str) -> tuple[float, list[str]]:
    """Charge complementarity against nearby partner side chains."""
    notes: list[str] = []
    q_mut, q_wt = CHARGE[mut], CHARGE[wt]
    partner = residue.partner_charge
    value = (q_mut - q_wt) * (-partner)

    if abs(q_mut) > 0.5 and abs(partner) < 0.1 and residue.buried_fraction > 0.6:
        value -= 0.5 * abs(q_mut)
        notes.append("buries a charge with no counter-charge nearby (desolvation cost)")
    if q_mut * partner > 0.3:
        notes.append("like-charge repulsion against the partner surface")
    elif q_mut * partner < -0.3:
        notes.append("complementary charge pairing with the partner surface")
    return value, notes


def _hbond_term(residue: Residue, wt: str, mut: str) -> tuple[float, list[str]]:
    notes: list[str] = []
    delta = HBOND_CAPACITY[mut] - HBOND_CAPACITY[wt]
    available = min(residue.partner_polar_atoms, 6) / 6.0
    value = 0.5 * delta * available
    if delta > 0 and available > 0.5:
        notes.append("new polar groups face partner donors/acceptors")
    return value, notes


def _conformational_term(residue: Residue, wt: str, mut: str) -> tuple[float, list[str]]:
    """Backbone-compatibility penalties for the special residues."""
    notes: list[str] = []
    value = 0.0

    if mut == "P":
        value -= 1.5
        notes.append("proline substitution restricts the loop backbone")
        if residue.backbone_nh_bonded:
            value -= 2.0
            notes.append("backbone NH is hydrogen bonded: proline would break it")
    if wt == "P":
        value -= 1.0
        notes.append("removing proline may change loop conformation")
    if mut == "G":
        value -= 0.5
        notes.append("glycine adds backbone entropy to the unbound loop")
    if wt == "G" and residue.phi is not None and residue.phi > 0:
        value -= 2.5
        notes.append(
            f"wild-type glycine sits in a left-handed backbone conformation "
            f"(phi = {residue.phi:+.0f} deg); non-glycine is strained there"
        )
    if mut == "C":
        value -= 1.0
        notes.append("free cysteine: mispairing and manufacturing risk")
    return value, notes


def _prior_term(wt: str, mut: str) -> float:
    """A mild preference for exchanges that occur in natural repertoires.

    Harshly scored substitutions are likelier to disturb folding, pairing or
    expression, none of which the structural terms can see.
    """
    try:
        score = float(BLOSUM62[wt, mut])
    except KeyError:
        score = 0.0
    return score / 4.0


def score_mutation(residue: Residue, mutant: str) -> Mutation:
    """Score a single substitution at a paratope position."""
    wt = residue.aa
    mutation = Mutation(
        chain=residue.chain,
        number=residue.number,
        wild_type=wt,
        mutant=mutant,
        loop=residue.loop,
        contacts_peptide=residue.contacts_peptide,
        contacts_mhc=residue.contacts_mhc,
        mhc_contact_fraction=residue.mhc_contact_fraction,
        buried_fraction=residue.buried_fraction,
        void_slots=residue.void_slots,
    )

    notes: list[str] = []
    mutation.packing, n = _packing_term(residue, wt, mutant)
    notes += n
    mutation.apolar_burial, n = _apolar_term(residue, wt, mutant)
    notes += n
    mutation.electrostatic, n = _electrostatic_term(residue, wt, mutant)
    notes += n
    mutation.hbond, n = _hbond_term(residue, wt, mutant)
    notes += n
    mutation.conformational, n = _conformational_term(residue, wt, mutant)
    notes += n
    mutation.substitution_prior = _prior_term(wt, mutant)

    mutation.triage_score = (
        WEIGHTS["packing"] * mutation.packing
        + WEIGHTS["apolar_burial"] * mutation.apolar_burial
        + WEIGHTS["electrostatic"] * mutation.electrostatic
        + WEIGHTS["hbond"] * mutation.hbond
        + WEIGHTS["conformational"] * mutation.conformational
        + WEIGHTS["substitution_prior"] * mutation.substitution_prior
    )
    mutation.notes = tuple(notes)

    # The specificity gate is deliberately kept out of the score so that it
    # cannot be averaged away by a large predicted affinity gain.
    fraction = residue.mhc_contact_fraction
    if residue.contacts_partner == 0:
        mutation.specificity_flag = "no-contact"
    elif fraction >= MHC_FRACTION_REJECT:
        mutation.specificity_flag = "mhc-dominated"
    elif fraction >= MHC_FRACTION_WATCH:
        mutation.specificity_flag = "mhc-leaning"
    else:
        mutation.specificity_flag = "peptide-focused"
    return mutation


def scan_complex(
    complex_obj: Complex,
    min_contacts: int = 3,
    loops: tuple[str, ...] | None = None,
    include_wild_type: bool = False,
    ddg_model=None,
) -> list[Mutation]:
    """Run saturation mutagenesis over the contacting CDR positions.

    Positions are restricted to the structurally defined paratope. Cysteines in
    a disulfide are skipped outright, since mutating them unfolds the domain.
    """
    mutations: list[Mutation] = []
    for residue in complex_obj.paratope(min_contacts=min_contacts):
        if loops and (residue.loop or "") not in loops:
            continue
        if residue.in_disulfide:
            continue
        for mutant in AA20:
            if mutant == residue.aa and not include_wild_type:
                continue
            mutation = score_mutation(residue, mutant)
            if ddg_model is not None:
                mutation.predicted_ddg = ddg_model.predict(residue, mutation)
            mutations.append(mutation)

    # Report each score relative to the scan as a whole. A z-score is the
    # honest way to read an uncalibrated function: it ranks, it does not
    # pretend to a physical unit.
    _add_z_scores(mutations, "z_score", mutations)

    # Rank the peptide-reading positions against each other as well. Across a
    # whole TCR-pMHC interface most positions contact MHC, so a global ranking
    # is dominated by substitutions that the specificity gate will reject
    # anyway. Scoring the positions that actually touch the peptide as their
    # own subset is what makes the surviving candidates comparable.
    peptide_readers = [m for m in mutations if m.contacts_peptide > 0]
    _add_z_scores(peptide_readers, "peptide_z", peptide_readers)
    return mutations


def _add_z_scores(targets: list[Mutation], field: str, reference: list[Mutation]) -> None:
    if not targets or not reference:
        return
    values = np.array([m.triage_score for m in reference])
    mean, sd = float(values.mean()), float(values.std())
    for mutation in targets:
        setattr(mutation, field, (mutation.triage_score - mean) / sd if sd > 1e-9 else 0.0)
