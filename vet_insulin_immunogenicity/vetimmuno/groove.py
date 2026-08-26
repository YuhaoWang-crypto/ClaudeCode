"""Binding-groove pseudosequences and the applicability-domain test.

Why this module exists
----------------------
NetMHCIIpan is a *pan-specific* model: it represents an MHC-II molecule by the
residues that line the peptide-binding groove, so it can score a molecule it has
never seen -- as long as that molecule's groove looks like something in its
training data. NetMHCIIpan-4.3 was trained on human HLA-DR/DQ/DP, mouse H-2 and
bovine BoLA-DRB3. Canine DLA and feline FLA are not in there.

So before trusting a single prediction for a dog or a cat, we measure *how far
outside the training space* each DLA/FLA molecule actually sits, using the same
representation the model itself uses: the groove-contact residues.

The threshold is not invented. We compute the leave-one-out nearest-neighbour
identity *within* the training set and ask where each DLA/FLA molecule falls in
that distribution. A molecule whose nearest training neighbour is further away
than essentially every training molecule's own nearest neighbour is, by
construction, extrapolation.

Honesty label
-------------
The contact-position list below is derived from published HLA-DR1/peptide
crystal structures (Brown et al. 1993; Stern et al. 1994) and is an
*approximation* of NetMHCIIpan's internal 34-residue pseudosequence, whose exact
position list ships inside the tool. Distances computed here are therefore
indicative of, not identical to, the model's own notion of allele similarity.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from Bio import Align
from Bio.Align import substitution_matrices

from . import data

# --- groove geometry ------------------------------------------------------
# Numbering is the mature class II domain (beta1 or alpha1), IMGT/Kabat style,
# which is shared across DRB / DQB and across species by homology.

BETA_CONTACT: Tuple[int, ...] = (
    9, 11, 13, 26, 28, 30, 37, 38, 47, 57, 60, 61, 64, 67, 70, 71, 74, 77, 78,
    81, 82, 85, 86, 89, 90,
)
ALPHA_CONTACT: Tuple[int, ...] = (
    9, 11, 22, 24, 31, 43, 51, 52, 53, 54, 55, 58, 61, 62, 65, 66, 68, 69, 72, 76,
)

# Pocket -> lining residues, used by the illustrative scorer in predict.py.
BETA_POCKETS: Dict[str, Tuple[int, ...]] = {
    "P1": (85, 86, 89, 90),
    "P4": (13, 26, 70, 71, 74, 78),
    "P6": (11, 30, 71),
    "P7": (28, 47, 61, 67, 71),
    "P9": (9, 37, 57, 60, 61),
}
ALPHA_POCKETS: Dict[str, Tuple[int, ...]] = {
    "P1": (31, 32, 43, 52, 53, 54, 55),
    "P4": (),
    "P6": (11, 62, 65, 66),
    "P7": (),
    "P9": (68, 69, 72, 76),
}

# --- reference frames -----------------------------------------------------
# Signal-peptide lengths and domain extents for the human loci used as the
# numbering frame. Verified against the IPD-IMGT/HLA records at load time.
LOCUS_REF = {
    "DRB": {"allele": "DRB1*01:01:01:01", "leader": 29, "domain": 95,
            "chain": "beta", "imgt": "DRB"},
    "DQB": {"allele": "DQB1*05:01:01:01", "leader": 32, "domain": 94,
            "chain": "beta", "imgt": "DQB1"},
    "DQA": {"allele": "DQA1*01:01:01:01", "leader": 23, "domain": 84,
            "chain": "alpha", "imgt": "DQA1"},
}

# Mature beta1/alpha1 domains expected after slicing; used purely as a tripwire
# so a silent change in an upstream release cannot shift the numbering frame.
REF_PREFIX = {"DRB": "GDTRPRFLWQLKFECHFFNGTERVR",
              "DQB": "RDSPEDFVYQFKGLCYFTNGTERVR",
              "DQA": "EDIVADHVASCGVNLYQFYGPSGQY"}

_MATRIX = substitution_matrices.load("BLOSUM62")


def _aligner() -> Align.PairwiseAligner:
    al = Align.PairwiseAligner()
    al.substitution_matrix = _MATRIX
    al.open_gap_score = -11
    al.extend_gap_score = -1
    al.mode = "global"
    # Free end gaps: candidates are often exon-2 fragments or carry a leader.
    try:                              # Biopython >= 1.85
        al.end_insertion_score = 0.0
        al.end_deletion_score = 0.0
    except AttributeError:            # older releases
        al.target_end_gap_score = 0.0
        al.query_end_gap_score = 0.0
    return al


@dataclass
class Molecule:
    """One MHC-II chain reduced to its groove-contact residues."""

    name: str
    species: str
    locus: str                 # DRB / DQB / DQA
    chain: str                 # beta / alpha
    sequence: str              # the raw record as fetched
    pseudoseq: str = ""        # contact residues, '-' where not covered
    coverage: float = 0.0      # fraction of contact positions resolved
    identity_to_ref: float = 0.0
    aligned: Dict[int, str] = field(default_factory=dict)
    source: str = ""

    @property
    def key(self) -> str:
        return f"{self.locus}:{self.pseudoseq}"


# --------------------------------------------------------------------------
# reference construction
# --------------------------------------------------------------------------

def reference_domain(locus: str) -> str:
    """Mature beta1/alpha1 domain of the human reference allele for ``locus``."""
    spec = LOCUS_REF[locus]
    table = data.imgt_hla(spec["imgt"])
    full = table.get(spec["allele"])
    if full is None:
        raise RuntimeError(f"reference allele {spec['allele']} missing from IPD-IMGT/HLA")
    mature = full[spec["leader"]:spec["leader"] + spec["domain"]]
    if not mature.startswith(REF_PREFIX[locus]):
        raise RuntimeError(
            f"{locus}: numbering frame moved -- {spec['allele']} mature domain starts "
            f"{mature[:25]!r}, expected {REF_PREFIX[locus]!r}"
        )
    return mature


def contact_positions(locus: str) -> Tuple[int, ...]:
    return BETA_CONTACT if LOCUS_REF[locus]["chain"] == "beta" else ALPHA_CONTACT


# --------------------------------------------------------------------------
# alignment -> pseudosequence
# --------------------------------------------------------------------------

def map_to_reference(seq: str, ref: str) -> Tuple[Dict[int, str], float]:
    """Map ``seq`` onto reference numbering; return ``{ref_pos: aa}`` and identity."""
    aln = _aligner().align(ref, seq)[0]
    mapping: Dict[int, str] = {}
    same = total = 0
    ref_blocks, qry_blocks = aln.aligned
    for (rs, re_), (qs, qe) in zip(ref_blocks, qry_blocks):
        for off in range(re_ - rs):
            rpos = rs + off + 1          # 1-based reference numbering
            aa = seq[qs + off]
            mapping[rpos] = aa
            total += 1
            same += aa == ref[rs + off]
    return mapping, (same / total if total else 0.0)


def build_molecule(name: str, species: str, locus: str, seq: str, source: str = "") -> Molecule:
    ref = reference_domain(locus)
    mapping, ident = map_to_reference(seq, ref)
    positions = contact_positions(locus)
    pseudo = "".join(mapping.get(p, "-") for p in positions)
    coverage = sum(c != "-" for c in pseudo) / len(positions)
    return Molecule(
        name=name, species=species, locus=locus, chain=LOCUS_REF[locus]["chain"],
        sequence=seq, pseudoseq=pseudo, coverage=coverage,
        identity_to_ref=ident, aligned=mapping, source=source,
    )


# --------------------------------------------------------------------------
# similarity + applicability domain
# --------------------------------------------------------------------------

def pseudo_identity(a: str, b: str) -> float:
    pairs = [(x, y) for x, y in zip(a, b) if x != "-" and y != "-"]
    if not pairs:
        return 0.0
    return sum(x == y for x, y in pairs) / len(pairs)


def pseudo_blosum(a: str, b: str) -> float:
    """BLOSUM62 similarity over contact positions, normalised to self-score = 1."""
    pairs = [(x, y) for x, y in zip(a, b) if x != "-" and y != "-"]
    if not pairs:
        return 0.0
    score = sum(float(_MATRIX[x, y]) for x, y in pairs)
    self_score = sum(float(_MATRIX[x, x]) for x, _ in pairs)
    return score / self_score if self_score else 0.0


@dataclass
class DomainCall:
    molecule: str
    locus: str
    nearest: str
    identity: float
    blosum: float
    percentile_vs_training: float   # where this NN identity sits in the training LOO distribution
    verdict: str


def _encode(pool: Sequence[Tuple[str, str]]) -> "np.ndarray":
    return np.frombuffer("".join(ps for _, ps in pool).encode(), dtype="S1").reshape(
        len(pool), -1)


def nearest_neighbour(query: str, pool: Sequence[Tuple[str, str]], exclude: Optional[str] = None
                      ) -> Tuple[str, float, float]:
    """Closest training molecule to ``query`` by contact-position identity."""
    mat = _pool_matrix(pool)
    q = np.frombuffer(query.encode(), dtype="S1")
    valid = (mat != b"-") & (q != b"-")
    matches = ((mat == q) & valid).sum(axis=1)
    denom = valid.sum(axis=1)
    ident = np.divide(matches, denom, out=np.zeros(len(pool)), where=denom > 0)
    if exclude is not None:
        for i, (name, _) in enumerate(pool):
            if name == exclude:
                ident[i] = -1.0
    best = int(ident.argmax())
    return pool[best][0], float(ident[best]), pseudo_blosum(query, pool[best][1])


_POOL_CACHE: Dict[int, "np.ndarray"] = {}


def _pool_matrix(pool: Sequence[Tuple[str, str]]) -> "np.ndarray":
    key = id(pool)
    mat = _POOL_CACHE.get(key)
    if mat is None or mat.shape[0] != len(pool):
        mat = _encode(pool)
        _POOL_CACHE[key] = mat
    return mat


def loo_identities(pool: Sequence[Tuple[str, str]]) -> List[float]:
    """Leave-one-out nearest-neighbour identity for every member of ``pool``.

    This is the calibration curve: how close a molecule *inside* the training
    space typically is to its nearest neighbour. Everything the applicability
    test says about DLA/FLA is relative to this distribution.
    """
    mat = _pool_matrix(pool)
    n = mat.shape[0]
    out = np.zeros(n)
    valid_self = mat != b"-"
    for i in range(n):
        q = mat[i]
        valid = valid_self & (q != b"-")
        matches = ((mat == q) & valid).sum(axis=1)
        denom = valid.sum(axis=1)
        ident = np.divide(matches, denom, out=np.zeros(n), where=denom > 0)
        ident[i] = -1.0
        out[i] = ident.max()
    return out.tolist()


def _percentile(value: float, sample: Sequence[float]) -> float:
    if not sample:
        return float("nan")
    return 100.0 * sum(1 for s in sample if s <= value) / len(sample)


def applicability(molecules: Sequence[Molecule],
                  training: Sequence[Tuple[str, str]],
                  loo: Sequence[float]) -> List[DomainCall]:
    """Rate every query molecule against the model's training space."""
    calls: List[DomainCall] = []
    p05 = sorted(loo)[max(0, int(0.05 * len(loo)) - 1)] if loo else 1.0
    for mol in molecules:
        name, ident, blos = nearest_neighbour(mol.pseudoseq, training)
        pct = _percentile(ident, loo)
        if ident >= 1.0:
            verdict = "in-domain (exact groove match in training set)"
        elif ident >= p05:
            verdict = "near-domain (as close as a typical training molecule)"
        elif pct >= 1.0:
            verdict = "marginal (closer than the 1st percentile of training molecules)"
        else:
            verdict = "OUT-OF-DOMAIN (further from training data than any training molecule)"
        calls.append(DomainCall(mol.name, mol.locus, name, ident, blos, pct, verdict))
    return calls


# --------------------------------------------------------------------------
# training-space proxy
# --------------------------------------------------------------------------

def conserved_positions(locus: str, threshold: float = 0.98,
                        sample: int = 400) -> Dict[int, str]:
    """Positions of the domain that are invariant *across species*.

    Derived from the data rather than asserted: a reference-frame position
    counts as invariant when at least ``threshold`` of the human + bovine +
    mouse records at this locus carry the same residue. Used as a structural
    sanity filter for panel members -- a pseudogene product or a mis-annotated
    chain breaks these positions, a real allele of another species does not.

    Two mistakes this design avoids:
      * hardcoding a position that only *looks* conserved. The alpha1 Cys of
        HLA-DQA1*01:01 is polymorphic (C/Y at alpha 11), so asserting it would
        have rejected every perfectly good DLA-DQA1 molecule.
      * calibrating on humans alone. Positions invariant within one species are
        routinely variable across species, which makes a human-only baseline
        reject real canine and feline alleles for being canine and feline.
    """
    path = data.CACHE / f"conserved_{locus}.json"
    if path.exists():
        return {int(k): v for k, v in json.loads(path.read_text()).items()}

    records = locus_records(locus)
    names = sorted(records)[:: max(1, len(records) // sample)][:sample]
    counts: Dict[int, Dict[str, int]] = {}
    n_used = 0
    ref = reference_domain(locus)
    for name in names:
        seq = records[name]
        if len(seq) < 60:
            continue
        mapping, ident = map_to_reference(seq, ref)
        if ident < 0.4:
            continue
        n_used += 1
        for pos, aa in mapping.items():
            tally = counts.setdefault(pos, {})
            tally[aa] = tally.get(aa, 0) + 1
    out: Dict[int, str] = {}
    for pos, tally in counts.items():
        aa, n = max(tally.items(), key=lambda kv: kv[1])
        if n_used and n / n_used >= threshold and sum(tally.values()) >= 0.5 * n_used:
            out[pos] = aa
    path.write_text(json.dumps({str(k): v for k, v in out.items()}))
    return out


def conservation_match(mol: Molecule, locus: Optional[str] = None) -> float:
    """Fraction of the locus's invariant positions that this molecule keeps."""
    conserved = conserved_positions(locus or mol.locus)
    covered = [(p, aa) for p, aa in conserved.items() if p in mol.aligned]
    if not covered:
        return 0.0
    return sum(mol.aligned[p] == aa for p, aa in covered) / len(covered)


def conservation_floor(locus: str, percentile: float = 1.0, sample: int = 300) -> float:
    """Empirical pass mark for ``conservation_match``.

    Set to the ``percentile``-th percentile of the same statistic computed over
    the training records themselves, so a panel molecule is judged against how
    well *real* class II alleles score -- not against an invented cut-off.
    """
    path = data.CACHE / f"conservation_floor_{locus}.json"
    if path.exists():
        return float(json.loads(path.read_text())["floor"])
    records = locus_records(locus)
    names = sorted(records)[:: max(1, len(records) // sample)][:sample]
    ref = reference_domain(locus)
    values: List[float] = []
    for name in names:
        seq = records[name]
        if len(seq) < 60:
            continue
        mapping, ident = map_to_reference(seq, ref)
        if ident < 0.4:
            continue
        mol = Molecule(name, "training", locus, LOCUS_REF[locus]["chain"], seq,
                       aligned=mapping)
        values.append(conservation_match(mol, locus))
    floor = float(np.percentile(values, percentile)) if values else 0.0
    path.write_text(json.dumps({"floor": floor, "n": len(values)}))
    return floor


def locus_records(locus: str) -> Dict[str, str]:
    """All training-species records for one locus: human + bovine (+ mouse)."""
    records: Dict[str, str] = {}
    if locus == "DRB":
        records.update({f"HLA-{k}": v for k, v in data.imgt_hla("DRB").items()})
        records.update(data.ipd_alleles(("BoLA-DRB3",)))
        records.update(data.mouse_mhc2("beta"))
    elif locus == "DQB":
        records.update({f"HLA-{k}": v for k, v in data.imgt_hla("DQB1").items()})
        records.update(data.ipd_alleles(("BoLA-DQB",)))
    elif locus == "DQA":
        records.update({f"HLA-{k}": v for k, v in data.imgt_hla("DQA1").items()})
        records.update(data.ipd_alleles(("BoLA-DQA",)))
        records.update(data.mouse_mhc2("alpha"))
    else:
        raise ValueError(locus)
    return records


def training_space(locus: str, cache_name: Optional[str] = None) -> List[Tuple[str, str]]:
    """Pseudosequences of the NetMHCIIpan-4.3 training species for one locus.

    Human (IPD-IMGT/HLA) + bovine (IPD-MHC BoLA) + mouse (UniProt H-2) for beta
    chains; human + bovine for the DQ alpha chain. Deduplicated at the
    pseudosequence level -- alleles that differ outside the groove are the same
    point in the model's input space.
    """
    cache_name = cache_name or f"training_{locus}.json"
    path = data.CACHE / cache_name
    if path.exists():
        return [tuple(x) for x in json.loads(path.read_text())]

    records = locus_records(locus)
    seen: Dict[str, str] = {}
    for name, seq in records.items():
        if len(seq) < 60:
            continue
        mol = build_molecule(name, "training", locus, seq)
        # Reject anything that is not really this locus's domain: a wrong chain
        # still aligns end-to-end under free end gaps, but not at this identity.
        if mol.coverage < 0.9 or mol.identity_to_ref < 0.30:
            continue
        seen.setdefault(mol.pseudoseq, name)
    pool = [(name, ps) for ps, name in seen.items()]
    path.write_text(json.dumps(pool))
    return pool
