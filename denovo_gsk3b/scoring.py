"""Multi-objective scoring for the GSK3-beta de novo design campaign.

The production objective is the geometric mean of three normalised terms:

    score = (activity * qed * sa_norm) ** (1/3)

multiplied by a hard *makeability / drug-likeness gate* that returns 0.0 for
molecules outside a physicochemical envelope. The gate exists because a genetic
algorithm optimising a surrogate activity model will otherwise happily drift
into high-scoring but unsynthesisable or non-drug-like chemical space.

Activity is the pretrained TDC GSK3B oracle (a random forest over ECFP
features, from Li et al. / the GuacaMol-MOSES lineage used by TDC). It is a
SURROGATE, not a measurement -- see the report's Limitations section.
"""

from __future__ import annotations

import functools
import warnings

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors

warnings.filterwarnings("ignore")
RDLogger.DisableLog("rdApp.*")

# --- physicochemical gate bounds -------------------------------------------
MW_MIN, MW_MAX = 150.0, 600.0
HEAVY_MIN, HEAVY_MAX = 10, 50
LOGP_MAX = 6.0
SA_HARD_MAX = 6.0          # above this, treat as unmakeable
MAX_RING_SIZE = 7
ALLOWED_ELEMENTS = {"C", "N", "O", "S", "F", "Cl", "Br", "H"}


@functools.lru_cache(maxsize=1)
def _oracles():
    """Load the TDC oracles once (they unpickle sklearn models from disk)."""
    from tdc import Oracle

    return {
        "gsk3b": Oracle(name="GSK3B"),
        "qed": Oracle(name="QED"),
        "sa": Oracle(name="SA"),
    }


def sa_to_norm(sa: float) -> float:
    """Map synthetic accessibility (1 easy .. 10 hard) onto (1 good .. 0 bad)."""
    return float(np.clip((10.0 - sa) / 9.0, 0.0, 1.0))


def passes_gate(mol: Chem.Mol) -> bool:
    """Hard makeability / drug-likeness gate. Cheap checks only (no oracles)."""
    if mol is None:
        return False
    if any(a.GetSymbol() not in ALLOWED_ELEMENTS for a in mol.GetAtoms()):
        return False
    if mol.GetNumAtoms() == 0:
        return False
    n_heavy = mol.GetNumHeavyAtoms()
    if not (HEAVY_MIN <= n_heavy <= HEAVY_MAX):
        return False
    mw = Descriptors.MolWt(mol)
    if not (MW_MIN <= mw <= MW_MAX):
        return False
    if Crippen.MolLogP(mol) > LOGP_MAX:
        return False
    # ring sanity: no macrocycles / exotic large rings
    ring_info = mol.GetRingInfo()
    if any(len(r) > MAX_RING_SIZE for r in ring_info.AtomRings()):
        return False
    # reject charged species and multi-fragment molecules
    if Chem.GetFormalCharge(mol) != 0:
        return False
    if len(Chem.GetMolFrags(mol)) != 1:
        return False
    return True


def canonical(smiles: str) -> str | None:
    """Sanitise and canonicalise; returns None for anything RDKit rejects."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        return None
    return Chem.MolToSmiles(mol)


def score_batch(smiles_list: list[str]) -> list[dict]:
    """Score a batch of SMILES. Returns one record per input (score 0.0 if gated).

    Batching matters: the TDC oracles are far faster called once on a list than
    once per molecule.
    """
    orc = _oracles()
    records: list[dict | None] = [None] * len(smiles_list)
    to_score: list[str] = []
    idx_map: list[int] = []

    for i, smi in enumerate(smiles_list):
        mol = Chem.MolFromSmiles(smi) if smi else None
        if mol is None or not passes_gate(mol):
            records[i] = {
                "smiles": smi,
                "score": 0.0,
                "activity": 0.0,
                "qed": 0.0,
                "sa": 10.0,
                "gated": True,
            }
        else:
            to_score.append(smi)
            idx_map.append(i)

    if to_score:
        act = orc["gsk3b"](to_score)
        qed = orc["qed"](to_score)
        sa = orc["sa"](to_score)
        act = [act] if np.isscalar(act) else list(act)
        qed = [qed] if np.isscalar(qed) else list(qed)
        sa = [sa] if np.isscalar(sa) else list(sa)

        for j, i in enumerate(idx_map):
            a, q, s = float(act[j]), float(qed[j]), float(sa[j])
            if s > SA_HARD_MAX:
                score = 0.0
                gated = True
            else:
                # geometric mean; a zero in any term zeroes the product, which
                # is the intended behaviour of a conjunctive objective.
                score = float((max(a, 0.0) * max(q, 0.0) * sa_to_norm(s)) ** (1.0 / 3.0))
                gated = False
            records[i] = {
                "smiles": to_score[j],
                "score": score,
                "activity": a,
                "qed": q,
                "sa": s,
                "gated": gated,
            }

    return [r for r in records if r is not None]


def score_one(smiles: str) -> dict:
    return score_batch([smiles])[0]
