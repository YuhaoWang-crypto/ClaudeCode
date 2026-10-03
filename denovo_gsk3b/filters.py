"""Filter cascade and novelty assessment for generated molecules.

The cascade is deliberately ordered cheap-to-expensive and each stage records
its attrition so the report can show where candidates were lost.
"""

from __future__ import annotations

import functools

import numpy as np
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import Crippen, Descriptors, rdFingerprintGenerator
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

RDLogger.DisableLog("rdApp.*")

NOVELTY_TANIMOTO_MAX = 0.40   # must be BELOW this vs every known active
ACTIVITY_MIN = 0.50
QED_MIN = 0.60
SA_MAX = 4.5
MAX_RINGS = 5

_MORGAN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


@functools.lru_cache(maxsize=1)
def _pains_catalog() -> FilterCatalog:
    params = FilterCatalogParams()
    for cat in (FilterCatalogParams.FilterCatalogs.PAINS_A,
                FilterCatalogParams.FilterCatalogs.PAINS_B,
                FilterCatalogParams.FilterCatalogs.PAINS_C):
        params.AddCatalog(cat)
    return FilterCatalog(params)


def fingerprint(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    return None if mol is None else _MORGAN.GetFingerprint(mol)


def reference_fingerprints(smiles_list: list[str]) -> list:
    return [fp for fp in (fingerprint(s) for s in smiles_list) if fp is not None]


def max_similarity(smiles: str, ref_fps: list) -> float:
    """Max Tanimoto of `smiles` against the reference set (0.0 if unparseable)."""
    fp = fingerprint(smiles)
    if fp is None or not ref_fps:
        return 1.0
    return float(max(DataStructs.BulkTanimotoSimilarity(fp, ref_fps)))


def nearest_reference(smiles: str, ref_smiles: list[str], ref_fps: list):
    """Return (max_similarity, smiles_of_nearest_reference)."""
    fp = fingerprint(smiles)
    if fp is None or not ref_fps:
        return 1.0, None
    sims = DataStructs.BulkTanimotoSimilarity(fp, ref_fps)
    i = int(np.argmax(sims))
    return float(sims[i]), ref_smiles[i]


def is_pains_free(smiles: str) -> bool:
    mol = Chem.MolFromSmiles(smiles)
    return mol is not None and not _pains_catalog().HasMatch(mol)


def pains_hits(smiles: str) -> list[str]:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ["unparseable"]
    return [m.GetDescription() for m in _pains_catalog().GetMatches(mol)]


def ring_sane(smiles: str) -> bool:
    """At least one ring, no macrocycles, not an improbable ring pile-up."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    ri = mol.GetRingInfo()
    rings = ri.AtomRings()
    if not rings or len(rings) > MAX_RINGS:
        return False
    if any(len(r) > 7 or len(r) < 3 for r in rings):
        return False
    return True


def properties(smiles: str) -> dict:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {}
    return {
        "mw": float(Descriptors.MolWt(mol)),
        "logp": float(Crippen.MolLogP(mol)),
        "tpsa": float(Descriptors.TPSA(mol)),
        "hbd": int(Descriptors.NumHDonors(mol)),
        "hba": int(Descriptors.NumHAcceptors(mol)),
        "rotb": int(Descriptors.NumRotatableBonds(mol)),
        "rings": int(mol.GetRingInfo().NumRings()),
        "heavy_atoms": int(mol.GetNumHeavyAtoms()),
        "arom_rings": int(Descriptors.NumAromaticRings(mol)),
        "fsp3": float(Descriptors.FractionCSP3(mol)),
    }


def run_cascade(records: list[dict], ref_smiles: list[str],
                activity_min: float = ACTIVITY_MIN, qed_min: float = QED_MIN,
                sa_max: float = SA_MAX,
                novelty_max: float = NOVELTY_TANIMOTO_MAX) -> tuple[list[dict], list[dict]]:
    """Apply the cascade. Returns (survivors, attrition_log).

    `records` are scoring dicts (smiles/score/activity/qed/sa).
    """
    ref_fps = reference_fingerprints(ref_smiles)
    log: list[dict] = []

    def note(stage: str, kept: list, before: int):
        log.append({"stage": stage, "in": before, "out": len(kept),
                    "removed": before - len(kept)})

    # 1. valid + unique + ungated
    seen, pool = set(), []
    before = len(records)
    for r in records:
        smi = r.get("smiles")
        if not smi or r.get("gated") or smi in seen:
            continue
        if Chem.MolFromSmiles(smi) is None:
            continue
        seen.add(smi)
        pool.append(dict(r))
    note("valid, unique, un-gated", pool, before)

    # 2. novelty vs known actives
    before = len(pool)
    for r in pool:
        sim, near = nearest_reference(r["smiles"], ref_smiles, ref_fps)
        r["max_sim_to_known"] = sim
        r["nearest_known"] = near
    pool = [r for r in pool if r["max_sim_to_known"] < novelty_max]
    note(f"novel (max Tanimoto < {novelty_max})", pool, before)

    # 3. activity + QED gates
    before = len(pool)
    pool = [r for r in pool if r["activity"] >= activity_min and r["qed"] >= qed_min]
    note(f"activity >= {activity_min}, QED >= {qed_min}", pool, before)

    # 4. PAINS
    before = len(pool)
    pool = [r for r in pool if is_pains_free(r["smiles"])]
    note("PAINS-free", pool, before)

    # 5. ring sanity
    before = len(pool)
    pool = [r for r in pool if ring_sane(r["smiles"])]
    note("ring sanity", pool, before)

    # 6. synthetic accessibility
    before = len(pool)
    pool = [r for r in pool if r["sa"] <= sa_max]
    note(f"SA_Score <= {sa_max}", pool, before)

    pool.sort(key=lambda r: r["score"], reverse=True)
    return pool, log
