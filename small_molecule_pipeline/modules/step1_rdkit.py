"""Step 1 — RDKit: structure standardization, deduplication, structural alerts.

Alerts checked:
  - PAINS (Pan-Assay Interference Compounds)
  - Brenk (undesirable substructures)
  - NIH compounds (known problematic scaffolds)
  - Lipinski/Veber rule-of-5 basic properties
"""

import logging
from typing import Dict, List, Optional

from rdkit import Chem
from rdkit.Chem import Descriptors, inchi
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

logger = logging.getLogger(__name__)


def _build_alert_catalogs() -> Dict[str, FilterCatalog]:
    catalogs: Dict[str, FilterCatalog] = {}

    pains_params = FilterCatalogParams()
    pains_params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    catalogs["PAINS"] = FilterCatalog(pains_params)

    brenk_params = FilterCatalogParams()
    brenk_params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
    catalogs["Brenk"] = FilterCatalog(brenk_params)

    nih_params = FilterCatalogParams()
    nih_params.AddCatalog(FilterCatalogParams.FilterCatalogs.NIH)
    catalogs["NIH"] = FilterCatalog(nih_params)

    return catalogs


_CATALOGS = _build_alert_catalogs()

_NORMALIZER = rdMolStandardize.Normalizer()
_UNCHARGER = rdMolStandardize.Uncharger()
_LARGEST_FRAG = rdMolStandardize.LargestFragmentChooser()
_TAUTOMER_PARENT = rdMolStandardize.TautomerEnumerator()


def standardize_mol(mol: Chem.Mol) -> Optional[Chem.Mol]:
    try:
        mol = _NORMALIZER.normalize(mol)
        mol = _LARGEST_FRAG.choose(mol)
        mol = _UNCHARGER.uncharge(mol)
        Chem.SanitizeMol(mol)
        return mol
    except Exception as e:
        logger.debug("Standardization failed: %s", e)
        return None


def compute_inchikey(mol: Chem.Mol) -> Optional[str]:
    try:
        return inchi.MolToInchiKey(mol)
    except Exception:
        return None


def check_alerts(mol: Chem.Mol) -> Dict[str, List[str]]:
    hits: Dict[str, List[str]] = {}
    for name, catalog in _CATALOGS.items():
        entries = catalog.GetMatches(mol)
        if entries:
            hits[name] = [e.GetDescription() for e in entries]
    return hits


def compute_properties(mol: Chem.Mol) -> Dict:
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = Descriptors.NumHDonors(mol)
    hba = Descriptors.NumHAcceptors(mol)
    tpsa = Descriptors.TPSA(mol)
    rotb = Descriptors.NumRotatableBonds(mol)
    rings = Descriptors.RingCount(mol)
    arom_rings = Descriptors.NumAromaticRings(mol)

    ro5_violations = sum([mw > 500, logp > 5, hbd > 5, hba > 10])
    veber_ok = tpsa <= 140 and rotb <= 10

    return {
        "MW": round(mw, 2),
        "LogP": round(logp, 2),
        "HBD": hbd,
        "HBA": hba,
        "TPSA": round(tpsa, 2),
        "RotBonds": rotb,
        "Rings": rings,
        "AromaticRings": arom_rings,
        "Ro5_violations": ro5_violations,
        "Veber_ok": veber_ok,
        "druglike": ro5_violations <= 1 and veber_ok,
    }


def run(smiles_list: List[str]) -> List[Dict]:
    results = []
    seen_inchikeys = set()

    for raw_smi in smiles_list:
        entry: Dict = {"input_smiles": raw_smi}

        mol = Chem.MolFromSmiles(raw_smi)
        if mol is None:
            entry["error"] = "invalid_smiles"
            entry["passed"] = False
            results.append(entry)
            continue

        std_mol = standardize_mol(mol)
        if std_mol is None:
            entry["error"] = "standardization_failed"
            entry["passed"] = False
            results.append(entry)
            continue

        std_smi = Chem.MolToSmiles(std_mol)
        inchikey = compute_inchikey(std_mol)

        if inchikey in seen_inchikeys:
            entry["error"] = "duplicate"
            entry["passed"] = False
            entry["std_smiles"] = std_smi
            entry["inchikey"] = inchikey
            results.append(entry)
            continue

        if inchikey:
            seen_inchikeys.add(inchikey)

        alerts = check_alerts(std_mol)
        props = compute_properties(std_mol)

        entry.update(
            {
                "std_smiles": std_smi,
                "inchikey": inchikey,
                "properties": props,
                "structural_alerts": alerts,
                "alert_flags": list(alerts.keys()),
                "n_alerts": sum(len(v) for v in alerts.values()),
                "passed": True,
                "error": None,
            }
        )
        results.append(entry)

    valid = sum(1 for r in results if r.get("passed"))
    logger.info(
        "Step 1 complete: %d input → %d valid unique molecules, %d with alerts",
        len(smiles_list),
        valid,
        sum(1 for r in results if r.get("passed") and r.get("n_alerts", 0) > 0),
    )
    return results
