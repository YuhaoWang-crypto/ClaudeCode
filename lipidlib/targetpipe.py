"""Generic target-mining + evidence-fusion engine (Problem B).

One place for the network mining and the enrichment maths so any target can be run
with one command (see ``analysis/run_target.py``). Three evidence streams:

  * ChEMBL    — measured actives           (mine_chembl)
  * DrugCLIP  — predicted virtual hits      (mine_drugclip)
  * humanPPI  — predicted protein partners  (mine_humanppi)

plus the analyses used on GLP1R:

  * molecular_enrichment  — nearest-neighbour Tanimoto of a foreground set to a
                            reference set vs a decoy set (whole-molecule)
  * fragment_enrichment   — substructure containment of fragments in actives vs decoy
  * ppi_surface_enrichment— cell-surface fraction of an interactome vs decoys

All functions are network- or CPU-only and return plain data (DataFrames / dicts);
plotting lives in the CLI so this module stays import-light.
"""

from __future__ import annotations

import time
from typing import Iterable

import numpy as np
import pandas as pd
import requests
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdFingerprintGenerator, DataStructs, rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize
from scipy.stats import mannwhitneyu, fisher_exact, binomtest

RDLogger.DisableLog("rdApp.*")

CHEMBL_API = "https://www.ebi.ac.uk/chembl/api/data"
DRUGCLIP_API = "https://dtwgapi.yanyanlan.com"
HUMANPPI_API = "http://prodata.swmed.edu/humanPPI"
SM_MW_CUTOFF = 900.0
GENOME_MEMBRANE_FRAC = 0.25

_GEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
_UNCHARGE = rdMolStandardize.Uncharger()


# ======================================================================
# mining
# ======================================================================

def resolve_chembl_target(uniprot: str) -> str | None:
    """UniProt accession -> ChEMBL target id (prefers human SINGLE PROTEIN)."""
    r = requests.get(f"{CHEMBL_API}/target.json",
                     params={"target_components__accession": uniprot, "limit": 50},
                     timeout=60)
    r.raise_for_status()
    targets = r.json().get("targets", [])
    single = [t for t in targets if t.get("target_type") == "SINGLE PROTEIN"]
    for t in (single or targets):
        if t.get("organism") in ("Homo sapiens", None):
            return t.get("target_chembl_id")
    return targets[0]["target_chembl_id"] if targets else None


def mine_chembl(target_chembl_id: str, min_pchembl: float = 6.0) -> pd.DataFrame:
    """All potent activities for a ChEMBL target, deduped to best per molecule."""
    rows, offset, limit = [], 0, 1000
    while True:
        r = requests.get(f"{CHEMBL_API}/activity.json",
                         params={"target_chembl_id": target_chembl_id,
                                 "pchembl_value__gte": min_pchembl,
                                 "limit": limit, "offset": offset}, timeout=90)
        r.raise_for_status()
        data = r.json()
        acts = data.get("activities", [])
        rows.extend(acts)
        total = data.get("page_meta", {}).get("total_count", len(rows))
        offset += limit
        if offset >= total or not acts:
            break
        time.sleep(0.2)

    best: dict[str, dict] = {}
    for a in rows:
        smi, cid, pch = a.get("canonical_smiles"), a.get("molecule_chembl_id"), a.get("pchembl_value")
        if not (smi and cid and pch):
            continue
        pch = float(pch)
        if cid not in best or pch > best[cid]["pchembl"]:
            best[cid] = {"molecule_chembl_id": cid, "smiles": smi,
                         "standard_type": a.get("standard_type"), "pchembl": pch}
    out = pd.DataFrame(best.values())
    if len(out):
        out["mw"] = out["smiles"].map(_mw)
        out["class"] = np.where(out["mw"] < SM_MW_CUTOFF, "small_molecule", "peptide_or_large")
        out = out.sort_values("pchembl", ascending=False).reset_index(drop=True)
    return out


def mine_drugclip(uniprot: str) -> pd.DataFrame:
    """DrugCLIP predicted pocket-ligand hits + resolved SMILES for a UniProt.

    Returns an empty DataFrame if DrugCLIP has no coverage for this protein
    (404) rather than raising, so the pipeline degrades gracefully.
    """
    r = requests.get(f"{DRUGCLIP_API}/complexes/{uniprot}", timeout=60)
    if r.status_code == 404:
        return pd.DataFrame()
    r.raise_for_status()
    comps = r.json()
    if not comps:
        return pd.DataFrame()
    sr = requests.post(f"{DRUGCLIP_API}/get_smiles",
                       json={"complex_list": [c["ligand_id"] for c in comps]}, timeout=120)
    sr.raise_for_status()
    blob = sr.json(); txt = blob.get("data", blob)
    txt = txt.get("smiles") if isinstance(txt, dict) else txt
    smap = {}
    for line in str(txt).splitlines():
        p = line.split("\t")
        if len(p) >= 2:
            smap[p[1].strip()] = p[0].strip()
    rows = [{"ligand_id": c.get("ligand_id_origin"),
             "smiles": smap.get(c.get("ligand_id_origin"), ""),
             "drugclip_score": round(c.get("drugclip_score", np.nan), 4),
             "docking_score": round(c.get("docking_score", np.nan), 3),
             "source": c.get("source"), "pocket_id": c.get("pocket_id"),
             "nearby_residues": ";".join(c.get("nearby_residues", []))}
            for c in comps]
    df = pd.DataFrame(rows).drop_duplicates("ligand_id").reset_index(drop=True)
    df["mw"] = df["smiles"].map(_mw)
    return df.sort_values("drugclip_score", ascending=False).reset_index(drop=True)


def mine_humanppi(uniprot: str, high_confidence_only: bool = False) -> pd.DataFrame:
    r = requests.get(f"{HUMANPPI_API}/data/{uniprot}",
                     params={"filter": str(high_confidence_only).lower(),
                             "length": 100000, "start": 0, "draw": 1}, timeout=90)
    if r.status_code == 404:
        return pd.DataFrame()
    r.raise_for_status()
    df = pd.DataFrame(r.json().get("data", []))
    if len(df):
        loc = df["Subcellular Locality"].astype(str)
        df["is_membrane"] = loc.str.contains("membrane", case=False)
        df["is_cell_surface"] = loc.str.contains("Cell membrane", case=False)
        df["AF"] = pd.to_numeric(df["AF_Score"], errors="coerce")
    return df


# ======================================================================
# chemistry helpers
# ======================================================================

def _mol(s):
    return Chem.MolFromSmiles(str(s))


def _mw(s):
    m = _mol(s); return Descriptors.MolWt(m) if m else np.nan


def _neutral(s):
    m = _mol(s)
    if m is None:
        return None
    try:
        return _UNCHARGE.uncharge(m)
    except Exception:
        return m


def fps(smiles: Iterable[str]):
    out = []
    for s in smiles:
        m = _mol(s)
        out.append(_GEN.GetFingerprint(m) if m else None)
    return out


def nearest_sim(query_fps, ref_fps) -> np.ndarray:
    ref = [f for f in ref_fps if f is not None]
    out = np.full(len(query_fps), np.nan)
    for i, q in enumerate(query_fps):
        if q is not None and ref:
            out[i] = max(DataStructs.BulkTanimotoSimilarity(q, ref))
    return out


# ======================================================================
# enrichment analyses
# ======================================================================

def molecular_enrichment(fg_smiles, ref_smiles, decoy_smiles, thr: float = 0.35) -> dict:
    """Nearest-neighbour Tanimoto of foreground vs reference, tested against decoy."""
    ref = fps(ref_smiles)
    fg = nearest_sim(fps(fg_smiles), ref)
    dc = nearest_sim(fps(decoy_smiles), ref)
    fg, dc = fg[~np.isnan(fg)], dc[~np.isnan(dc)]
    u, p = mannwhitneyu(fg, dc, alternative="greater") if len(fg) and len(dc) else (np.nan, np.nan)
    ef = ((fg >= thr).mean() / max((dc >= thr).mean(), 1e-9)) if len(dc) else np.nan
    return {"fg_sim": fg, "decoy_sim": dc, "fg_median": float(np.median(fg)) if len(fg) else np.nan,
            "decoy_median": float(np.median(dc)) if len(dc) else np.nan,
            "mwu_p": float(p), "enrichment_factor": float(ef), "threshold": thr,
            "n_fg": int(len(fg)), "n_decoy": int(len(dc))}


def fragment_enrichment(frag_smiles, decoy_smiles, active_smiles,
                        min_heavy: int = 10) -> dict:
    """Substructure containment of fragments within actives vs decoy fragments."""
    actives = [a for a in (_neutral(s) for s in active_smiles) if a is not None]

    def rate(frs):
        hit = tot = 0
        matched = []
        for s in frs:
            m = _neutral(s)
            if m is None or m.GetNumHeavyAtoms() < min_heavy or rdMolDescriptors.CalcNumRings(m) < 1:
                continue
            tot += 1
            n = sum(a.HasSubstructMatch(m) for a in actives)
            if n:
                hit += 1; matched.append((s, n))
        return hit, tot, matched

    gh, gt, gm = rate(frag_smiles)
    dh, dt, _ = rate(decoy_smiles)
    odds, p = fisher_exact([[gh, gt - gh], [dh, dt - dh]], alternative="greater") \
        if gt and dt else (np.nan, np.nan)
    return {"fg_hit": gh, "fg_tot": gt, "decoy_hit": dh, "decoy_tot": dt,
            "odds": float(odds), "fisher_p": float(p), "matched": gm}


def ppi_surface_enrichment(target_df: pd.DataFrame, decoys: dict[str, pd.DataFrame]) -> dict:
    """Cell-surface fraction of an interactome vs decoy interactomes + genome baseline."""
    tm, tt = int(target_df["is_cell_surface"].sum()), len(target_df)
    out = {"target_surface": tm, "target_total": tt,
           "target_frac": tm / tt if tt else np.nan, "vs": {}}
    for name, d in decoys.items():
        dm, dt = int(d["is_cell_surface"].sum()), len(d)
        odds, p = fisher_exact([[tm, tt - tm], [dm, dt - dm]], alternative="greater")
        out["vs"][name] = {"frac": dm / dt if dt else np.nan, "odds": float(odds), "p": float(p)}
    if tt:
        out["vs_genome"] = {"frac": GENOME_MEMBRANE_FRAC,
                            "p": binomtest(tm, tt, GENOME_MEMBRANE_FRAC, alternative="greater").pvalue}
    return out
