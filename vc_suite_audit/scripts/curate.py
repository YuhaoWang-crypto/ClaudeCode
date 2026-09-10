"""Curate ChEMBL activity dumps into modelling tables, following the report's stated protocol."""
import json, os, math
import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = HERE + "/data"


def largest_fragment(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    frags = Chem.GetMolFrags(m, asMols=True, sanitizeFrags=True)
    if not frags:
        return None
    return max(frags, key=lambda f: f.GetNumHeavyAtoms())


def load(name):
    rows = [json.loads(l) for l in open(f"{DATA}/{name}.jsonl")]
    return pd.DataFrame(rows)


def curate(names, kind, types=("IC50",), mw_max=None, prange=(2, 10)):
    """kind='MIC' -> ug/mL + nM accepted; kind='ACT' -> nM potency types."""
    df = pd.concat([load(n) for n in names], ignore_index=True)
    df = df[df.standard_relation == "="]
    df = df[df.standard_type.isin(types)]
    df = df[df.standard_value.notna() & df.canonical_smiles.notna()]
    df = df[df.data_validity_comment.isna()]
    df["standard_value"] = pd.to_numeric(df.standard_value, errors="coerce")
    df = df[df.standard_value > 0]

    recs = []
    cache = {}
    for r in df.itertuples():
        smi = r.canonical_smiles
        if smi not in cache:
            m = largest_fragment(smi)
            cache[smi] = (Chem.MolToSmiles(m), Descriptors.MolWt(m)) if m else (None, None)
        can, mw = cache[smi]
        if can is None:
            continue
        u = (r.standard_units or "").strip()
        if u == "nM":
            p = 9 - math.log10(r.standard_value)
        elif u in ("ug.mL-1", "ug ml-1", "ug/ml"):
            if not mw or mw <= 0:
                continue
            p = -math.log10(r.standard_value * 1e-3 / mw)
        else:
            continue
        if not (prange[0] <= p <= prange[1]):
            continue
        if mw_max and mw > mw_max:
            continue
        recs.append(dict(smiles=can, mol_id=r.molecule_chembl_id, p=p, mw=mw,
                         year=r.document_year, doc=r.document_chembl_id,
                         stype=r.standard_type))
    raw = pd.DataFrame(recs)
    if raw.empty:
        return raw
    # median aggregation per compound; keep earliest document year for temporal split
    g = raw.groupby("smiles")
    out = g.agg(p=("p", "median"), n_meas=("p", "size"), mw=("mw", "first"),
                year=("year", "min"), mol_id=("mol_id", "first")).reset_index()
    out["scaffold"] = [MurckoScaffold.MurckoScaffoldSmiles(smiles=s) or s for s in out.smiles]
    return out


CONFIGS = {
    "NA_influenza":  dict(names=["na_CHEMBL6135", "na_CHEMBL2051"], kind="ACT",
                          types=("IC50", "Ki", "EC50"), mw_max=650),
    "PA_influenza":  dict(names=["pa_CHEMBL1169598"], kind="ACT",
                          types=("IC50", "Ki", "EC50"), mw_max=650),
    "MIC_H_influenzae": dict(names=["mic_h_influenzae"], kind="MIC", types=("MIC",)),
    "MIC_E_cloacae":    dict(names=["mic_e_cloacae"], kind="MIC", types=("MIC",)),
    "hERG":          dict(names=["herg_CHEMBL240"], kind="ACT", types=("IC50",)),
    "CC50_pan":      dict(names=["cc50_all"], kind="ACT", types=("CC50",)),
}

if __name__ == "__main__":
    import sys
    os.makedirs(HERE + "/curated", exist_ok=True)
    want = sys.argv[1:] or list(CONFIGS)
    for k, cfg in {k: v for k, v in CONFIGS.items() if k in want}.items():
        try:
            d = curate(**cfg)
        except FileNotFoundError as e:
            print(f"{k}: missing raw file ({e})"); continue
        d.to_csv(f"{HERE}/curated/{k}.csv", index=False)
        print(f"{k}: {len(d)} compounds, {d.scaffold.nunique()} scaffolds, "
              f"p={d.p.mean():.2f}+-{d.p.std():.2f}, years {d.year.min()}-{d.year.max()}")
