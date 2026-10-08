"""Read the metadata parts of a remote .h5ad without downloading the matrix.

h5py opens an fsspec HTTP file object, so only the byte ranges of ``obs``,
``var`` and a few group headers are fetched; ``X`` / ``layers`` / ``obsm``
are never read.  Supports the anndata >= 0.8 on-disk encoding (all ArchMap
atlases use it).
"""
from __future__ import annotations

import re

import fsspec
import h5py
import numpy as np

ENSEMBL_RE = re.compile(r"^ENS[A-Z]*G\d{6,}")
ENSEMBL_COLS = ["ensembl", "gene_id", "gene_ids", "EID", "Accession", "GeneID",
                "feature_id"]
SYMBOL_COLS = ["feature_name", "gene_symbol", "gene_symbols", "GeneName", "Gene",
               "gene_name", "features", "original_gene_names", "gene_name-new"]
MAX_STORED_CATEGORIES = 20000


def _decode(arr):
    arr = np.asarray(arr)
    if arr.dtype.kind in "OS":
        return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in arr],
                        dtype=object)
    return arr


def _index_name(group):
    name = group.attrs.get("_index", "_index")
    return name.decode() if isinstance(name, bytes) else name


def read_column(group, name):
    """Full column as a numpy array (categoricals decoded, nullable -> NaN/None)."""
    item = group[name]
    if isinstance(item, h5py.Group):
        if "codes" in item:
            codes = item["codes"][()]
            cats = _decode(item["categories"][()])
            out = np.empty(len(codes), dtype=object)
            ok = codes >= 0
            out[ok] = cats[codes[ok]]
            out[~ok] = None
            return out
        if "values" in item and "mask" in item:
            vals = _decode(item["values"][()]).astype(object)
            vals[item["mask"][()].astype(bool)] = None
            return vals
        raise ValueError(f"unknown encoding for {name}: {dict(item.attrs)}")
    return _decode(item[()])


def _summarise_obs_column(item):
    """-> (column-info dict, {category: n} | None, codes | None, categories | None)."""
    info = {"kind": None, "dtype": None, "n_categories": None, "n_missing": 0,
            "min": None, "max": None, "mean": None, "median": None}
    if isinstance(item, h5py.Group) and "codes" in item:
        codes = item["codes"][()]
        cats = _decode(item["categories"][()])
        info.update(kind="categorical", dtype=str(cats.dtype), n_categories=len(cats),
                    n_missing=int((codes < 0).sum()))
        counts = np.bincount(codes[codes >= 0], minlength=len(cats))
        cc = None
        if len(cats) <= MAX_STORED_CATEGORIES:
            cc = {str(c): int(n) for c, n in zip(cats, counts)}
        return info, cc, codes, cats
    if isinstance(item, h5py.Group) and "values" in item:
        vals = item["values"][()]
        mask = item["mask"][()].astype(bool)
        info["n_missing"] = int(mask.sum())
        vals = vals[~mask]
        etype = item.attrs.get("encoding-type", "")
        etype = etype.decode() if isinstance(etype, bytes) else etype
        if vals.dtype == bool or "boolean" in etype:
            info.update(kind="boolean", dtype="bool", n_categories=2)
            return info, {"True": int(vals.sum()), "False": int((~vals.astype(bool)).sum())}, None, None
        item_vals = vals
    elif isinstance(item, h5py.Dataset):
        if item.dtype.kind in "OSU":
            info.update(kind="string", dtype="str")  # free text / barcodes: not read
            return info, None, None, None
        item_vals = item[()]
        if item_vals.dtype == bool:
            info.update(kind="boolean", dtype="bool", n_categories=2)
            n_true = int(item_vals.sum())
            return info, {"True": n_true, "False": len(item_vals) - n_true}, None, None
    else:
        info.update(kind="other")
        return info, None, None, None
    v = item_vals.astype(float)
    finite = v[np.isfinite(v)]
    info.update(kind="numeric", dtype=str(item_vals.dtype),
                n_missing=info["n_missing"] + int((~np.isfinite(v)).sum()))
    if finite.size:
        info.update(min=float(finite.min()), max=float(finite.max()),
                    mean=float(finite.mean()), median=float(np.median(finite)))
    return info, None, None, None


def _crosstab(codes_a, cats_a, codes_b, cats_b):
    ok = (codes_a >= 0) & (codes_b >= 0)
    a, b = codes_a[ok].astype(np.int64), codes_b[ok].astype(np.int64)
    flat = np.bincount(a * len(cats_b) + b, minlength=len(cats_a) * len(cats_b))
    idx = np.nonzero(flat)[0]
    return [(str(cats_a[i // len(cats_b)]), str(cats_b[i % len(cats_b)]), int(flat[i]))
            for i in idx]


def _norm_ensembl(x):
    if x is None:
        return None
    s = str(x)
    return s.split(".")[0] if ENSEMBL_RE.match(s) else None


def _resolve_key(candidates, available):
    """First candidate present in obs (exact match, then case-insensitive)."""
    lower = {c.lower(): c for c in available}
    for cand in candidates:
        if cand in available:
            return cand
        if cand and cand.lower() in lower:
            return lower[cand.lower()]
    return None


def summarise_h5ad(url, cell_type_keys=(), batch_keys=(), block_size=2 ** 22):
    """cell_type_keys / batch_keys: candidate obs columns in order of preference."""
    fs = fsspec.filesystem("http")
    res = {}
    with fs.open(url, "rb", block_size=block_size, cache_type="blockcache") as fh, \
            h5py.File(fh, "r") as h:
        X = h["X"]
        if isinstance(X, h5py.Group):
            shape = tuple(int(s) for s in X.attrs["shape"])
            enc = X.attrs.get("encoding-type", "sparse")
            enc = enc.decode() if isinstance(enc, bytes) else enc
            x_dtype = str(X["data"].dtype)
        else:
            shape, enc, x_dtype = tuple(int(s) for s in X.shape), "dense", str(X.dtype)
        res.update(n_obs=shape[0], n_vars=shape[1], x_encoding=enc, x_dtype=x_dtype)
        res["obsm_keys"] = sorted(h["obsm"].keys()) if "obsm" in h else []
        res["layers_keys"] = sorted(h["layers"].keys()) if "layers" in h else []
        res["uns_keys"] = sorted(h["uns"].keys()) if "uns" in h else []

        obs = h["obs"]
        idx_name = _index_name(obs)
        columns, cat_counts, coded = [], [], {}
        for name in obs.keys():
            if name == idx_name or name == "__categories":
                continue
            info, cc, codes, cats = _summarise_obs_column(obs[name])
            info["column"] = name
            columns.append(info)
            if cc:
                cat_counts += [(name, c, n) for c, n in cc.items()]
            if codes is not None:
                coded[name] = (codes, cats)
        res["obs_columns"], res["category_counts"] = columns, cat_counts

        # cell type x batch, and cell type x Cell Ontology (CELLxGENE schema)
        cell_type_key = _resolve_key(cell_type_keys, coded)
        batch_key = _resolve_key(batch_keys, coded)
        res["cell_type_key_used"], res["batch_key_used"] = cell_type_key, batch_key
        res["celltype_batch"], res["celltype_ontology"] = [], []
        if cell_type_key in coded and batch_key in coded and batch_key != cell_type_key:
            res["celltype_batch"] = _crosstab(*coded[cell_type_key], *coded[batch_key])
        if cell_type_key in coded and "cell_type_ontology_term_id" in coded:
            res["celltype_ontology"] = _crosstab(*coded[cell_type_key],
                                                 *coded["cell_type_ontology_term_id"])

        var = h["var"]
        vidx = read_column(var, _index_name(var))
        vcols = {c: read_column(var, c) for c in ENSEMBL_COLS + SYMBOL_COLS if c in var}
        genes = []
        for i, fid in enumerate(vidx):
            ens = _norm_ensembl(fid)
            for c in ENSEMBL_COLS:
                if ens is None and c in vcols:
                    ens = _norm_ensembl(vcols[c][i])
            sym = None
            for c in SYMBOL_COLS:
                if c in vcols and vcols[c][i] is not None and not _norm_ensembl(vcols[c][i]):
                    sym = str(vcols[c][i])
                    break
            if sym is None and not _norm_ensembl(fid):
                sym = str(fid)
            genes.append((i, str(fid), sym, ens))
        res["genes"] = genes
        res["var_columns"] = sorted(var.keys())
    return res
