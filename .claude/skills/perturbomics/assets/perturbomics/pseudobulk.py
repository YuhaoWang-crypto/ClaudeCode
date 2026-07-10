"""From raw single-cell perturbation counts to a signature — the RIGOROUS way.

Implements the sample-level (pseudobulk) differential-expression workflow that
the single-cell best-practices book recommends over naive per-cell tests
(Wilcoxon/t-test on cells inflates false positives because cells from one
sample are not independent replicates — the pseudoreplication trap).

    raw counts (cells x genes, with a `sample` and a `condition` label)
        --pseudobulk (sum counts within sample x group)-->  samples x genes
        --DGE (PyDESeq2 if available, else a filtered pseudobulk fallback)-->
        --Signature (signed Wald stat per gene)

If ``pydeseq2``/``decoupler`` are installed the real pipeline runs; otherwise a
transparent numpy fallback (CPM-logged pseudobulk + Welch t on the sample
means) is used so the demo works with only numpy/pandas/scipy. The fallback is
adequate for teaching / a synthetic demo but for real data install pydeseq2.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .signature import Signature


def pseudobulk(
    counts: pd.DataFrame, sample: pd.Series, min_cells: int = 10
) -> pd.DataFrame:
    """Sum counts within each sample -> a (samples x genes) matrix. RIGOROUS.

    Parameters
    ----------
    counts : DataFrame, cells x genes (raw integer counts)
    sample : Series indexed like ``counts.index`` giving each cell's sample id
        (a sample = one biological replicate x condition, the unit of
        independence). Aggregating to this level is the whole point.
    min_cells : drop pseudobulk samples backed by fewer than this many cells.
    """
    sample = sample.reindex(counts.index)
    n_cells = sample.value_counts()
    keep = n_cells[n_cells >= min_cells].index
    pb = counts.groupby(sample).sum()
    pb = pb.loc[pb.index.isin(keep)]
    pb.attrs["n_cells"] = n_cells.to_dict()
    return pb


def _fallback_dge(
    pb: pd.DataFrame, design: pd.Series, treatment: str, control: str,
    min_total: int = 10
) -> pd.DataFrame:
    """CPM-log pseudobulk + Welch t per gene. Fallback only (no pydeseq2)."""
    from scipy import stats

    # library-size normalise to counts-per-million, then log1p
    lib = pb.sum(axis=1).replace(0, np.nan)
    cpm = pb.div(lib, axis=0) * 1e6
    logcpm = np.log1p(cpm)

    grp = design.reindex(pb.index)
    a = logcpm[grp == treatment]
    b = logcpm[grp == control]
    # gene filter: expressed in enough samples
    expressed = (pb > 0).sum(axis=0) >= max(2, int(0.5 * len(pb)))
    expressed &= pb.sum(axis=0) >= min_total
    genes = pb.columns[expressed]

    t, p = stats.ttest_ind(a[genes], b[genes], equal_var=False, axis=0)
    lfc = a[genes].mean(axis=0) - b[genes].mean(axis=0)
    out = pd.DataFrame({
        "log2FoldChange": np.asarray(lfc) / np.log(2),
        "stat": np.nan_to_num(t),
        "pvalue": np.nan_to_num(p, nan=1.0),
    }, index=genes)
    # BH-adjust
    order = out["pvalue"].sort_values().index
    m = len(order)
    ranked = out.loc[order, "pvalue"].values
    adj = np.minimum.accumulate(
        (ranked * m / (np.arange(m, 0, -1)))[::-1])[::-1]
    out.loc[order, "padj"] = np.clip(adj, 0, 1)
    return out


def signature_from_pseudobulk(
    pb: pd.DataFrame, design: pd.Series, treatment: str, control: str,
    name: str = "dge", modality: str = "scrna_dge", meta: dict | None = None,
) -> Signature:
    """Run pseudobulk DGE and wrap the result as a Signature.

    Uses PyDESeq2 when importable (the recommended engine); otherwise the
    numpy/scipy fallback. Either way the returned Signature carries the signed
    Wald-like statistic per gene, ready for connectivity scoring.
    """
    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats

        meta_df = pd.DataFrame({"condition": design.reindex(pb.index)})
        dds = DeseqDataSet(counts=pb.astype(int), metadata=meta_df,
                           design_factors="condition", quiet=True)
        dds.deseq2()
        st = DeseqStats(dds, contrast=["condition", treatment, control],
                        quiet=True)
        st.summary()
        res = st.results_df
        engine = "pydeseq2"
    except Exception:
        res = _fallback_dge(pb, design, treatment, control)
        engine = "fallback_welch"

    md = {"engine": engine, **(meta or {})}
    return Signature.from_deseq2(res, name=name, modality=modality, meta=md)
