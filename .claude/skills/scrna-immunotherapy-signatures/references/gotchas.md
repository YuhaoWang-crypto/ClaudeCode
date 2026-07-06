# Debugging log — issues hit while building this pipeline (and the fixes)

Every one of these actually occurred during the verified end-to-end run. Consult this
before re-deriving a fix.

## 1. Out-of-memory kill during AnnData construction
**Symptom:** process dies silently (SIGKILL, no traceback), `free` shows memory freed.
**Cause:** `sc.AnnData(df.T)` on a 55,737 × 16,291 dense TPM matrix holds `df`, `df.T`,
the AnnData copy, and a `np.asarray(...float32)` copy simultaneously → ~4 × 3.6 GB > 15 GB.
**Fix:**
```python
from scipy import sparse
genes = df.index.astype(str).to_numpy(); cells = list(df.columns)
arr = df.to_numpy(dtype="float32"); del df
X = sparse.csr_matrix(arr.T); del arr           # TPM is mostly zeros
adata = sc.AnnData(X=X, obs=pd.DataFrame(index=cells), var=pd.DataFrame(index=genes))
```
Also: `sc.pp.filter_genes(adata, min_cells=3)` early; do **not** set `adata.raw` or keep
a dense `tpm` layer (each doubles memory). Derive `TCF7_pos` from `X > 0` (log1p keeps zeros).

## 2. `read_csv(dtype="float32")` fails on the gene-name column
**Symptom:** `ValueError: could not convert string to float: 'TSPAN6'`.
**Fix:** apply float32 per data column, keep column 0 as object:
```python
dtypes = {0: "object"}; dtypes.update({i: "float32" for i in range(1, len(cells)+1)})
pd.read_csv(path, sep="\t", skiprows=2, header=None, index_col=0, dtype=dtypes)
```

## 3. Trailing-tab / off-by-one columns
- **GSE120575 TPM:** gene rows have a trailing tab → one extra all-NaN column. Detect and
  drop: `if df.shape[1] == len(cells)+1: df = df.iloc[:, :len(cells)]`.
- **GSE123813 counts (10x, R-style):** header has N cell names but data rows have N+1
  fields (gene + N values), and the awk column count read off the *header* dropped the
  last cell. Always compute library sizes over the **data** row field count (`NF`), then
  assert `expr.shape[1] == len(libsize)`.

## 4. Text encoding
GSE120575 metadata contains a `µ` byte (0xB5) → `UnicodeDecodeError` under UTF-8.
Open with `encoding="latin-1"`.

## 5. NaN p-values from Mann-Whitney after h5ad round-trip
**Symptom:** identical log2FCs but `pval = NaN` on a re-run that loads the cached `.h5ad`.
**Cause:** obs string columns become **categorical**; `groupby([...]).unstack()` then
creates rows for empty category combinations → NaN fractions → `mannwhitneyu` returns NaN.
**Fix:** `groupby(..., observed=True)` everywhere, and `mannwhitneyu(a, b,
nan_policy="omit")`.

## 6. matplotlib API
`Axes.boxplot(labels=[...])` removed → use `tick_labels=[...]` or `ax.set_xticklabels`.

## 7. Bulk validation of a single-cell state signature (scientific, not a bug)
The memory−exhaustion difference gives AUC ≈ 0.4–0.54 on Riaz bulk; components
(CD8 abundance) give AUC ≈ 0.78 on-treatment. This is expected: bulk conflates T-cell
*quantity* with *state*. Report it honestly; do not tune the grouping to manufacture a
positive result. Validate state signatures on single-cell external data instead.

## Performance notes
- The pandas parse of the wide (16k-column) TSV is the slow step (~3–4 min). Cache the
  processed AnnData to `.h5ad` and gate the heavy compute so re-runs are seconds.
- For scoring only a handful of genes on a huge counts file, extract just those gene rows
  (+ per-cell library sizes) with a single `awk` pass instead of loading the full matrix.
