# Micro → macro: statistical physics of omics ensembles (D9)

`netjac/statmech.py` is the engine and `demos/d9_micro_to_macro.py` is the
worked study. Numbers for each result are in `NETJAC_REPORT.md` §10.

## Dictionary

| Physics | Omics | Function |
|---|---|---|
| micro-state | one sample / one cell (vector of 6k–10k proteins or 8.5k mRNAs) | — |
| ensemble | a cohort (tumours, normals), or the control cells of a Perturb-seq screen | — |
| order parameter | collective mode: eigenvector of the correlation matrix above the column-shuffle (or Marchenko–Pastur) null | `collective_modes`, `modes` |
| effective dimension | participation ratio PR = (Σλ)²/Σλ² | `participation_ratio` |
| couplings | partial correlations = −precision (Ledoit–Wolf + Woodbury, never forms p×p) | `ShrinkPrecision` |
| mean field | top-k modes; subtract them before you call anything an "edge" | `remove_modes` |
| external field | a knockdown / drug; response = −J⁻¹B h | — |
| fluctuation–dissipation | detailed balance ⇒ response ∝ covariance column | `fdt_score`, D9 part D |
| energy | ½ xᵀC_N⁻¹x relative to the normal ensemble; soft vs stiff split | D9 part C6 |
| universality | subspace overlap of top-k modes across indications / cell types | `subspace_overlap` |

## Exact results to remember (part A, all checked vs Lyapunov)

- 2-node loop: `C12 = (k12·C22 + k21·C11)/(a+b)`.
  - Positive feedback adds up: strong, sign-visible correlation (median |r| 0.65 over random stable loops).
  - Negative feedback cancels: median |r| 0.17, lower than a one-way edge (0.34).
  - So homeostatic loops hide themselves in static data.
- Symmetric J with isotropic noise gives −J⁻¹ = C exactly (FDT).
- A homeostat acting through a hidden sensed variable (SREBP–sterol logic) gives `corr(E,T) > 0` but `dT > 0` when E is knocked down. The response has the opposite sign to the fluctuation: **anti-FDT**.

## Workflow for a new indication / dataset

1. **Equal n.** Every spectral descriptor (PR, number of modes, mean |r|) depends on n. Subsample both ensembles to the same n, several times.
2. **Null.** Count modes against a column-shuffle null, which keeps the marginals. Do not rely on Marchenko–Pastur alone, because omics data are not Gaussian.
3. **Name the modes.** Run g:Profiler on the top / bottom 100 loadings, with the measured genes as background (`data/gprofiler.py`). In normal tissue mode 1 is usually composition (muscle / stroma / blood vs parenchyma).
4. **Subtract the macro fields before calling edges.** `rewiring_scan(..., k_remove=k)` with a label-permutation null.
5. **Prior signs.** Load the prior networks from `data/omnipath.py`: signed edges, 2-loops, complexes. Expect complexes to stay coupled.
   - Inhibition edges are not negatively correlated at the abundance level.
   - Negative loops are weaker than positive loops.
6. **Directed J needs perturbations.** Use Perturb-seq (`data/perturbseq.py`).
   - Split the cells into halves and accept only edges that replicate.
   - Sparse (lasso) MRA recovers direct edges. Naive inversion of R does not replicate (split-half r ≈ 0).
   - Row-wise regression MRA cannot see negative 2-loops: its signs follow the symmetric co-response. Read loop signs from the asymmetric **total** responses instead.
   - The FDT / anti-FDT ranking: for each knockdown j, rank j's own top-50 co-fluctuating partners against every other gene's partner set. Anti-FDT means homeostatic negative feedback.

## Pitfalls hit

- Building `set(genes)` inside a list comprehension over 50k edges is O(n·p) and cost minutes per cohort. Build the set once.
- A random-partner null for the FDT index is NOT flat. Each knockdown's global response moves everyone's partners, so rank against all partner sets.
- PR saturates at n−1. Ten-cell donor subsamples (islet scRNA) cannot resolve population order.
- PDC `clinicalPerStudy` returns empty fields. Take CPTAC-3 grade / stage from GDC `cases` (`tcga.clinical("CPTAC-3")`).
- HBV-HCC (PDC000198) has no biospecimen sample type. Aliquot labels T<n> / P<n> carry it.

## D10 additions

- **Phosphosites**: PDC returns `NP_xxx.v:s123`. Map with NCBI gene2refseq: stream the 2.4 GB file once (≈2 min) and keep the human rows. Ensembl BioMart was down; useast returned "Service unavailable".
- **The raw phospho layer has a huge global mode**: a random kinase's activity correlates 0.4–0.55 with any site. Protein-normalise each site, then remove modes.
- **No split halves in genome-wide bulk**: empirical non-targeting tails are heavy (|z|>6: 9e-5 vs Gaussian 2e-9). Validate in the independent screen instead. The anti-FDT count alone is at null level; the evidence is family enrichment plus cross-screen rank agreement.
- **gwps bulk contains ±inf / NaN entries (6,881)**: set them to 0 before any SVD (LAPACK DLASCL errors otherwise).
- **HBV-HCC pairing**: the proteome uses aliquots T471 / P472 and the phosphoproteome uses 471 / 472. Match on the aliquot number.
- **Patch-seq exocytosis is voltage-clamped**: it measures amplification and capacity, not the triggering threshold. Take triggering thresholds from Ca²⁺ recruitment data.
