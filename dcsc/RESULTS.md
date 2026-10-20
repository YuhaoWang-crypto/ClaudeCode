# DCSC feasibility run — results

Synthetic world: 2000 genes, 9 cell types in a 2-level tree, 32 programs, 4 memory and 4 environment factors, 6 batches, NB counts (~2500 UMI/cell). Encoder fitted from simulated reference + calibration contrasts only.

## 1. Factorized State Barcode (64 B)

Identity tree (single cells, n=4000): lineage acc = 1.000, leaf acc = 1.000

Recovery of each block (Pearson r with truth, mean over coordinates) and leakage:

| method | mod | mem | env | leak | null | tern |
|---|---|---|---|---|---|---|
| joint posterior (DCSC), single cell | 0.840 | 0.858 | 0.765 | 0.025 | 0.013 | 0.958 |
| naive gene-set score, single cell | 0.818 | 0.849 | 0.747 | 0.161 | 0.013 | 0.936 |
| joint posterior (DCSC), pseudobulk n=100 | 0.944 | 0.988 | 0.963 | 0.082 | 0.056 | 0.982 |
| naive gene-set score, pseudobulk n=100 | 0.936 | 0.988 | 0.942 | 0.149 | 0.056 | 0.985 |

`leak` = mean |r| between an estimated program and the TRUE artifact it shares 12 genes with by construction (TNF_NFkB/Heat_shock↔Dissociation, EMT/Prolif_G2M↔Culture_2D, Oxphos↔Mito); `null` = expected |r| from sampling noise alone. `tern` = ternary-state accuracy over 32 programs.

Size: 64 B/cell vs 8,000 B for a dense float32 profile of 2000 genes (×125; ×1250 for 20k genes). 100 M cells = 6.4 GB → fits in RAM.

Example barcodes:

- `Stromal>Fibroblast | Hypoxia+ Antigen_pres- AhR+ Autophagy+ | mem: TOX_exhaustion=1.0 | env: Culture_2D=1.2`
- `Lymphoid>CD8_T | IL6_JAK_STAT3- Oxphos+ Glycolysis- mTORC1- Estrogen- IL10_STAT3+ Apoptosis+ SASP+ | mem: TOX_exhaustion=1.2 Hypoxic_memory=-0.5 | env: -`
- `Stromal>Fibroblast | Prolif_G2M+ NRF2_oxidative+ | mem: LPS_tolerance=0.8 | env: Dissociation_stress=0.6 Culture_2D=1.4`

Query modes (top-10 neighbours of 200 query cells): true program-state distance (lower = same functional state), share from another lineage, share of same leaf:

| mode | prog_dist | other_lineage | same_leaf |
|---|---|---|---|
| mod | 4.86 | 0.68 | 0.10 |
| same_id | 5.63 | 0.00 | 1.00 |
| random pairs | 8.89 | 0.67 | 0.11 |

## 2. Search speed on 64-B codes (single core)

| engine | n | per_query_ms | cmp_per_s |
|---|---|---|---|
| numpy bitwise_count | 2,000,000 | 97.289 | 2.06e+07 |
| C POPCNT kernel | 2,000,000 | 12.568 | 1.59e+08 |
| faiss IndexBinaryFlat | 2,000,000 | 5.041 | 3.97e+08 |
| faiss IndexBinaryHNSW (build 29s, recall@10 1.00) | 300,000 | 0.256 | nan |

Extrapolation: brute-force scan of 100 M codes ≈ 0.63 s/query/core with the C kernel, i.e. ~0.1 s on 8 cores; HNSW gives ms-level latency at some recall cost. Ternary Hamming has only 0–64 distinct values, so ties are frequent: re-rank the top few thousand with the continuous posterior or Mem/Env blocks.

## 3. Δ-state signatures + residual pursuit (cocktail recipe)

Library: 34 stimuli measured in CD14_Mono (pseudobulk 300 vs 300 cells). Mutual coherence max cos = 0.97; near-duplicate pairs (cos>0.9): IFNa~IFNb, IL1b~LPS, IL4~IL13.

60 random queries per row; target and baseline = 300 cells each; doses U(0.6,1.4):

| setting | exact | equiv | top_k_cos | dose_mae | unexpl | flagged |
|---|---|---|---|---|---|---|
| 1 stimulus, same cell type | 1.00 | 1.00 | 1.00 | 0.05 | 0.13 | 0.00 |
| 2 stimuli, same cell type | 0.97 | 0.97 | 0.70 | 0.06 | 0.12 | 0.02 |
| 3 stimuli, same cell type | 0.90 | 0.92 | 0.33 | 0.06 | 0.12 | 0.02 |
| 2 stimuli, other lineage (CD4_T) | 0.97 | 0.98 | 0.53 | 0.20 | 0.18 | 0.00 |
| 2 stimuli + dissociation/2D artifact (naive-score Δ: 0.90) | 0.98 | 0.98 | 0.67 | 0.06 | 0.15 | 0.02 |
| 2 stimuli + synergy (non-additive) | 0.20 | 0.22 | 0.30 | 0.15 | 0.47 | 0.82 |

`exact` = recovered set equals truth; `equiv` = equal up to near-duplicate stimuli; `top_k_cos` = baseline that just takes the k best single matches (k given!); `unexpl` = |residual|/|ΔS|; `flagged` = share with unexplained > 0.35 (reported as 'missing stimulus or non-additive').

## 4. Asymmetric partial-view query

150 test samples (200-cell pseudobulk, CD14_Mono, random programs/memory/artifacts). Same loadings W for every gene panel; only rows are sliced. Error bars = linear-Gaussian posterior (prior sd 1, single-cell noise/n + between-batch floor) × conformal factor κ fitted on a separate 150-sample calibration plate:

| panel | n_obs | n_identifiable | rmse | rmse_targets | kappa | coverage90 | frac_determined | tern_acc_determined | tern_acc_all | r95_targets | cond |
|---|---|---|---|---|---|---|---|---|---|---|---|
| full transcriptome | 2000 | 32 | 0.14 | 0.13 | 4.87 | 0.89 | 0.95 | 0.996 | 0.983 | 0.32 | 3.91 |
| targeted 800 | 800 | 32 | 0.14 | 0.14 | 4.85 | 0.89 | 0.95 | 0.996 | 0.981 | 0.32 | 4.13 |
| qPCR 12 (top-loaded) | 12 | 6 | 0.42 | 0.17 | 1.00 | 0.93 | 0.14 | 0.988 | 0.855 | 0.31 | inf |
| qPCR 12 (A-optimal design) | 12 | 4 | 0.40 | 0.18 | 1.00 | 0.95 | 0.12 | 0.984 | 0.864 | 0.15 | inf |
| 4 secreted cytokines | 4 | 4 | 0.40 | 0.17 | 1.00 | 1.00 | 0.01 | 1.000 | 0.869 | 1.45 | inf |

`n_identifiable` = programs whose posterior sd < 0.5 × prior sd; `rmse_targets` = the 4 programs the panel was built for (cytokine row: IFNg/TNF/IL6/IL10); `kappa` = error-bar inflation needed on the calibration plate (1 = model noise already honest); `coverage90` = share of TEST errors inside the 90% interval (target ≈0.90); `frac_determined` = share of program bits whose 90% interval sits in one ternary bin (only those bits are searched); `r95_targets` = longest semi-axis of the 95% confidence ellipsoid over the 4 target programs (programs the panel does not see fall back to the prior: radius √χ²₃₂ ≈ 6.8); `cond` = condition number of the Fisher matrix (inf = fewer observations than latent coordinates — the prior regularises, the mask reports it).

Which of the 34 library stimuli produced this sample? Library = full-transcriptome codes; query = panel posterior, undetermined bits masked; 3 samples per stimulus. 'covered' = the stimulus moves a program the panel can see. 'ambiguous' = >2 library entries tie at the best masked distance, so the system answers 'cannot tell' instead of guessing:

| panel | stimuli | correct | ambiguous | wrong |
|---|---|---|---|---|
| full transcriptome | covered (34) | 0.98 | 0.02 | 0.00 |
| qPCR 12 (A-optimal design) | covered (8) | 0.50 | 0.50 | 0.00 |
| qPCR 12 (A-optimal design) | not covered (26) | 0.00 | 1.00 | 0.00 |
| 4 secreted cytokines | covered (7) | 0.00 | 1.00 | 0.00 |
| 4 secreted cytokines | not covered (27) | 0.00 | 1.00 | 0.00 |

_Run time 128 s._
