# Data access — the four perturbation sources (+ enrichment MCP servers)

How to actually obtain each dataset and coerce it into a `Signature`. All facts
below were checked against the sources; where a value is version-specific it is
flagged. Behind an agent proxy, prefer the MCP tools listed at the end for
programmatic access.

---

## 1. Drug perturbation — Broad Drug Repurposing Hub

**What it is.** A curated, annotated library of ~**7,400 unique compounds**
(FDA-approved + clinical-trial + pre-clinical tool compounds) hitting ~**2,250
protein targets**, annotated with **MOA, target, clinical phase, disease area,
and indication**. Metadata is CC-BY 4.0.

**Portal / download.** <https://repo-hub.broadinstitute.org/repurposing>
(→ `#download-data`). No API needed for the annotations — two flat files:

- **Drug information** (one row per compound): key columns
  `pert_iname`, `clinical_phase`, `moa`, `target`, `disease_area`, `indication`.
- **Sample information** (one row per physical well/sample): plate/well, dose,
  vendor, and the `broad_id` that joins to expression profiles.

**The expression signatures** come from the **Connectivity Map (CMap) L1000**
(the two are meant to be used together):

- **L1000 assay** measures **978 "landmark" genes** (Luminex beads) and infers
  ~**11,350** more → ~12,328-gene profiles.
- Scope (Subramanian et al., Cell 2017): ~**19,811 small-molecule** perturbagens
  **plus genetic** perturbations — **shRNA knockdown** (3 hairpins/gene) and
  **cDNA overexpression** — of ~**5,075 genes**, across many cell lines/doses/
  times; >1,000,000 profiles.
- Access processed signatures at **clue.io** (login) — level-5 **moderated
  z-score** matrices in **`.gctx`** with `sig_info` / `gene_info` / `pert_info`
  metadata. GEO mirrors: **GSE92742** (Phase I) and **GSE70138** (Phase II).
- Read `.gctx` with **`cmapPy`** (`from cmapPy.pandasGEXpress.parse import parse`).

**→ Signature.** Each column of the level-5 matrix is one perturbation's
z-scores over genes: `Signature.from_l1000(zscore_series, name=sig_id,
modality="drug", meta={...})`. clue.io's own **Query app** computes WTCS/τ
against the full reference; `perturbomics.connectivity` reproduces the local math.

---

## 2. CRISPR perturbation — Broad GPP (Genetic Perturbation Platform)

**Portal.** <https://portals.broadinstitute.org/gpp/public/>

**Design tools.**
- **CRISPick** — the sgRNA/library designer for **CRISPRko / CRISPRi / CRISPRa**
  (picks and ranks guides on-target/off-target).
- **Beagle** — base-editor sgRNA (tiling) designer.
- plus vector-design (Fragmid) and RNAi hairpin tools.

**Optimized genome-wide libraries** (Doench/Sanson et al., *Nat. Commun.* 2018;
distributed via Addgene):

| Library | Modality | Notes |
|---|---|---|
| **Brunello** | CRISPRko (human) | ~77,441 sgRNAs, ~4/gene + ~1,000 controls |
| **Brie** | CRISPRko (mouse) | mouse analog of Brunello |
| **Dolcetto** | CRISPRi | fewer guides/gene, matches KO at detecting essentials |
| **Calabrese** | CRISPRa | outperforms SAM at resistance-gene screens |

**Screen data → Signature.** A pooled screen yields per-gene guide-abundance
changes (treatment vs control), scored with **MAGeCK** or **DESeq2** →
`log2FC` and significance per gene. `Signature.from_deseq2(mageck_or_deseq_df,
stat_col="stat", modality="crispr_ko", ...)`. For **Perturb-seq** (single-cell
CRISPR + scRNA-seq) each guide's cells are one condition → run the pseudobulk
DGE of `reference/dge-signatures.md`, one signature per targeted gene.

**Orientation.** KO / CRISPRi = **loss** of function; CRISPRa / cDNA-OE =
**gain**. Keep the sign convention consistent with the drug signatures you
compare against (see SKILL "gotchas").

---

## 3. Single-cell perturbation DGE — sc-best-practices

<https://www.sc-best-practices.org/conditions/differential_gene_expression.html>
Full recipe in `reference/dge-signatures.md`. One line: **pseudobulk, then
DESeq2/edgeR** — never a per-cell Wilcoxon. `perturbomics.pseudobulk`
implements it.

---

## 4. Foundation model — Geneformer

- **Weights:** <https://huggingface.co/ctheodoris/Geneformer> (Apache-2.0).
  Checkpoints: **V1-10M**, **V2-104M**, **V2-104M_CLcancer**, **V2-316M**.
- **BioNeMo** packaging & sizes: <https://docs.nvidia.com/bionemo-framework/latest/models/geneformer/>
  (e.g. 10M = 6 layers/256-dim; 106M = 12 layers/768-dim; BERT encoder,
  rank-value input of the top ~1024–2048 expressed genes; pretrained on
  ~23–95M human single cells).
- **Load:** `AutoModelForMaskedLM.from_pretrained("ctheodoris/Geneformer")`;
  tokenizer maps Ensembl IDs → tokens using the packaged gene-median dictionary.
- **Use → Signature:** cell/gene embeddings for similarity, or **in-silico
  perturbation** (delete/activate a gene token, measure the embedding shift) →
  a per-gene shift vector. See `reference/geneformer.md`.

---

## 5. In-session MCP servers that enrich a hit (available in this environment)

Once connectivity nominates a compound or gene, these annotate it (all found via
`ToolSearch`; do not assume a repo is in scope without checking `list_repos`):

- **ChEMBL** (`mcp__ChEMBL__*`) — `compound_search`, `get_mechanism`,
  `get_bioactivity` (IC50/Ki), `target_search`, `get_admet`, `drug_search`.
  Turn a Repurposing-Hub `pert_iname` into quantitative target/activity data.
- **Clinical Trials** (`mcp__Clinical_Trials__*`) — is a nominated repurposing
  candidate already in trials for the indication? (`search_trials`,
  `search_by_sponsor`, `analyze_endpoints`).
- **PubMed** / **bioRxiv** — evidence for a drug↔gene↔disease link; check
  whether a predicted combination is already reported.
- Structure/PK servers (**Boltz**, **Inductive Bio**) — for a nominated small
  molecule + target, structure/binding and ADME estimates.

Pattern: **compute** the connectivity/combination hit locally (rigorous), then
**annotate & sanity-check** it with these servers before calling it a lead
(hypothesis). Keep the ✅/⚠️ split intact in the write-up.
