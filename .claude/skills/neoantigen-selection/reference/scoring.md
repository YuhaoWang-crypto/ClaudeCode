# The score: every feature, its formula, its source, and how it fails

Nothing here is a vendor formula. Each feature encodes one published
observation about what separates an immunogenic neoantigen from a peptide that
merely exists. Weights are defaults to argue with, not a fitted model — refit
them with `benchmark.py` or `tesla.py` if you have labelled data.

> **The weights changed once already, and the reason is on the record.**
> They started literature-balanced (presentation 0.30, with 0.40 spread across
> agretopicity / dissimilarity / TCR prior / hydrophobicity). Two independent
> benchmarks — the presentation-controlled IEDB one and the TESLA mirror, which
> has real T-cell assay labels and real negatives — put those four features at
> or below the random baseline and scored the diluted composite *below the
> binding predictor alone*. So presentation went to 0.45 and the four dropped.
> `config.LITERATURE_BALANCED` still holds the old set so the change is
> reproducible; `config.PRESENTATION_ONLY` holds the setting that actually
> scored highest on TESLA. See `reference/benchmark.md`.

Two layers, deliberately separate:

* **gates** — binary, biological, non-negotiable. A candidate that fails one is
  reported with the gate that killed it.
* **score** — a weighted preference among whatever survived.

Collapsing the two (e.g. giving expression a weight instead of a gate) is the
most common way a ranking ends up recommending an unexpressed gene.

---

## Gates (`config.Gates`)

| gate | default | why |
|---|---|---|
| `min_tpm` | 1.0 | An unexpressed gene produces no peptide, whatever NetMHCpan says. This is the step the public workflow description calls out explicitly ("tumor RNA expression data"). |
| `max_rank_mhc1` | 2.0 | NetMHCpan-4.1 weak-binder cutoff. Above it, presentation is not credible. |
| `max_rank_mhc2` | 10.0 | NetMHCIIpan weak-binder cutoff (class-II thresholds are looser by construction). |
| `min_dna_vaf` | 0.05 | Below ~5% the variant is as likely a caller artifact as a subclone. |
| `min_ccf` | 0.0 (off) | Turn on to require clonality. Subclonal neoantigens are present in only part of the tumor (McGranahan 2016, *Science* 351:1463). |
| `require_novel_vs_self` | True | The mutant k-mer must not occur anywhere in the reference proteome. Catches mutations that recreate a peptide already present in another human protein — those are self, and targeting them is an autoimmunity risk, not a neoantigen. |
| `drop_anchor_only` | False | Optional: drop peptides whose only change is at an MHC anchor. Off by default because anchor mutations are genuinely immunogenic (Duan 2014). |

---

## Features (`features.py`), all mapped to [0, 1], higher = better

### `presentation` — weight 0.45
`1 / (1 + (rank/0.5)^1.5)` on the NetMHCpan-4.1 **eluted-ligand** %rank.
0.5% → 0.5, 0.05% → 0.97, 2% → 0.11.

EL (mass-spec-trained) rather than BA (affinity-trained): EL models the whole
presentation pathway, and outperforms affinity for identifying real ligands
(Reynisson 2020, *NAR* 48:W449). Highest weight because it is the best-validated
link in the chain — and still not sufficient on its own, which is exactly what
the benchmark shows.

**Breaks when:** the HLA type is wrong. Every downstream number inherits that error.

### `agretopicity` — weight 0.08
`clip((log10(WT rank / MUT rank) + 2) / 4, 0, 1)`.

Ratio of self-peptide binding to mutant binding — "differential agretopicity"
(Duan 2014, *J Exp Med* 211:2231; Ghorani 2018, *Ann Oncol* 29:271). A mutation
that *creates* binding presents something the thymus never negatively selected
against. A mutant that binds no better than its wild-type counterpart is
competing against a tolerized repertoire.

Returns 0.5 (neutral) when either rank is missing, so a missing wild-type
never masquerades as evidence.

**Breaks when:** the peptide has no wild-type counterpart (neo-ORF, frameshift) —
handled as neutral rather than as a maximum.

### `expression` — weight 0.18
`log10(TPM+1) / log10(101)`, clipped. Saturates at 100 TPM: the difference
between 100 and 1000 TPM does not predict immunogenicity, the difference between
0.5 and 50 does.

Best used with transcript-level TPM from the patient's own tumor RNA-seq, and
ideally with the mutant-allele fraction *in RNA* (allele-specific expression),
which catches mutations silenced by nonsense-mediated decay or allelic
imbalance. Gene-level RSEM is a coarser stand-in.

### `clonality` — weight 0.12
CCF, clipped to [0,1], from
`VAF x (purity x CN_tumor + 2 x (1-purity)) / (purity x multiplicity)`.

Clonal neoantigens are in every tumor cell; subclonal ones are an escape route.
Clonal neoantigen burden, not total burden, tracks checkpoint-inhibitor benefit
(McGranahan 2016).

**Breaks when:** purity and copy number are guessed. With CN=2 and
multiplicity=1 it degenerates to `2·VAF/purity` — usable for ranking, not for
claiming a CCF value. The demo says so explicitly.

### `dissimilarity` — weight 0.03
From the BLOSUM62 score of the substituted position(s):
`clip((3 - mean_blosum) / 7, 0, 1)`.

A conservative substitution (I→V, BLOSUM +3) leaves a surface the repertoire is
tolerant to; a radical one (G→W, BLOSUM −2) does not. Dissimilarity-to-self
predicts neoantigen immunogenicity (Richman 2019, *Cell Systems* 9:375).

Neo-ORF peptides with no self counterpart get 0.75 — foreign by construction,
but unverified, so not the maximum.

### `tcr_prior` — weight 0.07
Łuksza-style `R = Z/(1+Z)`, `Z = Σ_j exp(-k(a - s_j))`, published constants
`a = 26`, `k = 4.87` (Łuksza 2017, *Nature* 551:517), over IEDB epitopes with a
positive human T-cell assay.

**Approximation, stated plainly:** `s_j` here is an *ungapped* BLOSUM62 score
against same-length reference epitopes, not the Smith-Waterman score of the
original. For equal-length peptides the two agree closely; the approximation
buys a ~100x speedup, which is what makes the feature affordable over
10^4 candidates.

**Behaves as a near-binary flag, by design.** With the published constants the
exponential is extremely sharp: a candidate that is essentially a known
immunogenic epitope scores ~1, everything else scores ~0. On the demo run only
2 of 103 ranked candidates scored above 0.01. That is the intended semantics —
"this peptide looks like something a human T-cell repertoire has demonstrably
responded to" is a rare and strong statement, not a graded one — but it means
the feature contributes nothing to most rankings. If you want a graded version,
lower `LUKSZA_A` and say in the report that you did, because it is no longer the
published score.

**Breaks when:** benchmarked against IEDB itself — that is circular. The
benchmark drops exact self-matches and reports the composite score with and
without this feature.

### `hydrophobicity` — weight 0.02
Mean Kyte-Doolittle over the TCR-facing residues (positions 3..L-1),
rescaled to [0,1].

Hydrophobicity at TCR contact residues associates with immunogenicity
(Chowell 2015, *PNAS* 112:E1754) and was one of the few features that survived
the prospective TESLA analysis (Wells 2020, *Cell* 183:818). Small weight
because the effect is real but modest.

### `mhc2_support` — weight 0.05
1.0 if the same variant also yields a class-II binder at ≤2% rank, 0.5 at ≤10%,
else 0.

CD4 help matters: the first mRNA neoantigen vaccine responses were
predominantly CD4 (Kreiter 2015, *Nature* 520:692; Sahin 2017, *Nature*
547:222), and the 25-mer minigene format exists partly to preserve class-II
epitopes around the mutation.

### `anchor_penalty` — modifier, not a weight
Score is multiplied by 0.85 when the mutation sits only at an MHC anchor
(P1/P2/PΩ). Such a mutation creates binding but leaves the TCR-facing surface
identical to self. Down-weighted, not excluded — Duan 2014 show both classes
are immunogenic.

---

## Composite

`score = Σ wᵢ · featureᵢ` over the normalized weights, times the anchor
modifier. Linear and additive on purpose: the point is that you can read the
per-feature columns in `ranked.csv` and see exactly why a candidate is where it
is. A gradient-boosted model would score better on a benchmark and be useless in
the conversation where a clinician asks why slot 7 is slot 7.

Then `best_per_variant` collapses peptide × allele to one representative
epitope per **mutation** — because the vaccine encodes a minigene per mutation,
not per epitope — carrying `n_epitopes` / `n_alleles` / `n_strong` as support.

---

## Selection is not sorting (`select.py`)

Greedy with a submodular coverage bonus:
`effective = score + λ·(new allele covered) + 0.05·(clonal)`, λ = 0.15,
subject to a per-gene cap, optional per-allele cap, epitope de-duplication and
forced driver inclusions. Deterministic, and every slot records
`why_selected`.

The allele-spread term is a hedge: a tumor that loses one HLA allele — a
documented immune-escape route (McGranahan 2017, *Cell* 171:1259) — should not
silence the entire payload.

---

## What would make these numbers better

1. Real four-digit HLA typing (this is the single largest error source).
2. Transcript-level TPM **and** RNA mutant-allele fraction from the patient's tumor.
3. Segment-level copy number and an ABSOLUTE/FACETS purity call.
4. MHC-II prediction actually run (needs class-II typing).
5. Weights refit on labelled outcome data instead of literature defaults.
6. Mass-spec immunopeptidomics on the patient's own tumor, which converts
   `[unverified]` presentation into `[computed]` presentation.
