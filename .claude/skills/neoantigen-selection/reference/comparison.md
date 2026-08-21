# Compared with `design-neoantigen-vaccines` (Neoantigen_Vaccine_Package v1)

A second neoantigen-vaccine skill was reviewed alongside this one. Both were
run, both were read, and where the other package was better this one changed.
This document records what was compared, what was adopted, and what does not
hold up — with the checks that back each claim.

## Verdict in one line

The other package has the better **audit and governance layer** and shipped the
**better ground-truth dataset**; this package has the better **computation** —
it computes the features the other one asks the user to supply, and it actually
performs the junction screening the other one only reports a status for.

---

## What was adopted from it

**1. The TESLA mirror.** `assets/demo/tesla_deepimmuno_public.csv` — 522
peptide-HLA pairs, 35 experimentally immunogenic, 6 patients, with real T-cell
assay labels. This is exactly what `benchmark.py` could not obtain: its decoys
are *unlabelled*, so its AUCs are lower bounds and the true negative rate is
unknown. Here a 0 means **assayed and negative**. Now shipped as
`data/tesla_deepimmuno_public.csv` with `tesla.py` to score it.

Cross-check: this package's `average_precision` reproduces the other package's
published `cnn_regress` value (0.131464) to four decimals on the same rows, so
the two implementations measure the same quantity.

**2. An audit manifest** (`provenance.py`). Their `audit_manifest.json` records
input SHA-256s, tool versions, a capabilities snapshot, the weights in force,
and the predictor actually used. This package now writes the same thing, plus
the git revision, content hashes of the variant / prediction / selection tables,
and the reference-proteome hash.

**3. Evidence levels on presentation numbers.** Their E0/E2 scheme is sharper
than a general "unverified" label for the specific question *what produced this
affinity*. Adopted as E0–E3 in `provenance.py`. This package has no E0 mode by
design: it ships no hash-based stand-in predictor, because a number that looks
like an affinity but is a hash eventually gets quoted as an affinity.

**4. Their guardrail list is good and mostly agreed with** — evaluate split by
patient, report PR-AUC and top-N rather than accuracy on imbalanced data, never
call a deterministic hash an AI prediction, never impute missing expression
silently, do not emit an mRNA as CMC-ready. `tesla.evaluate()` reports per
patient and by PR-AUC because of it.

---

## Where this package computes what the other one asks for

| | `design-neoantigen-vaccines` | this package |
|---|---|---|
| `recognition_score` | **required input column** | computed — Łuksza-style alignment to IEDB positive-assay epitopes |
| `normal_similarity` | **required input column** | computed — BLOSUM distance from the self peptide, plus an exact self-proteome k-mer novelty gate |
| agretopicity | not computed; `wt_context` is validated but the wild-type peptide is never predicted | computed exactly — every mutant window carries its positionally matched WT window, both are predicted |
| junction epitopes | counted, never screened | 34×34 cost matrix over all junction peptides, greedy + 2-opt ordering, then a rescan over lengths 8–11 |
| protein context | supplied per row in the input | reconstructed from the UniProt reference proteome, with a reference-mismatch check that reports isoform disagreements instead of applying them |

The first two matter most. `recognition_score` and `normal_similarity` are the
two hardest quantities in neoantigen selection, and a pipeline that takes them
as inputs has moved the difficulty to the user rather than solving it. Their
own demo fixture fills both columns with fixture values.

---

## Findings that did not hold up

Each was verified against the shipped code, not inferred from the docs.

**Junction screening is documented but not implemented.** `SKILL.md` step 7 says
"Screen linker junctions with the same real HLA predictor when possible;
otherwise mark junction safety as unresolved." `build_construct()` makes no
prediction call at all. It sets

```python
"junction_screen_status": "unresolved" if evidence == "E0"
                          else "predictor_available_not_clinically_validated"
```

so with a real predictor installed the report reads
`predictor_available_not_clinically_validated`, which sounds like a screen was
run. It means only that a predictor exists. In this package's own demo, 38 of
1,122 junction peptides bind at ≤0.5% rank *after* optimization — junctions are
not a formality.

**Class I and class II are merged, against the package's own guardrail.** The
guardrails say "Keep HLA-I and HLA-II evidence separate through scoring and
reporting" and "Never mix percentile rank, IC50, probability, and model score
without explicit normalization". The scorer does:

```python
class1 = 1.0 / (1.0 + affinity / 500.0)     # from MHCflurry IC50 in nM
class2 = 1.0 / (1.0 + c2rank / 2.0)         # from a class-II percentile rank
presentation = max(class1, class2)
```

An IC50-derived score and a rank-derived score are combined by `max`, and one
number then represents both classes.

**Class-I scoring uses raw nM, not %rank.** Affinity distributions differ
substantially between alleles, which is why percentile rank exists. Ranking
candidates across six alleles on raw nM systematically favours the alleles whose
binders are tighter in absolute terms. This package uses NetMHCpan-4.1 **EL**
%rank throughout.

**The benchmark cannot evaluate the package's own model.** `benchmark()` scores
only the pre-computed columns already present in the TESLA CSV
(`rf_classify`, `cnn_regress`, the DeepImmuno score, …). It never runs the
package's own composite on those peptides. So the shipped benchmark, however
well constructed, structurally cannot answer "does this package's ranking
work?". Running that missing evaluation is what `tesla.py` now does — for both
packages' scores, on identical rows.

**Codon "optimization" is one fixed codon per amino acid.** Every alanine is
`GCC`, every leucine `CTG`. The shipped construct comes out at **GC 0.6816**,
and the QC block reports that number without flagging it, has no repair loop,
and checks no restriction sites, cryptic polyadenylation signals or splice
motifs. This package's optimizer draws codons by human usage frequency and then
repairs GC windows, homopolymers, restriction sites and `AATAAA` until clean
(demo construct: GC 52.0%, QC pass).

**Track B runs on fabricated sequence.** The 48-row fixture is labelled
`synthetic_functional_fixture` and the patient JSON says "synthetic functional
fixture; not patient data", so this is disclosed, not hidden. But it is worth
being concrete about what it means: **0 of 48 `wt_context` values occur anywhere
in the reviewed human proteome**, including the one labelled `BRAF p.K101I`. The
functional track therefore exercises the plumbing, not the mutation-to-protein
step. This package's demo runs on real TCGA-SKCM variants against real UniProt
sequences, and reports the 4 variants whose reference residue disagreed with the
proteome rather than applying them.

---

## Head to head on the TESLA mirror

Both scoring approaches, plus the mirror's published columns, on identical rows.
Random baseline AP = 0.067.

| score | AP | AUC | positives in a 34-slot budget |
|---|---|---|---|
| **NetMHCpan-4.1 EL %rank alone** | **0.207** | 0.791 | **31 / 35** |
| this package's composite | 0.163 | 0.763 | 26 / 35 |
| `cnn_regress` (best published column) | 0.132 | 0.654 | 19 / 35 |
| `rf_regress` | 0.108 | 0.619 | 19 / 35 |
| DeepImmuno `immunogenic score` | 0.083 | 0.477 | 13 / 35 |
| `IEDB` immunogenicity score | 0.070 | 0.523 | 17 / 35 |

The other package's composite could not be scored here — it requires
`recognition_score` and `normal_similarity` per row, which the mirror does not
carry, and which is the point made above.

The result that matters is not which package wins. It is that **a current
presentation predictor beats every published immunogenicity model on this
dataset**, and beats this package's own composite. That is why the default
weights are now presentation-dominant.

---

## Where the two agree, independently

Convergent design is worth noting, since it suggests these choices are forced by
the problem rather than by taste:

* 25-mer minigenes with the mutation centred
* ≤34 antigens, per-gene caps, an HLA-coverage bonus, greedy portfolio selection
* expression, clonality and self-similarity as separate axes rather than one blob
* an explicit refusal to claim reproduction of the V940 algorithm or construct
* research-use framing, with the missing CMC elements named

## If you were merging them

Take the other package's manifest discipline, evidence levels, input-schema
validation and RUO/CMC guardrails (this package has now taken the first three);
take this package's computed features, agretopicity, self-proteome novelty gate,
real junction screening, %rank-based scoring, codon repair loop, and the TESLA
evaluation of the model actually being shipped.
