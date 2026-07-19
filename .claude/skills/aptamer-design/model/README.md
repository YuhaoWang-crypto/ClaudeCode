# Aptamer prediction model — data integration + runnable baseline

A **Stage-0 pre-filter** for the aptamer pipeline: fast, offline, ML-based scoring/generation
that narrows candidates *before* the expensive Boltz-2.1 structural step (Stage 1).

```
  known aptamers (UTexas + Apta-Index + literature)
            │  ingest.py  → unified schema (schema.md)
            ▼
   ┌─────────────────────────────┐     ┌──────────────────────────┐
   │ predict_baseline.py         │     │ generate_lm.py           │
   │ score (aptamer,target) pair │     │ Markov nucleotide LM     │
   │ = P(bind) / similarity      │     │ → novel candidates       │
   └─────────────┬───────────────┘     └────────────┬─────────────┘
                 └──────────────┬────────────────────┘
                                ▼   top-k survivors
                     ViennaRNA fold filter  →  Boltz-2.1 co-fold (Stage 1)
                                ▼
                     specificity counter-screen → experimental SELEX
```

## Files
| file | role |
|---|---|
| `schema.md` | unified dataset schema + source access notes |
| `ingest.py` | normalize a UTexas/Aptagen/literature CSV → schema; make synthetic negatives |
| `featurize.py` | pure-Python k-mer (aptamer) + AA-composition (protein) features |
| `predict_baseline.py` | interaction scorer: LogisticRegression (if sklearn) else kNN cosine |
| `generate_lm.py` | order-k Markov nucleotide **language model** → generate novel candidates |
| `seed_dataset.csv` | 3 illustrative real aptamers (TBA, AS1411, Sgc8c) — **replace with UTexas export** |
| `queries.csv` | our EGFR + GFRα1 designed candidates, ready to score |

## Quickstart (no dependencies required)
```bash
cd model
python3 predict_baseline.py --train seed_dataset.csv --query queries.csv
python3 generate_lm.py --train seed_dataset.csv --k 3 --n 10 --chem DNA
```
`predict_baseline.py` auto-upgrades to LogisticRegression if scikit-learn is installed,
else runs the pure-Python kNN fallback. Everything runs on the 3-row seed; results only
become meaningful once you load the real corpus.

## Get the real training data
1. Download the **UTexas Aptamer Database** public dataset (~1,475 rows):
   https://sites.utexas.edu/aptamerdatabase/
2. `python3 ingest.py utexas_export.csv --preset utexas -o unified.csv --negatives`
   (adjust `--map` to the real headers).
3. Optionally back-fill `target_seq` from UniProt for the protein targets.
4. Re-run the baseline / LM on `unified.csv`.
Apta-Index (~800) has no bulk export and is commercial — use it to spot-check, not scrape.

## Scale-up path to a real model (beyond this baseline)
The baseline is deliberately swappable at two seams:
1. **Embeddings** — replace `featurize.py` with pretrained encoders:
   **RNA-FM** (aptamer) + **ESM-2 / ProtBert** (protein). These carry the "language model"
   prior that 2k training rows cannot teach from scratch.
2. **Head** — replace LogisticRegression with a **two-tower + cross-attention** binder head
   (regress `kd_nM`, classify binding), and a **target-conditioned generator** (fine-tuned
   RNA-FM or a small GPT) for de novo design. This mirrors published systems:
   - AptaTrans (BMC Bioinformatics 2023) — transformer encoders for aptamer–protein interaction
   - AptaBLE (NeurIPS 2024) — LLM that both predicts interactions and generates aptamers
   - SelfTrans-Ensemble (2025) — ProtBert + RNA-FM, ~88% test accuracy

## Honesty
- Only ~1.5–2k public aptamers with affinity exist → **transfer learning, not from-scratch**.
  Models extrapolate poorly to unseen target families; always validate structurally (Boltz)
  and experimentally (SELEX).
- Baseline scores are **relative pre-filter signals, not affinities**. The Markov LM has
  **no target conditioning** — it captures aptamer sequence statistics only.
