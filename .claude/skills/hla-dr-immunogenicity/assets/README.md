# Installing the pipeline into a project

The scripts resolve every path relative to **the parent of the `scripts/`
directory**, so the layout matters:

```
<project>/
  scripts/        <- assets/scripts/*
  config/         <- assets/config/config.yaml
  data/           <- sequences.fasta, sequences_metadata.tsv, human_sprot.fasta
  results/        <- created by the modules
  figures/        <- created by make_figures.py
  package.json    <- assets/package.json  (deck only)
```

```bash
mkdir -p myproject && cd myproject
cp -r <skill>/assets/scripts <skill>/assets/config .
cp <skill>/assets/package.json .
pip install pyyaml biopython pandas numpy matplotlib python-pptx
npm install                                    # pptxgenjs, deck only
bash <skill>/assets/run_all.sh
```

## Running on your own sequence

`m0_fetch_sequences.py` builds the demo batch — a public stand-in test article
plus every benchmark and control — by fetching each sequence live from
RCSB/UniProt and locating regions by motif. For a proprietary ligand:

1. Run M0 once to get the benchmark and control set.
2. Append your de-identified sequence to `data/sequences.fasta` and add a row
   to `data/sequences_metadata.tsv` with `role = test`.
3. Point `benchmarks.anchor_low` in `config/config.yaml` at the comparator you
   want the fold-change expressed against (default `ProteinA_Z`).
4. Re-run M1 onward. Nothing else changes.

Do not run the ligand without the benchmarks and controls — pIRS is a relative
scale and M6's system-suitability verdict is what makes a batch reportable.

## What to re-run when

| Change | Re-run |
|---|---|
| New test article, same panel and thresholds | M1–M9, then figures/report/deck |
| Panel widened or a new allele added | M3 (resumable — only the new allele is fetched), M4–M9, and the M10–M13 calibration layer |
| A threshold, gate or weight changed in `config.yaml` | M3 onward *and* M10–M13, because the calibration is a property of the rule |
| Nothing but report wording | `make_figures.py`, `make_report.py`, `make_deck.py`, `check_deck.py` |

M10–M13 are the slow ones (~10 min, 1–2 h, minutes, ~2 h). All are resumable
and cache as they go; the caches are regenerable and belong in `.gitignore`:

```
results/m11_scores_partial.tsv
results/m13_panel_scores.tsv
```

## The human proteome

Fetched once, ~20,400 reviewed entries, needed by M4:

```bash
curl -L 'https://rest.uniprot.org/uniprotkb/stream?query=reviewed:true+AND+organism_id:9606&format=fasta&compressed=true' \
  | gunzip > data/human_sprot.fasta
```
