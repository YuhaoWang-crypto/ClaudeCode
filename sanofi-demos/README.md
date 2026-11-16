# Sanofi-Public open-source models — install & demo

Reproducible install + working demos for three [Sanofi-Public](https://github.com/Sanofi-Public)
repositories, all verified end-to-end on CPU (4 cores, no GPU):

| Library | What it does | Demo result |
|---|---|---|
| [**SNAC-DB**](https://github.com/Sanofi-Public/SNAC-DB) | Extract & curate antibody / nanobody–antigen complexes from PDB/CIF structures | unit tests 3/3; 7 real PDBs → 6 complexes; Foldseek clusters the 3 coronavirus-RBD binders together |
| [**RiboNN**](https://github.com/Sanofi-Public/RiboNN) | Predict mRNA translation efficiency (TE) from sequence ([Nat. Biotech. 2025](https://www.nature.com/articles/s41587-025-02712-x)) | 30 transcripts × 78 human cell/tissue TE tasks |
| [**CodonBERT**](https://github.com/Sanofi-Public/CodonBERT) | mRNA language model for codon optimization / expression ([preprint](https://www.biorxiv.org/content/10.1101/2023.09.09.556981v1)) | 87M-param model; mRFP expression Pearson **0.66** vs 0.38 codon-freq baseline |

See [`REPORT.md`](REPORT.md) for the full write-up, numbers, and install gotchas.

## Run it

Each demo is a self-contained script that clones the upstream repo, builds an
isolated environment, downloads the pretrained weights, and runs. They do **not**
assume the other two — each makes its own env (the three libraries pin
conflicting dependencies, all on Python 3.10).

```bash
# SNAC-DB  — needs micromamba/conda (ANARCI, HMMER, Foldseek are conda-only)
bash snac-db/run_snacdb_demo.sh

# RiboNN   — needs uv (or python3.10); downloads ~200 MB of weights from Zenodo
bash ribonn/run_ribonn_demo.sh

# CodonBERT — needs uv (or python3.10); downloads ~950 MB model from Sanofi CDN
bash codonbert/run_codonbert_demo.sh
```

By default each script works inside its own `work/` subdir (git-ignored), so the
large clones, environments, and model weights stay out of version control.

## What's committed here

```
sanofi-demos/
├── README.md, REPORT.md          # this overview + full report
├── snac-db/run_snacdb_demo.sh
├── ribonn/run_ribonn_demo.sh
├── codonbert/run_codonbert_demo.sh
├── codonbert/codonbert_mrfp_demo.py   # frozen-embedding benchmark
└── outputs/                       # small example outputs (evidence)
    ├── snac-db/pdbs_curation_summary.csv
    ├── snac-db/antigen_cluster.tsv
    └── ribonn/prediction_output_human.csv
```

Upstream repos, pretrained weights, conda/venv environments, and bulk outputs
are intentionally **not** committed — the scripts regenerate them.

## Notes / caveats

- CodonBERT demo is a lightweight **frozen-feature probe** (embeddings + ridge),
  not the paper's fine-tuning — the paper's numbers are higher.
- Each library pins Python 3.10 with mutually incompatible deps → three separate
  environments.
- Install fixes applied: `snacdb-patch-anarci` (tolerate `X` residues);
  RiboNN needs `uv --index-strategy unsafe-best-match` for `lightning-utilities`.
- All demos were run on CPU; a GPU makes CodonBERT embedding and RiboNN much faster.
