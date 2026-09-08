---
name: af2bind-binding-site
description: >-
  Predict small-molecule binding-site residues on a protein with AF2BIND, which
  reads the binding signal out of AlphaFold2's pair representation using 20
  "bait" amino acids. Use when asked where a ligand, drug or small molecule
  could bind a protein; to find or rank pockets on an apo, designed or
  AlphaFold-predicted structure with no known holo homolog and no MSA of
  binding data; to pick a docking box or a site to engineer; or to check
  whether a designed binder's intended pocket is the one the model sees.
  Runs the AlphaFold2 pass on a Modal GPU and scores locally with numpy.
  Not for predicting which ligand binds, binding affinity, protein-protein
  interfaces, or docking poses.
---

# AF2BIND: small-molecule binding sites from the AlphaFold2 pair representation

Implements Gazizov, Lian, Goverde, Mou, Ovchinnikov & Polizzi, Nature Methods
(2026), doi:10.1038/s41592-026-03011-2. Working package: `af2bind_pipeline/`.

## The question this answers

Given **one protein structure and nothing else** — no ligand, no holo homolog,
no MSA of known binders — which residues line a small-molecule binding site?

That "nothing else" is the point. Pocket finders that rely on a homologous
complex fail exactly where you most want an answer: de novo designs, orphan
folds, AlphaFold models of proteins nobody has crystallised with a ligand.

## How it works, in one paragraph

Append 20 bait residues to the target, one of each amino acid, each as its own
isolated single-residue chain. Run one AlphaFold2 forward pass. The pair
representation block coupling target residue *i* to bait *a* encodes how much the
network wants an amino acid of type *a* next to residue *i*; binding-site
residues want company. A logistic regression over those 5120 numbers per residue
gives p(bind). The representation does the work, the classifier is trivial. Full
detail in `reference/methodology.md`.

## Quick start

```bash
pip install numpy 'modal[api-proxy-support]'

# GPU pass on Modal, scoring + analysis locally
python -m af2bind_pipeline.run --target 6w70 --chain A --out results/6w70

# offline checks, no GPU needed
python -m af2bind_pipeline.selftest
```

`--target` takes a local PDB file, a 4-character PDB ID, or a UniProt accession
(pulled from the AlphaFold database). Setup, GPU sizing, cost and large-protein
handling are in `reference/running.md`.

## The two-stage split, and why it matters

```
structure --[GPU: one AlphaFold2 pass]--> features (L, 5120) --[numpy]--> p(bind)
```

The GPU stage is the only expensive part and it runs once. Everything after it is
a standardise-and-dot-product that needs neither jax nor a GPU. So `features.npz`
is written on every run, and re-scoring with a different seed, ensemble or
threshold is free:

```bash
python -m af2bind_pipeline.run --features results/6w70/features.npz \
    --target results/6w70/6W70.pdb --seeds 0,1,2,3,4,5,6,7,8,9 --out results/6w70_ens
```

Never re-run the GPU pass to change a downstream parameter.

## Reading the output

Read the **ranking**, not the absolute number. p(bind) is a logistic score fitted
on one training distribution; it is not calibrated across targets, and its
absolute value shifts with protein size, structure quality and which AlphaFold
parameter release you used. See `reference/interpretation.md` for what the
numbers do and do not support.

Operationally:

- **Top ~15 residues** is the useful unit, and is what the upstream notebook
  prints. On the paper's own example that set is essentially all true pocket.
- **Pockets, not residues**, is usually the decision you need. `report.json`
  clusters the top-scoring residues into spatially distinct sites, ranked, each
  with a PyMOL selection string. This clustering is added here, not part of the
  published method.
- **A flat profile means no answer.** If the top score sits near the bulk of the
  distribution, the honest report is "no confident site", not the top residue.

## Validate before you believe it

Every run on a structure containing a real ligand automatically scores itself
against the residues within 5 Å of that ligand, and writes ROC-AUC, average
precision and precision@N into `report.json`. Run a holo structure first
whenever you are working in a new setting. A pipeline that cannot recover a known
site is not one whose novel predictions are worth acting on.

This implementation was checked that way on two targets:

| target | ROC-AUC | precision@10 | mean p(bind) |
|---|---|---|---|
| 6W70 chain A, designed apixaban binder (the paper's example) | 0.978 | 1.00 | 0.198 |
| 1M17 chain A, EGFR kinase with erlotinib | 0.934 | 0.70 | 0.087 |

On EGFR the top residues are the canonical ATP pocket (L694, V702, K721, T766,
M769, L820, T830, D831), recovered with no ligand shown to the model. Note the
mean score more than halves between the two targets: that is the reason to read
ranks and not absolute values. Full numbers, the figures, and what they do and
do not license are in `reference/interpretation.md`.

Ground truth excludes waters and a list of crystallisation additives, and keeps
only ligands whose nearest protein chain is the chain being scored, so a second
copy in the asymmetric unit cannot manufacture false positives.

## Gotchas that silently produce wrong numbers

None of these raise an error. All of them invalidate the result.

- **The bait residue-index offset.** The 20 baits must be pushed 50 indices apart
  so AlphaFold2 treats them as 20 separate chains, not one 20-residue peptide.
  Handled in `modal_app.pair_features`; do not remove it.
- **Mask setting and head must match.** `mask_sidechains=True` goes with the
  `_nosc` head, `False` with the plain one. Mixing them yields plausible nonsense.
- **AlphaFold parameter release.** The head was fitted to features from one
  checkpoint. Upstream is inconsistent: the notebook downloads 2021-07-14, the
  README says the experiments used 2022-03-02. This pipeline defaults to
  2021-07-14, the release shipped with the trained heads. Changing
  `--af2-params` changes the numbers; do not compare across releases.
- **Disorder on predicted models.** Disordered AlphaFold tails are unpacked and
  exposed, which is what the score rewards. Trim with `--min-plddt 70`, and split
  large multi-domain proteins by domain rather than scoring them whole.
- **Modal behind a proxy.** "Could not connect to the Modal server" while curl
  works means a missing optional dependency, not a network block. Install
  `modal[api-proxy-support]`.

## What this does not do

It finds *where* a small molecule could bind. It does not say **which** ligand
binds, predict affinity, produce a pose, score druggability, or find
protein-protein interfaces. Chain it into docking or design; do not read it as an
answer to those questions.

## Where to look

- **The method, precisely** (shapes, the bait trick, the two heads, the 10
  folds) → `reference/methodology.md`
- **Setup, inputs, GPU sizing, cost, outputs** → `reference/running.md`
- **Thresholds, calibration, failure modes, measured numbers**
  → `reference/interpretation.md`
- **Scoring many targets in one session** → `assets/batch_targets.py`
- **Code** → `af2bind_pipeline/` (`core` scoring math, `modal_app` GPU pass,
  `structure` parsing and ground truth, `pockets` clustering, `validate`
  metrics, `run` CLI, `plot` optional figures, `selftest` offline checks)

```bash
# profile + bait-activation figure for a finished run (needs matplotlib)
python -m af2bind_pipeline.plot results/6w70 --out figures/af2bind_6w70.png
```
