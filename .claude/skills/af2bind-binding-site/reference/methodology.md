# AF2BIND: how the method actually works

Source: Gazizov, Lian, Goverde, Mou, Ovchinnikov & Polizzi, *AF2BIND: predicting
small-molecule binding sites using the pair representation of AlphaFold2*,
Nature Methods (2026), doi:10.1038/s41592-026-03011-2. Preprint:
doi:10.1101/2023.10.15.562410. Reference code (MIT):
https://github.com/sokrypton/af2bind

## The central idea

AlphaFold2's Evoformer maintains a **pair representation** — a tensor of shape
(L, L, 128) whose entry (i, j) encodes what the network believes about the
relationship between residues i and j. AF2BIND's claim is that this tensor
already contains a usable small-molecule binding signal, and that you can read it
out without ever showing the network a ligand.

The read-out trick is the **bait**. Append 20 extra residues to the target, one
of each amino acid type, each as its own isolated single-residue chain. Run one
AlphaFold2 forward pass. The pair-representation block coupling target residue i
to bait residue a says, loosely, "how much does the network want an amino acid of
type a to sit next to residue i". Binding-site residues want company. A small
logistic regression on those coupling vectors turns that into p(bind).

## The computation, precisely

For a target of length L:

| Step | Object | Shape |
|---|---|---|
| AF2 forward pass on target + 20 baits | pair representation | (L+20, L+20, 128) |
| `pair_A` = target rows, bait columns | target to bait | (L, 20, 128) |
| `pair_B` = bait rows, target columns, transposed | bait to target | (L, 20, 128) |
| flatten and concatenate | features `x` | (L, 5120) |
| standardise, then one linear layer | per-bait logits | (L, 20) |
| sum over baits, sigmoid | `p(bind)` | (L,) |

The head is genuinely just logistic regression: a 5120-vector of weights, a
scalar bias, and the per-feature mean/std used to standardise. That is the
paper's point — the representation does the work, not the classifier. It also
means the head costs nothing to evaluate, which is why this pipeline caches the
(L, 5120) features and treats scoring as a separate, free step.

## Why the residue-index offset matters

The reference implementation does this after `prep_inputs`:

```python
r_idx = af_model._inputs["residue_index"][-20] + (1 + np.arange(20)) * 50
af_model._inputs["residue_index"][-20:] = r_idx.flatten()
```

AlphaFold2 infers chain connectivity from residue-index adjacency. Without the
offset the 20 baits look like a single 20-residue peptide, and AF2 spends its
capacity folding that peptide instead of reporting each amino acid's independent
preference for each target position. Pushing the indices 50 apart makes them 20
mutually invisible single-residue chains. Omitting this step does not error; it
silently produces features the trained head was never fitted to.

## The two heads, and the mask that selects between them

`prep_inputs(rm_target_sc=...)` controls whether the target's side chains are
shown to AF2. Two heads were trained, one per setting:

| `mask_sidechains` | `rm_target_sc` | head |
|---|---|---|
| `True` (default) | side chains removed | `split_nosc_pair_A_split_nosc_pair_B_{seed}` |
| `False` | side chains kept | `split_pair_A_split_pair_B_{seed}` |

Pairing a head with the wrong mask setting is the single easiest way to get
plausible-looking nonsense: the features shift, no error is raised, and the
scores are simply wrong. `mask_sidechains=True` is the default because it is the
harder, more honest setting — the model cannot cheat off a pocket's existing
side-chain packing, which is what you want when scoring a designed or predicted
structure.

## The 10 folds

The released archive contains seeds 0 through 9 for every model type: ten
independently trained folds, not ten epochs. The upstream notebook uses seed 0
only. This pipeline defaults to seed 0 for exact reproduction and offers
`--seeds 0,1,...,9`, which averages logits before the sigmoid. Ensembling
narrows fold-to-fold variance; it does not make p(bind) better calibrated in an
absolute sense, and it is not what the paper reports.

## The bait-activation profile

Because the logit is a sum of 20 per-bait terms, it decomposes for free. The
paper shows that the pattern across those 20 terms correlates with the chemistry
of the ligand a pocket actually binds — a hydrophobic pocket lights up the
hydrophobic baits. `bait_activations.csv` writes this matrix out in the paper's
BLOSUM-style column order, which groups chemically similar amino acids together
so the pattern is legible by eye.

Treat this as a qualitative hint about pocket chemistry. It is a correlation
reported in the paper, not a calibrated ligand-property predictor.

## What the method does not do

- It does not identify **which** ligand binds, only where a small molecule could.
- It does not give a pocket volume, a druggability score, or a binding affinity.
- It does not model the ligand, so it cannot rank poses or guide chemistry.
- It scores residues independently; the grouping into distinct pockets in this
  pipeline is post-processing added here, not part of the published method.

## Other feature sets in the released archive

The weights zip also ships heads trained on ESM-2 embeddings, ESM-IF embeddings,
the AF2 MSA first row, a serial (non-split) bait arrangement, and combinations of
these. The paper's finding is that the AF2 pair representation beats the
alternatives on its own; the combination heads are there for reproducing that
comparison. This pipeline implements the pair-only heads, which are the method.
Running the others would require producing the matching ESM features, which the
GPU step here does not do.
