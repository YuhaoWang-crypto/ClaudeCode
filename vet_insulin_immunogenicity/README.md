# Veterinary insulin immunogenicity — canine (DLA) and feline (FLA) workflows

Two MHC-II immunogenicity assessment workflows for therapeutic insulins given to
dogs and cats, plus the validation harness that says what each layer is and is
not worth.

**Terminology first, because it changes what tooling applies.** HLA is the human
locus name only. The canine system is **DLA** (Dog Leukocyte Antigen) and the
feline one is **FLA** (Feline Leukocyte Antigen); class II means DRB / DQA /
DQB orthologues. IEDB's MHC-II endpoints accept HLA allele names only, so a
human pipeline does not transfer by swapping an allele string — it has to be
rebuilt around pan-specific prediction from MHC sequence.

```bash
pip install -r requirements.txt
python demo.py                       # both species + head-to-head comparison
python -m vetimmuno run --species dog
python -m vetimmuno run --species cat --backend netmhciipan   # needs a licensed install
python tests/test_vetimmuno.py
```

Outputs land in `results/<species>/`: `report.md`, `foreign_burden.csv`,
`panel_summary.csv`, `applicability_domain.csv`, `non_self_cores.csv`,
`validation.csv`, figures, and ready-to-run NetMHCIIpan inputs under
`netmhciipan/`. `results/comparison.md` is the dog-vs-cat readout.

## What the workflow does, layer by layer

| # | Layer | Standing |
|---|---|---|
| 1 | Mature A/B chains from UniProt; positional difference map vs the recipient's own insulin | **Rigorous.** Chain boundaries come from each entry's own `Peptide` features, never hardcoded. |
| 2 | 9-mer register enumeration; drop every core the recipient's own insulin already contains | **Rigorous.** Set arithmetic. Exactly as reliable for a cat as for a human. |
| 3 | Species MHC-II panel (IPD-MHC for dog, GenBank for cat) reduced to groove-contact pseudosequences | **Measured.** Panel size and curation status are reported, not assumed. |
| 4 | Applicability-domain test against NetMHCIIpan's training space | **Measured.** Calibrated on the training set's own leave-one-out nearest-neighbour distribution. |
| 5 | MHC-II binding prediction | **NetMHCIIpan-4.3** when a licensed binary is available; otherwise a clearly labelled **illustrative** surrogate that predicts nothing. |
| 6 | Validation harness: 16 known-answer tests and controls | **Measured.** Runs on every execution and is printed in the report. |

## The three things this workflow was built to get right

### 1. The tolerance filter, not the binding score, is the defensible output

A dog is centrally tolerant to canine insulin. The only cores that can present a
novel signal are those its own insulin does not contain, so the risk list is
`cores(drug) \ cores(self)` — no model, no training data, no extrapolation. The
harness asserts an exact invariant on it (KAT-4): the non-self core set equals
the set of cores spanning a sequence difference, with any mismatch named rather
than absorbed into a threshold.

This layer produces the results that survive review. In the demo run:

| recipient | product | foreign residues | non-self cores |
|---|---|---|---|
| dog | porcine insulin | 0 | 0 |
| dog | human insulin | 1 (B30) | 1 |
| dog | bovine insulin | 2 (A8, A10) | 10 |
| cat | bovine insulin | 1 (A18) | 4 |
| cat | porcine insulin | 3 (A8, A10, A18) | 13 |
| cat | human insulin | 4 (A8, A10, A18, B30) | 14 |

Note what the core counts say that the residue counts do not: the dog's single
difference from human insulin is the last residue of the B chain and spans one
register, while bovine insulin's two differences sit inside the A chain and span
ten. Position dominates count.

### 2. The `%Rank` gap for custom molecules is closed, not worked around

NetMHCIIpan emits no `%Rank` in `-mhcfsa` custom-molecule mode, because rank
reference distributions exist only for its built-in alleles. `predict.BackgroundRank`
rebuilds one the way the tool does: score 20 000 background peptides with the
same molecule, then express any peptide as its percentile against that
molecule's own distribution. Backend-agnostic, so it works identically for the
surrogate and for a licensed NetMHCIIpan run, and KAT-7 verifies the calibration
is actually uniform on a held-out draw.

### 3. Cross-species extrapolation is measured, not asserted

NetMHCIIpan-4.3 is a pan-specific model trained on human HLA-DR/DQ/DP, mouse
H-2 and bovine BoLA-DRB3. It will happily score a DLA or FLA molecule. The
question is whether it should.

`groove.py` answers it quantitatively: reduce every molecule to its
groove-contact residues (the representation the model itself uses), find each
panel molecule's nearest neighbour in the training space, and read that distance
against the training set's *own* leave-one-out nearest-neighbour distribution.
The threshold is therefore derived from the model's data, not invented.

Demo result, DRB locus:

| | training molecules | training NN identity (median / 5th pct) | panel NN identity (median) | out-of-domain |
|---|---|---|---|---|
| dog (DLA-DRB1) | 1718 | 0.96 / 0.88 | 0.72 | 19/24 |
| cat (FLA-DRB) | 1718 | 0.96 / 0.88 | 0.61 | 22/22 |

Every feline molecule sits further from the training data than essentially any
molecule the model was trained on. The canine panel is mostly outside it too,
but has members close enough to call marginal. This is the ceiling on how much
weight a binding score can carry, and the report states it before showing any
score.

## What the validation harness actually validates

It cannot validate immunogenicity prediction accuracy — no canine or feline
benchmark exists, which is itself the finding. What it does validate:

* **KAT-1** — five sequence facts recomputed from freshly fetched UniProt
  records (canine insulin differs from human only at B30; feline at A8/A10/A18/B30;
  porcine is identical to canine; bovine differs from feline at A18 only; bovine
  differs from canine at A8/A10). An upstream release change fails loudly instead
  of silently shifting every downstream number.
* **KAT-2 / KAT-3** — negative control (own insulin → zero risk) and identity
  control (porcine insulin in a dog → zero risk). KAT-3 *skips for the cat*,
  with the reason recorded: no natural insulin is identical to feline insulin.
  That asymmetry is the shortest statement of the dog/cat difference.
* **KAT-4** — the exact core-set invariant described above.
* **KAT-5** — composition-matched scramble must look almost entirely non-self,
  so the filter cannot be fooled by amino-acid composition.
* **KAT-6** — a single-residue probe must flag exactly the cores spanning it.
* **KAT-7** — background `%rank` calibration is uniform on held-out peptides.
* **KAT-8** — binding-register assignment survives a window shift.
* **KAT-9** — panel members are separable from decoys (shuffles and
  opposite-chain-class molecules) on cross-species-invariant positions.
* **KAT-10** — the out-of-domain guard-rail fires rather than quietly passing an
  extrapolated score through.
* **KAT-11** — surrogate direction check: the HLA-DRB1\*01:01 P1 pocket must
  prefer hydrophobic anchors.
* **KAT-12** — scoring is bit-deterministic.

## Honest limitations

* **No allele-frequency database exists for DLA or FLA**, and IEDB's
  population-coverage tool is human-only. Canine DLA haplotype frequencies are
  strongly breed-dependent and feline data is sparser still. Nothing here can be
  converted into "X% of the population covered". A panel is a diversity sample.
* **IPD-MHC has no feline section.** The dog panel comes from a curated,
  ISAG-named allele set; the cat panel is assembled from GenBank/RefSeq records
  with no allele names and no denominator. `panel_summary.csv` reports which is
  which on every run.
* **DQ alpha/beta pairs are combinatorial.** Real DLA class II haplotypes are
  linked, so most generated pairs exist in no dog. Supply real haplotypes in the
  species config to replace this; the generated script says so at the top.
* **The surrogate backend is not a predictor.** It reads only the four anchor
  positions, so cores differing at a non-anchor residue score identically. Every
  artefact it produces is labelled.
* **Acylated analogues are modelled as peptide backbones only.** Detemir's
  myristoyl and degludec's hexadecanedioyl-γ-Glu chains change processing and
  can create hapten-like determinants; a sequence-level workflow does not see
  them, and the configs say so per product.
* **TCR-level prediction is deliberately absent.** NetTCR-class models are
  trained on human (and a little mouse) data and generalise poorly to unseen
  peptides even within humans; canine and feline TCR data is essentially
  nonexistent. Screening MHC-II presentation plus self-similarity, and stopping
  there, is the standard approach for human programmes too.
* **Everything here is research-use.** For a species-matched product the real
  immunogenicity risk sits in aggregation, process impurities and PTMs, which
  this workflow does not address at all.

## Layout

```
vetimmuno/
  data.py       cached fetchers: UniProt, IPD-MHC, NCBI E-utilities, IPD-IMGT/HLA
  insulin.py    mature chains, analogue edit grammar, difference map
  epitope.py    tiling, 9-mer registers, self-tolerance filter
  groove.py     pseudosequences, cross-species invariants, applicability domain
  panel.py      species panel construction + NetMHCIIpan custom-molecule inputs
  predict.py    NetMHCIIpan backend, illustrative surrogate, background %rank
  validate.py   the 16 known-answer tests and controls
  report.py     figures and the labelled markdown report
  workflow.py   per-species driver
config/         dog.yaml, cat.yaml — products, panel sources, interpretation
tests/          unit tests
demo.py         both species + head-to-head comparison
data/cache/     downloaded sources (git-ignored; re-fetched on demand)
```

Retarget to another species by adding a config: the code has no dog- or
cat-specific branches.
