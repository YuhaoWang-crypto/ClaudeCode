# Running a de novo glycoprotein design campaign

The paper's exact parameters, in the order you would use them. None of the
generative models (RFdiffusionAA, LigandMPNN, RF3, AF3) run inside this
package - it supplies the **scoring and triage** stages, which are the ones
that decide what gets expressed.

## Stage 0 - choose the glycan and its input context

RFdiffusionAA needs a glycan to diffuse around. The paper harvested
glycan-bearing asparagines (or whole PNGS motifs) from a deliberately diverse
panel - viral, human and designed bacterial:

| Source | Sites used |
|---|---|
| PDB 6VXX (SARS-CoV-2 spike) | 8: N61, N164, N235, N331, N603, N709, N717, N1134 |
| PDB 7PFP | 3: N232, N396, N513 |
| I53-50A-3gly (designed) | 3: N45, N75, N96 |
| PDB 6S7O, isolated OST-A binding-pocket peptide | 1: N3 |

The OST-A peptide fragment is there to raise **occupancy**, not shielding -
for some designs its motif was point-mutated at the -2 position to `FANATAA`,
a known occupancy-enhancing context.

Glycans were modelled with **GlycoShape**, with **5 conformations** per site
to diversify the diffusion input. Model a single Man9 (or Man6, see stage 3).

## Stage 1 - backbone generation (RFdiffusionAA)

```
~150,000 backbones
300-500 amino acids
50-75 diffusion steps
guiding potential: ligand_ncontacts, weight 0.5-1, guide scale 2, quadratic decay
```

The contact potential is what biases folds toward wrapping the glycan. The
paper also tried RASA-based conditioning during generation via a pre-release
build and found **no measurable effect**, so those designs are treated as
unconditioned. Do not spend effort there.

## Stage 2 - sequence design (LigandMPNN)

- **LigandMPNN with the glycan as ligand context** - this is the part that
  makes the sequence aware of the sugar.
- Recycle through **SolubleMPNN** or **EnhancedMPNN** for solubility and
  sequence diversity (SolubleMPNN specifically to mutate exposed hydrophobics).
- Sampling temperature **0.1**, **5-8 sequences per backbone**.
- Fix an **NXT or NXS** motif; prefer NXT for occupancy.
- **Omit cysteines.**

For redesigning an existing scaffold instead (the THRgly route): run
ProteinMPNN with **500 sequences per target**, save per-position probabilities,
and pick positions that already favour Asn with a sequon-compatible S/T two
residues downstream. Then graft Man9 with GlycoShape and run glycan-aware
LigandMPNN on the surrounding side chains.

## Stage 3 - structure prediction and the confidence gate

**Use a glycan-aware predictor.** This is the single most consequential
lesson in the paper (see `g4_triage`):

| Generation | Predictor | Gate |
|---|---|---|
| first | AF2, single-sequence, **no glycan** | pLDDT > 80, pAE < 10, Ca-RMSD < 3 A |
| second | **RF3** (Man9) or **AF3 server** (Man6) | pLDDT > 75, pAE < 15, Ca-RMSD < 8 A |
| THRgly redesign | AF2 (tighter) | pLDDT > 90, pAE < 5, Ca-RMSD < 2 A |

AF3's public server caps glycan size, which is why Man6 appears there while
RF3 runs get Man9. Keep the glycan constant *within* a comparison.

Of the four oligomannose designs that worked, **only O2 would have passed the
AF2 gate**. O1 and O3 scored pLDDT 42.7 and 67.7 with RMSD 21.8 and 16.8 A -
AF2, not seeing a glycan, interprets a glycan-housing pocket as a large
unfilled cavity and reports a broken model. Both were subsequently confirmed
oligomannose, and O4's crystal structure matched its design model over
residues 51-356.

## Stage 4 - DeltaRASA filter (this package)

Computed **post-prediction**, on the predicted structure with its modelled
glycan:

1. Score every design with `g1_delta_rasa.site_delta_rasa(..., context="protein")`.
2. Rank by **mean per-atom DeltaRASA** and separately by **number of contact
   atoms**.
3. Take roughly the **top 5% by either** and inspect them by eye.
4. Check the glycan is in a pocket or crevice, not merely near a bulge.

The payoff is large and measured: **20/43** DeltaRASA-filtered designs exceeded
25% oligomannose, against **1/97** without the filter (Fisher exact
p = 2.0e-11, odds ratio 83; the paper reports p < 1e-10).

Keep the glycan model identical across all designs or the ranking is
size-confounded - see `reference/methodology.md` section 4.

## Stage 5 - expression and readout

- Genes as eBlocks or synthesised, into a modified **pCMV/R**.
- **N-terminal secretion signal** (modified bovine prolactin or mu-phosphatase)
  and a **C-terminal 6xHis** (optionally Myc, TEV or WELQut cleavable).
- **Expi293F**, 3 days, harvest supernatant, PDADMAC, filter, Ni-NTA, then SEC.
- **Kifunensine at 20 uM** for the forced-oligomannose control.
- Cross-check in a second host - the paper repeated in **CHO** and found
  comparable glycan compositions, showing local protein structure rather than
  the expression host is the determinant.

Then go to `reference/wet-validation.md`.

## Realistic expectations

From the paper's own numbers, so you can size a campaign:

- 257 first-generation designs passed the AF2 filter; 126 glycoproteins were
  successfully detected; **97 had >5% occupancy**; of those, **94% carried
  complex-type** glycans and exactly **one** exceeded 25% oligomannose.
- With glycan-aware prediction plus the DeltaRASA filter: 82 sequences across
  24 backbones expressed, 51 detected, **20 of 43** over 25% oligomannose.
- Four designs were taken forward (O1-O4). Only **O3** showed a cooperative
  thermal unfolding transition (79 C); the others showed none below 95 C.
  Deeply enclosing a glycan and getting a cleanly two-state protein are not
  the same goal.
