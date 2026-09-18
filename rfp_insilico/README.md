# Can the RFP's validation be done in silico?

**Short answer: a useful slice of it can, as a pre-screen. The core comparison
the study exists to make cannot, and the reason is specific rather than
hand-wavy.**

The RFP asks a CRO for peptide synthesis plus three wet assays: a mouse
MHC-II / DO11.10 IL-2 readout, cell-free and FACS-based HLA tetramer binding on
two human alleles, and an MLR with proliferation, cytokine and activation
readouts. Below, each line item is marked for what computation can and cannot
do.

`rfp_demo.py` runs the parts that can, on live predictions, against the RFP's
own alleles.

---

## Feasibility, line by line

| RFP item | In silico? | What is actually possible |
|---|---|---|
| Peptide synthesis, >95%, TFA→acetate, Ac-/-NH₂ | **No** | Chemistry. Computation can only flag synthesis liabilities: aggregation-prone stretches, Met/Cys oxidation, poor solubility. |
| **Step 1** binding to mouse I-A^d | **Yes** | NetMHCIIpan scores H2-IAd directly and recovers the OVA 323-339 register (below). |
| **Step 1** IL-2 production by DO11.10 | **No** | No method predicts cytokine output from a hybridoma. Binding is necessary, not sufficient. |
| **Step 2** peptide binding to HLA-A\*32:01 / B\*57:01 | **Yes** | Both alleles are supported; the demo ranks 9-mers on each. |
| **Step 2** cell-free affinity *and stability* screen | **Partly** | Predictors give an eluted-ligand likelihood and a %rank, not a measured off-rate or complex half-life. Stability is the part that does not transfer. |
| Tetramer/pentamer assembly | **No** | Reagent production. |
| FACS baseline **precursor frequency** of antigen-specific TCRs | **No** | Not predictable from sequence. This is the RFP's most information-rich readout and has no computational surrogate. |
| MLR: proliferation, cytokines, activation markers | **No** | Cellular phenotype. |
| HLA typing of donors | **No** | Genotyping. Computation can help *plan* recruitment: both alleles have known population frequencies, so the number of donors needed to find carriers is estimable. |
| **The all-L vs 4×D-amino-acid comparison** | **No — and this is the important one** | See below. |

## The blocker: the predictors cannot see D-amino acids

The study's central contrast is an all-L peptide against the same peptide
carrying four D-residues. Every sequence-based MHC predictor in routine use —
NetMHCpan, NetMHCIIpan, MHCflurry — consumes a 20-letter alphabet with no
stereochemistry channel. A D-residue and its L-enantiomer are the same
character, so the two molecules are **the same input string**.

The demo prints this rather than asserting it:

```
L-form, as the model sees it : ISPRTLNAW
D-form, as the model sees it : ISPRTLNAW
identical input strings      : True
```

The prediction for the D-variant is therefore not merely uncertain. It is
*identical to the L-variant by construction*, so it carries zero information
about the one variable the study exists to measure, while looking exactly like a
normal result. That is the dangerous failure mode: a confident number that
cannot be wrong because it is not about anything.

MHC class I and II both grip the peptide backbone through a stereochemically
specific groove, and D-substitution is chosen precisely because it changes
backbone presentation and protease resistance. The quantity the tool silently
assumes is unchanged is the quantity under test.

Also entirely outside the encoding: N-terminal acetylation, C-terminal
amidation, counter-ion form, and the uniblock/diblock assembly state.

**If someone offers to replace the D-peptide arm of this study with a binding
prediction, that is the claim to refuse.**

## What the demo establishes

```
python rfp_demo.py                          # ovalbumin as example source antigen
python rfp_demo.py --peptides PEP1 PEP2 ...  # your own 9-mers
```

**1. All three RFP alleles are supported** by the backend, checked against the
live allele lists: `HLA-A*32:01`, `HLA-B*57:01`, `H2-IAd`.

**2. Positive controls reproduce the RFP's own immunology.** This is what makes
the rest readable — a predictor that fails here is not worth running.

| control | result |
|---|---|
| OVA 323-339 on H2-IAd (the DO11.10 restriction element) | %rank **0.17**, strong binder, core `VHAAHAEIN` |
| Documented B\*57:01-restricted 9-mer `ISPRTLNAW` on B\*57:01 | %rank **0.09**, strong |
| the same peptide on A\*32:01 | %rank **0.71**, weak — correct allele discrimination |

**3. Shuffle control: the motif is doing the work, not the composition.**
Against 20 composition-matched anagrams of the positive control, **0 of 20**
scored better than the real peptide (real %rank 0.09, anagram median 31.0). A
predictor keying on amino-acid composition would not show that gap.

**4. A ranked 9-mer shortlist** per allele, which is the deliverable RFP Step 2
needs in order to choose which peptides to spend custom tetramer synthesis on.
On ovalbumin as an example source: 11 strong binders for A\*32:01, 5 for
B\*57:01.

## A numbering trap worth catching before ordering peptides

The immunology literature numbers ovalbumin on the **mature** protein, without
the initiator methionine. So the classic "OVA 323-339" is UniProt P01012
residues **324-340**:

```
UniProt 323-339 : KISQAVHAAHAEINEAG   <- wrong, shifted by one
UniProt 324-340 : ISQAVHAAHAEINEAGR   <- the canonical peptide
```

Ordering by position rather than by sequence gives a one-residue frameshift on
a control peptide. Cheap to check, expensive to discover in a failed assay.

## What this is worth

The honest value is **triage, not replacement**. Running this before the CRO
work costs minutes and can:

- confirm the five uniblock 9-mers actually bind the two chosen alleles, so
  custom tetramer synthesis is not spent on a peptide that will not load;
- predict the class II binding register, which tells you whether a truncation
  keeps the core intact;
- catch sequence and numbering errors before synthesis.

It cannot reduce the assay list, and it cannot address the L-versus-D question
at all. Every cell-based readout in the RFP still has to be run.

## Which skills apply

- **`immunogenicity-multimodel`** — used here. Wraps IEDB NetMHCpan/NetMHCIIpan
  cloud REST plus a BioLib MHC-II population-immunogenicity model. The right
  tool for the binding-prediction parts of Steps 1 and 2.
- **`neoantigen-selection`** — built for choosing tumour neoantigens for a
  personalised vaccine. Its HLA presentation machinery overlaps, but its scoring
  (agretopicity, foreignness, clonality, expression) is designed for mutant
  self-peptides in a vaccine payload and does not map onto a designed D-peptide
  study.
- **`contact-probability-specificity`** or a Boltz-2 co-fold would give a
  *structural* view of the pMHC complex. Worth knowing that this shares the same
  blind spot: standard structure predictors also assume L-amino acids, so they
  do not rescue the D-peptide question either.

## What the ProteinTalks virtual cell model can do here: nothing

Asked directly, because it is a reasonable thing to wonder given the rest of
this repository.

The ProteinTalks reproduction in `../proteintalks/` cannot perform any
validation in this RFP, and the mismatch is categorical rather than a matter of
accuracy:

- **Wrong biology.** It is trained on breast cancer cell lines perturbed with
  small-molecule drugs. There are no immune cells, no T cells, no MHC, no TCR
  and no peptides anywhere in its training corpus.
- **Wrong input space.** Its perturbation channel is target occupancy over a
  fixed list of 5,585 proteins. A synthetic designed peptide is not an element
  of that space, so there is no way to even pose the question to the model.
- **Wrong output.** It emits a proteome trajectory at 6/24/48 h plus a binary
  drug-efficacy label. Nothing in the RFP asks that question.
- **And on its own task it is modest.** Our reproduction found it beats a
  no-proteomics drug-mean control by 0.016 AUROC, and its own supplement reports
  no significant advantage over random forest on unseen drugs.

Using it here would be a category error. The MHC binding predictors above are
the correct tool, within the limits stated.
