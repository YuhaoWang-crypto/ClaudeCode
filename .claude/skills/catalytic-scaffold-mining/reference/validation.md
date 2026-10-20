# From candidate to claim: the controls that have to come first

A structural pipeline produces *nominations*. This file is about the distance
between a nomination and a defensible statement that a protein catalyses a
reaction — and about the specific ways a purified protein can appear to do
chemistry it does not do.

The precedent paper is unusually disciplined here. Its authors say plainly
that they detected the activity in April 2025 and then **spent the time to
prove it was intrinsic** before publishing. The controls below are theirs.

---

## Why the paranoia is justified

Claiming a protein polymerises nucleotides into RNA has four classic failure
modes:

| Failure mode | What it looks like | Control that kills it |
|---|---|---|
| Co-purifying enzyme | activity tracks the prep, not the protein | multi-step orthogonal purification; washed crystals |
| Carried-over nucleic acid | "template-independent" synthesis that is actually primed | protease digest + ion exchange; nuclease pre-treatment |
| Covalent cofactor / adduct | activity needs something you did not know was there | intact-protein mass spectrometry |
| Mis-assigned product chemistry | nuclease sensitivity read as proof of linkage | NMR, or equivalent direct bond-level evidence |

---

## The five orthogonal controls

### 1. Intact-protein mass spectrometry
Measure the mass of the purified protein and compare to the theoretical mass.
A match to within a few Daltons excludes a covalently bound nucleotide or
cofactor.

*Precedent:* OsNTP4_CD measured **37,681.0 Da** against a theoretical
**37,681.53 Da**. No adduct.

### 2. Exclude free nucleotide and oligonucleotide carryover
Digest the protein completely (Proteinase K), then run the digest over ion
exchange and look for released nucleotides, against standards processed
identically.

*Precedent:* 200 uM protein digested, analysed on sequential Mono S then
Capto Q. Each of the four nucleotide standards (200 uM) was detected cleanly;
**no nucleotide and no oligonucleotide was detected in the protein sample.**

The matched-standard design is the point: it proves the assay would have seen
the contaminant had it been there.

### 3. Orthogonal multi-step purification
Affinity -> ion exchange -> size exclusion, ideally with a tag cleavage and a
second pass over the affinity resin to remove uncleaved protein and protease.

*Precedent:* Ni-NTA -> Heparin -> Superdex 200 for the full-length protein;
for the catalytic domain, Ni-NTA -> Source 15Q -> 3C protease cleavage ->
Ni-NTA subtraction -> Source 15Q -> Superdex 200.

### 4. The washed-crystal assay — the decisive one
Grow crystals, wash them repeatedly in well buffer, dissolve, and assay. Run a
**matched non-crystalline control collected from the same drop** and processed
identically.

*Precedent:* crystals washed three times; dissolved crystals synthesised
high-molecular-weight product from UTP, the matched non-crystalline control
from the same drop did not. Detected by Cy5 fluorescence and independently by
SYBR Gold.

This is the strongest available argument short of reconstitution from
synthetic gene product: a trace contaminant essentially cannot co-crystallise
stoichiometrically and survive three washes.

### 5. Active-site mutants, with a surface-residue control
Each catalytic residue to Ala must abolish activity; a surface residue
mutation must not.

*Precedent:* D100A, D102A, D157A each abolished activity; **K194A (surface)
had little effect.** The negative control is what makes this informative — it
shows the protein tolerates mutation and did not simply misfold.

---

## Product chemistry needs its own proof

Nuclease sensitivity does not establish bond chemistry. RNase I cleaves
3'-5' phosphodiester bonds, so sensitivity is *consistent with* an RNA product
but does not demonstrate the linkage.

*Precedent:* `1H,1H` COSY and `1H,31P` HMBC on the homopolymeric product, with
synthetic poly-(U)/poly-(A)/poly-(C) references. The HMBC correlation of the
bridging phosphate to **H3' and H5'** is diagnostic of a canonical 3'-5'
linkage; a 2'-5' linkage would correlate to **H2'** instead.

Decide up front which observable distinguishes your two candidate chemistries,
and measure that.

---

## Honest labelling of a computational shortlist

Until the ladder above is climbed, a candidate is `⚠️ hypothesis`. Specific
things to **not** write:

- "X is an RNA polymerase" — write "X conserves the catalytic configuration of
  known NTases and is a candidate for nucleotidyl-transfer activity".
- "X synthesises RNA de novo" — that is a claim about substrates and products
  that no structure supports.
- "we identified N new enzymes" — you identified N candidates.

Things to **do** write:

- what was computed, with the numbers;
- what the thresholds were and how they were chosen;
- which taxa were searched, and (crucially) which were **not**;
- that hit caps make absence uninformative;
- the single experiment that would most cheaply falsify each nomination.

---

## Interpreting a negative result

If a candidate fails in vitro, that is a result about **the assay conditions**
as much as the protein. Before concluding "not an enzyme", check:

- metal dependence — the precedent reaction needs divalent metal, and the
  preference is substrate-dependent (Mg2+ or Mn2+ sufficed for poly-(U), while
  poly-(A) and poly-(G) largely required **Mn2+**);
- substrate identity — the same enzyme gave strong product with UTP and barely
  any with GTP;
- product detectability — G-rich products can form G-quadruplexes, which stain
  poorly with SYBR Gold and resist nuclease, so they can look like "no
  product" when they are simply invisible to the readout;
- whether the construct is the catalytic domain or the full-length protein.

Report a negative as a measured negative with its conditions, not as a
property of the protein.
