# Proving the glycan state at the bench

Four assays, in increasing cost and information. Run them in this order; each
one can stop the campaign.

## 1. Enzyme pair - PNGase F vs Endo H (SDS-PAGE, hours)

The cheapest decisive experiment. Both enzymes are read by gel shift.

| Enzyme | Cleaves | Leaves |
|---|---|---|
| **PNGase F** | all N-glycans, at the Asn-GlcNAc amide | nothing (Asn -> Asp) |
| **Endo H** | **oligomannose and hybrid only**, between the two core GlcNAc | one GlcNAc, **+204 Da** |

Decision table:

| PNGase F shift | Endo H shift | Conclusion |
|---|---|---|
| yes | yes | glycan present, **oligomannose/hybrid** |
| yes | no | glycan present, **complex-type** |
| no | no | site unoccupied, or no glycan |

**Partial Endo H digestion is a result, not a failure.** THRgly3, THRgly4 and
the HIV-1 Env precedent all digest incompletely because the enzyme cannot
fully reach a buried glycan. Incomplete digestion of a design you *predicted*
to be enclosed is positive evidence for enclosure. Always run enzyme-only
lanes.

## 2. Kifunensine control (the orthogonal axis)

Kifunensine inhibits ER mannosidase I, forcing oligomannose **regardless of
protein structure**, at 20 uM in the culture.

It separates two things a single gel cannot:

- **Can this site carry oligomannose at all?** Kifunensine-treated material
  says yes for every site, including the fully exposed THRgly1.
- **Does the structure do it on its own?** Only the untreated sample answers
  that. The paper's whole claim rests on untreated THRgly2-4 matching their
  kifunensine-treated counterparts.

Kifunensine is also the right positive control for lectin binding assays - it
gives you the "maximally oligomannose" version of your own protein rather than
an unrelated standard.

## 3. Intact mass spectrometry (when peptides fail)

Use when the protein resists denaturation badly enough that peptide-level
glycan analysis fails - which is exactly what happened to the THRgly repeat
proteins.

- Reduce free cysteines with TCEP, desalt, RP-LC/MS, deconvolute.
- Read glycan composition from the **mass shift against the predicted
  unglycosylated mass**.
- After Endo H, an oligomannose site collapses to **+204 Da**. Seeing that
  shift appear confirms the assignment from step 1 at mass resolution.

Intact MS gives composition without site assignment, so it only works cleanly
for single-site constructs. It cannot distinguish isomers.

## 4. LC-MS/MS glycopeptide analysis (site-specific, the real answer)

The only method that gives **per-site composition and occupancy**.

- Denature in 6 M urea, buffer exchange, digest overnight with
  **trypsin-LysC**, **chymotrypsin** or **alpha-lytic** protease at 1:30 - pick
  whichever gives a tractable glycopeptide around your sequon.
- Desalt, LC-ESI-MS on an Orbitrap with **stepped HCD** (NCE 15/30/45).
  Stepped energy is what yields oxonium ions and peptide b/y ions in the same
  spectrum.
- Search against an N-glycan library at tight tolerances (4 ppm precursor,
  10 ppm fragment, 1% FDR), then **validate manually**: require the expected
  b/y ions *together with* oxonium ions.
- Quantify from summed **XIC areas across all charge states** of a given
  glycopeptide. Report relative abundance within the occupied fraction, and
  occupancy separately.

Short gradients (30 min) first; re-run anything with thin glycopeptide
coverage on a long one (180 min).

Feed the resulting table into `g3_underprocessing` to replace
`ILLUSTRATIVE_PROFILES`.

## 5. Does the glycan still *work*?

Steric burial that also blocks recognition has achieved nothing. The paper
checks both directions, and they disagree in an informative way:

- **FimH** (uropathogenic *E. coli* adhesin) prefers **small** oligomannose,
  Man5. Assays: RNase B competition ELISA; T24 bladder-cell adhesion;
  haemagglutination inhibition. The Man5-rich designs O3 and O4 blocked
  adhesion by 56% and 75%; O3 inhibited haemagglutination with EC50 ~8 uM,
  100-1000x better than free mannose.
- **MBL** (innate immune lectin) prefers **large** mannose-rich glycans. By
  MBL-C ELISA the Man8-rich O2 bound best, reaching 79-99% of the
  kifunensine-treated control.

So pick the target lectin *first* and then decide which rung of the ladder you
are aiming at. Maximum burial is not the goal; the right glycoform is. Note
also that O3 and O4 carried similar Man5 but differed substantially in FimH
blocking - the surrounding scaffold shapes engagement too, so a glycan-only
design rationale will not fully predict activity.

## Pitfalls

- An **unmodelled** glycan in a structure is not an **unoccupied** site. Only
  LC-MS/MS occupancy settles that.
- Occupancy falls as enclosure rises, because OST also has to reach the
  sequon. THRgly2 and THRgly4 had low occupancy precisely because they worked
  sterically. Budget for it: bias +2 to Thr, and consider the OST-A
  binding-pocket motif.
- Glycan heterogeneity is the norm. "Programmed" here means a shifted
  distribution, not a single species - even the fully enclosed designs gave
  mixtures.
