# Methodology — why the constraints are what they are

Source: Kipouros & Chang, Nature 656, 763–770 (2026), doi:10.1038/s41586-026-10716-z.
Reference implementation: `github.com/yannikipouros/hal-discovery` (MIT).

## 1. The problem this solves

Sequence identity is a poor proxy for enzyme function. Sequence evolves under
selection, drift, redundancy, expression cost and horizontal transfer; catalytic
architecture does not have that freedom. Inside a promiscuous superfamily, two
enzymes at >50% identity can catalyse different reactions, while two enzymes
performing the same reaction can sit in unrelated clusters. Sequence methods
also scale as O(N²) in pairwise comparisons — hopeless when the family has
N > 10⁵ and the difference you are hunting is **one residue**.

The three residues of the halogenase motif are far apart in 1D sequence and
differently placed in different subfamilies, so no alignment column carries the
signal. In 3D they are one object.

## 2. Anchor + discriminator

**Anchor** = the geometrically stereotyped, superfamily-wide feature that locates
the site. Here: two His on adjacent β-strands of the cupin double-stranded
β-helix. Three constraints suffice, and each does distinct work:

| Constraint | Default | What it enforces |
|---|---|---|
| Nτ(NE2)–Nτ(NE2) distance | < 4.0 Å | the two side chains are *preorganized* to chelate one metal |
| backbone O(i)···N(j) | < 4.0 Å | strand-pairing H-bond, one direction |
| backbone N(i)···O(j) | < 4.0 Å | strand-pairing H-bond, the other direction |

Side-chain distance alone admits accidental proximity anywhere in the fold;
backbone pairing alone admits every β-strand pair. Requiring **both** backbone
distances makes the pair an antiparallel-strand His pair — the cupin metal site.

**Discriminator** = the feature that mechanism makes mandatory for the target
function. Here it is an *absence*: no Asp/Glu near the metal, because a halide
must occupy that coordination position. Absence-based discriminators are strong
precisely because they are hard to satisfy by accident — a chelating carboxylate
near an exposed metal site is otherwise the overwhelmingly favoured arrangement
(456,585 facial-triad sites vs 946 halogenase-like in the published run).

Published classification rules (all reproduced by `mcmine.py`, defaults in
`assets/motifs/fe_akg_radical_halogenase.json`):

```
closest ASP  > 5.5 Å  and  closest GLU  > 5.5 Å      # no carboxylate third ligand
closest ASN/GLN/CYS/MET/HIS/TRP/TYR > 4.0 Å          # no other ligand took its place
(closest ALA < 7.0 Å or closest GLY < 7.0 Å)         # the small residue that makes the room
residue at (HisA + 2) in {ALA, GLY}                  # sequence-context sanity check
protein length > 200 aa                              # not a fragment
```

Contrast class (the sibling you must separate from):
`min(closest ASP, closest GLU) ≤ 5.0 Å` → 2His-1Asp/Glu facial triad.

### The two atom sets — a real trap

Deciding what is **absent** and deciding what is **present** need different atom
sets:

- **`coordinating`** — only atoms that can ligate a metal (ASP: OD1/OD2; GLU:
  OE1/OE2; HIS: NE2/ND1; CYS: SG; …). Use for typing and first-shell membership.
- **`proximity`** — coordinating atoms **plus side-chain carbons** (ASP adds
  CB/CG; GLU adds CB/CG/CD). Use for exclusion rules: a carboxylate whose CB/CG
  sits near the probe can rotate its oxygens in, so excluding on OD/OE alone
  under-rejects and inflates the hit list.

The published halogenase filter uses `proximity`; the coordination-typing survey
uses `coordinating`. `mcmine.py` keeps both and the spec picks one.

## 3. Probe point

The metal is absent from AFDB models, so the site is represented by the
**midpoint of the two anchor side-chain atoms** — a virtual metal position. All
shell distances are measured from it. This is the step that makes apo predicted
structures usable: what is mined is the protein's *preorganization* for a
cofactor it does not currently hold.

## 4. Cost and scale

Per structure the work is: find anchor pairs (O(H²) over the handful of anchor
residues, H ≈ tens) then one distance pass over candidate atoms. Linear in
database size, no pairwise comparison, no alignment, no GPU. The published screen
of 530,814 AF2 models ran on ordinary hardware; the 200-structure demo takes
~10 s. Budget roughly 30–60 ms per model single-threaded, plus download time
(see `data-access.md`) which dominates.

## 5. Published results worth knowing

- 220M InterPro sequences → 1.8M cupin-domain → 530,814 AF2 models mined.
- 458,000 sites with 2His+Asp/Glu; **946** with 2His+Ala/Gly (the halogenase call).
- SSN at 30% identity: all 6 known halogenase families recovered (positive
  control) **plus 70 clusters with no prior halogenase annotation**.
- DAH (eukaryotic halogenase): BLAST+MSA found 1–3 homologues; this method built
  a 40-member cluster spanning fungi, plants and bacteria.
- Experimental validation: **AspX** (*Vibrio campbellii*) chlorinates
  L-aspartate — first negatively-charged amino-acid substrate — kcat 33.3 ± 0.5 min⁻¹,
  Km 0.64 ± 0.02 mM; **BtnX** (*Dinoroseobacter shibae* killer plasmid)
  chlorinates biotin and is unprecedentedly promiscuous (bile acids, dyes,
  peptides, a drug candidate). Both also transfer Br⁻ and N₃⁻.
- Crystal structure of BtnX at 1.20 Å (PDB 9PV1) confirms the predicted
  2His-1Gly site with Fe, Cl, αKG and biotin bound.

## 6. The engineering corollary

If the discriminator *causes* the function, then installing it should *install*
the function. In BtnX, **G117D and G117E abolish halogenation and increase
hydroxylation 40 ± 9-fold and 18 ± 6-fold** — the mined motif read backwards as a
design rule. Whenever a motif is validated, the reciprocal mutation is the
cheapest possible test of the mechanistic claim, and a reusable switch.

## 7. What generalizes

The paper's own extension: drop the halogenase rules and enumerate all **2His-Xₙ**
environments (n = 0–4; X = Asp, Glu, His, Asn, Gln, Cys, Met, Tyr, none) across
538,160 cupin sites — a census that surfaces both known and unknown biological
metal-site types. That is `mcmine.py type`.

Beyond metals, the requirement is only that the function depends on a **small set
of atoms in a defined geometry**: catalytic triads and dyads, oxyanion holes,
proton-relay pairs, recognition motifs, disulfide/metal-switch sites. The anchor
must be superfamily-conserved and geometrically distinctive; the discriminator
must be mechanistically necessary rather than statistically associated.

## 8. Honest limits

- **Prediction, not annotation.** 946 candidates, 2 assayed. Hits are hypotheses.
- **Model quality.** Low-pLDDT regions give unreliable side-chain geometry;
  AFDB models are apo and single-conformer, so alternate rotamers that would
  complete a coordination sphere are invisible.
- **Substrate scope is not predicted.** The motif says *what chemistry*, never
  *on what molecule*; substrate came from genome context plus assays.
- **Non-enzyme contamination.** The cupin fold hosts transcription regulators and
  storage proteins; a UniProt keyword filter after mining is still needed
  (`transcription`, `regulator`, `AraC`, `globin`, `AlkB`, `glutelin`, `TehB`,
  `tet`, `fragment`, `chemotaxis`, `helix-turn-helix`, `tellurite`, `adenosyl`,
  `SAM`, `ferredoxin`).
- **Threshold sensitivity.** Cutoffs are physically motivated but not sacred;
  always sweep and report the plateau.
- **Multi-chain sites.** Sites bridging chains are only visible if the model
  contains both chains — AFDB monomer models cannot show them.
