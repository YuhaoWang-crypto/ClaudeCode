# Epitope occupancy: definitions, verification, conventions

Everything here is computed deterministically from a predicted complex. That is
the point — it is the part of the pipeline that does not depend on a model's
opinion of its own output.

## Deriving the epitope

When an experimental complex of the target with its natural partner exists, the
epitope is the set of target residues in contact with that partner. Compute it
rather than transcribing it from a paper's figure:

- heavy-atom distance cutoff 4.5 Å, any atom to any atom
- take the residues on the *target* side
- record hotspots separately if the literature identifies them; they earn their
  own term in the blocker criterion

**Mapping across species or numbering schemes.** PDB entries frequently use
mature-protein numbering while UniProt includes the signal peptide. Derive the
offset by aligning the PDB chain's observed sequence against the reference and
taking the shift that maximises identity — then *assert* it reproduces the
residue identities at every epitope position. Do not hardcode an offset from
memory.

When mapping human→mouse (or any ortholog pair), check whether the alignment is
ungapped over the region before treating positions as interchangeable. If it is,
numbering carries across directly; if not, align properly.

Worth computing and reporting: how many epitope residues actually differ between
species. In the TGFβR2 case the ectodomain was 81.8% identical overall but the
ligand interface was 11/12 identical — which reverses the naive reading. The
divergence that justifies a species-specific reagent is *peripheral*, where
antibody epitopes often sit, not at the functional interface.

## Choosing the target model

Designing against a predicted model is fine if you show the model is right where
it matters:

- superpose the model's structured core onto the experimental structure (Kabsch
  on CA atoms) and report RMSD over the core **and** over the epitope alone
- report mean and minimum pLDDT across the epitope residues
- crop residues below ~70 pLDDT — usually disordered termini. Leaving them in
  invites designs that bind coordinates that do not exist

## The neighbourhood check

Scan the target sequence for N-glycosylation sequons (`N-X-[ST]`, X ≠ P) and
compute the 3D distance from each sequon's side chain to the nearest epitope
residue. A sequon within ~10 Å matters: a complex N-glycan spans well beyond
that, so it can shade part of an epitope that looks fully exposed in an
aglycosylated model.

Two consequences, both worth stating explicitly in a methods section:
- discourage binder contact near the sequon in the design spec
- the screening antigen must be expressed in a system that glycosylates
  (mammalian), or you will select binders that cannot reach the native target

## Metrics

Computed per design from the predicted complex, heavy-atom contacts at 4.5 Å:

| metric | definition | what it catches |
|---|---|---|
| `n_contact` | target residues contacted | interface size |
| `recall` | \|contacted ∩ epitope\| / \|epitope\| | how much of the blocking face is covered |
| `purity` | \|contacted ∩ epitope\| / n_contact | whether the binder is *focused* on the epitope |
| `atom_purity` | atom-weighted purity | a paratope offset onto flankers that residue-level purity misses |
| `hotspots` | which known hotspots are contacted | mechanistic blockade, not just proximity |
| `n_remote` | off-epitope contacts neither sequence-adjacent (±2) nor within 8 Å of any epitope residue | a genuinely separate binding patch, i.e. the wrong site |

The rim-vs-remote distinction matters. Off-epitope contacts immediately flanking
the epitope are the normal footprint of a real interface. A second patch
elsewhere on the surface is a different binding mode wearing the same score.

**Report the attainable ceiling, not the nominal one.** Some epitope residues are
sterically hard to reach — in our case one residue was engaged by 8/96 designs,
so the practical recall ceiling was 11/12, not 12/12. Quoting recall against 12
makes every design look worse than it is. Similarly, if one hotspot is contacted
by 100% of designs, hotspot coverage is effectively a test of the *other* one —
say so rather than reporting a two-hotspot score that has one degree of freedom.

## The blocker criterion

Define it from the observed distribution and justify each term:

```
blocker = recall ≥ <quartile of observed>            # covers the face
      AND all known hotspots contacted                # mechanistic
      AND atom_purity ≥ <just below observed Q1>      # focused, not sticky
      AND n_remote = 0                                # no second site
```

Report stability: how the count moves as you shift each threshold. A criterion
that swings wildly with a 0.05 change is not a criterion.

Do not import a threshold from a different target. Absolute co-folding cutoffs
are measurably non-portable — in published labelled data, a `score ≥ 0.90` gate
*lowered* hit rate versus a within-target top-20% rule and eliminated most
targets entirely.

## Per-file verification

Before scoring any structure, assert:

1. the expected chains are present (target and binder)
2. the target chain spans the expected residue range
3. **the target chain's sequence equals the reference exactly**

Do not score a file that fails; count failures and report them. This costs
almost nothing and is the only thing standing between you and an entire
analysis that is confidently off by a constant.

## CIF parsing

Boltz-2 predicted complexes are mmCIF with whitespace-splittable ATOM records.
Field indices observed (verify on the first file of any new run rather than
trusting this table):

| index | field |
|---|---|
| 5 | residue name (3-letter) |
| 6 | residue number |
| 9 | chain id |
| 10/11/12 | x / y / z |

Target chains keep native (UniProt) numbering when the design spec cropped by
residue index — verify by sequence match, which also confirms the crop.

`scripts/occupancy.py` in this skill implements all of the above including the
assertions. Prefer it over rewriting the contact loop each time.

## Working with large result sets

Design-API result payloads are dominated by presigned URL text and are enormous
(~0.5 MB per page of 100). Never load them into an agent's context:

- paginate, and parse the spilled response files with `jq`/Python
- extract only the metric fields and the artifact URLs you need
- presigned URLs typically expire ~30 min after issue — fetch in batches of ~20
  and download each batch immediately rather than collecting all URLs first
- cache downloaded structures on disk; later analyses will want them again
