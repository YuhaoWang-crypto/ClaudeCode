---
name: metal-coordination-mining
description: >-
  Mechanism-guided 3D active-site motif mining of predicted-structure databases
  (AlphaFold DB / ESM Atlas) for targeted enzyme discovery and functional
  classification inside large, sequence-diverse superfamilies. Use when the task
  is to find rare members of an enzyme family that sequence search misses,
  discriminate a catalytic subclass from its look-alike siblings (halogenase vs
  hydroxylase, one metal-site type vs another), survey the metal-site diversity
  of a fold, turn a mechanistic hypothesis into a database-scale screen, or
  design the residue swap that switches one activity into another. Encodes the
  metal-coordination mining method of Kipouros & Chang, Nature 656, 763 (2026),
  generalized to arbitrary motifs via a JSON spec, and enforces
  benchmark-before-scale plus ✅-computed vs ⚠️-predicted labeling.
---

# Metal-coordination mining — mechanism as a database query

The method: distil an enzyme's *mechanistic requirement* into atomic-level 3D
constraints, then use those constraints as the search query over a
predicted-structure database. Function comes from geometry, not from sequence
identity — so this finds members that BLAST/SSN clustering cannot reach, and it
scales as **O(N)** instead of the O(N²) of pairwise alignment.

Source: Kipouros & Chang, *Targeted enzyme discovery using metal-coordination
mining*, Nature **656**, 763–770 (2026), doi:10.1038/s41586-026-10716-z; their
code: `github.com/yannikipouros/hal-discovery` (MIT).

## The idea in one worked case

Fe(II)/αKG-dependent enzymes all chelate Fe with a **2-His pair on adjacent
β-strands** of the cupin fold. Most add a third ligand — Asp or Glu — forming the
canonical **2His-1Asp/Glu facial triad**, and hydroxylate/desaturate/cyclize.
Radical **halogenases** must instead park a halide on the metal, so they *cannot*
have that third carboxylate: they show **2His-1Gly/Ala**, an open coordination
site. One residue's *absence*, mechanistically required, is the entire classifier.

Result on the AF2 database: 530,814 cupin models → 458,000 facial-triad sites vs
**946** halogenase-like sites → all 6 known halogenase families recovered plus
**70 previously unrecognized families**; two predictions (AspX, BtnX) were
experimentally confirmed as halogenases. For comparison, BLAST found 1–3
homologues of the eukaryotic halogenase DAH where this method found a 40-member
cluster.

**The transferable pattern**: a *conserved anchor* that every superfamily member
shares (finds the site) + a *discriminating feature that mechanism makes
mandatory* (assigns the function). Both halves are required.

## Pipeline

| Step | What | Tool |
|---|---|---|
| 1 | Mechanism → motif spec (anchor + discriminator) | `reference/motif-spec.md` |
| 2 | **Benchmark on knowns before mining at scale** | `mcmine.py benchmark` |
| 3 | Sequence space → accessions → AFDB structures | `fetch_structures.py` |
| 4 | Mine + classify | `mcmine.py mine` |
| 5 | Survey what else is out there | `mcmine.py type` |
| 6 | Triage: SSN, genome neighbourhood, taxonomy, literature | `reference/downstream-triage.md` |

## Quick start

```bash
pip install biopython pandas numpy          # that is the whole dependency list

# 0. prove the toolchain reproduces the paper (clones a 200-structure demo set)
bash assets/validate_against_paper.sh

# 2. benchmark first — recall must be 1.0 on known positives
python assets/mcmine.py benchmark \
  --struct-dir demo_structures/ \
  --motif assets/motifs/fe_akg_radical_halogenase.json \
  --positives known_members.txt --outdir bench/ \
  --sweep "require_absent.0.min_distance:4.0:7.0:0.25"

# 3. build the structure set
python assets/fetch_structures.py accessions --families PF13640,IPR005123 --outdir acc/
python assets/fetch_structures.py download --accessions acc/accessions_af2_len_gt_150.csv --outdir structs/

# 4. mine
python assets/mcmine.py mine --struct-dir structs/ \
  --motif assets/motifs/fe_akg_radical_halogenase.json --outdir results/

# 5. survey coordination diversity (hypothesis generation, no rules applied)
python assets/mcmine.py type --struct-dir structs/ \
  --motif assets/motifs/two_his_xn_typing.json --outdir survey/
```

**Validation status of this skill (✅ reproduced here):** on the paper's 200-structure
demo set, `mcmine.py mine` returns 203 2-His sites, 192 facial-triad sites and
**10 halogenase hits whose accessions match the published set exactly**
(10/10 known positives, 0/190 negatives → recall 1.0, specificity 1.0), matching
the reference implementation's output. Freshly downloaded AFDB models of SyrB2,
CytC3 and OocP are also recovered, so the motif holds against current AFDB
releases. `assets/validate_against_paper.sh` re-checks all of this on demand.

## Non-negotiable: benchmark before scale

A motif that has not been tested against known members is a guess with a CSV
attached. Before any large run assemble **positives** (experimentally
characterized members of the target class) and **negatives** (characterized
members of the sibling class using the same fold), then require:

- **recall = 1.0** on positives — a miss means the geometry or the mechanism
  reasoning is wrong, and loosening thresholds afterwards is fitting to noise;
- **specificity reported**, not assumed;
- a **threshold sweep** (`--sweep`) showing whether you sit on a plateau or a
  knife edge. On the demo set the acid-exclusion cutoff is flat from 4.0–7.0 Å
  ✅; a real database set will not be that forgiving — report where the cliff is.

If a positive has no anchor site at all, `benchmark` says so separately: that is a
detection failure (geometry/model quality), not a classification failure.

## Designing a motif for a new mechanism

Ask, in order:

1. **What must be true of the atoms for the chemistry to happen?** An open
   coordination site, a specific ligand set, a proton-shuttle pair, a distance
   between two catalytic groups. That sentence is the motif.
2. **What is the anchor?** A geometrically stereotyped pair the whole superfamily
   shares. Pin it with *side-chain* distance (preorganization) **and** backbone
   cross-distances (secondary-structure adjacency) — either alone over-collects.
3. **What is the discriminator?** Prefer a *required absence* (as here) or a
   *required presence* that mechanism forces. Do not use residues that merely
   correlate with the class in the training examples.
4. **What is the negative class?** If you cannot name the look-alikes, you cannot
   claim specificity.

Then write the spec (`reference/motif-spec.md`) — no code changes needed. Apo
predicted models are fine: the point is that the *protein* is preorganized for
the cofactor, so metals and ligands being absent from AFDB does not matter.

## Honesty labeling

Every claim carries one of:

- **✅ computed** — a number this pipeline produced (site counts, hit lists,
  recall/specificity on a declared benchmark, coordination-type census).
- **⚠️ predicted** — biological function of a hit. A hit is a *candidate*, never
  an annotation. Even in the source paper, only 2 of 946 candidates were assayed.

Report false-positive modes explicitly: DNA/RNA-binding and structural His pairs
with accidental geometry, disordered low-pLDDT regions, and non-enzyme cupins
(transcription regulators, seed-storage proteins) — the source pipeline still
needed a UniProt keyword filter (`transcription`, `regulator`, `AraC`, `globin`,
`fragment`, …) after mining. Never state a hit count without stating the
benchmark it was validated on.

## Reference files

- `reference/methodology.md` — the method in full: why each constraint, published
  numbers, complexity, what generalizes beyond metal sites, known limits.
- `reference/motif-spec.md` — spec schema, every field, sweep addressing, worked
  new-motif examples.
- `reference/data-access.md` — InterPro/AFDB/UniProt endpoints, pagination and
  throttling behaviour, disk/time budgets, ESM Atlas and PDB alternatives.
- `reference/downstream-triage.md` — from hit list to experiment: SSN, genome
  neighbourhood, taxonomy, PaperBLAST, candidate selection, and the
  motif-directed protein-engineering corollary (swap the discriminator, swap the
  reaction).
