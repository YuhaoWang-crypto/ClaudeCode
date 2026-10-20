---
name: catalytic-scaffold-mining
description: >-
  Structure-first mining of predicted proteomes for unrecognised enzymes and
  functional proteins: use one known catalytic architecture as a 3D query
  (DALI / Foldseek against the AlphaFold Database), tier the hits by structural
  similarity, then keep only those conserving the catalytic residue
  constellation, verify the local active-site geometry, cluster the survivors
  into subfamilies, and place them on a phylogeny to nominate testable
  candidates. Use when asked to find enzymes/functional proteins homologous to
  a known catalytic core but invisible to sequence search (remote homology),
  to extend a published structure-guided discovery into new taxa, to triage
  AlphaFold hits into an experimental shortlist, or to audit such a claim.
  Enforces ✅-rigorous vs ⚠️-hypothesis labelling and a carryover/contaminant
  control checklist before any activity claim.
---

# Structure-first mining of catalytic scaffolds

A reusable pipeline that turns **one solved catalytic architecture** into a
**ranked, honestly-labelled shortlist of candidate enzymes** across whole
predicted proteomes — including organisms where no sequence homolog is
detectable.

Derived from Li, Yu, Li, Gu, Zhan & Shi, *Vita* **2026**, 10.0078, "*De novo*
RNA synthesis from free nucleotides by a family of conserved plant proteins",
which used human cGAS as a 3D query against AlphaFold-predicted plant
proteomes and recovered OsNTP4/OsNTP10/AtURT1 — proteins that share only
~15% sequence identity with the query.

## The core premise

> Remote evolutionary relationships are often invisible at the sequence level.
> Predicted structures now cover whole proteomes, so **catalytic architecture
> can be the search key instead of sequence.**

That single substitution — query by fold, filter by chemistry — is the whole
method. Everything else is making it not lie to you.

## When this applies (and when it does not)

| Use it when | Don't use it when |
|---|---|
| A catalytic core is structurally well-defined (NTase, kinase, protease, GT-A/B, Rossmann, metallohydrolase…) | The function is distributed/conformational with no compact catalytic motif |
| The activity you want is plausibly carried by a divergent scaffold | A plain BLAST/HMM search already finds the family cleanly — do that first |
| You can name the catalytic residues on the query | You cannot say which residues do the chemistry (go solve that first) |
| Predicted structures exist for the target taxa | Target taxa have no proteome in AFDB (then fold-predict first, or fall back to HMM) |

## The six stages

| Stage | Question | Method | Rigor |
|---|---|---|---|
| **0 Query definition** | what architecture, which residues do the chemistry? | pick an experimentally-solved structure + its catalytic set | ✅ from literature |
| **1 Structural sweep** | what in these proteomes looks like it? | DALI or Foldseek vs AFDB proteomes | ✅ deterministic scoring |
| **2 Similarity tiering** | which hits are signal, not fold-noise? | two thresholds: permissive (inventory) + strict (shortlist) | ✅ computed; ⚠️ threshold choice |
| **3 Chemistry filter** | does it keep the catalytic constellation? | map query catalytic residues through the structural alignment; require conservative identity | ✅ deterministic |
| **4 Geometry check** | is the active site actually superimposable? | local superposition over catalytic shell only; RMSD + carboxylate spacing | ✅ computed |
| **5 Family structure** | how many distinct subfamilies, who is a sibling of whom? | all-vs-all structure + sequence matrices → clusters; ML phylogeny → clades | ✅ computed |
| **6 Triage & labelling** | what is proven vs nominated? | rank, label every claim, specify the control ladder | ✅/⚠️ explicit |

```bash
# Stages 1-4, against the whole AlphaFold DB, no downloads:
python3 scripts/mine.py --query ref.pdb --chain A --catalytic 225,227,319 \
    --mode web --databases afdb-proteome afdb-swissprot --window 2 --out out/run

# Stages 1-4 against one downloaded proteome (adds side chains + pLDDT):
python3 scripts/mine.py --query ref.pdb --catalytic 225,227,319 \
    --mode local --target-dir proteomes/YEAST --out out/yeast

python3 scripts/annotate.py  out/run.candidates.tsv --out out/run.annotated.tsv
python3 scripts/calibrate.py out/run.annotated.tsv --positives ... --negatives ...

# Stage 5 (add --emit-fasta to the mine.py call above; --reuse makes it instant)
python3 scripts/cluster.py out/run.domains.fasta --table out/run.annotated.tsv \
    --references <validated accessions> --per-lineage 12 --out out/run
```

Stage-by-stage detail with the paper's exact numbers is in
`reference/methodology.md`; API recipes in `reference/data-access.md`; the
anti-self-deception controls in `reference/validation.md`.

## Stage 2 — the two-threshold rule

Never use one cutoff. Use a **permissive** threshold to size the landscape and
a **strict** threshold to build the shortlist.

The precedent run:

```
human cGAS (PDB 4KM5)  →  DALI vs AFDB (Zea mays, A. thaliana, O. sativa)
  Z ≥ 2.0   →  266 proteins     (inventory: "this fold exists here")
  Z ≥ 10.0  →   29 proteins     (shortlist: "this is the same architecture")
  + catalytic triad conserved → 19 proteins  (testable)
  → 3 assayed, 3 active
```

**Foldseek equivalence is a calibration, not an identity (⚠️).** DALI Z and
Foldseek TM-score are different statistics. If you substitute Foldseek for
DALI, *calibrate on known members of the target family first*
(`scripts/calibrate.py`) and report the mapping you derived rather than
importing numbers from this file.

## Stage 3 — the chemistry filter is what makes it work

Structural similarity alone is a fold statement, not a function statement.
The filter that converts one into the other:

1. Take the query's catalytic residues.
2. Map each through the **structural** alignment (not a sequence alignment) to
   its counterpart residue in the hit.
3. Require conservation **of chemistry, not identity** — in the precedent,
   Glu↔Asp substitution was accepted (both preserve an acidic environment),
   anything non-acidic was fatal.
4. Fail any one position ⇒ eliminate the protein.

In the precedent this removed 10 of 29 (34%). A filter that removes almost
nothing is mis-specified; a filter that removes almost everything means your
conservation class is too narrow.

**Expect misregistration.** Below ~20% identity a structural alignment
routinely places a catalytic motif one or two residues off — verified here:
aligning cGAS to AtURT1 puts Glu225 in a gap and lands Asp227 on the wrong
Asp. Use `--window 2` there, and keep the three safeguards the script
enforces: offsets are reported per position, exact (offset-0) matches are
counted separately, and the assignment is **unique and order-preserving** so
two catalytic positions cannot both claim the same nearby residue. That last
constraint is not optional; without it the filter manufactures conservation.

**Check your reference structure for engineered mutations** before trusting
its catalytic residues. Substrate-bound structures are routinely solved with
an active-site residue mutated to Ala — both AtURT1 crystals are D547A.

## Stage 4 — global RMSD will look bad; that is expected

The precedent's headline numbers:

- global superposition OsNTP4 vs cGAS: **RMSD 4.623 Å over 1,182 atoms** — poor
- catalytic region only: **RMSD 0.735 Å over 66 atoms** — excellent

**Judge the active site, not the protein.** A divergent scaffold with a
conserved catalytic core is exactly the thing this method exists to find, and
it is precisely the thing a global-RMSD cutoff throws away.

**Trim outliers, or a flexible loop ruins the number.** An activation loop or
lid beside the catalytic residues inflates the shell RMSD while the core
superimposes perfectly. On a verified self-match (cGAS crystal vs its own
AlphaFold model, all 362 positions mapped correctly) the whole shell gave
3.31 A over 53 pairs while the trimmed core gave **1.33 A over 44**. Rank by
`core_rmsd`, and always read `n_core` beside it — a tight fit over 15 residues
is not a tight fit over 50.

## Stage 5 — clade membership is the strongest prior you get

After clustering, ask: **do all experimentally-validated actives fall in one
clade?** If yes, uncharacterised members of that clade are your best
candidates, and members of sister clades are the controls that tell you
whether the activity is clade-restricted or family-wide. In the precedent all
three actives landed in Clade I of five.

This is the highest-yield output of the whole pipeline — rank by it.

## The non-negotiable discipline: honesty labelling

Every claim carries exactly one label:

- **✅ rigorous** — a deterministic computation, a real measured statistic, or
  a decimal-exact reproduction of published values.
- **⚠️ hypothesis** — a plausible but unvalidated claim: a nominated candidate,
  a threshold analogy, an inferred function, an assumed mechanism.

**A structural hit is never evidence of activity.** The precedent paper is the
standard to copy — about its own 16 remaining candidates it wrote:

> "The present data do not establish *de novo* RNA synthesis for these
> additional candidates: structural similarity and conservation of an acidic
> carboxylate configuration nominate them for functional analysis but do not
> define their substrates, products, or biological roles."

Copy that register. A computational shortlist is a hypothesis generator, and
saying so is the finding, not a disclaimer.

Report negative and partial results as results: folds that were present but
lost the catalytic residues, taxa that returned nothing, thresholds that
failed to separate signal from noise.

## Before anyone claims activity: the carryover ladder

A purified protein that appears to make a polymer is guilty until proven
innocent. The precedent ran five orthogonal controls — reproduce them, or
label the activity ⚠️ (details and rationale in `reference/validation.md`):

1. **Intact-protein MS** — exclude a covalently bound nucleotide/cofactor.
2. **Protease digestion + ion exchange** — exclude free NTP/oligonucleotide
   carryover in the prep.
3. **Multi-step orthogonal purification** — affinity → ion exchange → SEC, to
   strip co-purifying enzymes.
4. **Washed-crystal assay** — recover activity from thrice-washed crystals
   with a matched non-crystalline control from the same drop. This is the
   decisive one: it is very hard for a contaminant to co-crystallise.
5. **Active-site mutants** — each catalytic residue to Ala must abolish
   activity, while a surface residue must not.

Chemical identity of the product needs its own orthogonal proof (the
precedent used NMR COSY/HMBC to establish the 3′–5′ linkage; a nuclease
sensitivity assay alone does not establish bond chemistry).

## Output contract

Produce a table with one row per candidate and these columns, nothing hidden:

`accession | organism | lineage | structure score (+ method) | aligned coverage |
catalytic residues found (+ positions, + registration offsets) |
core active-site RMSD + n_core | clade | prior characterisation |
model confidence (pLDDT over catalytic shell) | label ✅/⚠️ | what would test it`

Candidates with low pLDDT **over the catalytic shell specifically** are
unreliable regardless of global model quality — report that column separately
and never average it away.
