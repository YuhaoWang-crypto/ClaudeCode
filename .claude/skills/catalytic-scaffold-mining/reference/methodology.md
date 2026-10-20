# Methodology, stage by stage

The worked precedent throughout is Li *et al.*, *Vita* 2026, 10.0078 — human
cGAS used as a 3D query against AlphaFold-predicted plant proteomes, yielding
OsNTP4 / OsNTP10 / AtURT1 as enzymes that synthesise RNA from free NTPs with
no primer and no template.

Their exact numbers are quoted as a **calibration reference**, not as
thresholds to copy blindly into a different scaffold or a different search
tool.

---

## Stage 0 — define the query

You need two things, and the second is the one people skip:

1. **A structure of the catalytic architecture.** Experimentally solved is
   preferable; a high-confidence predicted model is acceptable if you check
   pLDDT over the active site.
2. **The catalytic residues, as author residue numbers in that file.**

In the precedent: human cGAS, PDB **4KM5**, catalytic triad **Glu225, Asp227,
Asp319**. cGAS was chosen because "its nucleotidyltransferase core constitutes
a well-defined architecture for nucleotide binding and phosphodiester-bond
formation" — i.e. the query was picked for having a *crisply defined chemical
job*, not for being famous.

### Two traps at this stage

**Trap 1: your reference structure may be a catalytically dead mutant.**
Substrate-bound structures are routinely solved with an active-site residue
mutated to Ala so the enzyme cannot turn over. Both AtURT1 structures
(**6L3F**, **6L8K**) carry `SEQADV ... ALA A 547 ... ASP 547 ENGINEERED
MUTATION`. If you take catalytic residues from such a file without checking,
you will "discover" that the wild-type enzyme lacks its own catalytic residue.
Always read `SEQADV` / `_struct_ref_seq_dif`, or compare against the UniProt
sequence.

**Trap 2: author numbering is not sequence position.** Crystal structures
start at the first ordered residue (4KM5 chain A starts at 161) and contain
gaps. The pipeline converts author numbers to chain positions internally; do
not do that arithmetic by hand.

---

## Stage 1 — the structural sweep

The precedent used **DALI** against **AlphaFold DB** proteomes for *Zea mays*,
*Arabidopsis thaliana* and *Oryza sativa*.

Three practical options today:

| Route | Coverage | Cost | Gives you |
|---|---|---|---|
| DALI server | pick your own target set | slow, queued | Z-score (the published statistic) |
| Foldseek web (`search.foldseek.com`) | all of AFDB, server-side | minutes, no download | E-value, prob, inline target Calpha + sequence |
| Foldseek local binary | one proteome you downloaded | fast, needs disk | full side chains, pLDDT, TM-score |

`scripts/fsweb.py` drives the web route; `scripts/mine.py --mode local` drives
the local one.

### Choose the database deliberately

- `afdb50` — broadest (clustered whole AFDB), but capped at `--max-accept 1000`
  and the cap fills with the query's **closest** relatives. A cGAS query
  returns ~1000 metazoan cGAS/cGLR proteins and never reaches plants or fungi.
- `afdb-proteome`, `afdb-swissprot` — smaller, curated, **annotated**, far more
  taxonomically even. Much better for a cross-kingdom inventory.
- A locally built DB from one downloaded proteome — the only route with no cap
  at all, and the only faithful reproduction of a per-proteome DALI sweep.

**Breadth comes from several queries and several databases, never from one
big search.** This is the single most important operational point in Stage 1.

### Getting structures when bulk endpoints are blocked

`alphafold.ebi.ac.uk/files/...` may be unavailable from a sandbox while
`alphafold.ebi.ac.uk/api/...` works. Fallbacks, in order:

1. Foldseek web results already carry `tCa` (full Calpha trace) and `tSeq` for
   every hit — enough for the chemistry filter and the geometry check.
2. `https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/` serves per-proteome
   tar archives (`UP..._v6.tar`).
3. RCSB for anything experimentally solved.

Long transfers get cut; always resume (`curl -C -`) and verify against
`Content-Length` rather than trusting a single attempt.

---

## Stage 2 — tier the hits

Two thresholds, always.

Precedent (DALI Z, three plant proteomes):

| Threshold | Survivors | Reading |
|---|---|---|
| Z >= 2.0 | 266 (100 maize, 90 *A. thaliana*, 76 rice) | the fold exists here |
| Z >= 10.0 | 29 (13 At, 10 rice, 6 maize) | same architecture |
| + catalytic triad conserved | 19 | testable |
| assayed | 3 of 3 active | |

Observed sequence identity to the query across the survivors: **~15%**. That
is the whole justification for the method — no sequence search finds these.

### Foldseek thresholds are yours to calibrate

Foldseek reports no DALI Z. Do not invent a conversion. Calibrate against
known members of the family you are chasing (`scripts/calibrate.py`), then
report the mapping you derived.

Useful reference points measured on this pipeline (cGAS query vs
`afdb-proteome` + `afdb-swissprot`, catalytic-shell **core RMSD**, see Stage 4):

| Class of hit | core RMSD | n_core |
|---|---|---|
| mammalian cGAS | 1.19 - 1.36 A | 42 - 44 |
| invertebrate cGAS-like receptors (cGLRs) | 1.27 - 2.0 A | 25 - 39 |
| bacterial CD-NTases (cyclic di-nucleotide synthases) | 1.36 - 1.67 A | 23 - 27 |
| MB21D2 (same superfamily, different job) | 7.4 - 7.7 A | 36 - 46 |
| ribitol-5-P transferase FKTN (fold match, unrelated chemistry) | 6.6 - 6.8 A | 49 |

The gap between ~2 A and ~5.5 A is wide and clean. **core RMSD <= 2.5 A with
n_core >= 20** separated every genuine family member from every fold-only
coincidence in that run. Treat it as a starting point and re-derive it for
your own query.

---

## Stage 3 — the chemistry filter

Map the query's catalytic residues through the **structural** alignment and
require conservation of chemical class.

Precedent: cGAS Glu225/Asp227/Asp319 -> plant candidates must keep three
acidic residues. "Substitution between Glu and Asp preserves an acidic
environment while allowing divergence in local active-site geometry." Ten of
29 violated this and were dropped, leaving 19.

The surviving plant proteins showed two conserved sequence blocks:
a **D(phi)D(phi)** motif (phi = Ile/Leu/Val) matching cGAS Glu225/Asp227, and a
**D(phi)** motif matching cGAS Asp319.

### Registration: the failure mode nobody warns you about

At ~10-15% identity a structural alignment routinely places a catalytic motif
**one or two residues off**. Verified here: aligning cGAS (4KM5) to AtURT1
(6L8K) puts Glu225 in a gap column and maps Asp227 onto **D491**. The
biologically correct correspondence, from the DxD motif
`...FGFPKS D491 I D493 V...` and `...GISC D547 IC...`, is:

```
cGAS E225  ->  AtURT1 D491     (+2 registration offset)
cGAS D227  ->  AtURT1 D493     (+2 registration offset)
cGAS D319  ->  AtURT1 D547     (Ala in 6L3F/6L8K: engineered D547A)
```

So the pipeline supports `--window N`: look for an in-class residue within +-N
of the aligned position. Three rules make it safe rather than sloppy:

1. **Offsets are always reported** (`max_offset`, and `[+2]` in the residue
   column). An offset-0 match is stronger evidence than an offset-2 match, and
   `n_conserved_exact` is tracked separately.
2. **Assignment is unique and order-preserving.** Without this constraint two
   catalytic positions will both claim the same nearby Asp — which looks like
   conservation and is an artefact. `_assign()` searches all combinations for
   the one maximising matches, then minimising total offset, subject to
   distinct and strictly increasing target residues.
3. **Window 0 for a close query, 2 below ~20% identity.** Never more than 3.

---

## Stage 4 — active-site geometry

Superpose **the catalytic shell only**: aligned residue pairs whose query
Calpha lies within `--shell-radius` (default 12 A) of any catalytic Calpha.

Precedent, OsNTP4 (AlphaFold) vs human cGAS (4KM5):

- global: **RMSD 4.623 A over 1,182 atoms**
- catalytic region: **RMSD 0.735 A over 66 atoms**

A global cutoff would have discarded the discovery. Judge the site.

### Trim outliers or a flexible loop will ruin the number

A loop next to the active site — an activation loop, a lid, or just a
crystal-vs-model difference — drags the whole-shell RMSD up while the core
superimposes perfectly. Measured on the cGAS self-hit (4KM5 crystal vs its own
AlphaFold model, a mapping verified correct at all 362 positions):

```
whole shell : 3.31 A over 53 pairs
core (trimmed, 4 A cutoff) : 1.33 A over 44 pairs  (83% kept)
per-residue deviation: median 0.62 A, p90 2.26 A, max 11.63 A
```

The outliers were residues 212-218 and 364-369 — the cGAS activation-loop
region. `core_rmsd()` iteratively drops pairs beyond 4 A and re-fits; both
numbers are reported, and **`core_rmsd` is the one to rank by**, with `n_core`
alongside so a tight fit over 12 residues is not mistaken for a tight fit over
40.

### Also check, where the data allow

- **Carboxylate spacing** — distance between catalytic side-chain carboxylate
  carbons (Asp CG / Glu CD) in target vs query. A conserved two-metal site
  keeps these tight. Needs side chains, so local mode only.
- **pLDDT over the catalytic shell specifically.** A model can be 85 overall
  and 50 at the active site. Never average this away; report it as its own
  column. Not available in web mode — report as missing, do not guess.

---

## Stage 5 — family structure and phylogeny

The precedent did three things here:

1. **All-vs-all pairwise sequence identity** across the 19 survivors — within
   subgroups 40-80%, versus ~15% to the query.
2. **All-vs-all structural Z-scores** — broadly recapitulating the
   sequence-based grouping, which is itself the evidence for functional
   diversification. Four subgroups.
3. **Maximum-likelihood phylogeny** over all 47 NTP-family proteins from the
   three species: MEGA 7.0.26, JTT substitution model, gamma-distributed rates
   in 5 categories, partial deletion at 95% site coverage, 50 bootstrap
   replicates, nearest-neighbour-interchange search. Result: **five clades**.

**All three validated enzymes fall in Clade I.** That single fact is the
highest-value output of the whole pipeline: it converts "19 structural
candidates" into "prioritise Clade I, and use sister clades as controls".

Rank by clade membership before you rank by RMSD.

---

## Stage 6 — report

One row per candidate, every claim labelled. See the output contract in
`SKILL.md` and the control ladder in `validation.md`.

The precedent's own wording about its unassayed candidates is the register to
copy:

> "The present data do not establish *de novo* RNA synthesis for these
> additional candidates: structural similarity and conservation of an acidic
> carboxylate configuration nominate them for functional analysis but do not
> define their substrates, products, or biological roles."

---

## What the method cannot do

- It finds **architecture**, not activity. Every output is a hypothesis.
- It is blind to proteins with no AlphaFold model, and to activities that do
  not live in a compact catalytic motif.
- A conserved constellation can be structural rather than catalytic — some
  acidic triads bind metal for stability, not turnover.
- Hit caps bias results toward the query's closest relatives; absence of a
  taxon from your result list is **not** evidence of absence in that taxon.
  Say so explicitly rather than implying a negative result.
