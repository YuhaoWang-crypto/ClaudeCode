# 4Δ-Grafted IgG1 Variants — Computational Design and *In Silico* Evaluation

**Scope.** Structural modelling, candidate design and evaluation of IgG1 constructs
carrying the 4Δ integrin-binding peptide `FETLRGDLRILSILRHQNLLKEL` (23 aa) in the
CDRH3 loop and in the AB, CD and EF loops of CH3 — the four sites named in the
proposal — followed by ranking and down-selection.

**Status of this document.** Every number below was computed in this project from
experimental structures and new structure predictions. Scripts, raw metrics and
model files are in this repository; each claim is traceable to a `results/*.json`
file. Confidence levels are stated explicitly in §8.

---

## 1. Is the project feasible?

**Yes as a piece of computational work — it is executed in full here.** The
scientific verdict on the *design hypothesis* is more interesting than a yes/no,
and is the substance of this report:

> Grafting 4Δ into a **closed loop** of IgG1 is geometrically disfavoured, for two
> independent and quantified reasons: a **28 Å span mismatch** (§3) and a **6.1 %
> receptor escape cone** (§5). Both are properties of the peptide and of αvβ6, not
> of any particular design, so they constrain every one of the 16 candidates.
> Designs that resolve one tend to worsen the other (§6). A ranked set of 12 is
> nevertheless delivered (§7), together with the single topology that escapes both
> constraints and which we recommend be carried alongside (§7.3).

---

## 2. What 4Δ actually is, and why that matters

4Δ is not an arbitrary RGD peptide. Two facts, both established here from primary
structural data, determine the entire design problem.

### 2.1 Identity and the role of the E46L substitution

| | sequence | relationship |
|---|---|---|
| human chromogranin A 39–63 | `FETLRGDERILSILRHQNLLKELQD` | PDB **6R2X**, solution NMR, 15 models |
| **4Δ** | `FETLRGDLRILSILRHQNLLKEL` | CgA(39–61) with **E46L**, C-terminally truncated |

So 6R2X is a near-exact experimental template for 4Δ: it differs only by E46L and
the Q62/D63 truncation. This is the reference conformation used throughout.

The E46L substitution is not cosmetic. The αvβ6 co-crystal **4UM9** contains the
pro-TGF-β1 LAP peptide `GRGDLGRL` bound in the αv/β6 cleft — the canonical
**RGDLxxL** specificity motif of αvβ6/αvβ8. Native CgA reads `RGD-E-RIL`; 4Δ reads
**`RGD-L-RIL`**. E46L therefore converts the native chromogranin sequence into a
canonical RGDLxxL αvβ6-binding motif. The pharmacophore to be preserved by any
graft is accordingly **R43–G44–D45–L46–R47–I48–L49**, not merely "RGD".

### 2.2 The bioactive conformation is a linear rod

From the 6R2X ensemble (`results/01_baseline_4delta.json`, Figure 1):

| property | value |
|---|---|
| helical core | CgA **E46–N56**, helix fraction **1.00** across all 15 models |
| RGD mobility | backbone RMSF **1.8 / 0.7 / 0.5 Å** for R43 / G44 / D45 (helix-core aligned) |
| RGD exposure | mean relative SASA **0.577** |
| Arg-CZ···Asp-CG | **9.0 ± 2.7 Å** |
| **end-to-end Cα span** | **35.8 ± 2.1 Å** (range 29.9–38.2) |

4Δ is a **mobile RGD on the end of a rigid amphipathic helix**, with its two termini
~36 Å apart. That span is the number every graft has to reconcile.

### 2.3 Binding requires the RGD to unwind out of the helix

Comparing the free-state motif (6R2X, CgA 43–49) with the receptor-bound motif
(4UM9 LAP 241–247):

| | DSSP over the 7-residue motif | |
|---|---|---|
| free 4Δ | `SHHHHH` / `HHHHHH` | RGD is **inside** the helix |
| bound LAP | `STTHHH` | RGD is in a **turn**; helix starts at +4 |
| backbone RMSD, free vs bound | **1.77–2.43 Å** (mean **2.08 Å**) | a real local conformational change |

**Design consequence.** αvβ6 engagement requires the RGD to unwind locally out of
the 4Δ helix. A graft must therefore leave the RGD *conformationally permissive*
while holding the downstream helix intact. Rigidifying the RGD — by a short,
strained junction, or by packing it against the scaffold — removes binding
competence even when the helix looks perfectly preserved.

This value, **2.08 Å**, is also the calibration baseline for every "motif fit"
number in this report: the *experimental free peptide itself* sits 2.08 Å from the
bound pose. A candidate near 2.08 Å is behaving like free 4Δ. Zero is not the target.

---

## 3. The graft sites, and the span mismatch

Loop boundaries were taken from the Fcab literature (AB 355–362, CD 383–391,
EF 413–422, EU numbering) and then **validated on the human IgG1 Fc crystal
structure 3AVE** by DSSP, solvent accessibility and distance from CH2
(`results/02_graft_sites.json`). CDRH3 was taken from the trastuzumab Fab/HER2
complex **1N8Z**.

| site | EU range | loop | anchor Cα–Cα | loop SASA | distance from CH2 | at CH3/CH3 dimer interface |
|---|---|---|---|---|---|---|
| CDRH3 | H99–109 | `WGGDGFYAMDY` | **6.9 Å** | 278 Å² | — | — (contacts HER2 at H99/103/104/105) |
| CH3-AB | 355–362 | `RDELTKNQ` | **8.5 Å** | 702 Å² | 48.3 Å | **yes — EU 356, 357, 360** |
| CH3-CD | 383–391 | `SNGQPENNY` | **13.0 Å** | 712 Å² | 35.0 Å | marginal — EU 390 |
| CH3-EF | 413–422 | `DKSRWQQGNV` | **9.8 Å** | 674 Å² | 46.9 Å | **none** |

All three CH3 loops sit 35–48 Å from the CH2 centroid, confirming they are the
C-terminal-face loops of the domain. **EF is the only one entirely clear of the
CH3/CH3 dimer interface**; AB is the most compromised, placing three interface
residues inside the graft region.

### The mismatch (Figure 2)

| site | anchor span | slack a cassette must absorb |
|---|---|---|
| CDRH3 | 6.9 Å | **+28.9 Å** |
| CH3-AB | 8.5 Å | **+27.3 Å** |
| CH3-CD | 13.0 Å | **+22.8 Å** |
| CH3-EF | 9.8 Å | **+26.0 Å** |

A Gly/Ser linker reaches ~3 Å per residue, so closing a 23–29 Å gap needs ~8–10
flexible residues *in addition to* the 23-residue payload. This single calculation
drove the design matrix.

---

## 4. The 16 candidates

4 sites × 4 graft strategies, each strategy a different answer to the span problem
(`results/03_candidates.json`; sequences in `results/sequences/`).

| strategy | construction | hypothesis under test |
|---|---|---|
| **S1 minimal** | `GS` + 4Δ + `SG`, native loop kept | naive insertion; no slack budget (reach 6 Å vs 23–29 Å needed) |
| **S2 span-matched** | `GS` + 4Δ + `(GGGGS)×2`, native loop kept | C-return linker sized to the measured slack |
| **S3 replace** | `GD` + 4Δ + `GPG(GGGGS)`, native loop deleted | maximum room; helix N-capped, β-strands protected by a GPG break |
| **S4 core-truncated** | `GS` + 4Δ(1–18) + `(GGGGS)×2` | rigid helix core only — least strain, full RGDLxxL epitope retained |

**Controls and benchmark.** `WT_Fv` and `WT_CH3` (unmodified scaffolds predicted
identically), `REF_4delta_free` (free peptide — validates the pipeline against
6R2X), and `BM_CH3_Cterm_fusion` (4Δ appended to the free C-terminus K447 via
`GGGGS`). The benchmark is **not** one of the 16: it has no closed-loop span
penalty and shows the attainable ceiling.

### Modelling

Boltz-2.1 (`v2026-03-01`), **5 samples per construct**, automatic MSA, **no
structural templates** — the scaffold was deliberately left free so that
graft-induced distortion is measurable rather than suppressed. CDRH3 designs were
modelled as the Fv (VH+VL); CH3 designs as the **CH3 homodimer**, since Fc
function depends on that dimer. Job IDs: `results/04_job_ids.json`.

*Scaffold assumption:* client heavy/light chain sequences were not supplied, so
**trastuzumab** was used as the reference human IgG1. All designs are expressed as
**EU-numbered edits** (`results/sequences/edit_table.tsv`) and transfer unchanged
to the client's own IgG1.

---

## 5. Can an integrin physically reach the graft?

Confidence scores alone cannot answer this, so three structure-based tests were
used instead.

### 5.1 The escape cone — a project-level constraint

The αvβ6 RGD site is a deep cleft between the αv β-propeller and the β6 βI domain.
Sampling 2000 directions out of the bound RGD position in 4UM9 and asking along how
many a **domain-sized body (6 Å radius) can travel 25 Å without hitting the
receptor**:

> **Only 6.1 % of directions are open (122 / 2000).**

The pro-TGF-β1 LAP chain's own exit direction lies **1.8°** from an open direction —
nature threads a whole domain through that narrow cone. So presenting 4Δ on a large
scaffold is *possible but tightly constrained*: the graft's exit vector must be
aimed into a ~6 % solid-angle window. This is a property of αvβ6, and it applies to
any antibody format, not just these 16 designs.

### 5.2 Exit-vector compatibility and motif fit (Figure 3)

Each predicted model's RGDLRIL backbone was least-squares superposed onto the LAP
RGDLGRL backbone in the 4UM9 complex, placing the antibody module exactly where it
would sit if its RGD engaged the receptor. Two numbers follow: how well the grafted
motif matches the bound conformation (**motif fit**, calibrated against 2.08 Å), and
how far the scaffold points off the nearest open direction (**exit angle**).

*(Per-candidate values: §7 and `results/06_integrin_steric.json`,
`results/07_exit_vector.json`.)*

---

## 6. What the evaluation showed

*(This section is completed once all 20 constructs have been scored — see §7.)*

---

## 7. Ranking and down-selection

*(Completed below.)*

---

## 8. Confidence, limitations and what would de-risk this

*(Completed below.)*
