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

Crucially, the openings are **at both ends of the motif**: the LAP chain's own
N-side and C-side directions lie **1.4°** and **1.8°** from open directions. The
RGD site is a *groove open at both ends*, and nature threads the pro-TGF-β1 chain
straight through it.

Two consequences, and they point in opposite directions:

- **In favour of loop grafting.** Because both ends are open, tethering the peptide
  at both termini — which is what a loop graft does — is *not* topologically
  forbidden. A loop graft can in principle thread the groove the way LAP does.
- **Against.** The open window is only ~6 % of solid angle, and the scaffold body
  must be carried clear of the cleft (~25 Å) before it widens. Both junctions must
  therefore act as accurately aimed spacers simultaneously. Random loop insertion
  will miss this; it has to be designed for.

This constraint belongs to αvβ6, not to any particular construct, and applies to
every antibody format.

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

All 20 constructs were modelled (5 samples each) and scored on eight structural
axes grouped into four families (`results/05_evaluation.json`,
`06_integrin_steric.json`, `07_exit_vector.json`, `09_scaffold_integrity.json`).
Several results are worth stating plainly because they overturn the naive intuition.

**6.1 The scaffold fold always survives; the peptide is the casualty.**
Scaffold Cα-RMSD to the wild-type module stayed **< 1.5 Å for every candidate**
(WT-identical flanks only), so no graft unfolds the host domain. But the graft
region's own pLDDT fell to **46–68** against **86.6** for the free peptide: the
peptide, not the antibody, loses conformational definition on grafting — exactly
what "a 36 Å rod forced into a ~10 Å closed loop" predicts.

**6.2 The span collapses, as predicted.** End-to-end Cα span fell from the
free-peptide 33.4 Å to **5.7–23.5 Å** in the grafts. The peptide folds back on
itself inside the loop; the only question is whether it does so while keeping the
RGD exposed (good) or by burying it (fatal).

**6.3 RGD burial, not helix loss, is the main failure mode.** Helix core was
largely retained (0.71–0.85) except where a strategy deliberately or accidentally
destroyed it — **C08 (helix 0.00)**, **C04 (0.18)**, **C16 (0.38)**, **C15 (0.49)**.
The more discriminating axis was RGD accessibility to a *bulky* partner: span-matched
S2 designs (C02 0.18, C14 0.20, C10 0.29) buried the motif against the returning
linker, and a 5 Å-probe SASA confirmed the RGD was no longer reachable by anything
integrin-sized. **The long linker that fixes the span problem creates the burial
problem.**

**6.4 Exit vector and motif fit trade off against each other.** The candidates with
the best exit angle (C13 1.1°, C03 1.8°, C05 2.9°, C01 3.3°) had mediocre motif fit
(~2.8 Å), while those with the best motif fit (C06 1.82 Å, C10 1.85 Å, C02 1.90 Å)
pointed the scaffold 8–83° off any open direction (Figure 3). Only a handful sit
acceptably on both.

**6.5 The AB loop damages the Fc dimer, as flagged structurally.** CH3/CH3 interchain
contacts (WT = 316) ballooned for AB-loop grafts — **C06 702, C08 855** — confirming
the §3 warning that the AB loop overlaps the dimer interface (EU 356/357/360). CD
and EF grafts stayed much closer to WT.

**6.6 Confidence scores are useless for binding; chemistry is not (Figure 6).**
The Tier-2 integrin complexes (`results/08_binding_mode.json`) gave the cleanest
result of the whole study. Boltz reproduced the αvβ6/LAP co-crystal to **0.91 Å**
(method control: Asp→MIDAS metal 2.23 Å vs crystal 2.13 Å; Arg→αv Asp218 3.45 Å vs
3.61 Å). The dead RGE mutant lost metal coordination entirely (Glu→metal 10.4 Å).
Yet ipTM could not tell them apart: **0.943 (active 4Δ) vs 0.944 (dead mutant)**.
Only the measured chemistry separates binder from non-binder — which is why the
scoring here is built on geometry, not on confidence.

A corollary worth noting: free 4Δ itself engaged αvβ6 in **only 1 of 5 samples**,
with the Arg salt bridge never fully formed. The free peptide is conformationally
loose — which is precisely why the original authors *stapled* it. A rigid scaffold
could in principle supply the same pre-organisation, which is an argument *for* the
grafting concept, provided the geometry is right.

---

## 7. Ranking and down-selection

Two scores are reported (`results/10_ranking.json`, Figures 4–5): a pure
**structural score (0–100)** and a **total score** after programme-level deductions
(insert burden, loss of the parent paratope for CDRH3, dimer-interface risk for AB,
span-incompetent junctions for S1). A candidate must additionally **pass five gates**
— helix core kept, RGD not buried, RGD reachable by a bulky partner, no excess
internal clash, Fc interface preserved — to be called *advanceable*.

### 7.1 The down-selected 12 (ranked)

| rank | id | site | strategy | total | exit ° | motif fit Å | helix | RGD SASA | iface ×WT | gates |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | **C09** | CH3-CD | S1 minimal | **63.0** | 12.4 | 2.62 | 0.84 | 0.50 | 0.97 | **PASS** |
| 2 | **C11** | CH3-CD | S3 replace | **62.5** | 14.1 | 2.62 | 0.82 | 0.58 | 1.02 | **PASS** |
| 3 | **C01** | CDRH3 | S1 minimal | 57.9 | 3.3 | 2.86 | 0.84 | 0.57 | 1.39 | **PASS** |
| 4 | **C05** | CH3-AB | S1 minimal | 57.3 | 2.9 | 2.77 | 0.76 | 0.41 | 1.12 | **PASS** |
| 5 | C15 | CH3-EF | S3 replace | 55.5 | 15.5 | 2.79 | 0.49 | 0.65 | 0.95 | fail: helix |
| 6 | C06 | CH3-AB | S2 span | 55.3 | 8.2 | 1.82 | 0.85 | 0.38 | 2.23 | fail: interface |
| 7 | **C03** | CDRH3 | S3 replace | 53.6 | 1.8 | 2.83 | 0.76 | 0.49 | 1.32 | **PASS** |
| 8 | C13 | CH3-EF | S1 minimal | 52.4 | 1.1 | 2.77 | 0.82 | 0.28 | 1.40 | fail: RGD buried |
| 9 | C04 | CDRH3 | S4 core | 51.9 | 19.8 | 1.29 | 0.18 | 0.63 | 1.62 | fail: helix, interface |
| 10 | **C12** | CH3-CD | S4 core | 50.4 | 62.9 | 2.08 | 0.71 | 0.41 | 1.00 | **PASS** |
| 11 | C10 | CH3-CD | S2 span | 49.7 | 48.1 | 1.85 | 0.84 | 0.29 | 0.96 | fail: RGD buried |
| 12 | C16 | CH3-EF | S4 core | 45.8 | 35.3 | 2.27 | 0.38 | 0.48 | 1.50 | fail: helix |

**Excluded (ranks 13–16):** C07 (RGD unreachable), C02 (RGD buried + dimer damage),
C14 (RGD buried, unreachable, dimer damage), C08 (fails four of five gates — helix
destroyed, RGD buried, dimer wrecked).

### 7.2 The headline results

- **CH3-CD is the best site.** Three of the four CD designs made the top 12 and
  **C09, C11 and C12 pass every gate** — the CD loop is long enough and far enough
  from the dimer interface to host the peptide without collateral damage.
- **Minimal junctions (S1) outperform span-matched linkers (S2).** This is
  counter-intuitive but mechanistically clear: S1 does not bridge the 23–29 Å gap,
  but it also does not fold the peptide back over its own RGD. S2 closes the span
  and buries the motif. The span mismatch is better left unclosed than closed badly.
- **Six candidates are genuinely advanceable** (pass all gates): **C09, C11, C01,
  C05, C03, C12**. The remaining six in the "12" are delivered as ranked but are
  flagged with their specific liability.
- **The C-terminal fusion benchmark scored 45.3** — *below* the best loop grafts.
  Appending 4Δ to K447 removes the span problem but points the peptide 61.6° off any
  open direction and buries the RGD (SASA 0.18). Simple fusion is **not** a free win.

### 7.3 Recommendation

**Primary: advance C09 and C11 (CH3-CD).** They are the top two by total score, both
pass all gates, both keep the helix, expose the RGD, and preserve the Fc dimer.
They differ in strategy (minimal insertion vs loop replacement), giving two
mechanistically distinct shots on goal at the same site.

**Secondary: C01 and C03 (CDRH3) and C05 (CH3-AB)** as a format-diversity hedge —
all gate-passing, but CDRH3 sacrifices the parent paratope and AB carries residual
dimer risk, hence the deductions.

**Carry one bivalent control: C12 (CH3-CD, core-truncated).** It passes all gates
with the full RGDLxxL epitope and the cleanest interface ratio (1.00), and tests
whether the rigid-core-only design buys back exit-vector freedom.

**Run the C-terminal fusion (BM) in parallel as a format comparator**, not because
it scored well — it did not — but because it is the simplest thing to build and the
data say a well-placed loop graft should beat it, which is a hypothesis worth
settling experimentally.

> **Bottom line.** The concept is viable but geometrically demanding. The deliverable
> is the ranked 12 plus a mechanistic account of *why* CD-loop minimal/replacement
> grafts are the way to present 4Δ on IgG1, and *why* the obvious moves (long
> linkers, C-terminal fusion) underperform.

---

## 8. Confidence, limitations, and what would de-risk this

**What is rigorous here (✅):**
- 4Δ identity, the E46L→RGDLxxL rationale, and the bioactive conformation — all from
  experimental structures (6R2X, 4UM9, 5FFG, 3AVE, 1N8Z).
- The span mismatch and the 6.1 % escape cone — geometric facts, independent of any
  model.
- The binding-mode method validation — Boltz reproduced the αvβ6/LAP co-crystal to
  0.91 Å and correctly failed the RGE mutant.
- Relative ranking of the 16 designs on structure/sterics — internally consistent,
  calibrated against WT and free-peptide controls predicted identically.

**What is a hypothesis, not a measurement (⚠️):**
- **Absolute binding affinity is not computed and must not be inferred.** As the
  proposal states, reported values for parent Compound 4 do not transfer to 4Δ or to
  these constructs. ipTM is explicitly shown here to be blind to RGD binding (§6.6).
- Boltz is a single predictor; no orthogonal method (AF3, physics-based docking, MD)
  was run to cross-check individual poses.
- No full-length IgG1 was modelled — Fv and CH3 modules stand in, which is standard
  but omits long-range Fab/Fc context and avidity.
- Client heavy/light sequences were not provided; trastuzumab is a placeholder. The
  **EU-numbered edit table transfers the designs unchanged** once the real sequences
  arrive (`results/sequences/edit_table.tsv`).
- Developability (aggregation, expression, thermostability, immunogenicity of the
  junctions) was out of scope.

**What would de-risk this cheaply, in order:**
1. Re-run the top 6 as **4Δ-graft + αvβ6 headpiece complexes** (the Tier-2 pipeline
   is already built and validated) to test engagement directly rather than by
   superposition — ~$0.2/construct.
2. Orthogonal structure prediction (AF3) on the top 6 to confirm the graft poses.
3. Short MD on C09/C11 to test persistence of RGD exposure and helix.
4. Only then, express C09 + C11 (+ C01 as format control) and run an αvβ6/αvβ8
   binding assay.

---

## Appendix — methods and files

- **Structure prediction:** Boltz-2.1 (`v2026-03-01`), 5 samples/construct, automatic
  MSA, no templates. Pipeline cost ≈ \$4.3 total (20 single/dimer jobs ≈ \$3.3 + 3
  integrin complexes ≈ \$1.0), i.e. ~0.02 % of the proposal's compute line.
- **Experimental references:** 6R2X (4Δ conformation), 4UM9 + 5FFG (αvβ6 + RGDLxxL
  bound state), 3AVE (IgG1 Fc / CH3), 1N8Z (trastuzumab Fab / HER2).
- **Scripts:** `scripts/01`–`12` (baseline → sequences). **Results:** `results/*.json`.
  **Figures:** `figures/fig1`–`fig6`. **Models:** `models/<id>/sample_*.cif` (115
  structures; not committed — regenerable from `results/04_job_ids.json`).
- **Figures:** (1) 4Δ baseline conformation; (2) span reconciliation; (3) exit-vector
  vs motif-fit; (4) evaluation matrix; (5) ranking; (6) binding-mode validation.
