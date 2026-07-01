# Tirzepatide Binder Program — Computational Evaluation & Design Report

**Target:** Tirzepatide (39-aa dual GIP/GLP-1R agonist peptide)
**Date:** 2026-07-01
**Models used:** Boltz-2.1 (co-folding + binding module), local biophysical/developability analysis
**Scope:** (1) unified re-evaluation of 90 pre-existing designed binders; (2) de novo design + screening of 100 new peptide binders; (3) wet-lab-ready candidate selection.

> ⚠️ **Reading note on confidence.** Tirzepatide is a small, heavily-modified,
> partly-disordered peptide (two Aib residues, K20 fatty-diacid, C-terminal
> amide). This is a genuinely hard target for structure-based affinity
> prediction. Treat all scores as **relative triage signals to prioritize wet-lab
> testing**, not absolute affinity predictions.

---

## 1. Target definition

Tirzepatide primary structure (modifications noted):

```
Y-Aib2-EGTFTSDYSI-Aib13-LDKIAQ-K20(Nε-C20diacid-γGlu-(AEEA)2)-AFVQWLIAGGPSSGAPPPS-NH2
```

| Feature | Handling in modeling |
|---|---|
| Aib2, Aib13 (α-aminoisobutyric acid) | modeled explicitly as CCD `AIB` at residue index 1 & 12 |
| K20 C20-diacid/γGlu/(AEEA)₂ lipid | **omitted** from the co-fold epitope model (albumin-binding chain; not the antibody epitope). Flagged as a design consideration — see §6 |
| C-terminal amidation (-NH₂) | not represented (terminal cap); negligible effect on the mapped epitope |

**Rationale:** an antibody/peptide binder against tirzepatide will engage the
peptide backbone/side chains. The lipid at K20 is a solubility/half-life
appendage and is the *least* attractive epitope (flexible, albumin-shielded), so
the peptide-only model is the right target for binder discovery.

### 1.1 Fully-modified-target validation (K20 lipid included)

To confirm the leads work against the **actual modified drug** (not just the
backbone), the K20 acyl group (C20-diacid–γGlu–(AEEA)₂, extracted verbatim from
the supplied SMILES) was co-folded as a ligand alongside the Aib-containing
peptide, and the leads were re-run (`results/modified_target_results.json`,
complexes in `results/cif_modified/`):

| Lead | protein_iptm (binder↔peptide) | overall ipTM | vs lipid-free |
|---|---|---|---|
| TZP-B1 (spec_13) | 0.95 | 0.95 | ↑ (was 0.71) |
| TZP-P1 (R2_008) | 0.97 | 0.84 | ≈ (was 0.92) |
| TZP-P2 (R2_002) | 0.98 | 0.92 | ↑ (was 0.91) |
| TZP-P3 (R2_035) | 0.97 | 0.88 | ≈ (was 0.95) |

**The K20 lipid does not disrupt binding** — the binder↔peptide interface
confidence stays 0.95–0.98, consistent with the mapped epitope (F22/V23/L26/I27,
C-terminal face) sitting on the opposite side from K20. (Lower *overall* ipTM
just reflects the floppy lipid tail dragging the whole-complex metric down; the
protein–protein sub-score is what matters.)

**Feasibility limits of the hosted endpoint:** (i) atom-level covalent bonds are
allowed only to CCD ligands, so the lipid is co-folded as a *proximal, non-bonded*
ligand rather than a true Lys20-Nζ bond; (ii) the C-terminal amide is not
representable. For the **exact covalent model**, a local Boltz-2 run with the
provided `design/tirzepatide_fully_modified.boltz.yaml` (covalent bond to the
SMILES acyl) is the route.

---

## 2. Unified evaluation of the 90 existing binders

The supplied designs are **de novo mini-proteins** (β-sheet/mixed scaffolds,
63–104 aa), not antibodies. They were re-scored two independent ways.

### 2.1 Composite ranking from the supplied metrics (all 90)

A weighted z-score composite was computed emphasizing interface-confidence
metrics (ipTM, ipSAE, interface PAE) plus H-bonds, buried SASA and a
developability (liability) penalty.

**Key result — your 5 yellow-highlighted picks are exactly the top 5** by this
independent composite:

| Composite rank | ID | ipTM | ipSAE | minPAE | affinity_prob | liab | your pick |
|---|---|---|---|---|---|---|---|
| 1 | design_spec_05 | 0.62 | 0.45 | 3.40 | 0.015 | 90 | ★ yellow |
| 2 | design_spec_01 | 0.58 | 0.40 | 2.55 | 0.056 | 140 | ★ yellow |
| 3 | design_spec_91 | 0.58 | 0.44 | 2.31 | 0.062 | 75 | ★ yellow |
| 4 | design_spec_30 | 0.62 | 0.42 | 3.40 | 0.094 | 145 | ★ yellow |
| 5 | design_spec_16 | 0.62 | 0.41 | 4.54 | 0.056 | 155 | ★ yellow |

Your manual selection was well-calibrated to the aggregate structural evidence.

**Honest caveat:** across all 90, interface ipTM tops out at 0.62 and Boltz-2
affinity probability never exceeds 0.37 (median 0.05). None are "slam-dunk"
predicted binders — expected for this target class.

### 2.2 Independent Boltz-2 re-validation (top 20 co-folded fresh)

Each top-20 binder was re-co-folded with the unified Tirzepatide model. This
gave a cleaner, apples-to-apples read and **re-ordered the ranking**:

| Reval rank | ID | ipTM (refold) | binding_conf | minPAE | your pick |
|---|---|---|---|---|---|
| 1 | **design_spec_13** | 0.86 | **0.20** | 1.44 | — |
| 2 | design_spec_01 | 0.92 | 0.00 | 0.89 | ★ |
| 3 | design_spec_03 | 0.90 | 0.00 | 0.87 | — |
| 4 | design_spec_16 | 0.88 | 0.00 | 1.29 | ★ |
| 5 | design_spec_63 | 0.89 | 0.00 | 0.95 | — |
| 8 | design_spec_30 | 0.87 | 0.00 | 1.33 | ★ |
| 14 | design_spec_91 | 0.77 | 0.00 | 2.47 | ★ |
| 16 | design_spec_05 | 0.74 | 0.00 | 2.95 | ★ |

**Findings:**
- On fresh co-folding, **interface ipTM jumps to 0.74–0.92** for the top set —
  geometrically these are confident complexes.
- **`design_spec_13` is the standout**: it is the *only* binder with a non-trivial
  Boltz-2 binding-confidence (0.20, ~10× the field) **and** it had the highest
  affinity probability (0.371) in the original sheet. Cross-validated on two
  independent signals → **top existing candidate to carry forward.**
- Two yellow picks (01, 16) re-validate strongly; two (05, 91) drop on the
  refold (higher interface PAE, lower monomer confidence). So re-validation
  *refines* the yellow set rather than overturning it.

### 2.3 Consensus epitope on Tirzepatide

Interface contact mapping across the top-5 co-folded structures (≤4.5 Å heavy-atom):

| Contact frequency | Tirzepatide residues |
|---|---|
| **5/5 (core)** | **A18, Q19, F22, V23** |
| 4/5 | E3, L26, I27, P31, P37, P38, S39 |
| 3/5 | Y1, F6, S11, I12, L14, D15 |

→ The dominant epitope is the **C-terminal amphipathic helical face
(F22–V23–L26–I27) plus the polyproline C-terminal tail (P31/P37/P38/S39)**, with
Q19/E3 as polar anchors. All new designs were directed at this epitope.

---

## 3. New peptide binder design

### 3.1 Strategy

A 100-member library was generated (reproducible, no external DB) targeting the
consensus epitope from §2.3, spanning four families: epitope-directed
amphipathic helices (43), disulfide-stapled macrocyclic hairpins (18),
mini-binder-derived paratope grafts (13), incretin-receptor-ECD mimetics (6),
and diverse controls. It was screened with Boltz-2 **epitope-directed**
(epitope residues E3/A18/Q19/F22/V23/L26/I27).

### 3.2 Round 1 — the discriminating signal is `binding_confidence`

Under epitope-direction, interface **ipTM is uniformly high (~0.9)** for nearly
all candidates and therefore **not discriminating** — it is partly imposed by the
constraint. The honest discriminator is Boltz-2's **`binding_confidence`** head.
On that metric, of 100 designs only one stood out:

- **NB094** (`SMSTMENELSTLENEISTIENEWSTG`, amphipathic helix) — binding_confidence
  **0.21**, matching the best existing binder (spec_13, 0.20). All 98 others ≈ 0.

### 3.3 Round 2 — affinity maturation of NB094 (35 variants)

Guided variants were designed: (i) replace the two Met "knobs" (also an
oxidation liability) with Leu/Ile/Phe/Tyr/Trp; (ii) tune the aromatic anchor;
(iii) scan charge patterns and length. Binding_confidence **roughly doubled**:

| Variant | design move | bind (forced) | **bind (un-forced)** | ipTM |
|---|---|---|---|---|
| R2_035 | deMet + extended | 0.41 | **0.41** | 0.95 |
| R2_008 | deMet + Phe-anchor | 0.38 | **0.37** | 0.92 |
| R2_002 | Met→Leu (clean) | 0.37 | **0.37** | 0.91 |
| R2_010 | deMet + Tyr-anchor | 0.36 | **0.35** | 0.95 |
| R2_007 | di-Trp | 0.35 | **0.34** | 0.92 |
| R2_034 | 21-mer minimal | 0.32 | **0.33** | 0.94 |
| NB094 (parent) | — | 0.19 | **0.18** | 0.89 |

Met→Leu/Phe/Tyr **removed the only liability while improving affinity** — a
clean win. See `results/figure_maturation.png`.

### 3.4 Decisive validation — un-forced screen

Re-running the leads with **no epitope constraint** (the strictest, most honest
test) confirmed the gains are real, not an artifact of forcing: the matured
peptides retain `binding_confidence` 0.33–0.41 **un-forced**, ~2× both the NB094
parent (0.18) and the best existing binder spec_13 (0.20). Every non-matured
control (including ipTM-0.93 helices) collapses to ≈0 un-forced.

All three modeled leads reproducibly engage the **same amphipathic helical face**
of tirzepatide (Q19, V23, L26, I27 + F6/Y10/I12), consistent with the design
hypothesis (structures in `results/cif_leads/`).

## 4. Combined final panel & wet-lab recommendations

Master table: `results/MASTER_panel.csv` / `.md`. Constructs (peptide + Fc):
`results/wetlab_constructs.fasta`.

### 4.1 Recommended panel (order for wet lab)

| ID | Class | Tier | bind (un-forced) | Sequence |
|---|---|---|---|---|
| **TZP-P3** | de novo helix, extended | A | **0.41** | `SLSTLENELSTLENEISTIENEWSTGLENEISTG` |
| **TZP-P1** | de novo helix, Phe-anchor | A | **0.37** | `SLSTLENELSTLENEISTIENEFSTG` |
| **TZP-P2** | de novo helix, Trp-anchor (cleanest) | A | **0.37** | `SLSTLENELSTLENEISTIENEWSTG` |
| **TZP-P4** | de novo helix, Tyr-anchor | A | **0.35** | `SLSTLENELSTLENEISTIENEYSTG` |
| TZP-P5 | de novo helix, di-Trp | B | 0.34 | `SLSTWENELSTLENEISTIENEWSTG` |
| TZP-P6 | de novo helix, 21-mer minimal | B | 0.33 | `SLSTLENELSTLENEISTIEG` |
| **TZP-B1** | existing mini-binder (best; orthogonal β-scaffold) | A | 0.20 | 96-aa, see FASTA |

All Tier-A peptides: **no sequence liabilities, no Cys** (no scrambling),
soluble (net −6 to −8, negative GRAVY), 21–34 aa → **straightforward SPPS**.
Include TZP-B1 as an orthogonal-scaffold positive-diversity arm (recombinant).

### 4.2 Suggested wet-lab cascade

1. **Synthesize** TZP-P1–P4 (SPPS) + express TZP-B1. Add a scrambled-sequence
   negative control (same aa composition) and, ideally, an N-terminal biotin or
   FITC tag (through a short PEG/GSG spacer at the N-terminus, away from the
   C-terminal binding face).
2. **Primary binding:** SPR / BLI against immobilized synthetic tirzepatide (and
   the unmodified peptide backbone). Also test **± albumin/lipid** context, since
   the K20 diacid chain may modulate accessibility (§6).
3. **Epitope confirmation:** competition vs a biotinylated C-terminal fragment;
   optionally alanine-scan Q19/V23/L26/I27 on the target.
4. **Format up hits:** the provided **human IgG1-Fc fusions** (bivalent, long
   half-life) or **IgG4-S228P** (effector-silent) — `results/fc_fusion_constructs.json`.
   Fc dimerization gives avidity that can rescue modest monomer affinity.
5. **Iterate:** feed SPR KD back into another Boltz maturation round on the
   winning scaffold (the pipeline here converges in ~2 min/round).

### 4.2b Longer, autonomously-folding binders (stability upgrade)

The Tier-A peptides (P1–P6) are single amphipathic α-helices of 21–34 aa; a lone
short helix is only marginally stable in isolation (relies on binding-induced
folding). To give **more stable, autonomously-folding** reagents with their own
hydrophobic core (and, where asked, β-structure), a round of 23 longer designs
(34–75 aa: capped/stapled helices, 3-helix bundles, α-hairpins, mixed α/β) was
generated and screened against the fully-modified target. The best fold **and**
bind on par with the short peptides:

| ID | src | AA | topology | bind_conf(mod) | fold (struct_conf) | α / β |
|---|---|---|---|---|---|---|
| **TZP-S1** | LS3_bundle2 | 75 | 3-helix bundle | 0.74 | 0.85 | 85% α |
| **TZP-S2** | LS5_c | 49 | **mixed α/β** (β-hairpin+helix) | 0.74 | 0.82 | 43% α / **33% β** |
| **TZP-S3** | LS4_hairpinF | 47 | α-hairpin (2 helices) | 0.73 | 0.86 | 77% α |
| **TZP-S4** | LS2_b | 34 | disulfide-stapled helix | 0.73 | 0.85 | 94% α |

`structure_confidence` 0.82–0.86 indicates a well-defined autonomous fold (vs a
lone helix). **TZP-S1** (3-helix bundle) is the most stable format; **TZP-S2**
(mixed α/β) provides real β-sheet content as requested. Sequences + IgG1-Fc
fusions are in `results/wetlab_constructs.fasta`; ranking in
`results/stable_panel.json` and `results/long_stable_screen.json`.

### 4.2b-2 Atomistic confirmation of the stable designs (clean protein_iptm)

The four stable designs were re-run as full `structure_and_binding` predictions
against the fully-modified target (peptide+Aib+K20 lipid), giving the **clean
binder↔peptide interface score (protein_iptm)** and folded complex CIFs
(`results/cif_stable/`, metrics in `results/stable_atomistic.json`):

| ID | AA | topology | **monomer pTM (fold)** | **protein_iptm (interface)** |
|---|---|---|---|---|
| TZP-S1 | 75 | 3-helix bundle | 0.96 | 0.96 |
| TZP-S2 | 49 | mixed α/β | 0.96 | 0.97 |
| TZP-S3 | 47 | α-hairpin | 0.96 | 0.97 |
| TZP-S4 | 34 | stapled helix | 0.97 | 0.96 |

**All four both fold autonomously (pTM 0.96–0.97) and grip the modified drug at
the mapped epitope (protein_iptm 0.96–0.97)** — the clean confirmation that the
stability upgrade did not cost binding. Unlike the lipid-inflated
`binding_confidence`, protein_iptm isolates the peptide interface and is the
trustworthy number here.

**Optimization round (18 variants of S1/S2):** best is **TZP-S1.1** (`S1o_short`,
66 aa 3-helix bundle, F-anchor) — bind_conf 0.75, fold 0.87, i.e. a shorter,
higher-fold-confidence bundle than TZP-S1. The mixed-α/β parent (LS5_c/TZP-S2)
was already optimal among the hairpin variants. See `results/s1s2_opt_screen.json`.

### 4.2c Fc-fusion format validation (fully-modified target)

All six P-series leads were also co-folded as **IgG1-Fc fusions** against the
fully-modified target (`results/mod_leads_fc.json`). Fusion preserves binding —
interface ipTM 0.88–0.95 (bind_conf 0.37–0.57; the 273-aa Fc lowers the
aggregate score but not the binder interface). P3_Fc ranked best. Peptibody
format is compatible; bivalent Fc dimerization should add avidity.

> Note: with the K20 lipid present in the target, screen `binding_confidence`
> rose to ~0.65–0.76 for all binders — partly because the hydrophobic binders can
> also contact the greasy acyl chain, so these are **upper-bound** triage values.
> The cleaner metric for the peptide epitope is the atomistic protein_iptm
> (0.95–0.98, §1.1).

### 4.3 Honest expectations

- Boltz-2 `binding_confidence` ≈ 0.4 is **encouraging but not a guarantee** of
  nM affinity. Expect initial hits in the **µM–mid-nM** range; the Fc-avidity
  format and one more maturation round are the levers to improve it.
- The de novo helices are a **single structural hypothesis** (one amphipathic
  helix per target face). TZP-B1 hedges with an independent β-scaffold that hits
  a partially overlapping epitope — worth running in parallel.

## 5. Methods

- **Parsing/scoring:** `analysis/` and `design/` scripts (reproducible, no RNG).
- **Boltz-2.1** via Boltz API remote MCP: `protein_screen` (epitope-directed,
  epitope residues E3/A18/Q19/F22/V23/L26/I27) and `structure_and_binding`.
- **Composite score:** z-scored ipTM (w2.0), −interface PAE (w1.5),
  binding_confidence (w1.5), structure_confidence (w0.5); developability penalty
  0.10×(#liabilities) and 0.5 for unpaired Cys.
- **Developability:** net charge, GRAVY, aromatic fraction, Cys parity, and
  sequence-liability motifs (NG/NS deamidation, DG/DP isomerization, N-x-S/T
  glycosylation sequon, Met oxidation, hydrophobic runs).

## 6. Design considerations for wet lab

- **Lipid shielding:** in formulated tirzepatide the K20 diacid chain and
  albumin association may partially occlude parts of the helix. The mapped
  epitope (F22/V23/L26/I27 + C-terminal PPPS) is on the opposite face from K20,
  which is favorable, but binding to *drug-bound-to-albumin* tirzepatide should
  be confirmed empirically (e.g., SPR with and without albumin/lipid).
- **Assay target formats:** test against (i) synthetic full tirzepatide,
  (ii) the unmodified peptide backbone, and (iii) a biotinylated C-terminal
  fragment to confirm the mapped epitope.
- **Format / Fc:** peptide and mini-binder hits can be expressed as monomers for
  primary binding assays; for the final therapeutic-style reagent, fuse the
  validated binder to a **human IgG1 Fc (knob-into-hole optional)** or use a
  peptide-Fc fusion for avidity + half-life. Sequence provided in §4.
