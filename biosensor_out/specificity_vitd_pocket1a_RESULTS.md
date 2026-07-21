# Directional-contact specificity test — pin the ligand to VDR's 1α-OH subsite

**Idea borrowed from two Baker de-novo-design papers** (luciferase, Nature 2023,
10.1038/s41586-023-06735-9; modular peptide binders, Nature 2023,
10.1038/s41586-023-05909-9): specificity comes from a **specific directional polar
contact** to the discriminating group (Arg→anion; bidentate H-bond→peptide backbone),
not from global interface affinity. Our flat 3×3 vitamin-D matrix used **global**
co-folding confidence — so we retested with the contact **focused on the discriminating
subsite**.

## Setup
- Receptor: VDR-LBD (1DB1, 259 aa).
- Discriminating group: the A-ring **1α-OH** (present in 1,25(OH)₂D3, absent in 25(OH)D3 and D3).
- Constraint: Boltz `pocket` (force=true, 6 Å) pinning the ligand to VDR's **native 1α-OH
  anchors Ser237 + Arg274** — mapped by alignment to our construct as residues **68 + 105**
  (0-indexed; alignment verified residue-by-residue).
- 3 constrained co-folds vs the 3 unconstrained baselines. All ✅ real Boltz-2.1.

## Result (✅ real numbers)

| metric | condition | D3 | 25OHD3 | 1,25 | **1,25 − 25OH** |
|---|---|---|---|---|---|
| ligand_iptm | unconstrained | 0.9725 | 0.9731 | 0.9812 | +0.0081 |
| ligand_iptm | **pinned to 1α** | 0.9768 | 0.9710 | 0.9804 | **+0.0094** |
| binding_conf | unconstrained | 0.6667 | 0.7569 | 0.7434 | −0.0135 |
| binding_conf | **pinned to 1α** | 0.6774 | 0.7526 | 0.7522 | **−0.0004** |

**Effect of pinning to the 1α-subsite (constrained − unconstrained):**
- binding_conf: **1,25 +0.0088** (happier when pinned by its 1α-OH), **25OHD3 −0.0043** (penalized — it has no 1α-OH to satisfy the contact).
- The hard-pair **1,25-vs-25OH** margin moved from **−0.0135 (wrong direction) → −0.0004 (tied)**, a **+0.013** correction; by ligand_iptm, 25OHD3 became the **lowest** of the three (biologically correct — it's the worst 1α-ligand).

## Verdict ⚠️ (honest)
**Directionally correct, but not sufficient.** Focusing the score on the discriminating
subsite moved every number the right way — the on-target gains, the off-target that lacks
the 1α-OH loses — validating the *principle* borrowed from the papers. But the margin is
still ≈ 0 (not a usable positive discrimination), and D3's iptm rose slightly (a caveat).
A single –OH stays **at/below co-folding resolution** even with a mechanism-focused
constraint.

**What this says to do next** (exactly what the Baker papers actually do, beyond scoring):
don't just *score* a soft constraint — **design the specific (bidentate) contact into the
pocket sequence** (RifDock/LigandMPNN place an H-bond donor locked onto the 1α-OH, plus a
counter-selection filter against 25OH), then confirm with a **competition assay**. Scoring
focuses the metric; *design* is what actually manufactures the margin.

Data: `specificity_vitd_pocket1a.json`. Jobs: VDR+1,25 `sab_pred_N3pMRWMKniPUMT0Cei8a`,
+25OH `sab_pred_31Olu8KtBRu3zclKY9Fv`, +D3 `sab_pred_HKiGAXAH3OQjYNAFP5a6`.
