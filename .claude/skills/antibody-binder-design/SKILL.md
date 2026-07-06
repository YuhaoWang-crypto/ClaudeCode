---
name: antibody-binder-design
description: >
  Computationally design, mature, and evaluate antibody/nanobody binders against a
  peptide, protein, or modified-peptide (lipidated/PTM) target using structure co-folding.
  Use when the task involves antibody or binder design, affinity maturation, epitope/paratope
  analysis, humanized IgG assembly from VH/VL, Fc/Fab formatting, developability liability
  screening, or Boltz/Chai co-fold scoring and ranking of antibody candidates. Bundles tested
  scripts for IgG assembly, epitope + PTM-occlusion contact analysis, and liability motif
  scanning, plus a methodology playbook covering modified-target modeling and cross-model scoring.
---

# Antibody binder design & evaluation

A playbook + toolkit for going from a target (peptide / protein / **modified** peptide such as a
lipidated or PEGylated drug) to ranked, wet-lab-ready humanized IgG leads, validated by
structure co-folding. Distilled from a multi-round anti-Tirzepatide campaign.

## When to use
Antibody/nanobody/binder design; affinity maturation; "does the Fc change binding"; epitope or
paratope mapping from a predicted complex; assembling full IgG from VH/VL; ranking candidates by
co-fold; screening antibody sequences for developability liabilities; modeling a drug/peptide
target that carries a lipid or other post-translational modification.

## The 6 core moves (in order)

1. **Model the *real* target, not a convenient one.** If the target carries a modification
   (lipid, glycan, PEG, acyl), that group occupies space and can occlude the paratope. Modeling
   it as the bare residue gives falsely good, misranked binders. See
   `references/methodology.md` § "Modification-aware targeting". Represent non-standard residues
   with CCD codes (e.g. Aib → `AIB`); attach bulky PTMs as a co-folded ligand + bond, or at
   minimum **check the modification anchor is not buried** (move 4).

2. **Pick the epitope deliberately.** Optimize an existing binder toward a modification-tolerant
   pose, OR design *de novo* against a modification-distal epitope. Track which you're doing —
   the scoring differs and cross-epitope scores are not comparable.

3. **Co-fold and score with two models.** Run Boltz-2 (see `references/boltz_cofold.md` for the
   exact MCP call sequence) and, when available, Chai-1. Their **absolute interface ipTM does not
   correlate** (expect Pearson r≈0.2–0.3): Boltz is optimistic, Chai conservative. Rank on
   **structural signals that agree across models** — modification burial, hotspot coverage,
   epitope identity, anchor distance — not on a single ipTM.

4. **Analyze the interface + occlusion:**
   `python scripts/contacts.py complex.cif --ab H,L --antigen P --occlusion P:20:NZ`
   Emits the epitope residue list and whether the modification anchor atom (e.g. lipidated Lys
   Nζ) is buried (`buried: true` ⇒ the real modified drug would clash).

5. **Assemble the full humanized IgG and Fab-test the format:**
   `python scripts/assemble_igg.py --vh ... --vl ... --light lambda|kappa --fasta out.fasta`
   Then co-fold the **Fab** (Fd = VH-CH1, light = VL-CL) with the antigen. The Fc is >100 Å from
   the paratope and, co-folded whole, misfolds into an artifact — so score the Fab. A small ipTM
   drop vs the bare Fv is "metric dilution" (more non-contacting residues), not weaker binding;
   confirm the epitope/paratope and anchor distance are preserved.

6. **Screen developability before committing:**
   `python scripts/liabilities.py --vh ... --vl ... --parent-vh ... --parent-vl ...`
   Flags N-glyc sequons, deamidation (NG/NS), isomerization (DG/DS), unpaired Cys — and, with
   `--parent-*`, motifs **newly introduced by maturation** (a frequent trap: a matured CDR gains a
   fresh Asp-Gly). Fix in-CDR liabilities (e.g. D→E, or break the Gly) and re-score.

## Bundled scripts (all tested, stdlib-only Python 3)
- `scripts/assemble_igg.py` — VH/VL → full human IgG1 heavy + λ/κ light + Fab chains.
- `scripts/contacts.py` — epitope/paratope contacts + PTM-occlusion check from a co-fold CIF.
- `scripts/liabilities.py` — variable-domain liability motif scan, with parent-diff for maturation.

## References
- `references/methodology.md` — the full playbook: modification-aware targeting, optimize-vs-de-novo,
  maturation tool choice (lipid-aware LigandMPNN vs plain ProteinMPNN), cross-model scoring, the
  Fab Fc-interference test, ranking rubric, and the honest limitations.
- `references/boltz_cofold.md` — exact Boltz-2.1 MCP workflow (estimate → start → poll → download
  CIF → parse metrics), chain/entity conventions, and cost notes.

## Hard rule
Every co-fold number is a **prediction**, not measured affinity. The deliverable of this skill is a
*ranked, de-risked shortlist to test* — final calls require SPR/BLI against the **real (modified)**
target plus DSF/SEC.
