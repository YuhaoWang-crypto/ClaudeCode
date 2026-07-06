# Antibody binder design — methodology playbook

Distilled from a multi-round anti-Tirzepatide antibody campaign (target = a lipidated 39-mer
peptide drug). Every step below was used in practice; the failure modes are ones that actually bit.

---

## 1. Modification-aware targeting (the central lesson)

Many real targets are **modified**: a lipidated/PEGylated peptide drug, a glycoprotein, a
phospho-epitope. The modification occupies space.

**The trap.** If you optimize a binder against the *bare* target (e.g. lipidated Lys modeled as
plain Lys), the paratope happily collapses onto the now-empty attachment point. Co-fold scores look
great, but the real modified drug sterically clashes and won't bind. In the campaign, 14/15 prior
binders buried the K20 Nζ (median 3.2 Å) — exactly where the 762 Da lipid sits — and all failed.
An independent re-check on 5 unrelated Fab co-folds found **5/5** buried the Nζ at 3.2–4.0 Å: bare
modeling cannot discriminate modification-tolerant binders at all.

**What to do.**
- Represent non-standard residues by CCD code (α-aminoisobutyric acid → `AIB`, phosphoserine →
  `SEP`, etc.) as polymer `modifications` in the co-fold input.
- For a bulky PTM (lipid, glycan, PEG): attach it as a co-folded **ligand** (SMILES/CCD) bonded to
  the anchor atom, so it competes for space — this is the faithful model.
- If you must run a controlled bare-residue comparison, **always** also run the occlusion check
  (`contacts.py --occlusion CHAIN:RES:ATOM`) and treat any `buried: true` as disqualifying for the
  real drug.
- A "modification-tolerant recognition score" (interface × modification-avoidance × epitope
  robustness) can invert a naive co-fold ranking — the naive top binders are often the worst on the
  real target.

---

## 2. Two epitope strategies (keep them separate)

- **Epitope near the modification** — optimize an existing binder so the paratope *accommodates or
  side-steps* the modification. Use a **modification-aware** sequence designer (LigandMPNN with the
  lipid present) so new CDRs are designed knowing the bulk is there.
- **Epitope distal to the modification** — design *de novo* against a clean epitope (e.g. a
  C-terminal hotspot 12–14 Å from the lipid). Use **RFantibody/RFdiffusion** to generate CDRs on a
  humanized framework, then **ProteinMPNN** (NOT LigandMPNN — there is no lipid to condition on;
  conditioning on a distant lipid drags the design toward it, a real bug that happened) anchored on
  the docked pose.

Scoring formulas differ between the two (one rewards modification-coexistence, the other
modification-avoidance) → **within-epitope ranking only; never compare absolute scores across
epitopes.**

---

## 3. Cross-model scoring (Boltz vs Chai)

Run two independent structure models. Empirically their **interface ipTM absolute values do not
correlate** (Pearson r ≈ 0.26): Boltz optimistic, Chai conservative. A candidate with high Boltz but
low Chai interface (e.g. an un-matured parent at 0.909 / 0.253) is *unstable recognition*, not a
lead.

Rank on signals that **agree across both models**:
- modification burial (lower = better),
- hotspot residue coverage (more = better),
- epitope identity / Jaccard stability,
- anchor-atom distance (modification accommodated, not occluded).

Distrust any lead that only one model likes. De novo designs that Boltz scores 0.78–0.91 while Chai
scores 0.11–0.19 are **not confirmed** — they are seeds for larger-scale redesign, not leads.

---

## 4. Full humanized IgG assembly + the Fab Fc-interference test

Assemble VH/VL into a complete human IgG1 (heavy = VH-CH1-hinge-CH2-CH3; light = VL-Cλ or -Cκ) with
`scripts/assemble_igg.py`. Then answer "does the constant region / Fc change binding?" correctly:

- **Co-fold the Fab** (Fd = VH-CH1 + light = VL-CL) with the antigen — not the whole IgG. The Fc is
  >100 Å from the paratope; co-folding the entire IgG dilutes and artefacts the interface (a lone Fc
  misfolds into a ~0-confidence homodimer).
- A small interface-ipTM drop from the isolated Fv to the Fab (~0.05–0.07 seen) is **metric
  dilution** — ~200 extra non-contacting CH1+CL residues thin the global score — not weaker binding.
- Confirm the epitope (Jaccard ≈ 0.8–1.0), paratope, and modification-anchor distance are unchanged.
- Watch for format-dependent regressions: some variants that avoid the lipid as an Fv **re-bury** it
  in the Fab context (seen: 0→38, 1→28 buried atoms). The lead should stay clean across formats.

Only the Fc N297 glycan is introduced — a canonical, distal, non-paratope site.

---

## 5. Developability screening (do it on the *matured* sequence)

Run `scripts/liabilities.py` on every lead. Motif classes:
- **N-glycosylation** sequon `N-X-[S/T]`, X≠P — can be introduced by maturation in a framework loop
  (seen: a matured VL gained `N67-G-T`; fix N67Q or T69A).
- **Deamidation** `NG` (fast), `NS/NN/NH`.
- **Isomerization / succinimide** `DG` (fast), `DS/DD/DT`. Asp-Gly in a CDR is the worst case.
- **Unpaired Cys** — an odd cysteine count beyond the conserved intradomain pair.

Always pass `--parent-vh/--parent-vl` to surface **maturation-introduced** liabilities. Real example:
a maturation mutation `L-M47D` created a brand-new `D-G` in CDR-L2 (`…VLV-M-GE…` → `…VLVF-D-G-E…`)
that a first-pass QC missed because it only reported the pre-existing framework `D-G`. In-CDR
liabilities matter more than framework ones — fix (D→E, or break the Gly) and re-score before
committing to expression.

---

## 6. Ranking rubric

1. Within-epitope only; cross-epitope absolute scores are incomparable.
2. Structural cross-model agreement > single-model interface score.
3. Modification avoidance / accommodation is decisive for modified targets.
4. End-to-end robustness: lead must survive Fv → Fab → full-IgG without re-burying the modification.
5. Developability breaks ties (a slightly lower but liability-free lead beats a marginally higher
   one carrying a new CDR glycosite/isomerization site).

---

## 7. Honest limitations (state them in every deliverable)

- All co-fold numbers are predictions, not measured Kd/kon/koff.
- Boltz/Chai disagree on absolute interface confidence; only structural agreement is trustworthy.
- Custom recognition scores are heuristic weightings — calibrate against wet-lab data.
- De novo short CDR-H3 (e.g. 5 aa) gives a small paratope → likely moderate affinity; flag it.
- **Mandatory next step:** express (HEK293 transient), SPR/BLI vs the **real modified** target,
  DSF/nanoDSF + SEC-HPLC, and an epitope-binning experiment (e.g. truncated-peptide competition).
  Include the prior "failed" binders as negative controls to prove the new design recognizes the
  real drug.
