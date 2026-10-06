# Stage 0 — LYTAC / EndoTag receptor epitope specification

Reference: Designed endocytosis-inducing proteins degrade targets and amplify signals
([PMC11839401](https://pmc.ncbi.nlm.nih.gov/articles/PMC11839401/)).

Stage 0 is the receptor-side work. It is **target-independent**: the output depends
only on which lysosome-targeting receptor (LTR) you recruit, so one run serves any
protein-of-interest you later decide to degrade. No paid compute was used.

Reproduce with:

```bash
./fetch_structures.sh structures_raw
python3 ltr_epitope_spec.py --struct-dir structures_raw --out-dir specs
```

---

## What came out of it

Four receptors specified, 491 index assertions passing, and **three numbering traps
plus two pipeline bugs** found — each of which would have produced a confident design
run against the wrong surface.

| Receptor | Structure | Residues | Exposed | Avoided | Verdict |
|---|---|---|---|---|---|
| **Sortilin** (SORT1) | 3F6K A | 657 | 265 | 122 | **Best first choice** |
| ASGPR (ASGR1) | 5JQ1 A | 130 | 63 | 59 | Tight but workable |
| IGF2R D6+D11 | 6UM2 A (**bovine**) | 1638 | 837 | 82 | Species + competition risk |
| IGF2R D11 human | 1GP0 A | 133 | 70 | 0 | Correct-species D11 template |

---

## 1. Three numbering traps

All three were resolved from the authoritative `_struct_ref_seq` mapping in the mmCIF
header, not by assuming the published frame.

| Structure | UniProt | auth − UniProt | Consequence if ignored |
|---|---|---|---|
| 3F6K / 4PO7 / 6X3L sortilin | Q99523 | **−33** | 33-residue propeptide. Paper quotes Site 1 in mature numbering; PDB uses precursor. |
| 6UM2 IGF2R | **P08169 (Bos taurus)** | 0 | **Wrong species**, see below. |
| 5JQ1 ASGPR | P07306 | −1 | Off-by-one across the whole CRD. |
| 1GP0 IGF2R D11 | P11717 (human) | 0 | — |

**Sortilin Site 1 confirmed.** Published as F92/V93/T546/T559/T561. Scanning offsets,
only −33 reproduces the residues (4/5 name match → auth 59, 60, 513, 526, 528). The
next-best offset matches 2/5. The one residue that does not reconcile (published T546
→ auth 513 is VAL) should be re-derived before it is used as a design hotspot.

**The species trap.** The paper designed the D6 arm on 6UM2 and the D11 arm on 1GP0 —
these are **different organisms**. Pairwise alignment of the domains:

```
D6  bovine 773-932  vs human 765-924   83.1% identity, no gaps (27/160 residues differ)
D11 bovine 1523-1657 vs human 1514-1648 89.6% identity, no gaps (14/135 differ)
```

27 differences across a 160-residue domain is a lot to put under an engineered
interface. A D6 binder optimised on bovine coordinates has no guarantee of
transferring to human IGF2R. Design the human D6 arm on a human model.

---

## 2. The cross-domain clamp — and why the paper's 3-residue linker works

EndoTag's hardest-won IGF2R result was empirical: `D11mb–GGS–D6mb` drives uptake,
**longer linkers reduce it, a shorter Gly-Ser linker abolishes it**, and two copies of
one minibinder are not taken up at all.

6UM2 is the only structure holding D6 and D11 in one frame, so the geometry is
measurable. IGF2 bridges both domains — it *is* the natural bivalent ligand:

```
IGF2 footprint on D11 : 16 residues
IGF2 footprint on D6  : 11 residues
footprint centroid gap: 18.4 A  (raw 4.5 A contacts) / 22.1 A (8 A-dilated)
closest approach      : 2.6 A   -- the domains are in direct contact
=> ~3-7 Gly-Ser residues, and ~0 at the contact region
IGF2's own internal span: 11.4 A centroid; residues 11/14/18 touch BOTH domains
```

That reproduces the paper's linker window quantitatively. A GGS linker is right
because the two footprints are ~18 Å apart on domains that touch; longer linkers pay
entropy for span they do not need; a shorter one cannot reach.

**A naive anchor gets this wrong by 4×.** Scoring the best-exposed patch on each
domain independently puts them on opposite outer faces and reports a 72.4 Å span
(≥21 GS residues) — a construct that would almost certainly fail. The pipeline now
reports the native-ligand footprint span as the primary anchor and labels the
patch-based number as a caveat.

### The IGF2R design tension, quantified

The junction a bivalent binder must clamp largely *is* the IGF2 site:

```
D6  junction 34 residues, 15 also IGF2 footprint (44%)
D11 junction 17 residues, 14 also IGF2 footprint (82%)   <- limiting side
```

So an IGF2-orthogonal bivalent clamp has very little room on D11. Whether it has
*any* is **threshold-dependent** — across 24 combinations of footprint dilation
(4.5/6/8 Å), junction width (8/12/15/20 Å) and exposure cutoff (0.15/0.25), the
limiting side holds **0–22 residues, with 62% of settings leaving ≥4**.

That is not a clean verdict and is not reported as one. The honest reading: **an
orthogonal bivalent IGF2R EndoTag is geometrically tight-to-marginal.** If you pick
IGF2R with the proven bivalent architecture, plan for competition with circulating
IGF2 — which is abundant in serum, and which Science 2023 already flags as the same
class of problem as endogenous M6P-glycoprotein occupancy of surface CI-M6PR.

---

## 3. Why sortilin is the recommended starting receptor

- **No clamp requirement.** Single-domain engagement worked (21 nM, the best affinity
  in the paper). None of the bivalent geometry risk above applies.
- **Ample orthogonal surface.** 265 exposed residues against 122 avoided — and the
  avoid set is deliberately conservative: it merges **three native-ligand sites from
  three crystal forms**, the two neurotensin sites (3F6K chain N; 4PO7 chains N+P) and
  the progranulin-site compound UMJ (6X3L). Even so, 254 exposed residues stay free at
  an 8 Å dilation.
- **Robust to the cutoff.** Avoid radius 4→10 Å moves the free exposed count only
  262→246.

ASGPR is workable but cramped: a 130-residue CRD with three Ca²⁺ ions and the
galactosamine mimic leaves 44 free exposed residues at 8 Å. Its multivalency span is
measurable here — the two CRD copies in 5JQ1 sit **31.9 Å centroid-to-centroid**
(closest approach 2.8 Å), which is the span a 2C/3C construct has to cover. Note
ASGPR is a trimer in vivo, and the paper's gain was almost entirely valency:
monomer ≈ no uptake → 2C 2.5× → 3C 4.5×, with a binder of only 2.7 µM.

---

## 4. Two bugs the built-in validation caught

Both produced plausible-looking output and would not have raised anything at runtime.

**(a) Index desynchronisation across merged structures.** `union_auth` and
`union_0based` were built as two independently sorted lists. Merging sortilin's
native-ligand sites from 4PO7 and 6X3L introduced auth 106, which **3F6K does not
resolve** (it has a chain break at 101→107). One unresolved id made the lists differ
in length by one, and zipping them shifted every subsequent constraint by one residue
— relocating the entire non-binding set. Fixed by intersecting transferred ids with
the primary structure and deriving both lists from one filtered sequence. The dropped
ids are now recorded rather than silently absorbed.

**(b) Construct termini topping the epitope ranking.** Chain ends and chain-break
edges are maximally solvent-exposed, so they swept the top of the patch score.
Sortilin's winning patch was auth 708–716 — the **last nine residues of the
construct**, where the crystallographer cut, not a surface a binder can use (the real
receptor continues into a transmembrane region). ASGPR's and D11's top patches were
likewise part-terminal. Now 5 residues are trimmed from each end and 2 either side of
every chain break, and real internal surfaces rank first.

Validation is now part of the pipeline: it re-reads each written coordinate file and
asserts that every 0-based index resolves to its stated auth residue. An off-by-one
here does not crash a design run — it quietly scores designs against the wrong site.

---

## 5. Deliverables

```
lytac_design/
  ltr_epitope_spec.py      pipeline (numbering, SASA, footprints, patches, clamp, validation)
  fetch_structures.sh      pulls the six PDB entries
  specs/<receptor>.json    per-receptor spec
  specs/all_receptors.json combined
  structures/*.cif         cleaned single-chain files -- UPLOAD THESE
```

**The contract:** every 0-based index is relative to the matching file in
`structures/`. It is neither an auth number nor a UniProt number. Upload that exact
file or the indices are meaningless.

Top candidate epitopes, ready to paste into a design run:

| Receptor | Scope | n | mean rSASA | hydrophobic | `epitope_residues` (0-based) |
|---|---|---|---|---|---|
| sortilin | — | 10 | 0.53 | 0.50 | `[180,181,182,185,186,190,208,209,226,227]` |
| sortilin | — | 9 | 0.53 | 0.56 | `[357,358,375,377,378,380,407,408,409]` |
| asgpr | — | 10 | 0.45 | 0.30 | `[20,21,24,25,28,31,32,118,120,122]` |
| igf2r_d11_human | — | 11 | 0.52 | 0.55 | `[26,27,28,29,30,32,54,110,112,113,114]` |
| igf2r_bovine | D6 | 8 | 0.49 | 0.75 | `[295,297,299,300,301,302,385,451]` |
| igf2r_bovine | D11 | 10 | 0.49 | 0.40 | `[1035,1037,1038,1039,1040,1041,1042,1098,1100,1101]` |

`non_binding_residues` for each receptor is `native_ligand_footprint.union_0based`
in its spec file (sortilin 122, ASGPR 59, IGF2R 82, D11-human 0).

---

## Honest limits

- **Nomination is not validation.** These patches are hypotheses. Only a design pass
  plus co-folding — and ultimately experiment — says whether an interface is real.
  The paper's own hit rates were 0.15% (ASGPR, 2,689 → 4) through ~5% (IGF2R, 170
  tested → 8), and its ASGPR binder was 2.7 µM. Success there came from valency
  geometry, not affinity.
- **Patch scoring is a transparent heuristic**, not a trained model: exposure × log(size)
  × hydrophobic-content, greedily de-overlapped. It ranks; it does not predict.
- **SASA is computed on the isolated chain.** Correct for what a binder sees, but it
  ignores glycans. These are heavily glycosylated receptors and N-glycans will shield
  some nominated surface. Worth re-checking the shortlist against glycosylation sites.
- **6UM2 is 4.32 Å cryo-EM.** Side-chain positions there are not reliable at the
  precision an interface design wants; treat D6/D11 patch residue identities as
  approximate and prefer a higher-resolution or predicted human model for the final run.
- **The orthogonal-clamp question is unresolved**, deliberately. It needs a decision
  about IGF2 competition, not a better threshold.

## Suggested next step

Sortilin, single-domain, patches #0 and #1, against `structures/sortilin_3F6K_A.cif`
with the 122-residue non-binding set. That is the $50 / 500-design pilot — enough to
see whether the interface is designable at all before committing to scale.
