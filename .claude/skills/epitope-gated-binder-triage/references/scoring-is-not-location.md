# Co-folding confidence does not locate your epitope

The empirical result behind this skill, how it was established, and how to
re-run the test on a new campaign. Re-running it is worth the hour: the result
is target- and modality-specific in principle, even though the mechanism
suggests it generalises.

## What was measured

1,600 de novo VHH designs against one target (mouse TGFβR2 ectodomain,
epitope-constrained to the TGF-β-blocking face). Every design's epitope
occupancy computed from its predicted complex. Four confidence scores tested
against epitope recall.

Two independent generation modes: 1,000 undirected designs (curated nanobody
scaffolds) and 600 parent-seeded CDR3 redesigns. The second is a constrained
population drawn from good parents, so agreement across both is meaningful.

## Result

Full-range Spearman ρ vs epitope recall looks encouraging for every score:
+0.32 (ipSAE) to +0.48 (binding_confidence). Restricting to designs that are
actually docked collapses it:

| stratum | n | ipSAE | binding_confidence |
|---|---|---|---|
| all | 1600 | +0.318 | +0.430 |
| top ~33% (≈ bc > 1e-2) | 529 | +0.092 | −0.029 |
| top ~20% (≈ bc > 0.1) | 315 | +0.101 (p=0.073) | −0.043 |
| round 1 alone, top 33% | 293 | +0.078 (p=0.18) | +0.067 (p=0.27) |
| round 2 alone, top 33% | 236 | +0.071 (p=0.28) | −0.152 (p=0.015) |

`ipTM` and `min_interaction_pae` behave the same way. As a blocker classifier,
AUCs cluster around 0.63–0.69 — above chance for "is this a plausible complex",
but that is a different question from "is it on my epitope", and we compute the
latter directly.

The full-range correlation is generated entirely in the low-score stratum, where
mean recall is ~0.51 versus ~0.81 at the top. Those designs are not mis-located;
they are barely docked.

## Why — the mechanism, which is the part that generalises

Confidence scores track **how much** interface exists, not **where**:
ρ(score, contacting atom pairs) stays at +0.28 to +0.44 in exactly the strata
where ρ(score, recall) is zero.

ipSAE deserves its own note because it is the natural "try this instead"
candidate. It was designed to fix ipTM's confounding by chain length, and an
external meta-analysis of thousands of experimentally characterised binders
supports it as a predictor of *binding*. But in a VHH library every binder is
112–132 aa (sd ≈ 5). There is no length artifact left to correct, so ipSAE
degenerates toward ipTM — measured ρ(ipSAE, ipTM) = +0.75, and plain ipTM beat
ipSAE in every stratum tested.

This predicts where ipSAE *would* help: campaigns mixing binder sizes (peptides
vs minibinders vs scFvs). In a length-homogeneous library, expect it not to.

One further observation: inside the epitope-gated pool, ipSAE correlated
*negatively* with interface purity (ρ ≈ −0.38). It prefers larger, stickier
interfaces. "Bigger is better" is not a safe default when the objective is a
focused blocking footprint.

## The prefilter control — run this before trusting any score-based filter

The tempting shortcut is to score only high-confidence designs' structures.
Test it rather than assuming:

- **Set A**: every design above the filter threshold → 47.3% blockers
- **Set B**: a random, seeded sample from *below* it → 21.7% blockers
- Fisher OR 3.25, p = 4.7e-04 — the score is a real enrichment
- But the discarded stratum was 13.5× larger, so it held **69% of all blockers**

A full scan later confirmed 284/1000 blockers overall against a 264 estimate
from the two strata (the extrapolation was unbiased; the CI contained the truth).
No score floor was defensible — even the bottom decile held blockers.

**Conclusion: score everything.** A 2.2× enrichment is not worth discarding two
thirds of the hits, and the analysis is free once the designs are paid for.

## How to re-run this on a new campaign

1. Compute epitope occupancy for **every** design (not a score-selected subset —
   that is the thing being tested).
2. For each score, compute Spearman vs recall: full-range, then restricted to
   progressively higher strata, then within each generation round separately.
3. Report the stratified number as the result. If you report a positive
   correlation, demonstrate it survives restriction.
4. Head-to-head as a blocker classifier: ROC-AUC and average precision with
   bootstrap CIs, against an explicit chance baseline.
5. The decision-relevant form: *if we had ranked by score S and taken the top
   96, how many blockers would we have?* Compare rankers directly. Beware
   composition artifacts — a pooled win can come entirely from one ranker
   preferring the round with the higher base rate. Check each round separately.

## Computing ipSAE

Implement from the reference implementation (`DunbrackLab/IPSAE`, accompanying
Dunbrack 2025, bioRxiv 10.1101/2025.02.10.637595) rather than from memory, and
validate against its output on a held-out complex before trusting a batch.

```
ptm(x, d0)   = 1 / (1 + (x/d0)^2)
d0(L)        = max(1.0, 1.24 * (L - 15)^(1/3) - 1.8)
valid[i,j]   = chain[i]==c1 and chain[j]==c2 and pae[i,j] < 10
n0res[i]     = sum_j valid[i,j]
d0res[i]     = d0(max(26, n0res[i]))
pSAE[i]      = mean over {j : valid[i,j]} of ptm(pae[i,j], d0res[i])
ipSAE_asym(c1→c2) = max_i pSAE[i]
```

The reported value is the **maximum** of the two chain directions (not the
minimum — a plausible-sounding misremembering). Variants worth computing:
`d0chn` (d0 from summed chain lengths), `d0dom` (d0 from residues passing the
cutoff in both chains).

The full PAE matrix is needed, not the scalar `min_interaction_pae`. In Boltz-2
it lives in the job's `archive` artifact as `result/pae.npz`. **Check one
archive contains it before downloading a thousand** — a feasibility check that
fails in two minutes beats one that fails in an hour.

## What this licenses, and what it does not

Licensed: that on this target, this modality, at this n, these four scores carry
no epitope-location information in the informative range; that geometry must be
computed; that score-based prefiltering of what to analyse is unsafe.

Not licensed: that these scores are worthless (they are real signal for complex
plausibility), that this holds for every target or modality, or any statement
about experimental hit rates. Computational confidence is not a hit rate, and a
within-campaign null is not a universal null.
