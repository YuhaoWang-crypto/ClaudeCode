# DeltaRASA and the underprocessing score, in depth

## 1. Why a steric metric works at all

N-glycan processing is enzymatic, and the enzymes are large. OST transfers
Glc3Man9GlcNAc2 co-translationally onto N-X-S/T (X != Pro). Glucosidases I/II
and ER mannosidase I then trim to Man5GlcNAc2, the substrate for GlcNAc-TI,
after which the glycan commits to complex-type. Each of those steps needs a
folded enzyme to dock on the glycan.

So composition is set by **accessibility at the moment each enzyme acts**, not
by sequence. Wall the glycan in and processing stalls wherever access runs
out. That is why a *protein design* problem can control a *carbohydrate*
outcome, and why one geometric scalar can predict it.

The processing axis, least to most processed:

```
Glc1Man9 -> Man9 -> Man8 -> Man7 -> Man6 -> Man5 -> hybrid -> complex -> truncated core
   ER intermediate            |              commitment point (GlcNAc-TI)
```

## 2. DeltaRASA, exactly

Per atom, with a 1.4 A probe:

```
RASA(i)       = SASA(i) / (4 * pi * (r_i + r_probe)^2)
pass 1 (context)   glycan atoms, protein present
pass 2 (isolated)  same glycan, protein atoms deleted, glycan atoms kept
DeltaRASA(i)  = | RASA_isolated(i) - RASA_context(i) |
```

Reported as a percentage.

**Why pass 2 keeps the glycan's own atoms.** A glycan occludes itself heavily -
ring atoms bury each other, and bonded neighbours bury the linkage atoms. If
pass 2 used a free monosaccharide, DeltaRASA would be dominated by that
self-burial and would barely vary with the protein. Keeping the glycan intact
in both passes cancels self-occlusion exactly, leaving protein-dependent
shielding. This is the single easiest thing to get wrong.

**The ceiling.** Because self-occlusion already consumes part of the
isolated-atom reference, a fully buried atom cannot reach DeltaRASA = 100%.
The paper's per-atom colour bar stops at 50%. Treat DeltaRASA as a relative
shielding score only; "40%" is not "40% buried".

**Summary statistics.** The paper ranks on two, and inspects designs that are
high on *either*:
- mean per-atom DeltaRASA over the tree
- number of glycan atoms with non-zero DeltaRASA ("contact atoms")

They answer different questions: the mean asks *how hard* the protein presses
on the glycan, the count asks *how much* of the glycan the protein touches. A
deep narrow pocket scores high on the mean; a broad shallow shelf (the paper's
O2) scores high on the count. Both produced oligomannose.

## 3. Numerical practice

| Parameter | Value | Why |
|---|---|---|
| probe radius | 1.40 A | paper; water |
| atomic radii | Biopython `ATOMIC_RADII` | matches the tool the paper names |
| sphere points | 960 (default here) | 100 leaves ~1% relative error |
| hydrogens | excluded | crystal structures rarely have them; keeps heavy-atom counts comparable |
| altlocs | `' '` or `'A'` only | double-counting altlocs inflates occlusion |

Convergence on a 116-atom Man8 (5FYL G:N332), deviation from `n_points=1920`:

```
100   0.0475%      960   0.0049%
240   0.0022%     1920   0  (reference)
480   0.0148%
```

Golden-spiral quadrature is not monotone in `n_points`, so do not read the
240-point row as "better than 480"; read the envelope. Everything from 240 up
is within 0.015%.

**The real precision floor** comes from replicates, not quadrature. In a
C3-symmetric assembly each site appears three times and each copy is an
independent trip through parsing, tree building, neighbour search and
quadrature. On 5FYL assembly 1 those copies agree to **0.012% median, 0.082%
max**. Differences below that are not real.

## 4. The two confounds, measured

**Size.** Mean per-atom DeltaRASA averages over the whole tree, so it partly
measures how much tree there is. On 5FYL, Spearman(atoms modelled,
DeltaRASA) = **-0.38**; truncating every tree to the shared GlcNAc-GlcNAc-Man
core drops it to **-0.07**, and the standardised ranking agrees with the raw
one at only rho = 0.82. The fix is what the paper does: model **one identical
glycan** onto every design (GlycoShape Man9; Man6 when using the AF3 server,
which caps glycan size) so the comparison is like for like.

**Missing coordinates.** Glycan arms are mobile and often unresolved. Median
modelled size in 5FYL is **55%** of an intact Man9GlcNAc2 and 45% of sites are
below half. Consequently the neighbour-glycan shielding term is a lower bound,
zero at 59% of sites, and a deposited structure cannot demonstrate
crowding-driven oligomannose even where it genuinely occurs. Use predicted
structures with a standard modelled glycan for screening; use crystal
structures to check the fold, not the glycan.

## 5. The underprocessing score

Abundance-weighted sum over an ordinal ladder:

```
score = sum over species of ( relative abundance % )  x  ( ordinal rung )
```

Range 0-2000 (100% at rung 20). Paper's selection gate: **score > 1500 and
occupancy > 50%**.

**The ladder is reconstructed, not published.** The paper states only the
endpoints (20 = Glc1Man9; low = processed complex/core), the 0-2000 axis, and
the naming rules. `glycans.UNDERPROCESSING_LADDER` rebuilds the interior from
the categorical axis the paper uses for its own site-specific plots, which
expands to exactly twenty species, and which puts the 1500 cutoff at
"abundance-weighted mean rung >= 15", i.e. Man5-or-less-processed on average.
That is corroboration, not proof.

**What survives the uncertainty.** If profile A stochastically dominates
profile B along the processing axis (upper-tail cumulative abundance at least
as large at every rung), then score(A) > score(B) under *every* strictly
increasing ladder. `g3` finds 34 such pairs and verifies them against 500
random monotone ladders each: **0 violations in 17000 sweeps**. So rankings
between cleanly ordered profiles are safe.

**What does not.** The absolute 1500 threshold. Across 2000 random monotone
ladders sharing the published endpoints, all four O-designs flip pass/fail -
O1 passes under 66% of ladders, O4 72%, O3 81%, O2 98%. Complex controls stay
robustly below. So: trust "more underprocessed than", distrust "scores 1620".

## 6. Classification rules (these *are* published)

- `HexNAc(2)Hex(9-3)` -> Man9...Man3
- fucosylated oligomannose -> FM
- complex -> graded by HexNAc count and core fucosylation
- anything smaller than Man3 -> core / truncated
- isomers of equal composition are not distinguished by LC-MS/MS

Occupancy and composition are separate axes and must stay separate: abundances
are relative *within the occupied fraction*, so a site can be 95% Man9 and
still irrelevant at 4% occupancy.
