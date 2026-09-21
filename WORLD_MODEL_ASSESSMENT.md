# Is the VCWM buildable? An assessment against this repository's measurements

Xing & Song, *A world model of the virtual cell*, Cell 189, 5831–5843 (17 Sept
2026), doi 10.1016/j.cell.2026.08.042.

The paper is a Perspective: it proposes an architecture `(E, F, D)` over a
biological network `G = (V, E)`, a three-tier data framework, and evaluation
principles. It trains nothing and reports no results. So "can it be built?"
has to be split into three questions that have three different answers, and
this project happens to hold direct measurements bearing on the third.

Everything quoted below is reproducible from `results/*.json` in this repo.

---

## 1. As code: yes, and a minimal version is already here

    z   = E(x; G, e)          encoder onto a cellular manifold M
    z'  = F(z, a)             action-conditioned transition
    x'  = D(z; G)             structured decoder, cross-modal consistent

Nothing in this is architecturally novel or blocked. A GNN over a pathway graph
with modality-specific foundation-model embeddings attached to nodes, a
flow-matching or diffusion transition core, and a multi-head decoder tied by a
shared latent — every piece has a reference implementation. Weeks, not years.

`virtualcell/` already contains a degenerate instance of exactly this shape:
`E` is a rank-40–80 projection chosen by cross-validation, `F` is a linear
state-conditioned operator, `D` is the transpose back into gene space. That it
is small is the point of §3 below, not an accident of effort.

## 2. As a training problem: the paper's own data tier 3 has no supervision

The paper factors the data requirement into three tiers, and is candid about
the third:

| tier | what it supervises | data situation |
|---|---|---|
| multi-modal paired/partially-paired | the manifold `M` | scarce but real and growing (CITE-seq, multiome) |
| perturbation–response | one-step `F` | abundant — this is what the whole field has |
| longitudinal / multi-step | `F⁽ᵏ⁾`, the persistent state | **essentially absent** |

Tier 3 is not a side condition. *Persistent state* and *long-horizon rollout*
are the entire claimed differentiator from the predictive models the paper
positions against. The paper's own fallback is RNA velocity, and it then states
that velocity inference is "technically fragile" and that "the quality of the
learned world model may depend as much on the fidelity of these inferred
dynamical signals as on the architecture used to model them."

That is the correct diagnosis, and it is a measurement problem, not an
engineering one: scRNA-seq is destructive, so the same cell cannot be read
twice. No amount of compute produces `(z_t, a_t, z_{t+1})` triples from
snapshots. The work in this repo hit the same wall from the other side —
`virtualcell/state_space.py` documents which parts of a state-space framework
the VCC data can and cannot reach, and velocity fields, fiber variables, path
dependence and population flow are all in the second category for the same
reason.

## 3. As the claimed capability: five measurements that constrain it

### 3.1 The transcriptome channel is already at its noise ceiling

`results/replicate_ceiling.json` — two independent arms of the *same* K562
CRISPRi screen, 2,053 knockdowns over 8,202 genes:

    replicate vs replicate      median r = 0.319     mean r = 0.340

`results/context.json` — this project's model, four held-out cell lines:

    ContextTransfer             mean r = 0.328       per-fold 0.304 0.384 0.348 0.276

These are not the identical quantity (different gene space, cross-line rather
than same-line), so this is an order-of-magnitude statement, not an equality:
the headroom on this channel is a few percent, not a few fold.

This is the sharpest thing this repo can say to the paper's closing claim that
"coherence can emerge from scale." On the transcriptome projection, scale has
almost nothing left to buy, because the ground truth's own reproducibility is
the binding constraint. A model reported to beat this margin substantially is
fitting batch structure or has an evaluation defect, and should be read that
way first.

### 3.2 Existing foundation models as "building blocks": measured contribution ≈ 0

The paper's Table 2 lists modality-specific models as the components a VCWM
"builds upon." Four were tested here on the one task that could be scored:

| model | result |
|---|---|
| Arc **Stack**-Large-Aligned | PDS **0.505** = chance; ours 0.667, naive transfer 0.713, baseline 0.506 — same fold, same 100 knockdowns, all arms restricted to the 4,606 of 6,642 genes Stack covers |
| Arc **SE-600M** protein embeddings | nothing on cross-context (identical to 4 d.p. for `esm_mix` 0→1) and on *unseen target*; on the one axis where the source has no measurement at all, r 0.234 → 0.242, PDS 0.503 → 0.509, overlap@100 0.173 → 0.169 |
| **ST-SE-Tahoe** | structurally inapplicable — `pert_rep: onehot` over `drugname_drugconc` |
| **AlphaGenome** | null on cross-context difference, with a positive control that worked |

The Stack result was predicted before it was computed, by the magnitude-spread
diagnostic: L1 retrieval pays for across-perturbation spread, so any inference
procedure that pulls every condition toward a common profile scores chance by
construction, however good its biology. Stated boundary: Stack was post-trained
on chemical and cytokine perturbations and is being asked for in-context
transfer of genetic knockdowns.

The paper says a VCWM "cannot be built by merely calling a library of disparate
modality-specific models via an agentic system." Agreed — and these numbers
extend it. Welding the same models together with a GNN does not automatically
produce cross-context transferability either, because on this task their
marginal contribution over a tuned linear transfer model is not distinguishable
from zero. The integration has to earn the transfer; inheriting it from the
components is not supported here.

### 3.3 The paper's central premise about `F` is correct, and cheap to satisfy

> "an action does not have a fixed outcome, but acts conditionally on the
> current cellular state"

`results/connection.json` quantifies this. Take a knockdown's effect measured
in line A, predict it in line B, 12 ordered pairs over 4 lines, k = 30:

| arm | mean cosine |
|---|---|
| as-is (state-free transport) | **+0.163** |
| low-rank projection only | +0.022 |
| random rotation (control) | +0.005 |
| transported by `R_AB = UVᵀ` | **+0.295** |

State-conditioning is worth +81% over moving the effect unchanged, and both
controls sit at zero, so the gain is the *alignment* and not the rank
restriction or a generic rotation.

The cost of estimating it is the surprise. From the sample-size curve, with
only **25 measured knockdowns** in the new context the transported cosine is
already 0.230 against 0.175 for as-is; 1,000 knockdowns buys 0.309. Most of the
state-dependence of `F` is recoverable from a couple of dozen anchors and a
linear orthogonal operator.

### 3.4 The manifold assumption holds; the cone version of it does not

`results/state_space.json` — Replogle genome-wide K562, 9,866 perturbations
× 8,248 genes, against a null with per-perturbation magnitudes preserved and
directions randomised:

    participation ratio      25.4   (null 332)
    directions for 50% var     18   (null 141)
    directions for 90% var    267   (null 337)

So "low-dimensional cellular manifold" is measurable, not rhetorical. But the
sharper geometric claim fails: the reachable set is **not a cone**. One-sidedness
of the 20 leading directions has median 0.518 against a null floor of 0.503,
with only 2 of 20 above 0.65. Anyone hardening a VCWM with a conic constraint on
reachable displacements would be enforcing something the data does not show.

### 3.5 The connection between contexts is flat, which is good news

`results/connection.json` — holonomy around every 3-cycle of the four lines
lands in 4.374–4.693 against a null mean of 7.737, and path inconsistency for
A→B→C vs A→C is the same magnitude. The transport is close to path-independent,
so contexts can be aligned to a single common frame rather than requiring a
separately learned operator per ordered pair. That is a real simplification for
anyone building the paper's `F`.

### 3.6 A ceiling on "organised by biology, not by context"

`results/robustness_index.json` — of 1,600 measured effect profiles, the
nearest neighbour *in another cell line* is the same knockdown 26.1% of the
time, 313× chance; the same test on control profiles with no perturbation
signal gives 0.25%. This is a property of the measurements, not of any model,
so it bounds what any VCWM trained on these profiles can decode.

---

## 4. Two places the paper should be pressed

**It exempts itself from the only existing benchmark.** Arc's leaderboards
"capture only a narrow projection of cellular state and are therefore not
indicative of the quality of a virtual cell," and a not-yet-existing virtual
cell "arena" is proposed in their place. The first half is fair — one modality
is one projection. The consequence is not: as written, the framework currently
has no falsifiable evaluation at all. And §3.2 shows "narrow" does not mean
"easy" — a published perturbation foundation model scores at chance on that
narrow projection. Clearing it should be a precondition for claims about
multi-modal coherence, not something to be routed around. Both authors are
co-founders of GenBio AI (chief scientist and CTO), and the proposed arena
aligns with their own AIDO programme; that does not make the argument wrong,
but it belongs in the margin while reading it.

**"Coherence can emerge from scale" collides with §3.1.** On the one channel
where the claim is checkable here, the measurement's own reproducibility is
already the ceiling.

---

## 5. Verdict, and the buildable subset

**The architecture is implementable. It is not scalable to the advertised
capability without new measurement technology, and that is the paper's own
conclusion read strictly.** Its most valuable section is the call for data, not
the architecture diagram.

What is worth building *now*, and what it would honestly not do:

1. `z` — low-rank latent, rank 40–80, which cross-validation here chose
   independently of any geometry analysis (`gwps_single_source.json`).
2. `F` — state-conditioned orthogonal operator per context pair, calibrated
   from ~25 anchor knockdowns in a new context, composed through a single
   common frame given the flat holonomy in §3.5.
3. `D` — gene-space decoder, scored against the replicate ceiling as
   denominator rather than against zero.

It would deliver: cross-context one-step transfer at roughly the noise floor of
the assay. It would **not** deliver trajectories, cross-modal coherence, or
compositional actions — not for want of architecture, but because none of the
three has supervision here.

What would actually unlock the paper's version is three kinds of data, in this
order:

1. **non-destructive or live-cell molecular readout** — without it tier 3 stays
   empty and `F⁽ᵏ⁾` stays unsupervised and unevaluable;
2. **combinatorial perturbation × order × time**, which is what makes a
   sequential world model beat an additive one — and is the exact baseline the
   paper's own Box 1 names as the falsification test;
3. **the same perturbation across many contexts with paired multi-modal
   readout**, which is what §3.3 shows is the cheap, high-yield axis.

None of the three is an AI problem.
