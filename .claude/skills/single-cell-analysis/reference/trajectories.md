# Trajectories: pseudotime, RNA velocity, lineage tracing

Source: sc-best-practices *Inferring trajectories*; Bergen et al. 2021,
*"RNA velocity — current challenges and future perspectives"* (Mol Syst Biol).

The whole area is ⚠️ **convention-grade**: powerful, widely used, and far
easier to over-interpret than anything else in the workflow.

## Precondition: is there a trajectory at all?

A trajectory method will return a trajectory from any dataset, including PBMCs
in steady state and randomly shuffled counts. Before running one, establish
that a continuous process actually exists: a known developmental/differentiation
axis, cells forming a connected continuum rather than discrete islands, and
markers changing gradually rather than switching.

PAGA is the honest first step — it tests *connectivity* between clusters
instead of assuming a path:

```python
sc.tl.paga(adata, groups="celltype")
sc.pl.paga(adata, threshold=0.1)      # weak edges = no evidence of a transition
```

## Pseudotemporal ordering

| Method | Good for |
|---|---|
| **DPT** (`sc.tl.dpt`) | simple linear/branching, diffusion-based |
| **Palantir** | branching + terminal-state probabilities + entropy |
| **Slingshot** (R) | cluster-based curves, well-behaved on modest data |
| **Monocle 3** | large branching graphs |
| **CellRank 2** | the general framework — combines pseudotime, velocity, similarity or metabolic labelling into a Markov transition matrix |

```python
adata.uns["iroot"] = np.flatnonzero(adata.obs["celltype"] == "HSC")[0]
sc.tl.diffmap(adata); sc.tl.dpt(adata)
```

**Caveats**
- The **root cell must be chosen by biology**, not convenience — pseudotime is
  a distance from the root, so an arbitrary root gives an arbitrary ordering.
- Pseudotime is **monotone ordering, not real time**; the spacing carries no
  temporal information.
- Branch assignment is unstable near branch points. Report the uncertainty
  (Palantir/CellRank give per-cell fate probabilities — use them).
- Trajectories drawn on a UMAP are a projection of a projection. Compute in
  the PCA/diffusion/latent space.
- Testing genes "along pseudotime" reuses the data that built the pseudotime —
  the same circularity as marker p-values on Leiden clusters.

## RNA velocity ⚠️ — the strictest caveats in the book

Spliced/unspliced ratios estimate d(spliced)/dt, giving a *direction* rather
than just an ordering. The counting must be decided at quantification time
(velocyto or alevin-fry with velocity mode) — the layers cannot be recovered later.

```python
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
scv.tl.velocity(adata, mode="stochastic")          # or "dynamical" (EM, slower)
# scv.tl.recover_dynamics(adata) first if mode="dynamical"
scv.tl.velocity_graph(adata)
```

**When velocity should not be trusted — check all four:**

1. **Timescale mismatch.** The process must move on the timescale of mRNA
   half-lives (hours). Not applicable to steady-state systems (peripheral
   blood) or slow disease processes (Alzheimer's, Parkinson's) — the model has
   nothing to measure.
2. **Model assumptions.** Constant transcription/splicing/degradation rates per
   gene. Violated in e.g. erythroid maturation, where the EM model returns
   confidently wrong directions. **Look at the phase portraits** — you need the
   almond/hysteresis loop shape; multiple kinetic regimes mean stop.
3. **Projection artifacts.** "Projecting velocity vectors onto a low-dimensional
   representation may be misleading" — the arrows depend on gene count,
   neighbour count and plotting parameters, and degrade at embedding
   boundaries. Streamlines on a UMAP are the least reliable output of the
   entire pipeline.
4. **Ambient RNA and unspliced contamination** inflate unspliced counts and
   fabricate velocity.

**Best practice:** treat velocity as *one input* to CellRank rather than an
answer. Quantify with fate probabilities and driver genes; validate the
direction against an independent source (known markers, time points, lineage
barcodes). `scv.tl.velocity_confidence` and consistency scores are the minimum
reporting bar.

Newer: veloVI / DeepVelo (deep generative velocity with uncertainty),
`scv.tl.latent_time` for a global gene-shared time.

## Lineage tracing

CRISPR/barcode-based clonal tracking (LARRY, GESTALT, ScarTrace, CellTag) is
the ground truth that velocity and pseudotime approximate. Tools: Cassiopeia
(tree reconstruction), CoSpar (transition maps from clonal data),
LineageOT. If lineage barcodes exist, they arbitrate — do not defer to
velocity arrows over clonal data.

## Reporting standard

- State the root and why it is biologically justified.
- Show the branch/fate probabilities, not just the assigned branch.
- Show the same conclusion in ≥2 methods (e.g. DPT + Palantir, or
  velocity + CellRank fate probabilities), or say that you could not.
- Never present a velocity streamplot as the sole evidence for a direction.
