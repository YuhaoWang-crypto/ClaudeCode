# From hit list to experiment

Mining produces candidates. This file is how the source study turned 946
candidates into two characterized enzymes, and how to repeat that.

## 1. Sequence-similarity network (SSN) — organize the hits

Build an SSN over the hit sequences to see families rather than a flat list.
Tool: EFI-EST (<https://efi.igb.illinois.edu/efi-est/>), visualize in Cytoscape.

- **30% identity** for the global picture — this is the cutoff at which the
  published atlas resolved the 6 known families plus 70 unannotated clusters.
- **50% identity** on one cluster of interest, after enriching it with BLAST
  homologues (the paper used E-value 5 for maximum inclusiveness), to resolve
  subfamilies and to see the *non-halogenase counterparts* interleaved with the
  hits — that mixture is the signal that sequence clustering alone cannot make
  this call.
- Annotate every node with the mined class: hit vs contrast class. A cluster
  containing both is the most informative kind, and the best place to pick a
  test case.

Sanity check that must pass: **known members land where they should**. If the
SSN does not recover the known families as clusters, something upstream is wrong.

## 2. Genome neighbourhood — get the substrate hypothesis

The motif tells you the chemistry, never the substrate. Genomic context does.
Tool: EFI-GNT (<https://efi.igb.illinois.edu/efi-gnt/>).

What the paper read off the neighbourhoods:
- co-clustering with **acyl-carrier-protein / PKS / NRPS** genes → the enzyme
  acts on an ACP-tethered substrate;
- no such genes, but **amino-acid transporters and amino-acid-modifying enzymes
  (e.g. ATP-grasp ligases)** nearby → free-standing amino-acid substrate. This is
  exactly how AspX was predicted to act on a free amino acid — and it chlorinates
  L-aspartate.
- **biosynthetic gene clusters** with neither → a novel substrate class; BtnX sat
  at that branch point and turned out to chlorinate biotin.

## 3. Taxonomy and literature

- Taxonomic annotation of clusters surfaces under-explored space (the study found
  new **eukaryotic** halogenase families — fungal and plant — a niche where
  homology search had found almost nothing).
- **PaperBLAST** (<https://papers.genomics.lbl.gov/>) finds published work on
  homologues. This is how BtnX's context emerged: its gene sits on the "killer
  plasmid" of *Dinoroseobacter shibae*, already studied in a symbiosis context,
  which made biotin a natural substrate hypothesis.

## 4. Choosing what to test

Rank candidates by:

1. **Cluster novelty** — no member of the cluster has any prior annotation for
   the target activity.
2. **Substrate hypothesis strength** — genome context or literature points at a
   specific molecule you can buy or make.
3. **Model quality at the site** — high pLDDT over the anchor and shell residues.
4. **Expressibility** — bacterial origin, moderate length, no membrane anchor.
5. **Branch-point position** — enzymes between two known classes carry the most
   information per experiment (BtnX, at the ACP/free-standing junction, is where
   the unprecedented promiscuity was found).

Optional in-silico step before the bench: co-fold a candidate with its
hypothesized substrate and cofactor (Boltz-2 or equivalent) and check that the
target C–H sits next to the open coordination position. Treat as ⚠️ supporting
evidence, never as validation — it shares the failure modes of the model that
produced the structure being mined.

## 5. Experimental validation pattern

The published sequence, worth copying: heterologous expression → in-vitro assay
against a **panel** (all 20 amino acids for AspX; a broad substrate set for BtnX)
by LC–QTOF → product identification by NMR → steady-state kinetics (kcat, Km) →
alternative-anion scope (Br⁻, N₃⁻) → crystal structure with substrate + cofactor
to confirm the predicted site.

Two details worth stealing: assay a panel rather than one guess (that is how the
promiscuity of BtnX was found at all), and run the **anion-swap** experiment —
it is cheap and it tests the mechanistic claim behind the motif directly.

## 6. The engineering corollary — read the motif backwards

If the discriminator causes the function, installing it should install the
function. In BtnX, **G117D/G117E** — putting the carboxylate back into the open
coordination site — abolished halogenation and raised hydroxylation **40 ± 9-fold
and 18 ± 6-fold**.

So every validated motif yields two things: a screen, and a design rule. When a
hit is confirmed, propose the reciprocal mutation as the control experiment; when
someone wants an enzyme to switch reaction outcome, the motif says which residue
to change. The complementary structural lesson from BtnX: reaction outcome is set
by first-shell coordination, while *substrate scope* is set by the pocket —
a solvent-exposed channel plus a few specific H-bonds around the target C–H gave
promiscuity, whereas enclosed pockets (CurA, AspX, WelO5) gave narrow scope.

## 7. Reporting

Every reported number carries its label:

- ✅ computed — structures scanned, sites found, hits, contrast-class counts,
  benchmark recall/specificity, cluster counts, coordination-type census.
- ⚠️ predicted — "these are halogenases", substrate assignments, novelty claims
  ("70 new families" is ⚠️ until a member of each is assayed).

State the benchmark, the thresholds and the database versions next to the counts,
and say plainly how many candidates were actually tested. In the source study:
946 candidates, 2 characterized.
