# Failure modes

Every trap hit during a real campaign, with symptom, cause, and the check that
catches it. The ones marked **silent** produce output that looks correct.

Read this before a long run. Most cost minutes to guard against and hours to
discover afterwards.

---

## Silent: residue numbering shift

**Symptom** — none. Occupancy numbers are produced, formatted, plausible, wrong.

**Cause** — PDB mature-protein numbering vs UniProt (signal peptide included);
a cropped chain renumbered from 1; an assumed offset that is off by one.

**Check** — assert per file that the target chain spans the expected range *and*
that its sequence equals the reference exactly. Refuse to score files that fail;
count and report them. Across 1,600 files this cost nothing and made the whole
analysis falsifiable.

---

## Silent: trusting a reported position instead of relocating it

**Symptom** — an edit lands one residue away from the defect. The defect
survives; a healthy residue is mutated.

**Cause** — a report said "free cysteine at position 32". It was at 33.

**Check** — find anomalies by pattern in the actual sequence, then edit. The
edit is then correct regardless of whether the report was. Applies to any
position handed to you by an upstream summary.

---

## Silent: a repair that creates a new defect

**Symptom** — QC passes on the thing you fixed. A different liability appears
elsewhere and nobody looks.

**Cause** — reverting a free cysteine `C→S` turned `N-T-C` into `N-T-S`: a new
N-glycosylation sequon, next to the paratope.

**Check** — make every substitution sequon-aware (try S, then A, then T, …,
accepting the first that does not increase the sequon count), and **re-run the
entire QC suite on the final sequences**, not just the property you edited. This
one was caught only by the final batch check.

---

## Silent: full-range correlation on a score with an undocked tail

**Symptom** — a strong, highly significant positive correlation. It is an
artifact.

**Cause** — the bottom of the score distribution is designs with essentially no
interface. "Low score → low recall" is true and vacuous, and it generates the
whole correlation.

**Check** — stratify. Report the restricted number as the result. If a positive
correlation is claimed, show it survives restriction to the informative range.
Applies equally to any new score someone proposes.

---

## Silent: score-based prefiltering of what to analyse

**Symptom** — an analysis that looks complete and has quietly discarded most of
the hits.

**Cause** — "designs below threshold aren't worth scoring" is an assumption, not
a measurement. The discarded stratum held 69% of blockers.

**Check** — before trusting the filter, score a random seeded sample of what it
would discard and compare rates. If the discarded stratum is much larger, even a
real 2× enrichment loses more than it saves.

---

## Silent: an aglycosylated model hiding a glycan over the epitope

**Symptom** — perfect designs that fail against the native target.

**Cause** — the predicted structure has no glycans. A sequon 9.8 Å from an
epitope residue carries a glycan that spans far further than that.

**Check** — scan the target for `N-X-[ST]` (X ≠ P) and compute 3D distance to
the epitope. If close: discourage contact there in the spec, and require a
glycosylated (mammalian-expressed) screening antigen.

---

## Silent: rules that only apply to designed regions

**Symptom** — a library built with strict developability rules still carries
free cysteines and methionines.

**Cause** — `rules` constrain designed segments. Liabilities inside a curated
scaffold's *fixed* framework pass through untouched.

**Check** — audit the framework of every gate-passer independently of the spec.
Better: lock the framework as explicit fixed sequence, which removes the class
of problem rather than detecting it.

---

## Silent: seeds that collapse to identical specs

**Symptom** — an "iteration from six diverse parents" that is really three
identical jobs. Wasted spend, no maturation.

**Cause** — in a sequence-defined spec, children of one parent share a
framework, so the spec reduces to framework + CDR3 length range. The learned
CDR3s are discarded.

**Check** — print the distinct spec strings before submitting and count them. If
`distinct < seeds`, the design is wrong — use anchored CDR3 maturation instead
(`design-spec-recipes.md`).

---

## Silent: pooled statistics hiding per-group failure

**Symptom** — "blocker rate improved 28% → 44%, p = 2.5e-10". True, and
misleading.

**Cause** — one of eleven lineages produced 56.5% plate-grade designs; the other
ten produced 1.1%, *worse than the undirected baseline*. On the joint bar that
actually fills a plate, the round was not better per design (5.33% vs 6.40%,
p = 0.45).

**Check** — always report per-group rates alongside the pooled one, and evaluate
on the bar you would actually order at, not the loosest bar available. Ask
whether the ceiling moved, not just the mean.

---

## Silent: composition artifacts in ranker comparisons

**Symptom** — ranker B captures more blockers than ranker A in the pooled data.

**Cause** — B's top-96 was drawn mostly from a round with a higher base rate. It
reversed within each round separately.

**Check** — evaluate ranker comparisons within each population separately, and
compare against a chance baseline computed for that population.

---

## Loud but expensive: context blowout from API result payloads

**Symptom** — an agent's context is consumed by presigned URL text.

**Cause** — one page of 100 design results is ~0.5 MB, almost all URL.

**Check** — never load result payloads into context. Paginate, let oversized
responses spill to disk, parse with `jq`/Python, keep only the metric fields.
Delegate bulk scoring to a subagent with instructions to report statistics only.

---

## Loud: presigned URLs expiring mid-job

**Cause** — typically 30-minute expiry. Collecting 1,000 URLs before downloading
guarantees the tail expires.

**Check** — fetch in batches of ~20 and download each batch immediately. Cache
structures on disk; later analyses will want them.

---

## Loud: extrapolating from a ramping generation rate

**Symptom** — an ETA of 70 hours reported for a job that took 55 minutes.

**Cause** — 5 designs in the first 21 minutes, then ~4.5/min once workers spun up.

**Check** — wait for the rate to stabilise before quoting an ETA, and say the
early rate is a ramp when reporting it.

---

## Loud: assuming a job can be stopped

**Symptom** — promising the person they can abort, then discovering no
stop/cancel tool exists even though the API clearly supports it.

**Check** — verify the cancel path exists before starting, and size jobs so that
running to completion is acceptable.

---

## Process: a feasibility check that fails fast

Before a bulk download, verify the data you need is actually in the artifact.
"Does one archive contain a PAE matrix?" is a two-minute question whose answer
determines whether a one-hour job is worth starting. Build the check into the
task rather than discovering it at the end.

---

## Process: predictions stated before the round

State falsifiable predictions before running, and report the outcome either way.
Ours — "score-driven optimisation will push designs off the epitope" — was
falsified (2.5% showed the signature, sign test p = 0.20). Saying so plainly was
more useful than the prediction would have been if it had held, and it prevented
a wrong conclusion from being carried forward.

---

## Process: the conflation that wastes the most time

A speed result is not an accuracy result. When evaluating whether some new tool
or paper helps, ask first: does it address the bottleneck, or a different one?
A 4× inference speedup on a pipeline whose bottleneck is *which designs to
believe* changes the cost of being wrong, not the rate of being wrong.

Related: computational confidence is not a hit rate. Numbers sharing the
interval [0,1] are not thereby comparable.
