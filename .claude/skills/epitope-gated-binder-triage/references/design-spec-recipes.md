# Design spec recipes (Boltz-2 protein design API)

Patterns that worked, and the reasoning behind each, so they transfer to a
different API rather than being copied blindly.

## Always estimate before you start

`boltz_estimate_protein_design` validates the spec and returns a cost without
starting a paid job. Use it as a schema check too — a spec that estimates
cleanly will run. Confirm the actual number with the person paying before
submitting; an estimate is free, a run is not.

Pricing observed: a flat per-design unit cost, so total = N × unit. This means
**generation is cheap and synthesis is expensive** — the two "how many"
questions have opposite answers (generate widely, order narrowly).

## Epitope-constrained target

```json
{
  "type": "structure_template",
  "structure": {"type": "url", "url": "<model cif url>"},
  "chain_selection": {
    "A": {
      "chain_type": "polymer",
      "crop_residues": [<0-indexed, structured core only>],
      "epitope_residues": [<0-indexed epitope>],
      "non_binding_residues": [<0-indexed sites to discourage>]
    }
  }
}
```

- indices are **0-indexed** into the chain; for UniProt-numbered models that is
  `uniprot_number - 1`. Write the mapping down and assert it.
- `crop_residues` excludes low-pLDDT termini — see `epitope-occupancy.md`.
- `non_binding_residues` is where glycosylation sequons near the epitope go. It
  reduces, not eliminates, the risk; the screening antigen still has to be
  glycosylated.
- A public model URL (e.g. AlphaFold DB) avoids base64-uploading a large CIF.
  Inline base64 of a ~120 KB structure is ~160 KB of text per spec, which is
  impractical to emit for more than one or two parents.

## Curated scaffolds for round 1

```json
{"type": "boltz_curated", "binder": "boltz_nanobody", "rules": {...}}
```

Right for an undirected first round: the provider maintains the scaffold list,
and you get framework diversity for free.

Note a limit that bit us: `rules` apply only to **designed** regions. Liabilities
present in a curated scaffold's *fixed* framework (an extra cysteine in CDR1, a
methionine) pass straight through. Round-1 output needs a framework audit even
though the rules looked strict.

## Developability rules

```json
"rules": {
  "excluded_amino_acids": ["C", "M"],
  "excluded_sequence_motifs": ["NXS", "NXT", "DP", "DG", "NG"],
  "max_hydrophobic_fraction": 0.40
}
```

Reasoning, so these can be adjusted rather than cargo-culted:
- **C** — an unpaired cysteine in a CDR is an aggregation and
  disulfide-scrambling risk. Drop this exclusion if you specifically want long
  CDR3s with an extra CDR1–CDR3 disulfide.
- **M** — oxidation liability. A real constraint on diversity; drop it if
  diversity matters more than oxidation risk for your assay.
- **NXS/NXT** — N-glycosylation sequons. `X` is a single-residue wildcard.
- **DP/DG** — acid-labile; **NG** — deamidation hotspot.
- **max_hydrophobic_fraction** — a hydrophobic epitope pulls designs toward
  sticky solutions. Cap it, and expect to tune the cap by epitope character.

## Framework locking (sequence-defined designs)

For iteration rounds, define the binder as literal fixed residues around a
designed length range rather than uploading a structure:

```
<fixed framework + CDR1 + CDR2 + FR3>  <min>..<max>  <fixed FR4>
```

Example: `EVQLVESGGG...TAVYYC` `14..18` `WGQGTLVTVSS`

Three things this buys:
- **Framework anomalies go to zero.** Anything not designed cannot be corrupted.
  Measured: 11–14% of round-1 gate-passers carried framework defects; 0/600 in a
  CDR3-only round.
- **You can repair a parent in the spec.** A parent with a damaged FR4 can be
  written with the canonical sequence, so its children inherit the fix.
- **No large uploads.** The whole spec is a few KB.

Use `uniformly_sampled_specifications` to put several parents in one job; it
samples uniformly across them, so allocation is equal. Weighting requires
separate jobs.

## Anchored CDR3 maturation — the iteration recipe that actually carries information

The obvious iteration ("seed from the best children, redesign their CDR3") has a
trap: in a sequence-defined spec, children of the same parent share a framework,
so the spec encodes only *framework + CDR3 length range*. The children's CDR3
sequences — the thing that was actually learned — are discarded, and six seeds
collapse to three identical specs. That is not maturation; it is re-rolling.

Instead, keep the proven CDR3's anchors and diversify the middle:

```
<framework> <CDR3 N-anchor, ~5 aa> <min>..<max> <CDR3 C-anchor, ~3 aa> <FR4>
```

The N-terminal anchor sets the loop's exit vector and the C-terminal anchor its
re-entry; the designed middle is where contact chemistry varies. Six seeds now
give six distinct specs.

Bias the length range upward if longer CDR3s performed better in the previous
round (they did in ours: +1/+2 versus parent length gave 50–57% blockers against
34% at parent length).

## Choosing parents for an iteration round

Cluster the gate survivors and pick representatives of distinct modes, not the
top N by score — the top N are frequently near-duplicates of one lineage.
Cluster on CDR3 sequence distance, CDR3 length, and epitope footprint.

Then, having run the round, check per-lineage rates before planning the next
one. In our round 2, one lineage of eleven produced 56.5% plate-grade designs
and the other ten produced 1.1% — worse than the undirected baseline. Most
"seed diversity" was seeding failure. A round 3 should drop non-responders.

But a caution that pulls the other way: seeding only the one winning lineage
concentrates the final plate on a single binding mode. If that mode is
experimentally wrong, the plate fails as a unit. Cap any lineage's share of the
final plate (~30%) even when the computational yield argues otherwise.

## Idempotency and job hygiene

Pass a stable, descriptive `idempotency_key` so a retry cannot double-charge.
Include target, epitope, round and N in the key.

There may be **no stop/cancel tool** exposed even when the API supports stopping
(a `stopped_at` field in job status implies it does). Check before promising the
person they can abort mid-run — and size the job accordingly.

Generation rate ramps: one run produced 5 designs in the first 21 minutes and
then held ~4.5/min. An early-rate extrapolation will be wildly pessimistic; wait
for the rate to stabilise before reporting an ETA.
