"""
Module 1 — enzymatic label QC for SARS-CoV-2 Mpro (3CLpro).

Assay being modelled in silico: a purified-enzyme FRET protease assay of the
BPS Bioscience #79955 type — recombinant 3CLpro + internally quenched
DABCYL-KTSAVLQ|SGFRKME-EDANS, cleaved at Q|S, read Ex360/Em460, GC376 as the
control inhibitor. Buffer and purified protein only: NO cell, no membrane, no
metabolism. Everything downstream here is the enzymatic axis only.

Labels come from ChEMBL target CHEMBL4523582 (SARS-CoV-2 replicase
polyprotein 1ab, which carries nsp5/3CLpro), standard_type=IC50,
assay_type=B (binding/enzymatic), pchembl_value present, no
data_validity_comment.

  ! Cellular antiviral activity does NOT live under this target. It sits under
    the organism-level target CHEMBL4303835. Joining on the protein target
    alone finds 15 pairs; the cross-target join finds 615. Do not confuse the
    two axes — this module is the enzymatic axis and stops there.

--------------------------------------------------------------------------
WHAT THIS MODULE CORRECTS
--------------------------------------------------------------------------
An earlier draft of the protocol asserted: "stratify by assay family, do not
pool the literature IC50 — ebselen's 3 log spread is the cost of pooling."
Measured on the actual data, BOTH halves of that are wrong:

  (a) Stratification is not affordable. 3,629 usable records are spread over
      507 distinct assays, median 2 records each, 180 of them singletons; the
      single largest assay holds 77 records. Per-assay models are not
      trainable at that size, so "stratify" would shatter the dataset rather
      than clean it.

  (b) Pooling is mostly fine. Of 265 compounds measured in >=2 distinct
      assays, the median cross-assay range is 0.25 log and 59.2% agree within
      0.5 log. The damage is not spread evenly across the set — it is
      concentrated in a heavy tail: 26.8% exceed 1 log, 10.9% exceed 2 log,
      P99 = 3.00 log. ebselen (4.99-8.00) is that tail's worst case, not the
      typical compound, so generalising from it was the error.

The rule that survives contact with the data: POOL, then use cross-assay
concordance as a per-compound label-quality tier, and quarantine the
discordant tail. That keeps the data volume and removes the noise, which
stratification cannot do.

--------------------------------------------------------------------------
PREINCUBATION (time-dependence of covalent inhibitors)
--------------------------------------------------------------------------
Nirmatrelvir (nitrile), GC376/GC373 (aldehyde + its bisulfite prodrug) and
PF-00835231 (hydroxymethylketone) all modify Cys145 covalently, so their
apparent IC50 depends on how long enzyme and inhibitor were preincubated.
Two records for one covalent compound at different preincubation times are
not repeat measurements of one quantity.

The earlier draft put preincubation reporting at 64.6% of records. That used
a loose regex that also matched plain substrate "incubated for 30 mins".
Tightened to 'pre-?incubat', it is ~35% of assays. So enforcing it as a hard
filter costs about two thirds of the data, not one third — too expensive to
be a filter. It is used here as a TIER instead: never silently average a
tier-A against a tier-B value for the same compound.

Tiers emitted per compound:
  concordance  OK (<=0.5 log) / SUSPECT (>1 log) / UNUSABLE (>2 log) / SINGLE
  preincub     A (assay reports preincubation) / B (does not) / MIXED
"""
import collections
import json
import os
import statistics

DATA = os.path.join(os.path.dirname(__file__), "data", "enz_dataset.json")

CONCORDANT, SUSPECT, UNUSABLE = 0.5, 1.0, 2.0   # log-unit range thresholds


def load():
    with open(DATA) as fh:
        return json.load(fh)


def tier_compounds(rows):
    """Per-compound label-quality tiers. Returns {mol: {...}}."""
    by_mol = collections.defaultdict(list)
    for r in rows:
        by_mol[r["mol"]].append(r)

    out = {}
    for mol, recs in by_mol.items():
        ps = [r["pIC50"] for r in recs]
        assays = {r["assay"] for r in recs}
        pre = {bool(r["preincub"]) for r in recs}

        if len(assays) < 2:
            conc, rng = "SINGLE", 0.0
        else:
            rng = max(ps) - min(ps)
            conc = ("UNUSABLE" if rng > UNUSABLE else
                    "SUSPECT" if rng > SUSPECT else
                    "OK" if rng <= CONCORDANT else "BORDERLINE")

        out[mol] = dict(
            pIC50=statistics.median(ps),      # median, not mean: tail-robust
            n=len(recs), n_assays=len(assays), range=rng, concordance=conc,
            preincub=("A" if pre == {True} else "B" if pre == {False} else "MIXED"),
            smiles=next((r["smiles"] for r in recs if r["smiles"]), None),
        )
    return out


def trainable(tiers, drop=("UNUSABLE", "SUSPECT")):
    """The modelling set: pooled, minus the discordant tail."""
    return {m: v for m, v in tiers.items() if v["concordance"] not in drop}


def report():
    rows = load()
    tiers = tier_compounds(rows)

    print("=" * 70)
    print("MODULE 1 — enzymatic label QC (Mpro / 3CLpro FRET assay axis)")
    print("=" * 70)
    print(f"records {len(rows)}   compounds {len(tiers)}   "
          f"assays {len({r['assay'] for r in rows})}")

    # assay fragmentation -> why stratification is not an option
    by_assay = collections.Counter(r["assay"] for r in rows)
    sizes = sorted(by_assay.values(), reverse=True)
    print("\nassay fragmentation (rules out per-assay stratification):")
    print(f"  median records/assay {statistics.median(sizes)}   "
          f"largest assay {sizes[0]}   singletons {sum(1 for s in sizes if s == 1)}")
    for thr in (5, 20):
        keep = [s for s in sizes if s >= thr]
        print(f"  assays with >={thr:>2} records: {len(keep):>3}  "
              f"covering {100 * sum(keep) / len(rows):.0f}% of records")

    # concordance distribution -> the heavy tail, not a uniform problem
    multi = [v for v in tiers.values() if v["n_assays"] >= 2]
    rngs = sorted(v["range"] for v in multi)
    print(f"\ncross-assay range, {len(multi)} compounds in >=2 assays:")
    for q in (50, 75, 90, 99):
        i = min(len(rngs) - 1, int(q / 100 * len(rngs)))
        print(f"  P{q:<2} = {rngs[i]:.2f} log")

    conc = collections.Counter(v["concordance"] for v in tiers.values())
    print("\nconcordance tiers:")
    for k in ("OK", "BORDERLINE", "SUSPECT", "UNUSABLE", "SINGLE"):
        if conc[k]:
            print(f"  {k:<11} {conc[k]:>5}  ({100 * conc[k] / len(tiers):.1f}%)")

    pre = collections.Counter(v["preincub"] for v in tiers.values())
    print("\npreincubation tiers (never average A against B):")
    for k in ("A", "B", "MIXED"):
        if pre[k]:
            print(f"  tier {k:<5} {pre[k]:>5}  ({100 * pre[k] / len(tiers):.1f}%)")

    keep = trainable(tiers)
    print(f"\ntrainable set after quarantining the discordant tail: "
          f"{len(keep)} / {len(tiers)} compounds "
          f"({100 * len(keep) / len(tiers):.1f}% retained)")
    print("  (pooling retained; stratification would have left <=77 per assay)")

    worst = sorted(multi, key=lambda v: -v["range"])[:5]
    print("\nworst-discordance compounds (labels not trustworthy at any model quality):")
    for v in worst:
        mol = next(m for m, x in tiers.items() if x is v)
        print(f"  {mol:<16} range {v['range']:.2f} log over {v['n_assays']} assays "
              f"(median pIC50 {v['pIC50']:.2f})")

    return dict(tiers=tiers, trainable=keep)


if __name__ == "__main__":
    report()
