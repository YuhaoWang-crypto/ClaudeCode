"""Does AlphaGenome know how well CRISPRi will silence a gene in a new context?

AlphaGenome is a better candidate than a pure sequence model, and the reason is
worth being precise about.  A DNA language model has no cell-type input at all,
so its score for a gene is the same everywhere and it cannot express a
cross-context term.  AlphaGenome's outputs are *per biosample* -- and the
biosamples include K562, HepG2 and ``Jurkat, Clone E6-1``, three of the four
lines this project benchmarks on, plus ``keratinocyte`` and ``T-cell``, the
closest available matches to official contexts C and A.

What it still cannot do is the task.  AlphaGenome is a **cis** model: sequence
around a gene predicts *that gene's* expression and accessibility.  Perturbation
transfer is a **trans** problem -- knock down G, and thousands of other genes
move.  Nothing in a locus-local model expresses "G's protein is gone, therefore H
changed".

So most of what it could offer is already measured here and measured better:
context baseline expression comes from 18,400 real control cells per context,
and context similarity comes from those same profiles.  Exactly one quantity is
missing from the official contexts and present in AlphaGenome's output space:

    **how completely CRISPRi silences the target, in a context where no
    perturbation has ever been measured.**

The model currently transfers K562's measured on-target residual to every
context.  If promoter activity predicted by sequence tracked the residual across
contexts, that would be genuinely new information about A, B and C.

This module tests that on data where the answer is known: measured on-target
residual for the same genes in K562, HepG2 and Jurkat, against AlphaGenome's
predicted promoter output in the matching biosamples.

The prior is not encouraging and is stated up front so the result cannot be
rationalised afterwards: on-target efficiency is conserved *between* lines at
Spearman 0.38-0.54, which says it is mostly a property of the guide and the
gene rather than of the context, and target expression does not predict effect
magnitude across contexts at all (Spearman -0.014).  A weak or absent
correlation is the expected outcome.

**Result, and the cheaper way to have got it.**  All six correlations came out
indistinguishable from zero, with a positive control confirming the pipeline
works (predicted-vs-measured *control expression* ratios at Spearman +0.16 to
+0.29) and the target genuinely varying (median |log2 ratio| of the residual
between contexts, 0.83).

But the ``--gate-only`` track-resolution check answers it in seconds, without
spending a single prediction, and answers it more decisively::

    diagonal_mean 0.592, off_diagonal_mean 0.583, margin 0.009
    best_match: K562->HepG2, HepG2->HepG2, Jurkat->HepG2      usable: False

Every line's predicted track matches *HepG2's* measured profile best.  The head
does not resolve these three lines from each other at all, so no cross-context
quantity built on it could have worked -- the on-target null is a symptom, not
the finding.  Run this gate before building anything on a sequence model's
cell-type head.

    python -m virtualcell.alphagenome_prior --gate-only     # the gate alone
    python -m virtualcell.alphagenome_prior --n-genes 150   # the full test
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np

from .data import load_all

OUT = Path("results/alphagenome_ontarget.json")

# AlphaGenome biosample names for the lines this project measures.  RPE1 has no
# counterpart in the track set, so the fourth fold cannot be tested.
BIOSAMPLE = {"K562": "K562", "HepG2": "HepG2", "Jurkat": "Jurkat, Clone E6-1"}
# The API selects tracks by ontology term, not by biosample name.
ONTOLOGY = {"K562": "EFO:0002067", "HepG2": "EFO:0001187",
            "Jurkat": "CLO:0007045"}
WINDOW = 2 ** 17          # 131 kb around the TSS: the smallest supported window
PROMOTER = 2_000          # bp either side of the TSS to summarise over


def measured_residual(lines, symbols) -> tuple[dict[str, dict[str, float]], list]:
    """Measured on-target residual expression per gene per line.

    ``expm1(perturbed) / expm1(control)`` at the targeted gene's own column, so
    0.10 means CRISPRi left a tenth of the transcript behind.  Only genes that
    are both targeted and measured, and expressed well enough in every line for
    the ratio to mean anything, are returned.
    """
    out: dict[str, dict[str, float]] = {}
    for cl in lines:
        if cl.name not in BIOSAMPLE:
            continue
        col = {str(s): i for i, s in enumerate(cl.symbols)}
        idx = {n: i for i, n in enumerate(cl.names)}
        per: dict[str, float] = {}
        for name, r in idx.items():
            g = col.get(str(name))
            if g is None or cl.mu[g] < 0.5:
                continue
            before = np.expm1(cl.mu[g])
            after = np.expm1(cl.pert[r, g])
            if before > 1e-6:
                per[str(name)] = float(np.clip(after / before, 0.0, 2.0))
        out[cl.name] = per
    shared = sorted(set.intersection(*(set(v) for v in out.values())))
    return out, shared


def predict_promoter(client, chrom: str, tss: int, strand: str,
                     biosamples: list[str]) -> dict[str, float]:
    """Mean predicted RNA-seq and DNase over the promoter, per biosample."""
    from alphagenome.data import genome
    from alphagenome.models import dna_client

    half = WINDOW // 2
    start = max(tss - half, 0)
    interval = genome.Interval(chrom, start, start + WINDOW)
    out = client.predict_interval(
        interval=interval,
        organism=dna_client.Organism.HOMO_SAPIENS,
        requested_outputs=[dna_client.OutputType.RNA_SEQ,
                           dna_client.OutputType.DNASE],
        ontology_terms=list(ONTOLOGY.values()))
    centre = tss - start
    lo, hi = centre - PROMOTER, centre + PROMOTER
    res: dict[str, float] = {}
    for kind in ("rna_seq", "dnase"):
        track = getattr(out, kind, None)
        if track is None:
            continue
        md = track.metadata
        for bs in biosamples:
            sel = np.flatnonzero(md["biosample_name"].astype(str).values == bs)
            if sel.size == 0:
                continue
            res[f"{kind}:{bs}"] = float(
                np.mean(track.values[lo:hi, sel]))
    return res


def track_resolution(rows: list[dict]) -> None:
    """Does each cell line's track match *its own* measured profile best?

    The gate ``transferability-prior-eval`` says to run before building anything
    on a sequence model's cell-type head.  Every line shares one genome, so all
    cell-line specificity has to come from that head; if line X's predicted
    profile correlates better with line Y's measured expression than with X's
    own, the head does not resolve these lines and no downstream context feature
    built on it can either.
    """
    import sys
    sys.path.insert(0, "/root/.claude/skills/synced/transferability-prior-eval")
    from kernel import track_resolution_check                    # noqa: E402

    pred, obs = {}, {}
    for ln, bs in BIOSAMPLE.items():
        key = f"rna_seq:{bs}"
        keep = [r for r in rows
                if key in r["pred"] and np.isfinite(r["baseline"][ln])]
        pred[ln] = np.array([r["pred"][key] for r in keep])
        obs[ln] = np.array([r["baseline"][ln] for r in keep])
    n = min(len(v) for v in pred.values())
    pred = {k: v[:n] for k, v in pred.items()}
    obs = {k: v[:n] for k, v in obs.items()}

    res = track_resolution_check(pred, obs)
    print(f"\ntrack resolution gate ({n} genes) -- does each line's track match "
          f"its own profile?")
    for k, v in res.items():
        if isinstance(v, dict):
            print(f"  {k}: " + ", ".join(f"{a}->{b}" for a, b in v.items()))
        elif isinstance(v, float):
            print(f"  {k}: {v:.3f}")
        else:
            print(f"  {k}: {v}")
    Path("results/alphagenome_track_gate.json").write_text(
        json.dumps(res, indent=2, default=str))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gate-only", action="store_true",
                    help="just run the track-resolution gate on saved output")
    ap.add_argument("--n-genes", type=int, default=80)
    ap.add_argument("--tss", type=Path, default=Path("/home/user/vcc_data/tss.tsv"),
                    help="TSV: symbol, chrom, tss, strand (hg38)")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    if args.gate_only:
        rows = json.loads(OUT.read_text())
        track_resolution(rows)
        report(rows)
        return

    lines, _, symbols = load_all()
    residual, shared = measured_residual(lines, symbols)
    # Positive control: measured control expression per line, so the same
    # AlphaGenome numbers can be checked against a quantity they certainly do
    # encode.  Without it a null is indistinguishable from a broken pipeline.
    baseline = {}
    for cl in lines:
        if cl.name in BIOSAMPLE:
            col = {str(x): i for i, x in enumerate(cl.symbols)}
            baseline[cl.name] = {g: float(cl.mu[col[g]]) for g in shared
                                 if g in col}
    print(f"{len(shared):,} genes with a measured on-target residual in all of "
          f"{', '.join(BIOSAMPLE)}")

    if not args.tss.exists():
        raise SystemExit(f"need {args.tss}: symbol/chrom/tss/strand, hg38")
    tss = {}
    for line in args.tss.read_text().splitlines()[1:]:
        f = line.split("\t")
        if len(f) >= 4:
            tss[f[0]] = (f[1], int(f[2]), f[3])
    usable = [g for g in shared if g in tss]
    print(f"  {len(usable):,} of them have a TSS in {args.tss.name}")

    rng = np.random.default_rng(args.seed)
    pick = sorted(rng.choice(usable, min(args.n_genes, len(usable)),
                             replace=False))

    from alphagenome.models import dna_client
    client = dna_client.create(os.environ["ALPHAGENOME_API_KEY"])

    rows = []
    for i, g in enumerate(pick):
        chrom, pos, strand = tss[g]
        try:
            pred = predict_promoter(client, chrom, pos, strand,
                                    list(BIOSAMPLE.values()))
        except Exception as exc:                      # noqa: BLE001
            print(f"  {g}: {type(exc).__name__} {exc}")
            continue
        rows.append({"gene": g, "pred": pred,
                     "measured": {ln: residual[ln][g] for ln in BIOSAMPLE},
                     "baseline": {ln: baseline[ln].get(g, float("nan"))
                                  for ln in BIOSAMPLE}})
        if (i + 1) % 10 == 0:
            print(f"  {i + 1}/{len(pick)} genes", flush=True)

    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(rows, indent=2))
    print(f"\nwrote {OUT} ({len(rows)} genes)")
    if rows:
        track_resolution(rows)
        report(rows)


def report(rows: list[dict]) -> None:
    """Do predicted promoter differences track measured residual differences?"""
    from scipy.stats import spearmanr

    names = list(BIOSAMPLE)
    print("\ncontext pair        n   Spearman(pred ratio, measured ratio)")
    print("-" * 62)
    for kind in ("rna_seq", "dnase"):
        for a, b in ((names[i], names[j]) for i in range(len(names))
                     for j in range(i + 1, len(names))):
            pa, pb, ma, mb = [], [], [], []
            for r in rows:
                ka = f"{kind}:{BIOSAMPLE[a]}"
                kb = f"{kind}:{BIOSAMPLE[b]}"
                if ka not in r["pred"] or kb not in r["pred"]:
                    continue
                pa.append(r["pred"][ka]); pb.append(r["pred"][kb])
                ma.append(r["measured"][a]); mb.append(r["measured"][b])
            if len(pa) < 8:
                continue
            pr = np.log2((np.array(pa) + 1e-3) / (np.array(pb) + 1e-3))
            mr = np.log2((np.array(ma) + 1e-3) / (np.array(mb) + 1e-3))
            rho = spearmanr(pr, mr)
            print(f"{kind:>8} {a}/{b:<8}{len(pa):>4}   rho={rho.statistic:+.3f} "
                  f"p={rho.pvalue:.3f}")
    # Positive control, and the variation being predicted.
    if all("baseline" in r for r in rows):
        print("\npositive control -- the same predictions against measured "
              "control expression")
        print("-" * 62)
        for a, b in ((names[i], names[j]) for i in range(len(names))
                     for j in range(i + 1, len(names))):
            pa, pb, ba, bb = [], [], [], []
            for r in rows:
                ka, kb = f"rna_seq:{BIOSAMPLE[a]}", f"rna_seq:{BIOSAMPLE[b]}"
                if ka not in r["pred"] or kb not in r["pred"]:
                    continue
                if not np.isfinite([r["baseline"][a], r["baseline"][b]]).all():
                    continue
                pa.append(r["pred"][ka]); pb.append(r["pred"][kb])
                ba.append(r["baseline"][a]); bb.append(r["baseline"][b])
            if len(pa) < 8:
                continue
            pr = np.log2((np.array(pa) + 1e-3) / (np.array(pb) + 1e-3))
            br = np.array(ba) - np.array(bb)
            rho = spearmanr(pr, br)
            print(f"{'rna_seq':>8} {a}/{b:<8}{len(pa):>4}   rho={rho.statistic:+.3f} "
                  f"p={rho.pvalue:.4f}")

        spread = []
        for a, b in ((names[i], names[j]) for i in range(len(names))
                     for j in range(i + 1, len(names))):
            mr = [abs(np.log2((r["measured"][a] + 1e-3) /
                              (r["measured"][b] + 1e-3))) for r in rows]
            spread.append(float(np.median(mr)))
        print(f"\nmedian |log2 ratio| of the measured on-target residual "
              f"between contexts: {np.mean(spread):.2f}")
        print("-- i.e. how much real between-context variation there is to "
              "predict")

    print("\nA correlation indistinguishable from zero means promoter activity "
          "predicted\nfrom sequence carries no information about how completely "
          "CRISPRi silences a\ngene in one context versus another -- which is "
          "the only quantity here that is\nmissing from the official contexts "
          "and inside AlphaGenome's output space.")


if __name__ == "__main__":
    main()
