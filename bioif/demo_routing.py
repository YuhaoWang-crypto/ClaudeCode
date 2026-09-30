"""
bioif.demo_routing -- two models for one hop, and letting the registry choose.

Run:  python3 -m bioif.demo_routing

Until now every hop had exactly one implementation, so "which model" was
never a question the interface had to answer. It is the question that matters
most in practice. Here the assay-transfer hop carries up to twelve competing
adapters -- four source assays x three predictors (identity / shift / linear),
all fitted on real ChEMBL pairs, all wrapped in the same split-conformal
procedure -- and the registry picks per claim on declared properties:

    in-domain first, then narrowest calibrated interval, then kind.

Two results come out of it, and the second is the one that matters.

⚠️  Downstream of the transfer, the occupancy/activity/fitness adapters are
still illustrative stubs, so chain endpoints are hypotheses. The transfer
hop's numbers are measured.
"""
from __future__ import annotations

import random
import statistics

from .adapters_demo import (ActivityToFitness, OccupancyToResidualActivity,
                            ProteinToResidualActivity)
from .core import (Claim, Context, Entity, Estimate, Provenance, MEASURED,
                   REFUSE)
from .ensemble import compare_routes
from .quantities import (FITNESS, PIC50, PIC50_REF, REFERENCE_ASSAY,
                         STRICT_LOSSLESS)
from .real import conformal as C, transfer_models as T
from .real.sources import KRAS
from .real.transfer_adapter import build_all_transfers
from .registry import Registry

ALPHA = 0.2
LEVEL = 1.0 - ALPHA


def strict_registry():
    """
    A registry with no edge from a raw pIC50 straight into the downstream
    model. The only way across is an explicit, calibrated assay transfer --
    so 'which assay is this number from?' becomes a routing question rather
    than a footnote.
    """
    reg = Registry(include_lossless=False)
    for c in STRICT_LOSSLESS:
        reg.register(c)
    occ = OccupancyToResidualActivity()
    occ.target = KRAS
    for a in (occ, ProteinToResidualActivity(), ActivityToFitness()):
        reg.register(a)
    for t in build_all_transfers(alpha=ALPHA):
        reg.register(t)
    return reg


def source_claim(assay: str, x: float) -> Claim:
    return Claim(entity=KRAS, quantity=PIC50,
                 context=Context(assay=assay, dose_uM=0.30,
                                 system="depmap:ACH-000019", time_h=120.0),
                 estimate=Estimate.point(x, 4000), evidence=MEASURED,
                 provenance=(Provenance("ChEMBL", "REST/34", assay),))


def hr(t):
    print("\n" + "=" * 78 + f"\n{t}\n" + "=" * 78)


def main() -> None:
    pts = C.build_points()
    reg = strict_registry()

    # -- R1 what is competing ---------------------------------------------
    hr("R1  twelve adapters for one hop")
    cands = reg.candidates(PIC50, PIC50_REF)
    print(f"  {len(cands)} registered edges pIC50 -> pIC50(ref {REFERENCE_ASSAY})")
    print(f"\n  {'adapter':<40}{'width@0.2':>11}   model")
    for a in sorted(cands, key=lambda a: (a.pair.source, a.model.predictor.name)):
        print(f"  {a.name:<40}{a.calibrated_width(ALPHA):>11.2f}   "
              f"{a.model.predictor.describe()}")
    print("\n  Four source assays, three predictors each. Only one source is\n"
          "  ever in domain for a given claim; among that source's three, all\n"
          "  are conformally valid, so width is the only axis left to pick on.")
    print("\n  Look at the identity and shift rows for each source: identical\n"
          "  widths. That is not a bug -- see R2.")

    # -- R2 all three are valid; they differ only in width ----------------
    hr(f"R2  validity does not separate the models (alpha={ALPHA}, 200 splits)")
    cmp = T.compare_predictors(pts, repeats=200, alpha=ALPHA)
    agg = {k: {"cov": [], "w": []} for k in T.PREDICTORS}
    for d in cmp.values():
        for k in T.PREDICTORS:
            if d[k]["cov"] == d[k]["cov"]:
                agg[k]["cov"].append(d[k]["cov"])
                agg[k]["w"].append(d[k]["w"])
    print(f"  {'predictor':<12}{'mean coverage':>15}{'mean width':>13}"
          f"{'pairs valid':>14}")
    for k in T.PREDICTORS:
        n_ok = sum(1 for v in agg[k]["cov"] if v >= LEVEL - 0.02)
        print(f"  {k:<12}{statistics.fmean(agg[k]['cov']):>15.3f}"
              f"{statistics.fmean(agg[k]['w']):>13.2f}"
              f"{n_ok:>10}/{len(agg[k]['cov'])}")
    big = max(cmp.items(), key=lambda kv: kv[1]["n"])
    print(f"\n  On the well-powered pair ({big[0]}, n={big[1]['n']}):")
    for k in T.PREDICTORS:
        print(f"    {k:<10} coverage {big[1][k]['cov']:.3f}   "
              f"width {big[1][k]['w']:.2f}")
    same = abs(statistics.fmean(agg["identity"]["w"]) -
               statistics.fmean(agg["shift"]["w"])) < 0.02
    print(f"\n  identity and shift are the SAME conformalised model: {same}")
    print("  A conformalised adapter emits predict(x) + a resampled residual,\n"
          "  and the residuals carry any bias in predict with the opposite\n"
          "  sign. So 'reuse the number as measured' and 'reuse it with the\n"
          "  median offset applied' are one model once calibrated. Three\n"
          "  candidate predictors, two distinct models.")
    print("\n  ⚠️  An earlier version of this demo scored them on a symmetric\n"
          "  interval around predict(x) while emitting the de-biased\n"
          "  distribution. That made identity look 4x wider than shift and\n"
          "  inflated the model-disagreement figure in R5 by roughly 6x. The\n"
          "  demo caught it: identity and shift were returning byte-identical\n"
          "  endpoints while advertising very different widths.")
    print("\n  What remains is a real contest between shift and linear, and\n"
          "  linear wins only where there is data to fit its extra parameter:\n"
          f"  narrower on the n={big[1]['n']} pair, wider on all six small ones.")

    # -- R3 the registry picking, per claim -------------------------------
    hr("R3  the pick is per claim, and it changes with the claim")
    sources = sorted({a.pair.source for a in cands})
    print(f"  {'claim from assay':<18}{'chosen adapter':<40}{'width':>8}")
    for s in sources:
        claim = source_claim(s, 6.0)
        pick = reg.select(PIC50, PIC50_REF, claim, alpha=ALPHA)
        print(f"  {s:<18}{pick.name:<40}"
              f"{pick.calibrated_width(ALPHA):>8.2f}")
    unknown = source_claim("CHEMBL9999999", 6.0)
    print(f"  {'CHEMBL9999999':<18}{str(reg.select(PIC50, PIC50_REF, unknown, ALPHA)):<40}")
    print(f"\n  The domain test alone eliminates {len(cands) - 3} of {len(cands)}; "
          f"width settles the\n  rest, with ties broken by name (identity and "
          f"shift being one model,\n  a tie between them is real). No source "
          f"assay at all -> no route,\n  not a default.")

    counts = {k: sum(1 for s in sources
                     if reg.select(PIC50, PIC50_REF, source_claim(s, 6.0),
                                   ALPHA).model.predictor.name == k)
              for k in T.PREDICTORS}
    print(f"\n  winner by predictor across the four sources: {counts}")
    print("  -- no model wins everywhere, which is why this is a routing\n"
          "     decision and not a library-wide default.")

    # -- R4 end to end, and the disagreement ------------------------------
    hr("R4  run every in-domain route, not just the chosen one")
    print("  Two source assays, swept across the input range. The flag fires\n"
          "  when switching model moves the answer further than the selected\n"
          "  model's own interval admits.\n")
    print(f"  {'source assay':<16}{'n':>4}{'pIC50 in':>10}{'routes':>8}"
          f"{'endpoint':>10}{'sel.width':>11}{'spread':>9}{'flag':>6}")
    for src_assay in ("CHEMBL5737243", "CHEMBL4368373"):
        npair = next((a.pair.n for a in cands if a.pair.source == src_assay), 0)
        for x in (4.5, 6.0, 7.0, 7.5):
            claim = source_claim(src_assay, x)
            rc = compare_routes(reg, PIC50, FITNESS, claim, alpha=ALPHA)
            if not rc.ok:
                continue
            sel = rc.results[rc.selected]
            print(f"  {src_assay:<16}{npair:>4}{x:>10.2f}{len(rc.results):>8}"
                  f"{sel.final.estimate.mean:>+10.3f}"
                  f"{rc.selected_width(LEVEL):>11.3f}{rc.point_spread():>9.3f}"
                  f"{'YES' if rc.choice_exceeds_interval(LEVEL) else '-':>6}")
    print("\n  The well-calibrated source (n=226) never trips the flag; the\n"
          "  thin one (n=13) trips it across most of the range. Model choice\n"
          "  bites exactly where there is not enough data to tell the models\n"
          "  apart -- which is also where you least want to be guessing.")

    src_assay = "CHEMBL4368373"
    claim = source_claim(src_assay, 7.5)
    rc = compare_routes(reg, PIC50, FITNESS, claim, alpha=ALPHA)
    print(f"\n  at pIC50 7.50 from {src_assay}, route by route:")
    for label, res in sorted(rc.results.items(),
                             key=lambda kv: kv[1].final.estimate.mean):
        lo, hi = res.final.estimate.ci(LEVEL)
        mark = "  <- selected" if label == rc.selected else ""
        print(f"    {label:<40}{res.final.estimate.mean:>+9.3f}   "
              f"[{lo:+.3f}, {hi:+.3f}]{mark}")
    print(f"    routes refused by the domain test: {len(rc.refused)}")
    env = rc.envelope(LEVEL)
    print(f"    envelope across all routes: [{env[0]:+.3f}, {env[1]:+.3f}]"
          f"  (width {env[1] - env[0]:.3f} vs selected "
          f"{rc.selected_width(LEVEL):.3f})")
    f = rc.flag(LEVEL)
    print(f"\n  {f if f else 'model choice fits inside the selected interval here'}")

    # -- R5 how often does that happen? -----------------------------------
    dis = T.disagreement_vs_width(pts, alpha=ALPHA)
    hr(f"R5  measured on the transfer hop itself, across {len(dis)} pairs")
    print(f"  {'pair':<18}{'n':>4}{'selected':>10}{'width':>8}"
          f"{'med spread':>12}{'inside?':>9}")
    for k, v in sorted(dis.items(), key=lambda kv: -kv[1]["n"]):
        print(f"  {k:<18}{v['n']:>4}{v['selected']:>10}{v['width']:>8.2f}"
              f"{v['median_spread']:>12.2f}{v['frac_contained']:>9.2f}")
    fc = [v["frac_contained"] for v in dis.values()]
    ratio = statistics.median([v["median_spread"] / v["width"]
                               for v in dis.values() if v["width"] > 0])
    maxr = max(v["max_spread"] / v["width"] for v in dis.values())
    print(f"\n  pairs where some compound's rival prediction falls OUTSIDE the\n"
          f"  selected model's {LEVEL:.0%} interval: "
          f"{sum(1 for v in fc if v < 1.0)}/{len(fc)}")
    print(f"  median (model spread / selected interval width): {ratio:.2f}")
    print(f"  worst  (max spread / width) over the pairs:      {maxr:.2f}")
    print(f"\n  So: model choice is a SECOND-ORDER term here. Typically it moves\n"
          f"  the answer about {ratio:.0%} of the calibrated interval width -- much\n"
          f"  less than the interval itself, and far less than the cross-assay\n"
          f"  failure in INTEROP §6.3, where coverage collapsed to 0.00.")
    print(f"\n  It is not negligible at the edges: on {max(dis, key=lambda k: dis[k]['max_spread'] / dis[k]['width'])} the worst\n"
          f"  compound sees the models differ by more than the whole interval.\n"
          f"  That is what the R4 flag is for -- it fires per claim, not per\n"
          f"  pair, because whether the choice matters depends on where in the\n"
          f"  input range you are asking.")
    print("\n  A conformal interval covers the residuals of the model you\n"
          "  picked; it is silent on whether picking it was right. Running the\n"
          "  alternatives is how you price that silence, and it costs compute\n"
          "  rather than an experiment.")

    print("\n" + "=" * 78)
    print("✅ transfer-hop coverage, widths and disagreement are measured on\n"
          "   held-out compounds from the committed ChEMBL snapshot.\n"
          "⚠️  downstream adapters remain stubs; chain endpoints are\n"
          "   hypotheses, and the spread between them is the result.")
    print("=" * 78)


if __name__ == "__main__":
    main()
