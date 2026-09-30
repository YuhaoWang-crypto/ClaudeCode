"""
bioif.demo_conformal -- replacing an asserted error bar with a checked one.

Run:  python3 -m bioif.demo_conformal

Every interval elsewhere in this package is asserted: a stub declares sd=0.55
and the chain believes it. Here one adapter's interval is *calibrated* by
split conformal prediction and its coverage is *measured* on held-out data.

The hop is assay transfer -- may a potency measured in assay A be reused
where assay B is needed -- on real ChEMBL IC50 records for KRAS. All splits
are on compound, never on row.

⚠️  Coverage numbers below are computed from the committed snapshot. They are
about this target and these assays; they are not a general claim about
conformal prediction or about ChEMBL.
"""
from __future__ import annotations

import random
import statistics

from .chain import run_chain
from .core import Context
from .quantities import FITNESS, PIC50
from .real import conformal as C
from .real.transfer_adapter import build_transfer
from .demo_real import build_registry

ALPHA = 0.2
NOMINAL = 1.0 - ALPHA


def hr(t):
    print("\n" + "=" * 76 + f"\n{t}\n" + "=" * 76)


def main() -> None:
    pairs = C.discover_pairs()
    cont = {k: v for k, v in pairs.items() if v.status == "continuous"}
    binned = {k: v for k, v in pairs.items() if v.status == "binned"}
    pts = C.build_points(pairs)

    # -- C0 the gate that had to exist before any of this was meaningful ---
    hr("C0  not every float column is a measurement")
    print(f"  candidate assay pairs with >=10 shared compounds: {len(pairs)}")
    print(f"    continuous: {len(cont)}      refused as binned: {len(binned)}")
    worst = max(binned.values(), key=lambda p: p.tie_x)
    print(f"\n  worst offender: {worst.key}  n={worst.n}")
    print(f"    {worst.reason}")
    print(f"    source levels: {worst.n_distinct_x}, target levels: "
          f"{worst.n_distinct_y}")
    print("\n  These are patent potency bands at one-log spacing arriving as a\n"
          "  float at the band midpoint. Fitting a regression and a conformal\n"
          "  interval to them is a category error, and nothing in the column\n"
          "  types says so -- only the tie fraction does.")

    # -- C1 does the interval cover, inside one assay pair? ---------------
    hr(f"C1  within-pair coverage (nominal {NOMINAL:.2f}), 300 compound splits")
    e1 = C.experiment_within_pair(pts, repeats=300)
    print(f"  {'assay pair':<18}{'n':>4}   {'conformal':>18}   {'gaussian +-1.28sd':>20}")
    print(f"  {'':<18}{'':>4}   {'cov':>8}{'width':>10}   {'cov':>9}{'width':>11}")
    cc, gg = [], []
    for name, d in sorted(e1.items(), key=lambda kv: -kv[1]["n"]):
        a = d[ALPHA]
        cc.append(a["cov"])
        gg.append(a["gcov"])
        print(f"  {name:<18}{d['n']:>4}   {a['cov']:>8.3f}{a['w']:>10.2f}   "
              f"{a['gcov']:>9.3f}{a['gw']:>11.2f}")
    print(f"\n  conformal: mean {statistics.fmean(cc):.3f}, "
          f"{sum(1 for v in cc if v >= NOMINAL - 0.02)}/{len(cc)} pairs at or "
          f"above nominal")
    print(f"  gaussian : mean {statistics.fmean(gg):.3f}, "
          f"{sum(1 for v in gg if v >= NOMINAL - 0.02)}/{len(gg)} pairs at or "
          f"above nominal")
    print("\n  The Gaussian interval is narrower and wrong. Conformal pays for\n"
          "  the guarantee in width -- that is the trade, stated honestly.")

    big = max(e1.items(), key=lambda kv: kv[1]["n"])
    print(f"\n  On the one well-powered pair ({big[0]}, n={big[1]['n']}):")
    for a in (0.2, 0.1):
        print(f"    alpha={a}: coverage {big[1][a]['cov']:.3f} "
              f"(nominal {1 - a:.2f}), width {big[1][a]['w']:.2f}")

    # -- C2 what a small calibration set cannot buy -----------------------
    hr("C2  when the guarantee is not available, the method says so")
    unattainable = [(n, d["n"]) for n, d in e1.items()
                    if d[0.1]["unattainable"] > 0.5]
    print(f"  at alpha=0.1, {len(unattainable)}/{len(e1)} pairs cannot support "
          f"the interval at all:")
    print("    ceil((n_cal+1)*0.9) > n_cal, so the conformal quantile is "
          "+infinity.")
    print(f"    e.g. {', '.join(f'{n} (n={k})' for n, k in unattainable[:4])} ...")
    print("\n  It returns an infinite interval rather than the widest observed\n"
          "  residual. A method that quietly returned the latter would be\n"
          "  claiming a guarantee it does not have -- the same failure as an\n"
          "  adapter that cannot refuse.")

    # -- C3 the guarantee does not transfer across assay pairs ------------
    hr(f"C3  calibrate on OTHER assay pairs, test on a held-out one")
    e2 = C.experiment_cross_pair(pts, repeats=100)
    covs = [d[ALPHA] for d in e2.values() if d[ALPHA] == d[ALPHA]]
    print(f"  {'held-out pair':<18}{'n':>4}{'coverage':>10}")
    for name, d in sorted(e2.items(), key=lambda kv: kv[1][ALPHA]):
        mark = "  <-- " if d[ALPHA] < 0.5 else ""
        print(f"  {name:<18}{d['n']:>4}{d[ALPHA]:>10.3f}{mark}")
    print(f"\n  mean coverage {statistics.fmean(covs):.3f} against a nominal "
          f"{NOMINAL:.2f}")
    print(f"  {sum(1 for v in covs if v < NOMINAL)}/{len(covs)} pairs "
          f"under-covered; worst {min(covs):.3f}")
    print("\n  Conformal guarantees coverage under exchangeability. Two assay\n"
          "  pairs are not exchangeable -- their systematic shifts differ by\n"
          "  up to 3 log units -- so the guarantee is void, and the failure is\n"
          "  total, not graceful. This is why the adapter REFUSES an\n"
          "  uncalibrated pair instead of widening for it.")

    # -- C4 marginal coverage can be met while every class fails ----------
    hr("C4  marginal coverage is met; conditional coverage is not")
    e3 = C.experiment_mondrian(pts, alpha=ALPHA, repeats=300)
    marg = e3.pop("_pooled_marginal")
    pl = [d["pooled"] for d in e3.values()]
    mo = [d["mondrian"] for d in e3.values()]
    print(f"  pooled calibration, MARGINAL coverage over all test points: "
          f"{marg:.3f}   (nominal {NOMINAL:.2f})")
    print(f"  ... and yet, per assay pair:\n")
    print(f"  {'pair':<18}{'pooled':>9}{'mondrian':>11}")
    for name, d in sorted(e3.items(), key=lambda kv: kv[1]["pooled"]):
        print(f"  {name:<18}{d['pooled']:>9.3f}{d['mondrian']:>11.3f}")
    print(f"\n  pooled   : {sum(1 for v in pl if v < NOMINAL)}/{len(pl)} pairs "
          f"under-covered, worst {min(pl):.3f}")
    print(f"  mondrian : {sum(1 for v in mo if v < NOMINAL)}/{len(mo)} pairs "
          f"under-covered, worst {min(mo):.3f}")
    print("\n  The pooled interval keeps its promise on average and is useless\n"
          "  for most individual pairs. The class that repairs it is the assay\n"
          "  -- the context field the contract already carries. Calibrate per\n"
          "  class and conditional coverage comes back.")

    # -- C5 the adapter, in the chain -------------------------------------
    hr("C5  the calibrated adapter, used in a chain")
    key = big[0]
    tx = build_transfer(key, alpha=ALPHA)
    info = tx.pair
    print(f"  {tx.name}   [{tx.kind}]  {tx.version}")
    print(f"    source {info.source}  ->  target {info.target}")
    print(f"    fitted shift {tx.cal.shift:+.3f} log units  "
          f"(n_train={tx.cal.n_train}, n_cal={tx.cal.n_cal})")
    print(f"    calibrated source range [{tx.cal.x_range[0]:.2f}, "
          f"{tx.cal.x_range[1]:.2f}]")

    x = 6.00
    lo, hi = tx.conformal_interval(x)
    print(f"\n  a compound reading pIC50 {x:.2f} in {info.source[-6:]}:")
    print(f"    guaranteed {NOMINAL:.0%} conformal interval in "
          f"{info.target[-6:]}: [{lo:.2f}, {hi:.2f}]")

    from .core import Claim, Entity, Estimate, Provenance, MEASURED
    from .real.sources import KRAS
    src = Claim(entity=KRAS, quantity=PIC50,
                context=Context(assay=info.source, dose_uM=0.30,
                                system="depmap:ACH-000019", time_h=120.0),
                estimate=Estimate.point(x, 4000), evidence=MEASURED,
                provenance=(Provenance("ChEMBL", "REST/34", info.source),))
    out = tx.apply(src, random.Random(3))
    s_lo, s_hi = out.estimate.ci(NOMINAL)
    print(f"    sampled  {NOMINAL:.0%} interval from the same adapter:        "
          f"[{s_lo:.2f}, {s_hi:.2f}]")
    print(f"    evidence: {src.evidence} -> {out.evidence}")
    print(f"    context.assay: {src.context.assay} -> {out.context.assay}")

    print("\n  and the refusal, offered a claim from an assay it never saw:")
    wrong = src.derive(context=src.context.with_(assay="CHEMBL9999999"))
    v = tx.domain(wrong)
    print(f"    {v.status.upper()}: {v.reason}")

    reg = build_registry()
    path = [tx] + reg.route(PIC50, FITNESS)
    res = run_chain(path, src, seed=7)
    print(f"\n  full chain: " + " -> ".join(a.name for a in path))
    print(f"    endpoint {res.final.estimate.summary()}  [{res.final.evidence}]")
    print(f"    transfer link owns {res.traces[0].var_share * 100:.1f}% of the "
          f"endpoint variance,")
    print(f"    and it is the only link in this chain whose share is backed by\n"
          f"    a measured coverage number rather than an assumed sd.")

    print("\n" + "=" * 76)
    print("✅ coverage figures are measured on held-out compounds from the\n"
          "   committed ChEMBL snapshot.\n"
          "⚠️  the downstream occupancy/activity/fitness adapters are still\n"
          "   stubs, so the endpoint remains a hypothesis, not a prediction.")
    print("=" * 76)


if __name__ == "__main__":
    main()
