"""
bioif.real.heterogeneity -- how much does assay identity actually cost?

The contract insists that a measurement's assay is part of its context, and
that pIC50 / pKi / pKd are different quantities. That insistence is only
worth its cost if ignoring it changes numbers materially. This module
measures the effect on real data.

Run: python3 -m bioif.real.heterogeneity

✅ Computed from the committed ChEMBL snapshot (measured data, not model
output). The headline result is reproduced in INTEROP.md.
"""
from __future__ import annotations

import collections
import statistics

from .chembl import load_snapshot, snapshot_provenance

AFFINITY_TYPES = ("IC50", "Ki", "Kd", "EC50")


def _pct(sorted_vals, q):
    if not sorted_vals:
        return float("nan")
    return sorted_vals[min(int(q * len(sorted_vals)), len(sorted_vals) - 1)]


def analyse(rows=None) -> dict:
    rows = rows if rows is not None else load_snapshot()
    rows = [r for r in rows if r["standard_type"] in AFFINITY_TYPES
            and r["pchembl_value"]]

    by_type = collections.Counter(r["standard_type"] for r in rows)

    # (a) same compound, measured under more than one readout TYPE
    per_cmp_type = collections.defaultdict(lambda: collections.defaultdict(list))
    # (b) same compound, same readout type, more than one ASSAY
    per_cmp_assay = collections.defaultdict(lambda: collections.defaultdict(list))
    for r in rows:
        v = float(r["pchembl_value"])
        per_cmp_type[r["molecule_chembl_id"]][r["standard_type"]].append(v)
        per_cmp_assay[(r["molecule_chembl_id"], r["standard_type"])][
            r["assay_chembl_id"]].append(v)

    across_types = sorted(
        max(statistics.fmean(v) for v in d.values()) -
        min(statistics.fmean(v) for v in d.values())
        for d in per_cmp_type.values() if len(d) > 1)

    within_type = sorted(
        max(statistics.fmean(v) for v in d.values()) -
        min(statistics.fmean(v) for v in d.values())
        for d in per_cmp_assay.values() if len(d) > 1)

    return {
        "provenance": snapshot_provenance(),
        "n_records": len(rows),
        "n_compounds": len({r["molecule_chembl_id"] for r in rows}),
        "n_assays": len({r["assay_chembl_id"] for r in rows}),
        "by_type": dict(by_type),
        "across_readout_types": {
            "n_compounds": len(across_types),
            "median": statistics.median(across_types) if across_types else None,
            "p90": _pct(across_types, 0.9),
            "max": max(across_types) if across_types else None,
        },
        "within_type_across_assays": {
            "n_compound_type_pairs": len(within_type),
            "median": statistics.median(within_type) if within_type else None,
            "p90": _pct(within_type, 0.9),
            "max": max(within_type) if within_type else None,
        },
    }


def report(res: dict | None = None) -> str:
    r = res or analyse()
    p = r["provenance"]
    L = [
        "ChEMBL affinity heterogeneity for a single target",
        f"  source   : {p.get('source', '?')}  target {p.get('target', '?')}",
        f"  fetched  : {p.get('fetched_utc', '?')}",
        f"  records  : {r['n_records']} pChEMBL-bearing affinity measurements",
        f"  compounds: {r['n_compounds']}   distinct assays: {r['n_assays']}",
        f"  readouts : {r['by_type']}",
        "",
        "Same compound, DIFFERENT readout type (IC50 vs Ki vs Kd vs EC50):",
    ]
    a = r["across_readout_types"]
    L.append(f"  n={a['n_compounds']} compounds   spread in pChEMBL  "
             f"median {a['median']:.2f}  p90 {a['p90']:.2f}  max {a['max']:.2f} log units")
    L += ["", "Same compound, SAME readout type, DIFFERENT assay:"]
    w = r["within_type_across_assays"]
    L.append(f"  n={w['n_compound_type_pairs']} compound-type pairs   spread  "
             f"median {w['median']:.2f}  p90 {w['p90']:.2f}  max {w['max']:.2f} log units")
    L += [
        "",
        "Reading:",
        "  The rule everyone knows -- do not mix IC50 with Ki -- is not the",
        f"  dominant term here ({a['median']:.2f} log median). Holding the readout type",
        f"  FIXED and only changing assay is worse ({w['median']:.2f} log median, i.e. a",
        f"  {10 ** w['median']:.0f}x difference in apparent potency for the same compound on",
        "  the same target). So assay identity has to be part of the context,",
        "  not a column you drop on the way into a model.",
    ]
    return "\n".join(L)


if __name__ == "__main__":
    print(report())
