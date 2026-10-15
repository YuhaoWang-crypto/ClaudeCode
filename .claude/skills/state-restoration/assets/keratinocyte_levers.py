"""
Template: rank levers for a SKIN indication with the default keratinocyte
cytokine library (step 5 of the pipeline).

The library is built by `netjac.data.keratinocyte.lever_library()` and
contains nine measured cytokine responses in primary keratinocytes:

    IFN-alpha  IFN-gamma  IL-4  IL-13  IL-17A  TNF      (GSE36287, 3 donors)
    IL-9       IL-22      IL-26                         (GSE305982, 1 donor)

IL-4 is in the default on purpose: with IL-13 as the only type-2 lever the
atopic-dermatitis question cannot be answered (IL-13 crossed specificity
+0.004 n.s. vs IL-4 +0.103 [0.074, 0.129]).

Run:  python3 keratinocyte_levers.py          (from the repository root)
Check the library itself:  python3 -m netjac.data.keratinocyte
"""
import sys

import numpy as np

sys.path.insert(0, ".")
from netjac.data import keratinocyte as kc          # noqa: E402
from netjac.demos import d17_keratinocyte_lever as d17   # noqa: E402

# ---------------------------------------------------------------- CONFIG
# Disease displacement(s): indication -> {"genes": array, "d": array}
# d = cohort-mean (lesional - own non-lesional) on a log scale.
# The default uses the two cohorts already wired up; replace with your own.
DISEASE = None            # None -> d17.disease_vectors() (psoriasis + atopic dermatitis)
# Which levers are approved drug targets, for the crossed check.
VALIDATED = kc.VALIDATED  # {"psoriasis": {"IL-17A", "TNF"}, ...}
TIME = "24 hours"         # only affects the GSE305982 columns


if __name__ == "__main__":
    lib = kc.lever_library(time=TIME)
    dis = DISEASE or d17.disease_vectors()
    print(f"library: {len(lib['names'])} levers, {len(lib['genes'])} genes "
          f"({', '.join(f'{n} [{s}]' for n, s in zip(lib['names'], lib['source']))})")

    for ind, dv in dis.items():
        r = d17.rank(dv["d"], dv["genes"], lib["R"], lib["genes"], lib["names"])
        val = VALIDATED.get(ind, set())
        print(f"\n{ind}: top lever {r['top']}"
              + (f"  [{'HIT' if r['top'] in val else 'miss'}]" if val else ""))
        for j in np.argsort(-r["energy"]):
            n = r["names"][j]
            # coef > 0 means the lesion lies ALONG this input, so restoring the
            # state means removing it - what a blocking antibody does.
            action = "block" if r["coef"][j] > 0 else "supply"
            mark = " <- validated target" if n in val else ""
            print(f"   {n:10s} cos {r['cos'][j]:+.3f}  coef {r['coef'][j]:+.3f}"
                  f"  p={r['p'][j]:.4f}  -> {action}{mark}")

    if len(dis) >= 2:
        print("\ncrossed specificity (more aligned with its own indication than the other):")
        for ind, levers in VALIDATED.items():
            if ind not in dis:
                continue
            for lev in sorted(levers):
                if lev not in lib["names"]:
                    print(f"   {lev}: NOT IN LIBRARY - add it before concluding anything")
                    continue
                d, lo, hi, p = kc.crossed_specificity(lib, dis, ind, lev)
                print(f"   {lev:8s} for {ind:19s} {d:+.3f} [{lo:+.3f}, {hi:+.3f}] "
                      f"p={p:.4f} [{'OK' if lo > 0 else 'n.s.'}]")
