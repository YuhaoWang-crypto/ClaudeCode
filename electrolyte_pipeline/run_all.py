"""
Run every analysis module on whatever trajectories are present and print one
consolidated, honesty-labelled summary.

    python -m electrolyte_pipeline.run_all            # analyses only (needs elyte_work/)
    python -m electrolyte_pipeline.run_all --qc       # also (re)run the CPU quantum-chemistry parts

Trajectories are produced by
    modal run electrolyte_pipeline/modal_run.py                      # E1 bulk series
    modal run electrolyte_pipeline/modal_run.py --stage nemd         # E4 viscosity
    modal run electrolyte_pipeline/modal_run.py --stage interface    # E7
    modal volume get elyte-work / elyte_work/
"""
from __future__ import annotations

import json
import os
import sys

from electrolyte_pipeline.e0_systems import WORK, FIG
from electrolyte_pipeline import e2_rdf_cn, e3_clusters, e4_properties, e5_qc_clusters, e6_desolvation, \
    e7_interface, e8_ml, e9_reactive


def main(argv=None):
    argv = argv if argv is not None else sys.argv[1:]
    print("#" * 72 + "\n#  ELECTROLYTE COMPUTATION PIPELINE — Yao et al. Chem. Rev. 2022 digest, reproduced\n" + "#" * 72 + "\n")
    summary = {}
    for name, fn in (("E2", e2_rdf_cn.report), ("E3", e3_clusters.report), ("E4", e4_properties.report)):
        try:
            summary[name] = fn(); print()
        except Exception as e:
            print(f"{name} skipped ({type(e).__name__}: {e})\n")
    for name, fn in (("E5", e5_qc_clusters.report), ("E6", e6_desolvation.report)):
        try:
            summary[name] = fn(force="--qc" in argv); print()
        except Exception as e:
            print(f"{name} skipped ({type(e).__name__}: {e})\n")
    for name, fn in (("E7", e7_interface.report), ("E8", e8_ml.report), ("E9", e9_reactive.report)):
        try:
            summary[name] = fn(); print()
        except Exception as e:
            print(f"{name} skipped ({type(e).__name__}: {e})\n")
    os.makedirs(WORK, exist_ok=True)
    json.dump(summary, open(os.path.join(WORK, "summary.json"), "w"), default=str)
    print(f"figures -> {FIG}")


if __name__ == "__main__":
    main()
