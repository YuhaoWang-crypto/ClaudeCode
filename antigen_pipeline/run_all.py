"""Run the whole cell-surface antigen discovery pipeline, step 1 to step 9.

    python3 -m antigen_pipeline.run_all

Each step writes its results and a provenance sidecar to results/ before the
next one starts, so any step can also be re-run on its own against the files
already on disk.
"""

from __future__ import annotations

import time

from . import (s1_surfaceome, s2_census, s3_topology, s4_druggability,
               s5_safety, s6_score, s7_validate, s8_literature,
               s9_export, s9_figures, s9_report)


def main() -> None:
    t0 = time.time()
    surfaceome = s1_surfaceome.run()
    s2_census.run(surfaceome)
    accessible = s3_topology.run(surfaceome)
    annotated = s4_druggability.run(accessible)
    safe = s5_safety.run(annotated)
    ranked = s6_score.run(safe)
    s7_validate.run(ranked)
    s8_literature.run(ranked)
    s9_figures.run()
    s9_export.run()
    s9_report.run()
    print(f"\npipeline complete in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
