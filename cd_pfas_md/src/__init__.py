"""β-cyclodextrin / PFAS host-guest binding free energy pipeline.

Two deliverables:
  * APR (attach-pull-release) absolute binding free energy  -> calibration
    against a measured dye·CD Ka   (parameterize -> build_apr -> run_apr -> analyze_apr)
  * FEP/TI relative binding free energy ΔΔG for CD modifications
    (fep_ti_ddg)

Everything is config-driven from ../config/*.yaml.
"""

__all__ = [
    "utils",
    "parameterize",
    "build_apr",
    "run_apr",
    "analyze_apr",
    "fep_ti_ddg",
]
