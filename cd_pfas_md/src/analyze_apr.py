"""Analyze APR windows -> absolute binding free energy, and calibrate vs experiment.

Uses pAPRika's fe_calc (TI over the attach leg + MBAR over the pull leg, plus the
analytic release / reference-state correction) to produce:

    ΔG_bind = -(ΔG_attach + ΔG_pull + ΔG_release)      (kcal/mol)

then compares against the experimental ΔG derived from the measured dye·CD Ka in
system.yaml. The calibration report tells you whether your force field + protocol
reproduce the known number before you trust it on new modifications.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from . import utils

log = utils.get_logger("cd_pfas.analyze_apr")


def compute_binding_free_energy(work: Path, cfg: dict) -> dict:
    """Run pAPRika fe_calc over the finished windows. Returns a result dict."""
    manifest = json.loads((work / "apr_manifest.json").read_text())
    try:
        from paprika.analysis import fe_calc
    except ImportError as exc:
        raise ImportError("paprika required for analyze_apr") from exc

    T = manifest["temperature_K"]
    fecalc = fe_calc()
    fecalc.temperature = T
    fecalc.prmtop = manifest["prmtop"]
    fecalc.trajectory = "prod.dcd"
    fecalc.path = str(work / "windows")
    fecalc.restraint_list = []          # populated from the pAPRika restraint objects
    fecalc.methods = ["ti-block", "mbar-block"]
    fecalc.boot_cycles = 1000
    fecalc.compute_free_energy()
    fecalc.compute_ref_state_work([])   # analytic release correction

    # pAPRika stores per-leg results in fecalc.results; we summarize defensively.
    res = getattr(fecalc, "results", {})

    def leg(name, method="mbar-block"):
        try:
            return (res[name][method]["fe"], res[name][method]["sem"])
        except Exception:  # noqa: BLE001
            return (float("nan"), float("nan"))

    a_fe, a_sem = leg("attach")
    p_fe, p_sem = leg("pull")
    r_fe = res.get("ref_state_work", float("nan"))

    dG_bind = -(a_fe + p_fe + r_fe)
    sem = (a_sem ** 2 + p_sem ** 2) ** 0.5 if a_sem == a_sem else float("nan")
    return {
        "attach_fe": a_fe, "attach_sem": a_sem,
        "pull_fe": p_fe, "pull_sem": p_sem,
        "ref_state_work": r_fe,
        "dG_bind_kcal": dG_bind, "dG_bind_sem": sem,
        "temperature_K": T,
    }


def calibrate(result: dict, guest: utils.Guest) -> dict:
    """Attach an experimental comparison if the guest has a measured value."""
    out = dict(result)
    if guest.exp_dg is not None:
        err = result["dG_bind_kcal"] - guest.exp_dg
        out["exp_dG_kcal"] = guest.exp_dg
        out["error_kcal"] = err
        out["calc_Ka_M_inv"] = utils.dg_to_ka(result["dG_bind_kcal"], result["temperature_K"])
        out["exp_Ka_M_inv"] = utils.dg_to_ka(guest.exp_dg, result["temperature_K"])
        verdict = "GOOD (<1 kcal/mol)" if abs(err) < 1.0 else (
            "FAIR (1-2 kcal/mol)" if abs(err) < 2.0 else "POOR (>2 kcal/mol) — revisit FF/protocol")
        out["calibration_verdict"] = verdict
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description="Analyze APR + calibrate vs experiment.")
    ap.add_argument("--config", default="config/system.yaml")
    ap.add_argument("--work", required=True, help="apr_<guest> dir")
    ap.add_argument("--guest", required=True)
    args = ap.parse_args()

    cfg = utils.load_config(args.config)
    work = utils.resolve(args.work)
    guest = utils.Guest.from_config(args.guest, cfg["guests"][args.guest], cfg["thermo"]["temperature_K"])

    result = compute_binding_free_energy(work, cfg)
    report = calibrate(result, guest)
    (work / "binding_result.json").write_text(json.dumps(report, indent=2))

    log.info("=== %s binding free energy ===", guest.name)
    log.info("ΔG_bind = %.2f ± %.2f kcal/mol", report["dG_bind_kcal"], report.get("dG_bind_sem", float("nan")))
    if "exp_dG_kcal" in report:
        log.info("experiment = %.2f kcal/mol | error = %+.2f | %s",
                 report["exp_dG_kcal"], report["error_kcal"], report["calibration_verdict"])
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
