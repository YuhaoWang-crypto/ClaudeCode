"""Heuristic binding + displacement prefilter  (deliverables b + c, CPU-only).

WHAT THIS IS
------------
A cheap, transparent, physics-motivated scoring function to TRIAGE a large CD
modification library BEFORE spending FEP/APR GPU cycles, and to combine the
per-guest numbers into a single "displacement favorability" figure of merit.

WHAT THIS IS NOT
----------------
It is NOT converged molecular dynamics. It is an order-of-magnitude estimator
built from three interpretable terms, every constant exposed in system.yaml:

    ΔG_bind(host, guest) ≈  cavity(guest)                      # hydrophobic burial
                          +  electrostatics(host, guest)        # screened Coulomb
                          +  fluorophilic(host, guest)          # fluorous complementarity

Use it to decide WHICH modifications are worth the real fep_ti_ddg run — not to
report final affinities. Anything it flags should be confirmed by APR/FEP.

Displacement favorability (the sensor figure of merit):

    ΔG_displace(host, pfas) = ΔG_bind(host, pfas) - ΔG_bind(host, dye)

    negative  => PFAS out-competes the dye => displacement is downhill => good sensor
"""
from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass, asdict
from pathlib import Path

from . import utils

log = utils.get_logger("cd_pfas.prefilter")


@dataclass
class Binding:
    host: str
    guest: str
    cavity: float
    electrostatic: float
    fluorophilic: float
    dG_bind: float


@dataclass
class HostScore:
    host: str
    rim_charge: int
    fluorophilic: bool
    dG_dye: float
    dG_pfas: dict          # guest -> ΔG_bind
    dG_displace: dict      # guest -> ΔG_displace (neg = favorable)
    favorability: float    # mean (-ΔG_displace) over PFAS guests; higher = better
    verdict: str
    note: str = ""


def electrostatic_term(rim_charge: int, head_charge: int, pf: dict) -> float:
    """Screened Coulomb between the (effective) rim charge and the guest head.

    Only ~`max_effective_rim_charges` rim charges can simultaneously contact a
    single monovalent head group, so we clamp the interacting charge — a linear
    sum over all 7 rim charges would badly overcount.
    """
    q_eff = max(-pf["max_effective_rim_charges"],
                min(pf["max_effective_rim_charges"], rim_charge))
    r = pf["contact_distance_ang"]
    coulomb = pf["coulomb_prefactor"] * q_eff * head_charge / (pf["eps_effective"] * r)
    screening = math.exp(-r / pf["debye_length_ang"])
    return coulomb * screening


def fluorophilic_term(host_fluorophilic: bool, guest_fluorinated: bool, pf: dict) -> float:
    return pf["fluorophilic_bonus_kcal"] if (host_fluorophilic and guest_fluorinated) else 0.0


def bind(host_name: str, rim_charge: int, host_fluoro: bool,
         guest_key: str, gcfg: dict, pf: dict) -> Binding:
    cavity = float(gcfg["cavity_dG_kcal"])
    elec = electrostatic_term(rim_charge, int(gcfg.get("head_charge", gcfg["net_charge"])), pf)
    fluo = fluorophilic_term(host_fluoro, bool(gcfg.get("fluorinated", False)), pf)
    return Binding(host_name, guest_key, cavity, elec, fluo, cavity + elec + fluo)


def score_host(host_name: str, rim_charge: int, host_fluoro: bool,
               guests: dict, dye_key: str, pfas_keys: list[str], pf: dict) -> HostScore:
    b = {gk: bind(host_name, rim_charge, host_fluoro, gk, guests[gk], pf)
         for gk in [dye_key, *pfas_keys]}
    dG_dye = b[dye_key].dG_bind
    dG_pfas = {gk: b[gk].dG_bind for gk in pfas_keys}
    dG_disp = {gk: dG_pfas[gk] - dG_dye for gk in pfas_keys}
    favorability = sum(-v for v in dG_disp.values()) / len(dG_disp)

    # verdict combines: is displacement favorable, and is the probe usable?
    notes = []
    if dG_dye > pf["baseline_dye_min_kcal"]:
        notes.append("dye too weak to load (low baseline signal)")
    if dG_dye < pf["dye_too_tight_kcal"]:
        notes.append("dye bound very tightly (may resist displacement)")
    if favorability > 0:
        verdict = "FAVORABLE — PFAS out-competes dye"
    elif favorability > -0.5:
        verdict = "MARGINAL"
    else:
        verdict = "UNFAVORABLE — dye holds the cavity"
    return HostScore(host_name, rim_charge, host_fluoro, dG_dye, dG_pfas, dG_disp,
                     favorability, verdict, "; ".join(notes))


def run(config_path: str, mods_path: str) -> dict:
    cfg = utils.load_config(config_path)
    mods = utils.load_config(mods_path)
    pf = cfg["prefilter"]
    guests = cfg["guests"]
    dye_key = "dye"
    pfas_keys = [k for k in guests if k != dye_key]

    hosts = [("beta_cyclodextrin (reference)", 0, False, "unmodified β-CD")]
    for m in mods["modifications"]:
        hosts.append((
            m["id"],
            int(m.get("net_charge_delta", 0)),
            bool(m.get("fluorophilic", False)),
            m.get("description", ""),
        ))

    scores = []
    detail = []
    for name, q, fluo, desc in hosts:
        hs = score_host(name, q, fluo, guests, dye_key, pfas_keys, pf)
        hs.note = (hs.note + (" | " if hs.note and desc else "") + desc).strip(" |")
        scores.append(hs)
        for gk in [dye_key, *pfas_keys]:
            detail.append(asdict(bind(name, q, fluo, gk, guests[gk], pf)))

    scores.sort(key=lambda s: -s.favorability)  # best sensor first
    return {
        "dye": dye_key,
        "pfas": pfas_keys,
        "ranked_hosts": [asdict(s) for s in scores],
        "binding_detail": detail,
        "disclaimer": "HEURISTIC prefilter — order-of-magnitude, not converged FEP/APR.",
    }


def print_report(res: dict) -> None:
    dye = res["dye"]
    pfas = res["pfas"]
    print("\n" + "=" * 78)
    print("  HEURISTIC DISPLACEMENT PREFILTER  (NOT converged MD — triage only)")
    print("=" * 78)
    print(f"  reporter dye = {dye} (TNS) | PFAS = {', '.join(pfas)}")
    print(f"  ΔG_displace = ΔG_bind(PFAS) − ΔG_bind(dye);  negative = PFAS wins = good\n")
    hdr = f"{'host':<34}{'q':>3}{'Fφ':>4}{'ΔG_dye':>8}"
    for gk in pfas:
        hdr += f"{'disp['+gk+']':>13}"
    hdr += f"{'FOM':>7}  verdict"
    print(hdr)
    print("-" * len(hdr))
    for s in res["ranked_hosts"]:
        row = f"{s['host'][:33]:<34}{s['rim_charge']:>+3}{'Y' if s['fluorophilic'] else '·':>4}"
        row += f"{s['dG_dye']:>8.2f}"
        for gk in pfas:
            row += f"{s['dG_displace'][gk]:>+13.2f}"
        row += f"{s['favorability']:>+7.2f}  {s['verdict']}"
        print(row)
    print("-" * len(hdr))
    print("  q = net rim charge · Fφ = fluorophilic cavity · FOM = mean(−ΔG_displace)\n")


def main() -> None:
    ap = argparse.ArgumentParser(description="Heuristic binding/displacement prefilter.")
    ap.add_argument("--config", default="config/system.yaml")
    ap.add_argument("--modifications", default="config/modifications.yaml")
    ap.add_argument("--out", default="work/prefilter/prefilter_report.json")
    args = ap.parse_args()

    res = run(args.config, args.modifications)
    out = utils.resolve(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(res, indent=2))
    print_report(res)
    print(f"  full report -> {out}")


if __name__ == "__main__":
    main()
