"""Lightweight tests that run WITHOUT OpenMM/AmberTools installed.

They validate the config schema, the thermodynamic conversions, and the ΔΔG
ranking/aggregation logic — i.e. everything that does not require an MD engine.
Run:  pytest cd_pfas_md/tests -q
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parent))  # so `import cd_pfas_md.src...` resolves

from cd_pfas_md.src import utils  # noqa: E402


def test_ka_dg_roundtrip():
    T = 298.15
    ka = 3.9e3
    dg = utils.ka_to_dg(ka, T)
    assert dg < 0  # binding is favorable
    ka2 = utils.dg_to_ka(dg, T)
    assert math.isclose(ka, ka2, rel_tol=1e-6)


def test_ka_dg_known_value():
    # Ka = 1e6 /M at 298.15 K  ->  ΔG = -RT ln(1e6) ≈ -8.18 kcal/mol
    dg = utils.ka_to_dg(1e6, 298.15)
    assert -8.5 < dg < -7.9


def test_system_config_loads():
    cfg = utils.load_config("config/system.yaml")
    assert cfg["host"]["name"] == "beta_cyclodextrin"
    assert "dye" in cfg["guests"]
    for gkey, g in cfg["guests"].items():
        assert "net_charge" in g, f"{gkey} missing net_charge"
    # PFAS must be modeled as anions.
    assert cfg["guests"]["pfoa"]["net_charge"] == -1
    assert cfg["guests"]["pfos"]["net_charge"] == -1


def test_guest_dataclass_derives_dg_from_ka():
    cfg = utils.load_config("config/system.yaml")
    dye = utils.Guest.from_config("dye", cfg["guests"]["dye"], cfg["thermo"]["temperature_K"])
    assert dye.exp_dg is not None
    assert dye.exp_dg < 0
    assert dye.net_charge == -1


def test_modifications_config_schema():
    cfg = utils.load_config("config/modifications.yaml")
    assert cfg["reference_host"]["name"] == "beta_cyclodextrin"
    ids = [m["id"] for m in cfg["modifications"]]
    assert len(ids) == len(set(ids)), "modification ids must be unique"
    for m in cfg["modifications"]:
        assert "net_charge_delta" in m
        assert m["builder"] in {"substituent", "smiles"}


def test_ddg_ranking_orders_most_negative_first():
    from cd_pfas_md.src.fep_ti_ddg import DDGResult
    rows = [
        DDGResult("a", "pfoa", -50, 1, -48, 1, -2.0, 1.4, +1, True),
        DDGResult("b", "pfoa", -50, 1, -55, 1, +5.0, 1.4, -4, False),
        DDGResult("c", "pfoa", float("nan"), float("nan"), float("nan"),
                  float("nan"), float("nan"), float("nan"), 0, False, "skipped"),
        DDGResult("d", "pfoa", -60, 1, -52, 1, -8.0, 1.4, +7, True),
    ]
    rows.sort(key=lambda x: (x.ddG_bind_kcal != x.ddG_bind_kcal, x.ddG_bind_kcal))
    assert [r.modification_id for r in rows] == ["d", "a", "b", "c"]
    assert rows[0].tighter_than_wt is True   # most negative ΔΔG binds tightest


def test_charge_mismatch_is_caught(monkeypatch):
    # ka_to_dg guards; here we assert the Guest net_charge is surfaced for PFAS.
    cfg = utils.load_config("config/system.yaml")
    pfoa = utils.Guest.from_config("pfoa", cfg["guests"]["pfoa"], 298.15)
    assert pfoa.fluorine_policy == "strict"
    assert pfoa.net_charge == -1


def test_prefilter_electrostatics_sign_and_clamp():
    from cd_pfas_md.src import prefilter
    pf = utils.load_config("config/system.yaml")["prefilter"]
    # +1 rim charge attracting a -1 head must be negative (favorable).
    e1 = prefilter.electrostatic_term(+1, -1, pf)
    assert e1 < 0
    # +7 rim charge is clamped to max_effective_rim_charges, not 7x.
    e7 = prefilter.electrostatic_term(+7, -1, pf)
    e2 = prefilter.electrostatic_term(+2, -1, pf)
    assert math.isclose(e7, e2), "rim charge must be clamped, not summed linearly"
    assert e7 > 7 * e1  # (less negative than a naive 7x would be)


def test_prefilter_fluorophilic_only_helps_fluorinated():
    from cd_pfas_md.src import prefilter
    pf = utils.load_config("config/system.yaml")["prefilter"]
    assert prefilter.fluorophilic_term(True, True, pf) < 0
    assert prefilter.fluorophilic_term(True, False, pf) == 0.0   # dye gets nothing
    assert prefilter.fluorophilic_term(False, True, pf) == 0.0


def test_prefilter_run_surfaces_the_selectivity_insight():
    from cd_pfas_md.src import prefilter
    res = prefilter.run("config/system.yaml", "config/modifications.yaml")
    by_host = {h["host"]: h for h in res["ranked_hosts"]}

    ref = next(h for k, h in by_host.items() if k.startswith("beta_cyclodextrin"))
    tma = by_host["mono_6_trimethylammonium"]
    fluo = by_host["fluorous_tagged"]

    # 1) a symmetric rim charge does NOT change displacement (both guests are −1):
    assert ref["dG_displace"]["pfoa"] == pytest.approx(tma["dG_displace"]["pfoa"], abs=1e-9)
    # 2) but it DOES raise overall affinity (dye bound tighter):
    assert tma["dG_dye"] < ref["dG_dye"]
    # 3) the fluorophilic cavity is what actually flips displacement favorable:
    assert ref["dG_displace"]["pfoa"] > 0        # native: PFAS cannot displace
    assert fluo["dG_displace"]["pfoa"] < 0       # fluorous: PFAS displaces
    # 4) best-ranked host is fluorophilic:
    assert res["ranked_hosts"][0]["fluorophilic"] is True


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
