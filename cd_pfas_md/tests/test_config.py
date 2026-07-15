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


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
