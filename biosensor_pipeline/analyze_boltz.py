"""
analyze_boltz.py -- turn Boltz predictions into the reported validation.

Reads biosensor_out/boltz_results.json (the compact metric table captured from
the Boltz MCP tool) and any downloaded chimera CIF models, then:

  * ✅ measures the TEM-1 catalytic constellation (Cα–Cα distances of
    S70/K73/S130/D131/N132/E166) in each predicted chimera model and in native
    TEM-1 (3GMW) as a reference;
  * ⚠️ computes the illustrative switch proxy per construct;
  * prints the apo-vs-holo and chimera-vs-native comparisons.

    python3 -m biosensor_pipeline.analyze_boltz
"""

from __future__ import annotations
import json
import os
import warnings
warnings.filterwarnings("ignore")

from .systems import TEM1, SYSTEMS, get_system
from .screen import build_library
from .scoring import SwitchFeatures, switch_proxy

OUT = os.path.join(os.getcwd(), "biosensor_out")
CATALYTIC_AMBLER = {"S70": 70, "K73": 73, "S130": 130, "D131": 131, "N132": 132, "E166": 166}


def _ca_coords_by_resid(cif_path, wanted_resids, chain=None):
    import biotite.structure.io.pdbx as pdbx
    import biotite.structure as struc
    import numpy as np
    cif = pdbx.CIFFile.read(cif_path)
    arr = pdbx.get_structure(cif, model=1)
    arr = arr[struc.filter_amino_acids(arr)]
    if chain is None:
        chains = {c: int((arr.chain_id == c).sum()) for c in set(arr.chain_id)}
        chain = max(chains, key=chains.get)
    arr = arr[arr.chain_id == chain]
    ca = arr[arr.atom_name == "CA"]
    out = {}
    for label, rid in wanted_resids.items():
        sel = ca[ca.res_id == rid]
        if len(sel) >= 1:
            out[label] = sel.coord[0]
    return out


def _pairwise(coords):
    import numpy as np
    labels = list(coords)
    d = {}
    for i in range(len(labels)):
        for j in range(i + 1, len(labels)):
            a, b = labels[i], labels[j]
            d[f"{a}-{b}"] = round(float(np.linalg.norm(coords[a] - coords[b])), 2)
    vals = list(d.values())
    return {"pairs": d, "n": len(coords),
            "max": round(max(vals), 2) if vals else None,
            "mean": round(sum(vals) / len(vals), 2) if vals else None}


def native_reference():
    """Catalytic constellation of native TEM-1 from the deposited 3GMW model."""
    import biotite.database.rcsb as rcsb
    import tempfile
    path = rcsb.fetch("3GMW", "cif", tempfile.mkdtemp())   # text CIF (CIFFile.read)
    coords = _ca_coords_by_resid(path, CATALYTIC_AMBLER, chain="A")   # 3GMW author numbering == Ambler
    return _pairwise(coords)


def chimera_constellation(key, kind):
    """Constellation in a predicted chimera CIF (chain renumbered 1..N)."""
    sysm = get_system(key)
    lib = build_library(sysm)
    ch = lib[len(lib) // 2]
    # reporter catalytic residues are before the insertion -> chimera res_id = reporter_idx + 1
    wanted = {lab: TEM1.catalytic[lab] + 1 for lab in CATALYTIC_AMBLER}
    cif = os.path.join(OUT, f"{key}_{kind}.cif")
    if not os.path.exists(cif):
        return None
    coords = _ca_coords_by_resid(cif, wanted)
    return _pairwise(coords)


def main():
    with open(os.path.join(OUT, "boltz_results.json")) as fh:
        res = json.load(fh)["results"]

    print("=" * 74)
    print("IN-SILICO VALIDATION (Boltz-2.1)")
    print("=" * 74)

    # --- ✅ catalytic constellation ---
    print("\n[✅ geometric] TEM-1 catalytic constellation (Cα–Cα, Å)")
    try:
        ref = native_reference()
        print(f"  native TEM-1 (3GMW)      : n={ref['n']} mean={ref['mean']} max={ref['max']}")
    except Exception as e:
        ref = None
        print(f"  native reference unavailable ({e})")
    for key in sorted({r["system"] for r in res.values()}):
        for kind in ("holo", "apo"):
            c = chimera_constellation(key, kind)
            if c:
                tag = ""
                if ref and c["max"] is not None:
                    tag = f"  (Δmax vs native {c['max']-ref['max']:+.1f} Å)"
                print(f"  {key:5s} {kind:4s} chimera    : n={c['n']} mean={c['mean']} max={c['max']}{tag}")

    # --- ⚠️ switch proxy + comparisons ---
    print("\n[metrics] fold + binding (as returned by Boltz)")
    hdr = f"  {'construct':13s} {'kind':7s} {'sconf':>6s} {'ptm':>6s} {'cplddt':>7s} {'lig_iptm':>8s} {'bind':>6s}"
    print(hdr)
    for k, r in res.items():
        print(f"  {k:13s} {r['kind']:7s} {r['structure_confidence']:6.3f} {r['ptm']:6.3f} "
              f"{r['complex_plddt']:7.3f} "
              f"{(r['ligand_iptm'] if r['ligand_iptm'] is not None else 0):8.3f} "
              f"{(r['binding_confidence'] if r['binding_confidence'] is not None else 0):6.3f}")

    # --- calibrated triage ranking of holo candidates (grounded on wet-lab labels) ---
    # reference/boltz-score-calibration.md: Boltz score triages WITHIN a panel by
    # percentile; it predicts WHETHER it binds, not affinity. Annotate accordingly.
    from .calibration import calibrated_rank, BASELINE_HIT_RATE
    holo_cands = [{"construct": k, "ligand_iptm": r["ligand_iptm"]}
                  for k, r in res.items() if r["kind"] == "holo" and r["ligand_iptm"] is not None]
    ranked = calibrated_rank(holo_cands, score_key="ligand_iptm")
    if ranked and ranked[0].get("rank"):
        print("\n[calibrated triage] holo candidates by Boltz ligand_iptm "
              f"(baseline hit-rate {BASELINE_HIT_RATE}%; within-panel percentile, NOT affinity)")
        print(f"  {'rank':>4s} {'construct':13s} {'lig_iptm':>8s} {'pctile':>7s} "
              f"{'exp_hit%':>8s} {'enrich':>6s}")
        for it in ranked:
            if it.get("rank"):
                print(f"  {it['rank']:>4d} {it['construct']:13s} {it['ligand_iptm']:8.3f} "
                      f"{it['within_panel_percentile']:7.2f} {it['expected_hit_rate_pct']:8.1f} "
                      f"{it['enrichment_vs_baseline']:5.2f}x")
        print("  note: " + ranked[0]["triage_note"])

    print("\n[⚠️ hypothesis] switch proxy per system (illustrative, NOT a predicted DR)")
    summary = {"_calibrated_triage": [
        {kk: it[kk] for kk in ("rank", "construct", "ligand_iptm", "within_panel_percentile",
                               "expected_hit_rate_pct", "enrichment_vs_baseline")}
        for it in ranked if it.get("rank")]}
    keys = sorted({r["system"] for r in res.values()})
    for key in keys:
        holo = res[f"{key}_holo"]; apo = res[f"{key}_apo"]; ctrl = res[f"{key}_ctrl"] if f"{key}_ctrl" in res else res[f"ctrl_{key}"]
        c = chimera_constellation(key, "holo")
        active_ok = None
        if c and c["max"] is not None:
            try:
                active_ok = c["max"] <= (native_reference()["max"] + 6.0)
            except Exception:
                active_ok = None
        feat = SwitchFeatures(
            fold_confidence=holo["structure_confidence"],
            receptor_plddt=holo["complex_plddt"] * 100.0,
            ligand_binding=holo["ligand_iptm"],
            active_site_ok=active_ok,
            apo_holo_shift=abs(holo["complex_plddt"] - apo["complex_plddt"]),
        )
        proxy = switch_proxy(feat)
        # chimera-vs-native binding retention (a real, comparable metric)
        retention = round(holo["ligand_iptm"] / ctrl["ligand_iptm"], 2) if ctrl["ligand_iptm"] else None
        summary[key] = {
            "role": get_system(key).role,
            "analyte": get_system(key).receptor.ligand_name,
            "chimera_ligand_iptm": holo["ligand_iptm"],
            "native_ligand_iptm": ctrl["ligand_iptm"],
            "binding_retention_vs_native": retention,
            "apo_to_holo_plddt_gain": round(holo["complex_plddt"] - apo["complex_plddt"], 3),
            "active_site_intact": active_ok,
            "switch_proxy": proxy["proxy"],
        }
        print(f"  {key:5s} ({summary[key]['role']:12s}, {summary[key]['analyte']}): "
              f"proxy={proxy['proxy']}  binding_retention_vs_native={retention}  "
              f"apo->holo_plddt_gain={summary[key]['apo_to_holo_plddt_gain']:+.3f}  "
              f"active_site_intact={active_ok}")

    with open(os.path.join(OUT, "analysis_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    print(f"\nwrote {OUT}/analysis_summary.json")
    print("\nReminder: proxy is ⚠️ illustrative. Ground-truth switchability is a "
          "wet-lab kobs(+L)/kobs(−L) titration.")


if __name__ == "__main__":
    main()
