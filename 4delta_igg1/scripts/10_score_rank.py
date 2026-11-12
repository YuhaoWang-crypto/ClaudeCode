#!/usr/bin/env python3
"""Task 3d/4: combined evaluation matrix, ranking and down-selection to 12.

Every anchor below is calibrated on a measured reference, not invented:
  exit angle      0 deg = straight down an open direction out of the alpha-v/beta-6
                  cleft; the natural pro-TGF-beta1 exit is 1.8 deg (script 07)
  motif fit       2.08 A = where the EXPERIMENTAL free peptide (6R2X) sits
                  relative to the receptor-bound LAP pose (script 01/06)
  helix core      1.00 = the 6R2X ensemble
  RGD rel SASA    0.577 = the 6R2X ensemble
  interface       WT_CH3 / WT_Fv predicted the same way

Two scores are reported separately and on purpose:
  score_structural  pure structure/steric quality, 0-100
  score_total       after programme-level deductions (insert burden, site risk,
                    loss of the parent antibody's own paratope)
so the client can re-weight the deductions without touching the structural work.
"""
import json

import numpy as np

REF_HELIX, REF_RGD_SASA, REF_MOTIF_FIT = 1.00, 0.577, 2.08


def clamp(x):
    return float(max(0.0, min(1.0, x)))


def lin(x, good, bad):
    """1.0 at `good`, 0.0 at `bad`, linear in between (either direction)."""
    if x is None:
        return 0.0
    return clamp((x - bad) / (good - bad))


ev = json.load(open("results/05_evaluation.json"))["results"]
st = json.load(open("results/06_integrin_steric.json"))["results"]
xv = json.load(open("results/07_exit_vector.json"))["candidates"]
sc = json.load(open("results/09_scaffold_integrity.json"))
cand = json.load(open("results/03_candidates.json"))
conf = json.load(open("results/04_confidence.json"))
by_id = {c["id"]: c for c in cand["candidates"] + cand["controls"]}

WT = {"CH3": sc.get("WT_CH3", {}).get("metrics", {}),
      "Fv": sc.get("WT_Fv", {}).get("metrics", {})}
wt_iface = {k: (v.get("interchain_contacts") or {}).get("mean") for k, v in WT.items()}
print(f"WT interchain contacts: {wt_iface}")

# calibrate bulky-probe RGD accessibility on the FREE peptide control: that is the
# physical ceiling, since nothing is occluding the motif there
free = ev.get("REF_4delta_free", {}).get("metrics", {})
BULKY_FREE = ((free.get("rgd_bulky_probe_sasa_A2") or {}).get("mean")) or 400.0
# an internal-clash gate relative to the WT scaffolds predicted the same way
wt_clash = []
for w in ("WT_CH3", "WT_Fv"):
    v = (ev.get(w, {}).get("metrics", {}).get("graft_scaffold_clashes") or {}).get("mean")
    if v is not None:
        wt_clash.append(v)
CLASH_GATE = (max(wt_clash) + 3.0) if wt_clash else 6.0
print(f"free-peptide bulky-probe RGD SASA (ceiling) = {BULKY_FREE:.0f} A^2; "
      f"internal-clash gate = {CLASH_GATE:.1f}")

rows = []
for cid, c in by_id.items():
    if cid not in ev or cid not in xv or cid not in st:
        continue
    m = ev[cid]["metrics"]
    s = sc.get(cid, {}).get("metrics", {})
    kind = sc.get(cid, {}).get("module")

    def g(d, key, stat="mean"):
        v = (d or {}).get(key)
        return None if v is None else v.get(stat)

    helix = g(m, "helix_frac_core")
    rmsd6 = g(m, "rmsd_to_6R2X_A")
    rgd_sasa = g(m, "rgd_rel_sasa")
    rgd_near = g(m, "rgd_min_dist_to_scaffold_A")
    rgd_bulky = g(m, "rgd_bulky_probe_sasa_A2")
    clash_int = g(m, "graft_scaffold_clashes")
    e2e = g(m, "end_to_end_A")
    fit_best = st[cid]["motif_fit_rmsd_A"]["min"]
    ang_best = xv[cid]["angle_to_open_deg"]["min"]
    ang_mean = xv[cid]["angle_to_open_deg"]["mean"]
    angs = [r["angle_to_nearest_open_direction_deg"] for r in xv[cid]["per_sample"]]
    ang_sd = float(np.std(angs)) if angs else None
    scaf_rmsd = g(s, "scaffold_rmsd_vs_wt_A")
    iface = g(s, "interchain_contacts")
    graft_plddt = g(s, "graft_plddt")
    core_clash = st[cid]["core_receptor_clashes"]["min"]

    iface_ratio = None
    if iface and kind in wt_iface and wt_iface[kind]:
        iface_ratio = iface / wt_iface[kind]

    # ----------------------------------------------------------- component scores
    c1a = 22 * lin(ang_best, 0.0, 60.0)
    c1b = 13 * lin(fit_best, 1.5, 3.0)
    c2a = 14 * clamp(helix if helix is not None else 0.0)
    c2b = 8 * lin(rmsd6, 1.5, 5.0)
    c3a = 8 * clamp((rgd_sasa or 0.0) / REF_RGD_SASA)
    # can a BULKY partner reach the motif? 5 A probe, scaled to the free peptide
    c3b = 5 * clamp((rgd_bulky or 0.0) / BULKY_FREE)
    c4a = 8 * lin(scaf_rmsd, 0.5, 3.0)
    c4b = 7 * (clamp(1 - abs(iface_ratio - 1) / 0.6) if iface_ratio else 0.0)
    c4c = 5 * lin(clash_int, 0.0, CLASH_GATE)
    c5a = 6 * lin(graft_plddt, 80.0, 40.0)
    c5b = 4 * lin(ang_sd, 0.0, 45.0)

    structural = sum([c1a, c1b, c2a, c2b, c3a, c3b, c4a, c4b, c4c, c5a, c5b])

    # ----------------------------------------------------------- gates
    gates = {
        "no_excess_internal_clash": bool(clash_int is not None and clash_int <= CLASH_GATE),
        "helix_core_kept": bool(helix is not None and helix >= 0.50),
        "rgd_not_buried": bool(rgd_sasa is not None and rgd_sasa >= 0.30),
        "rgd_reachable_by_bulky_partner": bool(
            rgd_bulky is not None and rgd_bulky >= 0.35 * BULKY_FREE),
        "interface_preserved": bool(iface_ratio is None or 0.6 <= iface_ratio <= 1.5),
    }
    passed = all(gates.values())

    # ----------------------------------------------------------- deductions
    ded, notes = 0.0, []
    il = c.get("insert_len")
    if il and il > 34:
        ded += 4; notes.append(f"insert {il} aa: large foreign surface")
    elif il and il > 30:
        ded += 2; notes.append(f"insert {il} aa")
    if c.get("site") == "CH3-AB":
        ded += 4; notes.append("AB loop overlaps the CH3/CH3 dimer interface (EU 356/357/360)")
    if c.get("site") == "CDRH3":
        ded += 6; notes.append("CDRH3 graft destroys the parent antibody's own paratope "
                               "(trastuzumab CDRH3 contacts HER2 at H99/H103/H104/H105)")
    if c.get("strategy") == "S1_minimal":
        ded += 3; notes.append("C-junction cannot span the slack (reach 6 A vs 23-30 A needed)")

    rows.append({
        "id": cid, "site": c.get("site"), "strategy": c.get("strategy"),
        "module": kind, "insert_len": il,
        "score_structural": round(structural, 1),
        "deductions": round(ded, 1),
        "score_total": round(structural - ded, 1),
        "gates": gates, "gates_passed": passed,
        "deduction_notes": notes,
        "metrics": {
            "exit_angle_best_deg": ang_best, "exit_angle_mean_deg": ang_mean,
            "exit_angle_sd_deg": None if ang_sd is None else round(ang_sd, 1),
            "motif_fit_best_A": fit_best,
            "helix_frac_core": helix, "rmsd_to_6R2X_A": rmsd6,
            "rgd_rel_sasa": rgd_sasa, "rgd_min_dist_to_scaffold_A": rgd_near,
            "rgd_bulky_probe_sasa_A2": rgd_bulky,
            "payload_end_to_end_A": e2e,
            "internal_graft_clashes": clash_int,
            "scaffold_rmsd_vs_wt_A": scaf_rmsd,
            "interchain_contacts": iface,
            "interface_ratio_vs_wt": None if iface_ratio is None else round(iface_ratio, 2),
            "graft_plddt": graft_plddt,
            "min_core_receptor_clashes": core_clash,
        },
        "component_scores": {
            "exit_geometry_22": round(c1a, 1), "motif_fit_13": round(c1b, 1),
            "helix_core_14": round(c2a, 1), "conformation_rmsd_8": round(c2b, 1),
            "rgd_exposure_8": round(c3a, 1), "rgd_bulky_access_5": round(c3b, 1),
            "scaffold_rmsd_8": round(c4a, 1), "interface_7": round(c4b, 1),
            "internal_clash_5": round(c4c, 1),
            "graft_confidence_6": round(c5a, 1), "sample_consistency_4": round(c5b, 1),
        },
    })

designs = [r for r in rows if r["id"].startswith("C")]
controls = [r for r in rows if not r["id"].startswith("C")]
designs.sort(key=lambda r: -r["score_total"])
for k, r in enumerate(designs, 1):
    r["rank"] = k
    r["selected_top12"] = k <= 12

hdr = (f"{'rank':<5}{'id':<22}{'site':<9}{'strat':<6}{'struct':>7}{'ded':>5}{'total':>7}"
       f"{'exit°':>7}{'fit':>6}{'helix':>7}{'RGDsasa':>8}{'ifaceR':>7}{'pLDDT':>7}  gates")
print("\n" + hdr); print("-" * len(hdr))
for r in designs:
    m = r["metrics"]
    print(f"{r['rank']:<5}{r['id']:<22}{r['site']:<9}{r['strategy'].split('_')[0]:<6}"
          f"{r['score_structural']:>7.1f}{r['deductions']:>5.1f}{r['score_total']:>7.1f}"
          f"{m['exit_angle_best_deg']:>7.1f}{m['motif_fit_best_A']:>6.2f}"
          f"{(m['helix_frac_core'] or 0):>7.2f}{(m['rgd_rel_sasa'] or 0):>8.2f}"
          f"{str(m['interface_ratio_vs_wt']):>7}{str(m['graft_plddt']):>7}  "
          f"{'PASS' if r['gates_passed'] else 'FAIL:' + ','.join(k for k,v in r['gates'].items() if not v)}")

print("\ncontrols / benchmark:")
for r in controls:
    m = r["metrics"]
    print(f"     {r['id']:<22}{str(r['site']):<9}{'-':<6}{r['score_structural']:>7.1f}"
          f"{r['deductions']:>5.1f}{r['score_total']:>7.1f}"
          f"{m['exit_angle_best_deg']:>7.1f}{m['motif_fit_best_A']:>6.2f}"
          f"{(m['helix_frac_core'] or 0):>7.2f}{(m['rgd_rel_sasa'] or 0):>8.2f}"
          f"{str(m['interface_ratio_vs_wt']):>7}{str(m['graft_plddt']):>7}")

json.dump({
    "scoring_rubric": {
        "calibration": {
            "exit_angle_deg": "0 = straight out of an open direction; natural pro-TGF-beta1 exit = 1.8 deg",
            "motif_fit_A": f"experimental free peptide (6R2X) baseline = {REF_MOTIF_FIT} A vs bound LAP pose",
            "helix_frac_core": f"6R2X ensemble = {REF_HELIX}",
            "rgd_rel_sasa": f"6R2X ensemble = {REF_RGD_SASA}",
            "interface_ratio": "1.00 = same interchain contact count as the WT module predicted identically",
        },
        "weights_sum_to_100": {
            "exit_geometry": 22, "motif_fit": 13, "helix_core": 14,
            "conformation_rmsd": 8, "rgd_exposure": 8, "rgd_bulky_access": 5,
            "scaffold_rmsd": 8, "interface": 7, "internal_clash": 5,
            "graft_confidence": 6, "sample_consistency": 4,
        },
        "gates": "a candidate must pass all four to be called advanceable",
        "deductions": "programme-level, reported separately from score_structural",
    },
    "wt_reference_interchain_contacts": wt_iface,
    "designs_ranked": designs,
    "controls": controls,
}, open("results/10_ranking.json", "w"), indent=2)
print("\nwrote results/10_ranking.json")
