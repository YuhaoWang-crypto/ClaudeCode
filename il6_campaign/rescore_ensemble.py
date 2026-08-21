"""Combine the Boltz-2 and Chai-1 ipSAE_min scores into the paper-style ranking, and
gate the result on self-consistency DockQ.

Ranking score  = mean of per-run z-scores of ipSAE_min from each predictor
Gate           = DockQ(design model, independent Chai-1 prediction) >= threshold

Chai-1 never saw the design run, so it is the orthogonal judge; DockQ answers whether
its prediction is the *designed* pose rather than some other way of sticking the two
chains together.

Usage:
    python rescore_ensemble.py --boltz ipsae_scores.csv --chai chai_scores.json \
        --interface interface_all.csv --design-root pae/ --out rescored.csv
"""

import argparse
import csv
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from dockq import dockq, load  # noqa: E402


def ipsae_from_chai(rec, ipsae_script, workdir: Path, pae_cut=10, dist_cut=15):
    """Write Chai-1 PAE + structure to disk and score with the reference ipsae.py."""
    d = workdir / rec["design_id"].replace("|", "__") / f"seed{rec['seed']}"
    d.mkdir(parents=True, exist_ok=True)
    cif = d / "pred.cif"
    cif.write_text(rec["cif"])
    npz = d / "pae.npz"
    np.savez(npz, pae=np.array(rec["pae"], dtype=np.float32))

    tag = f"{int(pae_cut):02d}_{int(dist_cut):02d}"
    out = d / f"pred_{tag}.txt"
    if not out.exists():
        subprocess.run([sys.executable, str(ipsae_script), str(npz), str(cif),
                        str(pae_cut), str(dist_cut)], capture_output=True, text=True, timeout=600)
    if not out.exists():
        return None, cif
    vals = []
    with out.open() as fh:
        header = []
        for line in fh:
            header = line.split()
            if header:
                break
        for line in fh:
            p = line.split()
            if len(p) >= len(header):
                row = dict(zip(header, p))
                if row.get("Type") == "asym":
                    vals.append(float(row["ipSAE"]))
    return (min(vals) if vals else None), cif


def zscore(values):
    v = np.asarray(values, dtype=float)
    s = v.std(ddof=0)
    return (v - v.mean()) / s if s > 0 else np.zeros_like(v)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--boltz", required=True, type=Path)
    ap.add_argument("--chai", required=True, type=Path)
    ap.add_argument("--interface", required=True, type=Path)
    ap.add_argument("--design-root", required=True, type=Path)
    ap.add_argument("--ipsae", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--dockq-gate", type=float, default=0.23)
    args = ap.parse_args()

    boltz = {r["design_id"]: r for r in csv.DictReader(args.boltz.open())}
    iface = {r["design_id"]: r for r in csv.DictReader(args.interface.open())}
    chai_raw = json.loads(args.chai.read_text())

    work = Path(tempfile.mkdtemp(prefix="chai_ipsae_"))
    rows = []
    for rec in chai_raw:
        run, did = rec["design_id"].split("|")
        ips_chai, chai_cif = ipsae_from_chai(rec, args.ipsae, work)
        design_cif = args.design_root / run / did / "predicted_structure.cif"
        gate = {}
        if design_cif.exists():
            try:
                gate = dockq(load(design_cif), load(chai_cif), "B", "X", "A", "B", -2, 0)
            except Exception as e:
                gate = {"error": str(e)[:80]}
        b, f = boltz.get(did, {}), iface.get(did, {})
        rows.append({
            "run": run, "design_id": did,
            "iptm_boltz": float(b.get("iptm", "nan")),
            "ipsae_boltz": float(b.get("ipsae_min", "nan")),
            "iptm_chai": rec["iptm"], "plddt_chai": rec["plddt"],
            "ipsae_chai": ips_chai if ips_chai is not None else float("nan"),
            "chai_clash": rec["has_inter_chain_clashes"],
            "dockq": gate.get("dockq", float("nan")), "fnat": gate.get("fnat", float("nan")),
            "irmsd_A": gate.get("irmsd_A", float("nan")),
            "dockq_quality": gate.get("quality", "n/a"),
            "bsa_A2": float(f.get("bsa_A2", "nan")),
            "siteI_recall": float(f.get("siteI_recall", "nan")),
            "siteII_recall": float(f.get("siteII_recall", "nan")),
            "sequence": f.get("sequence", ""),
        })

    # ensemble: per-run z-score of each predictor's ipSAE_min, averaged
    for run in {r["run"] for r in rows}:
        idx = [i for i, r in enumerate(rows) if r["run"] == run]
        zb = zscore([rows[i]["ipsae_boltz"] for i in idx])
        zc = zscore([np.nan_to_num(rows[i]["ipsae_chai"], nan=0.0) for i in idx])
        zt = zscore([rows[i]["iptm_boltz"] for i in idx])
        for k, i in enumerate(idx):
            rows[i]["z_ipsae_boltz"] = round(float(zb[k]), 3)
            rows[i]["z_ipsae_chai"] = round(float(zc[k]), 3)
            rows[i]["ensemble_score"] = round(float((zb[k] + zc[k]) / 2), 3)
            rows[i]["z_iptm_boltz_only"] = round(float(zt[k]), 3)
        # ranks under the old and new schemes, plus the gated new scheme
        for key, name in (("z_iptm_boltz_only", "rank_old_iptm"),
                          ("ensemble_score", "rank_new_ensemble")):
            order = sorted(idx, key=lambda i: -rows[i][key])
            for pos, i in enumerate(order, 1):
                rows[i][name] = pos
        passed = [i for i in idx if (rows[i]["dockq"] == rows[i]["dockq"])
                  and rows[i]["dockq"] >= args.dockq_gate]
        order = sorted(passed, key=lambda i: -rows[i]["ensemble_score"])
        for i in idx:
            rows[i]["rank_gated"] = ""
        for pos, i in enumerate(order, 1):
            rows[i]["rank_gated"] = pos

    fields = ["run", "design_id", "rank_gated", "rank_new_ensemble", "rank_old_iptm",
              "ensemble_score", "ipsae_boltz", "ipsae_chai", "z_ipsae_boltz", "z_ipsae_chai",
              "iptm_boltz", "iptm_chai", "plddt_chai", "chai_clash", "dockq", "fnat",
              "irmsd_A", "dockq_quality", "bsa_A2", "siteI_recall", "siteII_recall", "sequence"]
    with args.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(sorted(rows, key=lambda r: (r["run"], -r["ensemble_score"])))
    print("wrote", args.out, f"({len(rows)} designs)")

    for run in sorted({r["run"] for r in rows}):
        sub = [r for r in rows if r["run"] == run]
        print(f"\n== {run}")
        print(f"{'design':<22}{'ens':>7}{'ipsB':>7}{'ipsC':>7}{'iptmB':>7}{'iptmC':>7}"
              f"{'DockQ':>7}  {'gate':<11}{'old→new rank'}")
        for r in sorted(sub, key=lambda x: -x["ensemble_score"])[:10]:
            print(f"{r['design_id'].replace('pres_','')[:20]:<22}{r['ensemble_score']:>7.2f}"
                  f"{r['ipsae_boltz']:>7.3f}{r['ipsae_chai']:>7.3f}{r['iptm_boltz']:>7.3f}"
                  f"{r['iptm_chai']:>7.3f}{r['dockq']:>7.3f}  {r['dockq_quality']:<11}"
                  f"{r['rank_old_iptm']}→{r['rank_new_ensemble']}")


if __name__ == "__main__":
    main()
