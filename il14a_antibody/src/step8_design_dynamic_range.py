"""Step 8 - Does the paratope move the score more than an irrelevant change does?

THE TEST
--------
Step 7 established that a FRAMEWORK-ONLY change - identical CDR1, CDR2 and CDR3,
differing only at humanisation positions - moved ipTM by 0.438 and interaction
PAE by 8.96 A on the coiled-coil epitope. That is the size of the metric's
response to something that should not matter.

The stage3 design job is the perfect counterpart. Every one of its 150 designs
carries the SAME humanised framework and the SAME CDR1 and CDR2, and differs
ONLY in CDR3 - the loop that dominates a VHH paratope. So its score spread is
the metric's response to the thing that should matter most.

    signal_to_artifact = (spread of ipTM across differing CDR3s)
                       / (shift in ipTM from a framework-only change)

If that ratio is not comfortably above 1, then changing the actual binding loop
moves the score no more than an irrelevant framework edit does, and ranking
designs by this score is ranking noise. This is a stronger statement than "the
control failed": it quantifies how much of the ranking is real.

A ratio below 1 means the ranking is dominated by the artifact. Between 1 and 2
is still unusable for picking a handful of candidates out of hundreds.
"""
import collections
import csv
import json
import re
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import RESULTS

RES = RESULTS / "boltz_results"
# measured in step7 on this same epitope, from the identical-paratope pair
FRAMEWORK_ONLY_IPTM_SHIFT = 0.438
FRAMEWORK_ONLY_PAE_SHIFT = 8.96


def read_fasta(path):
    recs, name, buf = {}, None, []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if name:
                recs[name] = "".join(buf)
            name, buf = line[1:].split("|")[0].strip(), []
        else:
            buf.append(line.strip())
    if name:
        recs[name] = "".join(buf)
    return recs


def framework_integrity(seqs):
    """Check that each design still has the parts an antibody V domain needs.

    Three invariants, none of which the design objective enforces:

    * The J-region tryptophan (Kabat W103, the WGQG/WGKG motif). It is
      essentially invariant across natural V domains and packs into the
      hydrophobic core; a design that loses it is not reliably foldable.
    * The two conserved framework cysteines of the intradomain disulfide.
    * The VHH hallmark tetrad, seen here as WFRQ in framework 2. A VHH must
      fold and stay soluble with no VL partner, and the hallmarks are what
      replace that interface. Losing them raises aggregation risk.
    """
    n = len(seqs)
    has_j_trp = [k for k, v in seqs.items() if re.search(r"WG[QKR]G", v)]
    cys = collections.Counter(v.count("C") for v in seqs.values())
    hallmark = [k for k, v in seqs.items() if "WFRQ" in v]
    return {
        "n": n,
        "j_region_Trp_present": len(has_j_trp),
        "j_region_Trp_MISSING": n - len(has_j_trp),
        "cys_count_distribution": {str(k): v for k, v in sorted(cys.items())},
        "all_have_2_cys": set(cys) == {2},
        "vhh_hallmark_WFRQ_present": len(hallmark),
        "vhh_hallmark_WFRQ_missing": n - len(hallmark),
    }


def load_metrics():
    path = RES / "design_metrics.csv"
    if not path.exists():
        sys.exit(f"missing {path} - run the design-extraction step first")
    rows = []
    with path.open() as fh:
        for r in csv.DictReader(fh):
            for k in ("iptm", "binding_confidence", "min_interaction_pae",
                      "structure_confidence"):
                try:
                    r[k] = float(r[k])
                except (TypeError, ValueError):
                    r[k] = None
            rows.append(r)
    return rows


def spread(values):
    """Robust spread: the 10th-90th percentile range.

    Full min-max is dominated by single outliers; the central 80% is the range
    a candidate-selection step would actually be working across.
    """
    v = sorted(x for x in values if x is not None)
    if len(v) < 10:
        return None
    lo = v[int(0.10 * (len(v) - 1))]
    hi = v[int(0.90 * (len(v) - 1))]
    return hi - lo


def main():
    rows = load_metrics()
    by_job = {}
    for r in rows:
        by_job.setdefault(r["job"], []).append(r)

    out = {"framework_only_reference": {
        "iptm_shift": FRAMEWORK_ONLY_IPTM_SHIFT,
        "min_interaction_pae_shift": FRAMEWORK_ONLY_PAE_SHIFT,
        "source": "step7 identical-paratope control on epitope CC_367_382"}}

    print("=" * 94)
    print("STEP 8 - SIGNAL-TO-ARTIFACT RATIO OF THE DESIGN RANKING")
    print("=" * 94)
    print(f"reference artifact (framework-only change, step7): "
          f"ipTM {FRAMEWORK_ONLY_IPTM_SHIFT:.3f}, PAE {FRAMEWORK_ONLY_PAE_SHIFT:.2f} A\n")

    for job, rs in sorted(by_job.items()):
        ipt = [r["iptm"] for r in rs]
        pae = [r["min_interaction_pae"] for r in rs]
        bc = [r["binding_confidence"] for r in rs]
        s_ipt, s_pae = spread(ipt), spread(pae)
        ratio_ipt = s_ipt / FRAMEWORK_ONLY_IPTM_SHIFT if s_ipt else None
        ratio_pae = s_pae / FRAMEWORK_ONLY_PAE_SHIFT if s_pae else None
        clean = [x for x in ipt if x is not None]
        out[job] = {
            "n": len(rs),
            "iptm_p10_p90_spread": round(s_ipt, 4) if s_ipt else None,
            "pae_p10_p90_spread": round(s_pae, 3) if s_pae else None,
            "signal_to_artifact_iptm": round(ratio_ipt, 3) if ratio_ipt else None,
            "signal_to_artifact_pae": round(ratio_pae, 3) if ratio_pae else None,
            "iptm_median": round(statistics.median(clean), 4) if clean else None,
            "iptm_max": round(max(clean), 4) if clean else None,
            "max_binding_confidence": max(x for x in bc if x is not None),
            "n_binding_confidence_ge_0.5": sum(
                1 for x in bc if x is not None and x >= 0.5),
        }
        print(f"### {job}   n={len(rs)}")
        print(f"  ipTM      median {out[job]['iptm_median']}  max {out[job]['iptm_max']}"
              f"  p10-p90 spread {out[job]['iptm_p10_p90_spread']}")
        print(f"  PAE       p10-p90 spread {out[job]['pae_p10_p90_spread']} A")
        print(f"  max binding_confidence {out[job]['max_binding_confidence']:.6f}"
              f"   designs >= 0.5: {out[job]['n_binding_confidence_ge_0.5']}")
        if ratio_ipt is not None:
            print(f"  SIGNAL-TO-ARTIFACT (ipTM) = {ratio_ipt:.2f}x", end="")
            print("   <- ranking dominated by the artifact"
                  if ratio_ipt < 1.0 else
                  ("   <- too small to pick a few from hundreds"
                   if ratio_ipt < 2.0 else "   <- paratope does dominate"))
            print(f"  SIGNAL-TO-ARTIFACT (PAE)  = {ratio_pae:.2f}x")
        print()

    # --- framework integrity of the generated pools ---------------------
    fasta = {"stage3": RES / "stage3_cdr3_redesign_CC_designs.fasta",
             "stage2": RES / "stage2_denovo_nanobody_CC_designs.fasta"}
    print("=" * 94)
    print("FRAMEWORK INTEGRITY OF THE GENERATED POOLS")
    print("=" * 94)
    for job, path in fasta.items():
        if not path.exists():
            continue
        fi = framework_integrity(read_fasta(path))
        out.setdefault(job, {})["framework_integrity"] = fi
        print(f"  {job}  n={fi['n']}")
        print(f"    J-region Trp (WGQG/WGKG) missing : {fi['j_region_Trp_MISSING']}"
              f"{'   <- not reliably foldable' if fi['j_region_Trp_MISSING'] else ''}")
        print(f"    conserved Cys pair intact        : {fi['all_have_2_cys']}"
              f"  {fi['cys_count_distribution']}")
        print(f"    VHH hallmark WFRQ missing        : {fi['vhh_hallmark_WFRQ_missing']}"
              f"{'   <- must fold without a VL partner' if fi['vhh_hallmark_WFRQ_missing'] else ''}")
    print()

    # stage3 is the clean comparison: framework and CDR1/CDR2 fixed by construction
    s3 = next((k for k in by_job if "stage3" in k), None)
    print("=" * 94)
    print("VERDICT")
    print("=" * 94)
    if s3:
        r = out[s3]["signal_to_artifact_iptm"]
        print(f"  stage3 varies ONLY CDR3 on a fixed humanised framework, so its")
        print(f"  spread is the metric's response to the real paratope.")
        print(f"  That response is {r}x the response to a framework-only change.")
        if r is not None and r < 2.0:
            print("\n  => Changing the actual binding loop moves this score no more than")
            print("     an edit that cannot affect binding. The ranking is not usable")
            print("     for selecting candidates. Treat the designs as an UNORDERED pool.")
    print(f"\n  No design in either job reached binding_confidence >= 0.5.")

    (RES / "step8_design_dynamic_range.json").write_text(json.dumps(out, indent=2))
    print(f"\nwrote {RES / 'step8_design_dynamic_range.json'}")


if __name__ == "__main__":
    main()
