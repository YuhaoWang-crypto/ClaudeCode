#!/usr/bin/env python3
"""
seed_from_db.py -- Stage-0 SEED lookup for the aptamer pipeline.

Instead of cold-starting design from a random ViennaRNA pool, look up KNOWN aptamers
against the target (or a homolog family) and use them as starting sequences ("seeds").
A seed that survives the specificity counter-screen becomes a maturation parent; a seed
that FAILS becomes a negative anchor that calibrates the score scale. Either way the
literature buys you a warmer start than random.

Sources (see schema.md for access notes, verified 2026):
  - bundled  : model/seed_index.csv       (curated, offline, always available)
  - utexas   : UTexas Aptamer DB public download  (--live; best-effort HTTP)
  - apta-index / aptamer-base : web-search / MediaWiki only, NO bulk API -> spot-check by hand
The bundled corpus is the reliable default. --live only augments it.

Matching tiers (most to least specific):
  1. DIRECT   : same target (name substring or exact UniProt accession)
  2. FAMILY   : homolog / same target_family keyword  -> related scaffold, still likely de novo
  3. NONE     : no seed  -> fall back to de novo design (report this honestly)

Usage:
  python3 seed_from_db.py --target "Thrombin"
  python3 seed_from_db.py --target "alpha-1-acid glycoprotein" --uniprot P02763 \
                          --family "lipocalin,orosomucoid,ORM1,ORM2,acute-phase" --json
  python3 seed_from_db.py --calibrate
"""
import argparse, csv, json, os, sys, urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
SEED_INDEX = os.path.join(HERE, "seed_index.csv")
CAL_ANCHORS = os.path.join(HERE, "calibration_anchors.csv")
UTEXAS_HINT = "https://sites.utexas.edu/aptamerdatabase/"


def load_csv(path):
    if not os.path.exists(path):
        return []
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def norm(s):
    return (s or "").strip().lower()


def keywords(csv_str):
    return [norm(k) for k in (csv_str or "").split(",") if k.strip()]


def match_seeds(rows, target, uniprot, family_kw):
    """Return (direct, family) lists; a row lands in the most specific tier it qualifies for."""
    t = norm(target)
    up = norm(uniprot)
    fam = family_kw or []
    direct, family = [], []
    for r in rows:
        tname, tfam = norm(r.get("target_name")), norm(r.get("target_family"))
        tup, note = norm(r.get("target_uniprot")), norm(r.get("note"))
        is_direct = False
        if up and tup and up == tup:
            is_direct = True
        elif t and tname and (t in tname or tname in t):
            is_direct = True
        if is_direct:
            direct.append(r)
            continue
        hay = " ".join([tname, tfam, note])
        if any(k and k in hay for k in fam):
            family.append(r)
    return direct, family


def try_live_utexas():
    """Best-effort: the UTexas DB is a public *download*, not a REST API. We attempt a
    fetch so a provisioned env can augment the bundled corpus; failure is non-fatal."""
    try:
        req = urllib.request.Request(UTEXAS_HINT, headers={"User-Agent": "aptamer-seed/1.0"})
        with urllib.request.urlopen(req, timeout=15) as resp:
            html = resp.read(4096).decode("utf-8", "replace")
        return True, f"reachable ({len(html)}+ bytes); locate the CSV/XLSX export link and " \
                     f"ingest with ingest.py --preset utexas"
    except Exception as e:  # noqa: BLE001 - network is optional
        return False, f"unreachable ({type(e).__name__}: {e}). Use the bundled corpus."


def flags_for(r):
    f = []
    if norm(r.get("seq_verified")) == "no" or not r.get("sequence"):
        f.append("METADATA-ONLY (no usable sequence - verify from DOI before seeding)")
    elif norm(r.get("seq_verified")) == "partial":
        f.append("SEQUENCE-PARTIAL (core only; mods/full construct from source)")
    if "promiscuous" in norm(r.get("note")) or "negative anchor" in norm(r.get("note")):
        f.append("NEGATIVE-ANCHOR (known non-specific - use to calibrate, do NOT seed as a hit)")
    if r.get("pdb_anchor"):
        f.append(f"PDB-ANCHOR ({r['pdb_anchor']}) - usable for absolute ipTM calibration")
    return f


def emit_calibration():
    anchors = load_csv(CAL_ANCHORS)
    print("=" * 74)
    print("CALIBRATION ANCHORS  (turn Boltz ipTM into a meaningful scale)")
    print("=" * 74)
    print(
        "Absolute ipTM is not comparable across targets. To read a new candidate's\n"
        "decoy gap, first run a VALIDATED binder through the SAME co-fold pipeline and\n"
        "record what it scores. That reference gap is your yardstick.\n"
    )
    for a in anchors:
        print(f"- {a['name']}  [{a['chem']}]  Kd~{a['kd_nM']} nM  PDB:{a['pdb'] or 'n/a'}")
        print(f"    seq: {a['sequence'] or '(from source)'}")
        print(f"    {a['note']}\n")
    print("PROCEDURE")
    print("  1. Co-fold anchor x its target  -> reference on-target ipTM  (ref_on)")
    print("  2. Co-fold scramble(anchor) x same target -> reference decoy ipTM  (ref_decoy)")
    print("  3. ref_gap = ref_on - ref_decoy  is the gap a REAL ~100 nM binder produces here.")
    print("  4. Judge each design's (on - decoy) gap against ref_gap, not against 0.")
    print("  NOTE: high-res aptamer x protein co-crystals are rare; TBA x thrombin is the")
    print("        most solid. Treat the anchor as an order-of-magnitude yardstick.")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--target", default="", help="target protein name (substring match)")
    ap.add_argument("--uniprot", default="", help="target UniProt accession (exact match)")
    ap.add_argument("--family", default="",
                    help="comma-separated homolog/family keywords for tier-2 fallback")
    ap.add_argument("--live", action="store_true", help="also probe the UTexas public download")
    ap.add_argument("--calibrate", action="store_true", help="print calibration anchors and exit")
    ap.add_argument("--json", action="store_true", help="machine-readable output")
    args = ap.parse_args()

    if args.calibrate:
        emit_calibration()
        return

    if not (args.target or args.uniprot):
        ap.error("provide --target and/or --uniprot (or use --calibrate)")

    rows = load_csv(SEED_INDEX)
    direct, family = match_seeds(rows, args.target, args.uniprot, keywords(args.family))

    live_msg = ""
    if args.live:
        ok, live_msg = try_live_utexas()
        live_msg = f"[live UTexas] {live_msg}"

    # verdict
    seedable = [r for r in direct if r.get("sequence") and norm(r.get("seq_verified")) != "no"
                and "NEGATIVE-ANCHOR" not in " ".join(flags_for(r))]
    if direct and seedable:
        verdict = "DIRECT_SEED"
    elif direct:
        verdict = "DIRECT_METADATA_ONLY"   # target has known aptamers but no usable seq / negative anchor
    elif family:
        verdict = "FAMILY_ONLY"
    else:
        verdict = "NO_SEED"

    if args.json:
        out = {
            "query": {"target": args.target, "uniprot": args.uniprot,
                      "family": keywords(args.family)},
            "verdict": verdict,
            "direct": [dict(r, flags=flags_for(r)) for r in direct],
            "family": [dict(r, flags=flags_for(r)) for r in family],
            "live": live_msg,
            "recommendation": RECS[verdict],
        }
        print(json.dumps(out, indent=2))
        return

    print("=" * 74)
    print(f"SEED LOOKUP  target={args.target or '-'}  uniprot={args.uniprot or '-'}")
    print("=" * 74)
    if live_msg:
        print(live_msg + "\n")

    def show(tier, items):
        if not items:
            return
        print(f"\n[{tier}]  {len(items)} hit(s)")
        for r in items:
            print(f"  * {r['name']}  ({r['chem']}, {r['target_name']})")
            if r.get("sequence"):
                print(f"      seq: {r['sequence']}  ({r.get('length','?')} nt)")
            if r.get("doi"):
                print(f"      doi: {r['doi']}")
            for fl in flags_for(r):
                print(f"      ! {fl}")

    show("DIRECT (same target)", direct)
    show("FAMILY (homolog / same fold family)", family)

    print("\n" + "-" * 74)
    print(f"VERDICT: {verdict}")
    print(RECS[verdict])


RECS = {
    "DIRECT_SEED":
        "Usable literature seed(s) for this exact target. -> Seed the design pool with\n"
        "  these sequences, run each through the specificity counter-screen FIRST (a known\n"
        "  binder can still be promiscuous), then dope-mutate survivors. Skip cold-start.",
    "DIRECT_METADATA_ONLY":
        "Known aptamer(s) exist for this target but none is a usable seed here (no\n"
        "  verified sequence, or flagged as a NEGATIVE anchor). -> Retrieve the sequence\n"
        "  from the DOI, OR treat the negative anchor as a score-scale calibrator and\n"
        "  proceed with de novo design + the specificity gate.",
    "FAMILY_ONLY":
        "No aptamer for your exact target, but the fold/family has aptamers. -> Their\n"
        "  scaffolds are weak priors at best; run de novo design and use the family\n"
        "  aptamer only as a sanity anchor. Do NOT present it as a hit for your target.",
    "NO_SEED":
        "No literature aptamer for this target OR its family in the corpus. -> This is a\n"
        "  genuine cold-start: fall back to de novo (ViennaRNA pool + Boltz gate). Report\n"
        "  the absence honestly - it means there is no validated prior to lean on.",
}


if __name__ == "__main__":
    main()
