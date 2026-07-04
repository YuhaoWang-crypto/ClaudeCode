#!/usr/bin/env python3
"""Finish the two pending guide-QC steps: Rule Set 3 on-target scoring and
off-target specificity, then merge into a single ranked table.

This runs anywhere the prerequisites exist (a Latch task, or any machine with
an hg38 index). It is deliberately tool-agnostic at the boundaries:

  * on-target  -> the `rs3` package (Broad Rule Set 3), scoring the 30-mer
                  contexts prepared in data/rs3_context.tsv
  * off-target -> a specificity table produced by GuideScan2 (preferred) or
                  CRISPOR, keyed by spacer; this script only *parses & merges*
                  it, so we never hand-roll a CFD matrix.

Inputs (already generated, committed under data/):
  data/rs3_context.tsv   id, gene, spacer, pam, context_30mer, source
  data/all_spacers.txt   one 20-nt spacer per line (for the off-target tool)

Usage:
  # 1. score on-target (needs: pip install rs3)
  python3 scripts/score_guides.py ontarget \
      --context data/rs3_context.tsv --out data/ontarget_rs3.tsv

  # 2. after running GuideScan2/CRISPOR on data/all_spacers.txt, merge:
  python3 scripts/score_guides.py merge \
      --context data/rs3_context.tsv \
      --ontarget data/ontarget_rs3.tsv \
      --offtarget data/offtarget_guidescan.tsv \
      --out data/guides_scored.tsv

  # verify the plumbing without any genome/tool:
  python3 scripts/score_guides.py selftest
"""
import argparse, csv, sys


def read_context(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def score_ontarget(context_path, out_path):
    """Rule Set 3 sequence-based scoring of the 30-mer contexts."""
    rows = read_context(context_path)
    contexts = [r["context_30mer"] for r in rows]
    bad = [r["id"] for r in rows if len(r["context_30mer"]) != 30]
    if bad:
        sys.exit(f"{len(bad)} guides lack a 30-mer context (e.g. {bad[:3]}) — regenerate inputs")
    try:
        from rs3.seq import predict_seq
    except ImportError:
        sys.exit("rs3 not installed. Run:  pip install rs3   (Broad Rule Set 3)")
    scores = predict_seq(contexts, sequence_tracr="Hsu2013")
    with open(out_path, "w", newline="") as fo:
        w = csv.writer(fo, delimiter="\t")
        w.writerow(["id", "spacer", "rs3_score"])
        for r, s in zip(rows, scores):
            w.writerow([r["id"], r["spacer"], f"{s:.4f}"])
    print(f"wrote {out_path}: {len(rows)} Rule Set 3 scores "
          f"(median {sorted(scores)[len(scores)//2]:.3f})")


def parse_offtarget(path):
    """Parse a spacer-keyed off-target/specificity table (GuideScan2 or CRISPOR).

    Accepts any TSV/CSV with a spacer/sequence column and a specificity column;
    common header names are auto-detected. Returns spacer -> (specificity, n_offtargets).
    """
    delim = "\t" if path.endswith((".tsv", ".txt")) else ","
    out = {}
    with open(path) as f:
        r = csv.DictReader(f, delimiter=delim)
        cols = {c.lower(): c for c in (r.fieldnames or [])}
        sp_c = next((cols[k] for k in ("spacer", "sequence", "grna", "sgrna", "gRNA".lower()) if k in cols), None)
        spec_c = next((cols[k] for k in ("specificity", "cfd_specificity", "specificity_score", "guidescan2_specificity") if k in cols), None)
        n_c = next((cols[k] for k in ("n_offtargets", "offtargets", "num_offtargets", "n_off") if k in cols), None)
        if not sp_c or not spec_c:
            sys.exit(f"off-target file needs a spacer column and a specificity column; saw {r.fieldnames}")
        for row in r:
            sp = row[sp_c].strip().upper()[:20]
            try:
                spec = float(row[spec_c])
            except (ValueError, TypeError):
                spec = float("nan")
            n = row.get(n_c, "") if n_c else ""
            out[sp] = (spec, n)
    return out


def merge(context_path, ontarget_path, offtarget_path, out_path, cfd_min=0.2):
    rows = read_context(context_path)
    on = {}
    if ontarget_path:
        with open(ontarget_path) as f:
            for r in csv.DictReader(f, delimiter="\t"):
                on[r["id"]] = r["rs3_score"]
    off = parse_offtarget(offtarget_path) if offtarget_path else {}
    kept = dropped = 0
    with open(out_path, "w", newline="") as fo:
        w = csv.writer(fo, delimiter="\t")
        w.writerow(["id", "gene", "spacer", "pam", "source",
                    "rs3_score", "cfd_specificity", "n_offtargets", "pass_filter"])
        for r in rows:
            sp = r["spacer"]
            spec, n = off.get(sp, (float("nan"), ""))
            passed = not (spec == spec) or spec >= cfd_min   # keep if unknown or >= threshold
            kept += int(passed); dropped += int(not passed)
            w.writerow([r["id"], r["gene"], sp, r["pam"], r["source"],
                        on.get(r["id"], "NA"),
                        "" if spec != spec else f"{spec:.3f}", n,
                        "PASS" if passed else "FAIL"])
    print(f"wrote {out_path}: {len(rows)} guides | pass {kept} / fail {dropped} "
          f"(CFD specificity >= {cfd_min})")


def selftest():
    """Exercise parsing/merging on tiny synthetic data — no genome or rs3 needed."""
    import tempfile, os
    d = tempfile.mkdtemp()
    ctx = os.path.join(d, "ctx.tsv")
    with open(ctx, "w") as f:
        f.write("id\tgene\tspacer\tpam\tcontext_30mer\tsource\n")
        f.write("G1_1\tGENE1\tGAATACGAGGGCTGCAAAGT\tTGG\tACAGTGGCAAAATCATCACTACGAAGGCAG\tBrunello\n")
        f.write("G1_2\tGENE1\tACCTCAGGAGCAATGTAGTC\tAGG\tCATGAACACCTCAGGAGCAATGTAGTCAGG\tde-novo\n")
    ont = os.path.join(d, "on.tsv")
    with open(ont, "w") as f:
        f.write("id\tspacer\trs3_score\nG1_1\tGAATACGAGGGCTGCAAAGT\t0.71\nG1_2\tACCTCAGGAGCAATGTAGTC\t0.55\n")
    offt = os.path.join(d, "off.tsv")
    with open(offt, "w") as f:
        f.write("spacer\tspecificity\tn_offtargets\n")
        f.write("GAATACGAGGGCTGCAAAGT\t0.82\t3\n")
        f.write("ACCTCAGGAGCAATGTAGTC\t0.11\t altrue\n".replace(" altrue", "17"))
    out = os.path.join(d, "scored.tsv")
    merge(ctx, ont, offt, out)
    print("--- merged output ---")
    print(open(out).read().rstrip())
    # assert the low-specificity guide is failed
    body = [l.split("\t") for l in open(out).read().splitlines()[1:]]
    assert body[0][-1] == "PASS" and body[1][-1] == "FAIL", "filter logic broken"
    print("\nselftest OK")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("ontarget"); a.add_argument("--context", required=True); a.add_argument("--out", required=True)
    m = sub.add_parser("merge")
    m.add_argument("--context", required=True); m.add_argument("--ontarget", default=None)
    m.add_argument("--offtarget", required=True); m.add_argument("--out", required=True)
    m.add_argument("--cfd-min", type=float, default=0.2)
    sub.add_parser("selftest")
    args = p.parse_args()
    if args.cmd == "ontarget": score_ontarget(args.context, args.out)
    elif args.cmd == "merge": merge(args.context, args.ontarget, args.offtarget, args.out, args.cfd_min)
    elif args.cmd == "selftest": selftest()


if __name__ == "__main__":
    main()
