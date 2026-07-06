#!/usr/bin/env python3
"""Off-target + on-target scoring for a fixed guide list, self-contained via FlashFry.

FlashFry (a single JAR, needs only Java + a genome FASTA) builds its own index and
computes, per guide: Doench 2016 on-target (Rule Set 2), CFD off-target specificity,
and off-target counts. This wraps it around OUR fixed kinome library so the two
pending QC fields in data/*guides*.tsv get filled — the exact step the Latch
workflow runs.

Subcommands
-----------
  prep   build a FlashFry-ready FASTA of target windows from rs3_context.tsv
         (testable without any genome/JAR)
  run    full pipeline: index -> discover -> score -> merge -> guides_scored.tsv
         (needs java, FlashFry.jar, and a genome FASTA)

`run` shells out to FlashFry; `prep` and `merge` are pure Python.
"""
import argparse, csv, os, subprocess, sys, tempfile


def read_context(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_targets_fasta(context_path, fasta_path):
    """One FASTA record per guide: the 30-mer context. FlashFry rediscovers the
    guide inside and scores it; we later match back by spacer."""
    rows = read_context(context_path)
    n = 0
    with open(fasta_path, "w") as fo:
        for r in rows:
            ctx = r["context_30mer"].strip().upper()
            if len(ctx) == 30:
                fo.write(f">{r['id']}|{r['gene']}|{r['spacer']}\n{ctx}\n")
                n += 1
    print(f"prep: wrote {n} target windows -> {fasta_path}")
    return n


def _f(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return float("nan")


def parse_flashfry(scored_path):
    """Map spacer(20nt) -> dict(on_target, cfd_specificity, offtargets).

    FlashFry `score` output is a TSV; column names differ slightly across
    versions, so detect by fuzzy header match. `target` holds the 23-mer
    (20 protospacer + 3 PAM)."""
    out = {}
    with open(scored_path) as f:
        r = csv.DictReader(f, delimiter="\t")
        cols = {c.lower(): c for c in (r.fieldnames or [])}

        def pick(*subs):
            for want in subs:
                for lc, orig in cols.items():
                    if want in lc:
                        return orig
            return None

        c_target = pick("target")
        c_on = pick("doench2016", "doenchon", "ruleset")
        c_spec = pick("cfd_specificity", "specificity", "hsu2013")
        c_ot = pick("otcount", "offtarget", "ot_count")
        if not c_target:
            sys.exit(f"FlashFry output missing a target column; saw {r.fieldnames}")
        for row in r:
            tgt = (row[c_target] or "").strip().upper()
            spacer = tgt[:20]
            out[spacer] = dict(
                on_target=_f(row.get(c_on)) if c_on else float("nan"),
                cfd_specificity=_f(row.get(c_spec)) if c_spec else float("nan"),
                offtargets=(row.get(c_ot, "") if c_ot else ""),
            )
    return out


def merge(context_path, flashfry_scored, out_path, cfd_min=0.2):
    rows = read_context(context_path)
    sc = parse_flashfry(flashfry_scored)
    kept = dropped = missing = 0
    with open(out_path, "w", newline="") as fo:
        w = csv.writer(fo, delimiter="\t")
        w.writerow(["id", "gene", "spacer", "pam", "source",
                    "ontarget_doench2016", "cfd_specificity", "offtargets", "pass_filter"])
        for r in rows:
            sp = r["spacer"]
            s = sc.get(sp)
            if not s:
                missing += 1
                spec = float("nan"); on = float("nan"); ot = ""
            else:
                spec, on, ot = s["cfd_specificity"], s["on_target"], s["offtargets"]
            passed = not (spec == spec) or spec >= cfd_min
            kept += int(passed); dropped += int(not passed)
            w.writerow([r["id"], r["gene"], sp, r["pam"], r["source"],
                        "" if on != on else f"{on:.3f}",
                        "" if spec != spec else f"{spec:.3f}", ot,
                        "PASS" if passed else "FAIL"])
    print(f"merge: {len(rows)} guides | pass {kept} fail {dropped} | "
          f"{missing} not found in FlashFry output (CFD >= {cfd_min})")


def run(genome, context_path, out_path, ff_jar, cfd_min, workdir):
    os.makedirs(workdir, exist_ok=True)
    db = os.path.join(workdir, "ff_db")
    fasta = os.path.join(workdir, "targets.fa")
    disc = os.path.join(workdir, "discover.tsv")
    scored = os.path.join(workdir, "scored.tsv")
    write_targets_fasta(context_path, fasta)
    if not os.path.exists(db):
        subprocess.run(["java", "-Xmx8g", "-jar", ff_jar, "index",
                        "--tmpLocation", os.path.join(workdir, "tmp"),
                        "--database", db, "--reference", genome,
                        "--enzyme", "spcas9ngg"], check=True)
    subprocess.run(["java", "-Xmx8g", "-jar", ff_jar, "discover",
                    "--database", db, "--fasta", fasta, "--output", disc], check=True)
    subprocess.run(["java", "-Xmx8g", "-jar", ff_jar, "score",
                    "--input", disc, "--output", scored, "--database", db,
                    "--scoringMetrics", "doench2016cfd,hsu2013,dangerous"], check=True)
    merge(context_path, scored, out_path, cfd_min)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("prep"); a.add_argument("--context", required=True); a.add_argument("--fasta", required=True)
    b = sub.add_parser("run")
    b.add_argument("--genome", required=True); b.add_argument("--context", required=True)
    b.add_argument("--out", required=True); b.add_argument("--ff", default="/opt/FlashFry.jar")
    b.add_argument("--cfd-min", type=float, default=0.2); b.add_argument("--workdir", default="/tmp/ff")
    m = sub.add_parser("merge")
    m.add_argument("--context", required=True); m.add_argument("--scored", required=True)
    m.add_argument("--out", required=True); m.add_argument("--cfd-min", type=float, default=0.2)
    args = p.parse_args()
    if args.cmd == "prep":
        write_targets_fasta(args.context, args.fasta)
    elif args.cmd == "run":
        run(args.genome, args.context, args.out, args.ff, args.cfd_min, args.workdir)
    elif args.cmd == "merge":
        merge(args.context, args.scored, args.out, args.cfd_min)


if __name__ == "__main__":
    main()
