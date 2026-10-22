"""CLI:  python3 -m gismo <command>

  lookup GENE [GENE ...]   scores + interpretation (offline; LOEUF if fetched)
  fetch                    download Zenodo inputs + reference data to $GISMO_DATA
  build [N_ITER]           recompute gismo/data/gismo_scores.tsv.gz (default 2000 sims)
  benchmark                reproduce the paper's numbers -> gismo/data/benchmark.json
"""
import json
import sys


def main(argv):
    if not argv or argv[0] in ("-h", "--help"):
        print(__doc__)
        return
    cmd, rest = argv[0], argv[1:]
    if cmd == "lookup":
        from gismo.interpret import fmt, lookup
        print(fmt(lookup(rest)))
    elif cmd == "fetch":
        from gismo.paths import fetch_all
        print(fetch_all())
    elif cmd == "build":
        from gismo.core import build
        build(n_iter=int(rest[0]) if rest else 2000)
    elif cmd == "benchmark":
        from gismo.benchmark import run, table
        r = run()
        print(table(r).to_string(index=False))
        print(json.dumps({k: r[k] for k in ["recessive_rule", "x_chrom"]}, default=float))
    else:
        print(__doc__)
        sys.exit(1)


if __name__ == "__main__":
    main(sys.argv[1:])
