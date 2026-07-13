"""
cli.py — one-line entry points for the common omics/pathway tasks.

    python -m ppi_pathway.cli download --source all
    python -m ppi_pathway.cli enrich   genes.txt            # STRING API, online
    python -m ppi_pathway.cli pathways genes.txt --fdr 0.05
    python -m ppi_pathway.cli subnet   genes.txt --expand 1 --min-score 700 \
                                        --out module          # offline, needs cache

``genes.txt`` is one gene symbol per line (blank lines / '#' comments ignored).
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from . import HUMAN_TAXON


def _read_genes(path: str) -> list[str]:
    out = []
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            out.append(line.split()[0])
    return out


def _cmd_download(a):
    from . import download
    if a.source in ("string", "all"):
        download.string_files(taxon=a.taxon, force=a.force)
    if a.source in ("biogrid", "all"):
        download.biogrid_file(taxon=a.taxon, force=a.force)


def _cmd_enrich(a):
    from . import string_api
    genes = _read_genes(a.genes)
    print(f"# {len(genes)} genes -> STRING functional enrichment (taxon {a.taxon})",
          file=sys.stderr)
    ppi = string_api.ppi_enrichment(genes, taxon=a.taxon)
    if ppi:
        print(f"# PPI enrichment p-value = {ppi.get('p_value')} "
              f"(expected edges {ppi.get('expected_number_of_edges')}, "
              f"observed {ppi.get('number_of_edges')})", file=sys.stderr)
    rows = string_api.enrichment(genes, taxon=a.taxon)
    rows = [r for r in rows if float(r.get("fdr", 1)) <= a.fdr]
    if a.json:
        print(json.dumps(rows, indent=2))
    else:
        print("category\tfdr\tgenes\tdescription")
        for r in sorted(rows, key=lambda r: float(r["fdr"])):
            print(f"{r['category']}\t{r['fdr']:.2e}\t{r['number_of_genes']}\t{r['description']}")


def _cmd_pathways(a):
    from . import string_api
    genes = _read_genes(a.genes)
    rows = string_api.top_pathways(genes, taxon=a.taxon, fdr_max=a.fdr, n=a.n)
    print("category\tfdr\tgenes\tdescription")
    for r in rows:
        print(f"{r['category']}\t{r['fdr']:.2e}\t{r['number_of_genes']}\t{r['description']}")


def _cmd_annotate(a):
    from . import interpro
    accs = _read_genes(a.accessions)  # UniProt accessions, one per line
    print(f"# InterPro annotation for {len(accs)} proteins", file=sys.stderr)
    print("uniprot\tinterpro\ttype\tname\tgo_terms")
    for acc in accs:
        try:
            for e in interpro.entries_for_protein(acc):
                gos = ",".join(g["id"] for g in e["go_terms"][:4])
                print(f"{acc}\t{e['accession']}\t{e['type']}\t{e['name']}\t{gos}")
        except Exception as e:
            print(f"{acc}\tERROR\t\t{e}\t", file=sys.stderr)


def _cmd_reactome(a):
    from . import reactome, string_api
    genes = _read_genes(a.genes)
    try:
        rows = reactome.top_pathways(genes, fdr_max=a.fdr, n=a.n)
        print("fdr\tfound/total\tstId\tname")
        for r in rows:
            print(f"{r['fdr']:.2e}\t{r['found']}/{r['total']}\t{r['stId']}\t{r['name']}")
    except reactome.ReactomeBlocked as e:
        print(f"# {e}", file=sys.stderr)
        print("# falling back to STRING's Reactome-category enrichment:", file=sys.stderr)
        for r in string_api.top_pathways(genes, categories=("RCTM",),
                                         fdr_max=a.fdr, n=a.n):
            print(f"{r['fdr']:.2e}\t{r['number_of_genes']}\t{r['term']}\t{r['description']}")


def _cmd_subnet(a):
    from . import network, subnetwork
    genes = _read_genes(a.genes)
    print(f"# loading STRING (min_score={a.min_score}) ...", file=sys.stderr)
    g = network.load_string(taxon=a.taxon, min_score=a.min_score)
    print(f"# background graph: {g.number_of_nodes():,} nodes / "
          f"{g.number_of_edges():,} edges", file=sys.stderr)
    sub = subnetwork.induced(g, genes, expand=a.expand)
    print("# " + json.dumps(subnetwork.summary(sub)), file=sys.stderr)
    if a.out:
        paths = subnetwork.write_cytoscape(sub, a.out)
        print(f"# wrote {paths['edges']} and {paths['nodes']}", file=sys.stderr)
    print("gene\tis_seed\tdegree\tbetweenness\tclustering")
    for r in subnetwork.rank_genes(sub)[:a.top]:
        print(f"{r['gene']}\t{r['is_seed']}\t{r['degree']}\t{r['betweenness']}\t{r['clustering']}")


def main(argv=None):
    ap = argparse.ArgumentParser(prog="ppi_pathway", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--taxon", type=int, default=HUMAN_TAXON)
    sub = ap.add_subparsers(dest="cmd", required=True)

    d = sub.add_parser("download", help="cache STRING/BioGRID bulk files")
    d.add_argument("--source", choices=["string", "biogrid", "all"], default="all")
    d.add_argument("--force", action="store_true")
    d.set_defaults(func=_cmd_download)

    e = sub.add_parser("enrich", help="STRING functional enrichment (online)")
    e.add_argument("genes")
    e.add_argument("--fdr", type=float, default=0.05)
    e.add_argument("--json", action="store_true")
    e.set_defaults(func=_cmd_enrich)

    p = sub.add_parser("pathways", help="pathway-only enrichment (KEGG/Reactome/WikiPathways)")
    p.add_argument("genes")
    p.add_argument("--fdr", type=float, default=0.05)
    p.add_argument("--n", type=int, default=25)
    p.set_defaults(func=_cmd_pathways)

    an = sub.add_parser("annotate", help="InterPro domains/families/GO for UniProt accessions")
    an.add_argument("accessions", help="file of UniProt accessions, one per line")
    an.set_defaults(func=_cmd_annotate)

    rx = sub.add_parser("reactome", help="Reactome pathway enrichment (falls back to STRING if blocked)")
    rx.add_argument("genes")
    rx.add_argument("--fdr", type=float, default=0.05)
    rx.add_argument("--n", type=int, default=25)
    rx.set_defaults(func=_cmd_reactome)

    s = sub.add_parser("subnet", help="offline PPI subnetwork + hub ranking")
    s.add_argument("genes")
    s.add_argument("--min-score", type=int, default=700)
    s.add_argument("--expand", type=int, default=0)
    s.add_argument("--out", default=None, help="write <out>.edges.tsv / .nodes.tsv")
    s.add_argument("--top", type=int, default=30)
    s.set_defaults(func=_cmd_subnet)

    args = ap.parse_args(argv)
    return args.func(args) or 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except BrokenPipeError:
        # downstream (e.g. `| head`) closed the pipe — exit quietly
        sys.exit(0)
