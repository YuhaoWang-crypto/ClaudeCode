"""
T2 - one-to-one orthologues and coding sequences across a divergence ladder.

Why this module exists
----------------------
Ensembl STOPPED publishing pairwise dN/dS: the REST `/homology/` endpoint now
returns `dn_ds: None` and no `*_homolog_dn` / `*_homolog_ds` attribute exists in
BioMart (release 116, checked). So the pipeline cannot look evolutionary rates
up - it has to compute them (T3). This module assembles the inputs:

  1. high-confidence ONE-TO-ONE orthologue pairs, human vs each ladder species
     (BioMart `homologs` attribute page, queried per chromosome so a single
     timeout cannot lose the whole species);
  2. the CDS FASTA for every species (Ensembl FTP).

Transcript choice: for each gene we keep the LONGEST protein-coding CDS whose
length is a multiple of 3. This is symmetric between the two species and needs
no extra queries; it is not always Ensembl's canonical transcript, which is why
T3 re-derives percent identity rather than trusting Compara's.

Only 1:1 orthologues are kept. One-to-many/many-to-many pairs would let gene
duplication - itself strongly tissue-biased - masquerade as a change in
substitution rate, which is exactly the confound this pipeline is built to
avoid.
"""
import os
import csv

from .common import (RAW, DERIVED, SPECIES, HUMAN_CHROMS, ENSEMBL_FTP,
                     biomart_query, download, http_get, log, read_fasta)

ORTHO_COLS = ["human_gene", "ortholog_gene", "orthology_type",
              "perc_id", "perc_id_r1", "confidence"]


def fetch_orthologs(sp, force=False):
    """BioMart 1:1 orthologue table for human vs `sp`, cached as TSV.

    Cached PER CHROMOSOME: BioMart is flaky enough that a single failure two
    hours into a four-species run must not cost the completed chromosomes.
    """
    out = os.path.join(RAW, f"orthologs_human_{sp['key']}.tsv")
    if os.path.exists(out) and not force:
        return out
    p = sp["biomart"]
    attrs = ["ensembl_gene_id",
             f"{p}_homolog_ensembl_gene",
             f"{p}_homolog_orthology_type",
             f"{p}_homolog_perc_id",
             f"{p}_homolog_perc_id_r1",
             f"{p}_homolog_orthology_confidence"]
    part_dir = os.path.join(RAW, "biomart_parts")
    os.makedirs(part_dir, exist_ok=True)

    rows = []
    for chrom in HUMAN_CHROMS:
        part = os.path.join(part_dir, f"{sp['key']}_chr{chrom}.tsv")
        if os.path.exists(part) and os.path.getsize(part) > 0:
            with open(part) as fh:
                txt = fh.read()
        else:
            txt = biomart_query("hsapiens_gene_ensembl", attrs,
                                filters=[("chromosome_name", chrom)])
            with open(part, "w") as fh:
                fh.write(txt)
        n_before = len(rows)
        for line in txt.rstrip("\n").split("\n")[1:]:      # drop header
            f = line.split("\t")
            if len(f) < 6 or not f[1]:
                continue
            rows.append(f[:6])
        log(f"  {sp['key']}: chr{chrom} +{len(rows) - n_before} "
            f"-> {len(rows)} cumulative homology rows")
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(ORTHO_COLS)
        w.writerows(rows)
    log(f"  {sp['key']}: wrote {len(rows)} rows -> {os.path.basename(out)}")
    return out


def _cds_url(ftp_name):
    """Resolve the `*.cds.all.fa.gz` filename by listing the FTP directory."""
    base = f"{ENSEMBL_FTP}/{ftp_name}/cds/"
    listing = http_get(base, timeout=300)
    for token in listing.split('"'):
        if token.endswith(".cds.all.fa.gz"):
            return base + token.split("/")[-1]
    raise RuntimeError(f"no cds.all.fa.gz found under {base}")


def fetch_cds(ftp_name):
    dest = os.path.join(RAW, f"cds_{ftp_name}.fa.gz")
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        return dest
    url = _cds_url(ftp_name)
    log(f"  downloading CDS {os.path.basename(url)}")
    return download(url, dest)


def index_longest_cds(fa_path):
    """gene id -> longest in-frame protein-coding CDS string."""
    best = {}
    for header, seq in read_fasta(fa_path):
        fields = header.split()
        biotype, gene = None, None
        for f in fields[1:]:
            if f.startswith("gene:"):
                gene = f.split(":", 1)[1].split(".")[0]
            elif f.startswith("transcript_biotype:"):
                biotype = f.split(":", 1)[1]
        if gene is None or biotype != "protein_coding":
            continue
        if len(seq) < 60 or len(seq) % 3 != 0:
            continue
        if gene not in best or len(seq) > len(best[gene]):
            best[gene] = seq
    return best


def build_pairs(sp):
    """Write a TSV of 1:1 orthologue pairs with both CDS sequences attached."""
    out = os.path.join(DERIVED, f"pairs_human_{sp['key']}.tsv")
    if os.path.exists(out) and os.path.getsize(out) > 0:
        log(f"  {sp['key']}: pairs already built")
        return out

    ortho = fetch_orthologs(sp)
    hs_cds = index_longest_cds(fetch_cds("homo_sapiens"))
    tg_cds = index_longest_cds(fetch_cds(sp["ftp"]))
    log(f"  {sp['key']}: CDS indexed human={len(hs_cds)} target={len(tg_cds)}")

    n_all = n_one2one = n_withseq = 0
    with open(ortho) as fh, open(out, "w", newline="") as ofh:
        rd = csv.DictReader(fh, delimiter="\t")
        w = csv.writer(ofh, delimiter="\t")
        w.writerow(["human_gene", "ortholog_gene", "perc_id", "confidence",
                    "human_cds", "ortholog_cds"])
        for r in rd:
            n_all += 1
            if r["orthology_type"] != "ortholog_one2one":
                continue
            n_one2one += 1
            h, t = hs_cds.get(r["human_gene"]), tg_cds.get(r["ortholog_gene"])
            if not h or not t:
                continue
            n_withseq += 1
            w.writerow([r["human_gene"], r["ortholog_gene"], r["perc_id"],
                        r["confidence"], h, t])
    log(f"  {sp['key']}: homology rows={n_all} 1:1={n_one2one} "
        f"with both CDS={n_withseq} -> {os.path.basename(out)}")
    return out


def main():
    print("=" * 72)
    print("T2 - orthologue pairs + CDS across the divergence ladder")
    print("=" * 72)
    for sp in SPECIES:
        log(f"{sp['key']} ({sp['latin']}, ~{sp['divergence_mya']} Mya)")
        build_pairs(sp)
    print("done")


if __name__ == "__main__":
    main()
