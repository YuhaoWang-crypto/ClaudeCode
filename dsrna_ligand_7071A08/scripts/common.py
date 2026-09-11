#!/usr/bin/env python3
"""Shared paths, config loading, FASTA IO and the IEDB REST clients."""
import os
import time
import urllib.parse
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _p(*parts):
    path = os.path.join(ROOT, *parts)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    return path


def data_path(*p):
    return _p("data", *p)


def results_path(*p):
    return _p("results", *p)


def figures_path(*p):
    return _p("figures", *p)


def load_config(path=None):
    import yaml
    with open(path or os.path.join(ROOT, "config", "config.yaml")) as f:
        return yaml.safe_load(f)


def read_fasta(path):
    """Ordered {id: sequence}."""
    seqs, name, buf = {}, None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line)
    if name:
        seqs[name] = "".join(buf)
    return seqs


def read_metadata():
    """sequences_metadata.tsv -> {id: {col: value}}"""
    path = data_path("sequences_metadata.tsv")
    rows, header = {}, None
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if header is None:
                header = parts
                continue
            rows[parts[0]] = dict(zip(header, parts))
    return rows


# --------------------------------------------------------------------------
# Population coverage
# --------------------------------------------------------------------------
class CoverageModel:
    """
    Phenotypic coverage of an HLA-DR panel, matching the IEDB population
    coverage tool exactly.

    Single-locus Hardy-Weinberg over the IEDB DRB1 allele-frequency tables:

        cov = 1 - (1 - f)^2,  f = sum(panel allele frequencies) / N

    N is the population's *total* DRB1 allele frequency when that total
    exceeds 1, and 1 otherwise. The renormalisation matters: several curated
    populations sum above 1 (Europe = 1.040), and skipping it overstates
    coverage - Europe by ~0.8 points on a 20-allele panel. Verified against
    `calculate_population_coverage.py -c II` to two decimals.

    DRB3/4/5 carry no frequencies in these tables and therefore contribute
    nothing here; they are reported separately as presentation breadth.
    """

    def __init__(self, tables, weights):
        self.tables = tables
        self.weights = weights
        self.norm = {p: max(sum(t.values()), 1.0) for p, t in tables.items()}

    def population(self, panel, pop):
        drb1 = [a for a in panel if a.startswith("HLA-DRB1")]
        f = min(sum(self.tables[pop].get(a, 0.0) for a in drb1) / self.norm[pop], 1.0)
        return 1.0 - (1.0 - f) ** 2

    def weighted(self, panel):
        return sum(w * self.population(panel, p) for p, w in self.weights.items())


def population_weights(cfg):
    """US composite (census-weighted) + Europe, combined per config."""
    us_w, eu_w = cfg["panel"]["us_eu_split"]
    tot = sum(cfg["panel"]["us_weights"].values())
    w = {p: us_w * v / tot for p, v in cfg["panel"]["us_weights"].items()}
    eu = cfg["panel"]["eu_population"]
    w[eu] = w.get(eu, 0.0) + eu_w
    return w


# --------------------------------------------------------------------------
# IEDB REST clients
# --------------------------------------------------------------------------
def _post(url, fields, timeout=600, tries=4):
    body = urllib.parse.urlencode(fields).encode()
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(url, data=body,
                                         headers={"Content-Type":
                                                  "application/x-www-form-urlencoded"})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read().decode("utf-8")
        except Exception as e:
            last = e
            time.sleep(3 * (i + 1))
    raise RuntimeError(f"IEDB request failed after {tries} tries: {last}")


def iedb_mhcii(endpoint, method, sequences, alleles, length=15, timeout=600):
    """
    sequences: {id: seq}. alleles: list of HLA-DR molecule names.
    Returns a list of dicts with the raw IEDB columns plus 'id'.

    The IEDB endpoint takes parallel comma-separated allele/length lists and a
    multi-record FASTA, so one call covers every sequence for a block of
    alleles. seq_num in the response is the 1-based FASTA record index.
    """
    ids = list(sequences)
    fasta = "".join(f">{i}\n{sequences[i]}\n" for i in ids)
    txt = _post(endpoint, {
        "method": method,
        "sequence_text": fasta,
        "allele": ",".join(alleles),
        "length": ",".join([str(length)] * len(alleles)),
    }, timeout=timeout)

    lines = [l for l in txt.strip().split("\n") if l.strip()]
    if not lines or "\t" not in lines[0]:
        raise RuntimeError(f"unexpected IEDB response: {txt[:400]}")
    header = lines[0].split("\t")
    out = []
    for line in lines[1:]:
        parts = line.split("\t")
        if len(parts) != len(header) or parts[0] == "allele":
            continue
        rec = dict(zip(header, parts))
        rec["id"] = ids[int(rec["seq_num"]) - 1]
        out.append(rec)
    return out


def iedb_bcell(endpoint, method, sequence, timeout=300):
    """Per-residue B-cell epitope propensity. Returns [(pos, residue, score)]."""
    txt = _post(endpoint, {"method": method, "sequence_text": sequence},
                timeout=timeout)
    rows = []
    for line in txt.strip().split("\n")[1:]:
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        try:
            rows.append((int(parts[0]), parts[1], float(parts[2])))
        except ValueError:
            continue
    return rows
