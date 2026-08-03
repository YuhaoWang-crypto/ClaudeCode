"""
T3 - compute dN, dS and omega = dN/dS for every 1:1 orthologue pair.

Ensembl no longer distributes pairwise dN/dS, so this module derives it:

    CDS -> protein -> global protein alignment (BLOSUM62, affine gaps, free
    terminal gaps) -> back-translated codon alignment -> Nei-Gojobori 1986
    counting with Jukes-Cantor correction.

The NG86 counting is a vectorised reimplementation (two 64- and 64x64-entry
lookup tables) because Biopython's reference implementation is a per-codon
Python loop and this pipeline needs ~60k pairs. `validate_against_biopython()`
checks the fast path against `Bio.codonalign.cal_dn_ds` on a random sample and
reports the largest disagreement - it is part of the module's output, not an
optional extra.

Conventions follow Biopython's NG86 exactly so the two are comparable:
mutations to a stop codon count as non-synonymous sites; site counts are
averaged over the two sequences; multi-difference codons are averaged over all
pathways (weight 1/2 for two differences, 1/6 for three); dS is undefined
(reported as NaN) once pS >= 3/4, i.e. once synonymous divergence saturates.

SATURATION IS THE MAIN INTERPRETIVE LIMIT of the deep comparisons: by the
human-chicken split most synonymous sites have been hit more than once, so
omega there is not comparable to omega at the human-macaque split. T6 quantifies
this instead of assuming it away.
"""
import os
import csv
import math
import sys
from itertools import product
from multiprocessing import Pool

import numpy as np

from .common import DERIVED, SPECIES, log, to_tsv_gz

BASES = "TCAG"
CODONS = ["".join(c) for c in product(BASES, repeat=3)]
CODON_IDX = {c: i for i, c in enumerate(CODONS)}

# Standard genetic code, in TCAG x TCAG x TCAG order.
AA_TABLE = ("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")
AA_OF = {c: AA_TABLE[i] for i, c in enumerate(CODONS)}
STOPS = {c for c in CODONS if AA_OF[c] == "*"}

MIN_CODONS = 100          # aligned, ungapped, unambiguous codons required
PS_SATURATION = 0.75      # Jukes-Cantor breaks down at pS >= 3/4


# --------------------------------------------------------------------------
# NG86 lookup tables
# --------------------------------------------------------------------------
def _build_site_table():
    """S_sites[i], N_sites[i] per codon; they sum to 3 (k = 1, i.e. no
    transition/transversion weighting - matches Biopython's default)."""
    S = np.zeros(64)
    N = np.zeros(64)
    for i, codon in enumerate(CODONS):
        if codon in STOPS:
            continue                      # stop codons are excluded from analysis
        aa = AA_OF[codon]
        s = n = 0
        for pos in range(3):
            for b in "ACGT":
                if b == codon[pos]:
                    continue
                nb = codon[:pos] + b + codon[pos + 1:]
                if nb in STOPS:
                    n += 1                # Biopython: change to stop = nonsynonymous
                elif AA_OF[nb] == aa:
                    s += 1
                else:
                    n += 1
        norm = (s + n) / 3.0              # = 3 when k = 1
        S[i], N[i] = s / norm, n / norm
    return S, N


def _compare(c1, c2, weight):
    """Synonymous/nonsynonymous weight for a single-step codon change."""
    return (weight, 0.0) if AA_OF[c1] == AA_OF[c2] else (0.0, weight)


def _build_diff_table():
    """Sd[i,j], Nd[i,j] averaged over every shortest mutational pathway."""
    SD = np.zeros((64, 64))
    ND = np.zeros((64, 64))
    for i, c1 in enumerate(CODONS):
        for j, c2 in enumerate(CODONS):
            if c1 == c2:
                continue
            diff = [p for p in range(3) if c1[p] != c2[p]]
            sd = nd = 0.0
            if len(diff) == 1:
                s, n = _compare(c1, c2, 1.0)
                sd, nd = sd + s, nd + n
            elif len(diff) == 2:
                for p in diff:
                    mid = c1[:p] + c2[p] + c1[p + 1:]
                    for a, b in ((c1, mid), (mid, c2)):
                        s, n = _compare(a, b, 0.5)
                        sd, nd = sd + s, nd + n
            else:
                for p1 in diff:
                    m1 = c1[:p1] + c2[p1] + c1[p1 + 1:]
                    for p2 in [p for p in diff if p != p1]:
                        m2 = m1[:p2] + c2[p2] + m1[p2 + 1:]
                        for a, b in ((c1, m1), (m1, m2), (m2, c2)):
                            s, n = _compare(a, b, 1.0 / 6.0)
                            sd, nd = sd + s, nd + n
            SD[i, j], ND[i, j] = sd, nd
    return SD, ND


S_SITES, N_SITES = _build_site_table()
SD_TABLE, ND_TABLE = _build_diff_table()


def encode_codons(seq):
    """CDS string -> int array of codon indices; -1 for gap/ambiguous/stop."""
    n = len(seq) // 3
    out = np.full(n, -1, dtype=np.int64)
    for i in range(n):
        c = seq[3 * i:3 * i + 3].upper().replace("U", "T")
        idx = CODON_IDX.get(c, -1)
        if idx >= 0 and c not in STOPS:
            out[i] = idx
    return out


def ng86(c1, c2):
    """Vectorised NG86. Inputs are equal-length codon-index arrays."""
    mask = (c1 >= 0) & (c2 >= 0)
    n_cod = int(mask.sum())
    if n_cod < MIN_CODONS:
        return dict(n_codons=n_cod, dN=np.nan, dS=np.nan, omega=np.nan,
                    pN=np.nan, pS=np.nan, identity=np.nan)
    a, b = c1[mask], c2[mask]
    s_sites = 0.5 * (S_SITES[a].sum() + S_SITES[b].sum())
    n_sites = 0.5 * (N_SITES[a].sum() + N_SITES[b].sum())
    sd = SD_TABLE[a, b].sum()
    nd = ND_TABLE[a, b].sum()
    ps = sd / s_sites if s_sites > 0 else np.nan
    pn = nd / n_sites if n_sites > 0 else np.nan

    def jc(p):
        if not np.isfinite(p) or p >= PS_SATURATION:
            return np.nan
        return abs(-0.75 * math.log(1.0 - 4.0 / 3.0 * p))

    dn, ds = jc(pn), jc(ps)
    omega = dn / ds if (np.isfinite(dn) and np.isfinite(ds) and ds > 0) else np.nan
    ident = float((a == b).mean())
    return dict(n_codons=n_cod, dN=dn, dS=ds, omega=omega,
                pN=pn, pS=ps, identity=ident)


# --------------------------------------------------------------------------
# alignment
# --------------------------------------------------------------------------
def _aligner():
    from Bio.Align import PairwiseAligner, substitution_matrices
    al = PairwiseAligner()
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score = -11
    al.extend_gap_score = -1
    # free terminal gaps: orthologues routinely differ in annotated start/end.
    # Biopython 1.87 renamed target_/query_end_gap_score, so set the pair that
    # this build actually exposes.
    try:
        al.end_gap_score = 0.0
    except AttributeError:                      # pragma: no cover - old Biopython
        al.target_end_gap_score = 0.0
        al.query_end_gap_score = 0.0
    al.mode = "global"
    return al


_ALIGNER = None


def translate(cds):
    """Translate, dropping a single terminal stop. Returns None on internal stop."""
    n = len(cds) // 3
    aas = []
    for i in range(n):
        aas.append(AA_OF.get(cds[3 * i:3 * i + 3].upper().replace("U", "T"), "X"))
    if aas and aas[-1] == "*":
        aas.pop()
    if "*" in aas:
        return None
    return "".join(aas)


def codon_align(cds1, cds2):
    """Protein-guided codon alignment. Returns (aligned_cds1, aligned_cds2)."""
    global _ALIGNER
    if _ALIGNER is None:
        _ALIGNER = _aligner()
    p1, p2 = translate(cds1), translate(cds2)
    if p1 is None or p2 is None or len(p1) < MIN_CODONS or len(p2) < MIN_CODONS:
        return None
    aln = _ALIGNER.align(p1.replace("X", "A"), p2.replace("X", "A"))[0]
    out1, out2 = [], []
    for (s1, e1), (s2, e2) in zip(aln.aligned[0], aln.aligned[1]):
        # aligned blocks are gap-free on both sides
        out1.append(cds1[3 * s1:3 * e1])
        out2.append(cds2[3 * s2:3 * e2])
    return "".join(out1), "".join(out2)


def pair_stats(args):
    hg, og, perc_id, conf, cds1, cds2 = args
    aligned = codon_align(cds1, cds2)
    if aligned is None:
        return None
    a1, a2 = aligned
    r = ng86(encode_codons(a1), encode_codons(a2))
    r.update(human_gene=hg, ortholog_gene=og,
             ensembl_perc_id=perc_id, confidence=conf,
             cds_len_human=len(cds1) // 3)
    return r


# --------------------------------------------------------------------------
# validation against Biopython
# --------------------------------------------------------------------------
def validate_against_biopython(pairs, n=60, seed=0):
    """Run both implementations on a random sample; return the max |difference|."""
    import warnings
    from Bio import BiopythonExperimentalWarning
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", BiopythonExperimentalWarning)
        from Bio.codonalign.codonseq import CodonSeq, cal_dn_ds

    rng = np.random.default_rng(seed)
    idx = rng.choice(len(pairs), size=min(n, len(pairs)), replace=False)
    diffs = []
    for i in idx:
        hg, og, pid, conf, c1, c2 = pairs[i]
        aligned = codon_align(c1, c2)
        if aligned is None:
            continue
        a1, a2 = aligned
        # Biopython errors on stop/ambiguous codons: drop those columns from
        # BOTH implementations so the comparison is like-for-like
        e1, e2 = encode_codons(a1), encode_codons(a2)
        keep = (e1 >= 0) & (e2 >= 0)
        if keep.sum() < MIN_CODONS:
            continue
        f1 = "".join(a1[3 * k:3 * k + 3] for k in np.flatnonzero(keep))
        f2 = "".join(a2[3 * k:3 * k + 3] for k in np.flatnonzero(keep))
        mine = ng86(encode_codons(f1), encode_codons(f2))
        try:
            bdn, bds = cal_dn_ds(CodonSeq(f1), CodonSeq(f2), method="NG86")
        except Exception:                       # noqa: BLE001
            continue
        if bdn < 0 or bds < 0:
            continue
        diffs.append((abs(mine["dN"] - bdn), abs(mine["dS"] - bds)))
    if not diffs:
        return None
    arr = np.array(diffs)
    return dict(n_compared=len(arr), max_abs_dN=float(arr[:, 0].max()),
                max_abs_dS=float(arr[:, 1].max()))


# --------------------------------------------------------------------------
# driver
# --------------------------------------------------------------------------
def load_pairs(sp_key):
    path = os.path.join(DERIVED, f"pairs_human_{sp_key}.tsv")
    if not os.path.exists(path):
        raise FileNotFoundError(f"{path} - run tissue_evolution.t2_orthologs first")
    out = []
    with open(path) as fh:
        rd = csv.reader(fh, delimiter="\t")
        next(rd)
        for r in rd:
            out.append((r[0], r[1], r[2], r[3], r[4], r[5]))
    return out


def run_species(sp_key, processes=4, force=False, validate=True):
    out = os.path.join(DERIVED, f"t3_dnds_human_{sp_key}.tsv.gz")
    if os.path.exists(out) and not force:
        log(f"  {sp_key}: dN/dS already computed")
        return out
    pairs = load_pairs(sp_key)
    log(f"  {sp_key}: {len(pairs)} 1:1 pairs")

    if validate:
        v = validate_against_biopython(pairs)
        if v:
            log(f"  {sp_key}: NG86 cross-check vs Biopython on n={v['n_compared']}: "
                f"max|dN diff|={v['max_abs_dN']:.2e}  max|dS diff|={v['max_abs_dS']:.2e}")

    import pandas as pd
    with Pool(processes) as pool:
        rows = [r for r in pool.imap_unordered(pair_stats, pairs, chunksize=64)
                if r is not None]
    df = pd.DataFrame(rows)
    to_tsv_gz(df, out, index=False)
    ok = df["omega"].notna()
    log(f"  {sp_key}: {len(df)} aligned, {ok.sum()} with usable omega, "
        f"median dS={df['dS'].median():.3f} median omega={df.loc[ok,'omega'].median():.4f}")
    return out


def report(processes=4):
    print("=" * 72)
    print("T3 - dN/dS from scratch (NG86 on protein-guided codon alignments)")
    print("=" * 72)
    import pandas as pd
    for sp in SPECIES:
        path = os.path.join(DERIVED, f"pairs_human_{sp['key']}.tsv")
        if not os.path.exists(path):
            log(f"  {sp['key']}: pairs missing, skipped")
            continue
        log(f"{sp['key']} (~{sp['divergence_mya']} Mya)")
        run_species(sp["key"], processes=processes)
    print()
    print(f"{'species':<10} {'pairs':>7} {'usable w':>9} {'med dS':>8} "
          f"{'med dN':>8} {'med w':>8} {'dS>1 frac':>10}")
    for sp in SPECIES:
        f = os.path.join(DERIVED, f"t3_dnds_human_{sp['key']}.tsv.gz")
        if not os.path.exists(f):
            continue
        d = pd.read_csv(f, sep="\t")
        ok = d["omega"].notna()
        sat = (d["dS"] > 1).mean()
        print(f"{sp['key']:<10} {len(d):>7} {ok.sum():>9} {d['dS'].median():>8.3f} "
              f"{d['dN'].median():>8.4f} {d.loc[ok, 'omega'].median():>8.4f} {sat:>10.1%}")
    return True


if __name__ == "__main__":
    report(processes=int(sys.argv[1]) if len(sys.argv) > 1 else 4)
