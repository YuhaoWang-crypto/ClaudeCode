"""
T24 - the subunit-swap matrix: is there evidence that partners coevolve?

The proposal being tested is a 2x2 (here NxN) design:

                  partner from species A   partner from species B
    protein A          native                    mismatched
    protein B        mismatched                    native

If both native combinations score better than both mismatched ones, the two
proteins carry information about each other - their interface residues are
matched to a specific partner rather than to a generic one. That is the
definition of coevolution, and unlike a comparison of two RER values it cannot
be produced by the two proteins merely sharing a rate.

WHAT IS SCORED. Not expression, not abundance - the score is built only from
which amino acids sit at structurally contacting positions, so "the difference
is really expression" is not available as an explanation. For a contact (i in
protein A, j in protein B) taken from a solved human complex, and a panel of
~180 vertebrate orthologues on human coordinates (T22):

    S(s, t) = sum over contacts of  log [ P(a_i, b_j) / P(a_i) P(b_j) ]

with a_i from species s and b_j from species t. This is pointwise mutual
information summed over the interface: a combination that the panel has
actually explored scores positive, one it never explores scores negative.

THREE WAYS THIS COULD BE WRONG, AND WHAT IS DONE ABOUT EACH.

1. Circularity. The native combination (s, s) is itself one of the
   observations used to estimate P(a_i, b_j), so it would be scored against a
   distribution it helped build. Every score here is therefore LEAVE-TWO-OUT:
   species s and t are both removed from the frequency estimate before
   S(s, t) is computed. Without this the diagonal wins automatically and the
   result is meaningless.

2. Phylogeny. Close relatives share amino acids at every position for reasons
   that have nothing to do with the interface, so a native pair is favoured
   simply because the two chains come from the same branch of the tree. The
   control is internal: the identical computation is run on NON-contacting
   position pairs drawn from the same two proteins and the same species panel,
   matched to the contacts on the marginal entropy of both columns. Phylogeny
   inflates contacts and non-contacts alike; coevolution does not.

3. Conservation. A perfectly invariant interface gives every species the same
   residues, so native and mismatched become identical and the test has no
   power - it cannot show coevolution and it cannot show its absence. Columns
   with no variation are therefore reported as excluded, not silently dropped,
   because for this gene set (OXPHOS cores are among the most conserved
   proteins in the genome) that exclusion may well be most of the interface.

Complexes used, all human, all experimentally solved:
    8GS8  complex II   SDHA/B/C/D - the group-3 candidate
    5Z62  complex IV   mtDNA-encoded CO1/2/3 against nuclear-encoded subunits
    5XTE  complex III  mtDNA-encoded CYB against nuclear subunits
    5XTD  complex I    mtDNA-encoded ND subunits against nuclear subunits
    2FP4  succinyl-CoA ligase  SUCLG1 + SUCLG2 - the group-1 candidate

The mitochondrial complexes are the sharpest test available: mtDNA evolves
about an order of magnitude faster than the nuclear genome, so if any protein
pair in a genome is under pressure to track a partner, it is a nuclear subunit
bolted onto an mtDNA-encoded core.
"""
import gzip
import io
import os
import urllib.request
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy import stats

from .common import RAW, DERIVED, FIGDIR, apply_figure_style, PALETTE, log
from . import t22_orthologs

PDBDIR = os.path.join(RAW, "pdb")
CONTACT_A = 5.0            # heavy-atom distance defining a contact
N_NULL = 8                 # independent matched non-contact sets per interface
AA = "ACDEFGHIKLMNPQRSTVWY"

COMPLEXES = {
    "8GS8": "complex II (SDH)",
    "5Z62": "complex IV (cytochrome c oxidase)",
    "5XTE": "complex III (cytochrome bc1)",
    "5XTD": "complex I",
    "2FP4": "succinyl-CoA ligase",
}


# ------------------------------------------------------------ structures

def fetch_cif(pdb_id):
    os.makedirs(PDBDIR, exist_ok=True)
    p = os.path.join(PDBDIR, f"{pdb_id}.cif.gz")
    if not os.path.exists(p):
        url = f"https://files.rcsb.org/download/{pdb_id}.cif.gz"
        with urllib.request.urlopen(url, timeout=300) as fh:
            data = fh.read()
        with open(p, "wb") as fh:
            fh.write(data)
    with gzip.open(p, "rt") as fh:
        return fh.read()


def parse_cif(text):
    """Entity sequences, chain->entity map, and heavy-atom coordinates.

    label_seq_id is used rather than auth_seq_id: it indexes directly into the
    entity's canonical sequence, so no author-numbering gymnastics are needed
    to line coordinates up with residues.
    """
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    d = MMCIF2Dict(io.StringIO(text))

    def col(k):
        v = d.get(k, [])
        return v if isinstance(v, list) else [v]

    seqs = {}
    for eid, s in zip(col("_entity_poly.entity_id"),
                      col("_entity_poly.pdbx_seq_one_letter_code_can")):
        seqs[eid] = "".join(s.split())
    chain_entity = {}
    for eid, chains in zip(col("_entity_poly.entity_id"),
                           col("_entity_poly.pdbx_strand_id")):
        for c in chains.split(","):
            chain_entity[c.strip()] = eid

    coords = defaultdict(list)
    for grp, ch, lseq, el, x, y, z in zip(
            col("_atom_site.group_PDB"), col("_atom_site.auth_asym_id"),
            col("_atom_site.label_seq_id"), col("_atom_site.type_symbol"),
            col("_atom_site.Cartn_x"), col("_atom_site.Cartn_y"),
            col("_atom_site.Cartn_z")):
        if grp != "ATOM" or el == "H" or lseq in (".", "?"):
            continue
        coords[(ch, int(lseq))].append((float(x), float(y), float(z)))
    return seqs, chain_entity, coords


def chain_contacts(coords, chain_a, chain_b, cutoff=CONTACT_A):
    """Residue pairs (1-based label_seq_id) within cutoff, via a grid index."""
    A = {k[1]: np.array(v) for k, v in coords.items() if k[0] == chain_a}
    B = {k[1]: np.array(v) for k, v in coords.items() if k[0] == chain_b}
    if not A or not B:
        return []
    bkeys = list(B)
    ball = np.vstack([B[k] for k in bkeys])
    bidx = np.concatenate([[k] * len(B[k]) for k in bkeys])
    # coarse grid so this stays O(atoms) rather than O(atoms^2)
    cell = cutoff
    grid = defaultdict(list)
    for n, pt in enumerate(ball):
        grid[tuple((pt // cell).astype(int))].append(n)
    out = set()
    for ra, pts in A.items():
        for pt in pts:
            base = (pt // cell).astype(int)
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    for dz in (-1, 0, 1):
                        for n in grid.get((base[0] + dx, base[1] + dy,
                                           base[2] + dz), ()):
                            if np.linalg.norm(pt - ball[n]) <= cutoff:
                                out.add((ra, int(bidx[n])))
    return sorted(out)


# ------------------------------------------------- structure <-> panel map

def align_index(pdb_seq, human_seq):
    """pdb position (1-based) -> human panel position (0-based)."""
    from Bio import Align
    al = Align.PairwiseAligner(mode="global", open_gap_score=-11,
                               extend_gap_score=-1, substitution_matrix=None)
    al.match_score, al.mismatch_score = 2, -1
    try:
        al.target_end_gap_score = 0.0
        al.query_end_gap_score = 0.0
    except AttributeError:
        al.end_gap_score = 0.0
    aln = al.align(pdb_seq, human_seq)[0]
    m = {}
    ident = 0
    for (ps, pe), (hs, he) in zip(aln.aligned[0], aln.aligned[1]):
        for k in range(pe - ps):
            m[ps + k + 1] = hs + k
            ident += pdb_seq[ps + k] == human_seq[hs + k]
    denom = min(len(pdb_seq), len(human_seq))
    return m, ident / max(denom, 1)


def identify_chains(seqs, chain_entity, panel, min_ident=0.90):
    """Assign each chain to a panel gene by sequence identity."""
    out = {}
    for chain, eid in chain_entity.items():
        s = seqs.get(eid, "")
        if len(s) < 40:
            continue
        best = (None, 0.0, None)
        for gene, d in panel.items():
            if abs(len(d["human"]) - len(s)) > 0.6 * max(len(s), 1):
                continue
            m, ident = align_index(s, d["human"])
            if ident > best[1]:
                best = (gene, ident, m)
        if best[0] and best[1] >= min_ident:
            out[chain] = dict(gene=best[0], ident=best[1], map=best[2],
                              entity=eid)
    return out


# --------------------------------------------------------- coevolution

def columns(panel, gene, species):
    """species x position character matrix for one gene."""
    cols = panel[gene]["cols"]
    return np.array([list(cols.get(s, "-" * len(panel[gene]["human"])))
                     for s in species])


def _mi(a, b):
    """Mutual information (nats) of two aligned character vectors."""
    ok = np.array([x in AA and y in AA for x, y in zip(a, b)])
    if ok.sum() < 20:
        return np.nan
    a, b = a[ok], b[ok]
    ua, ia = np.unique(a, return_inverse=True)
    ub, ib = np.unique(b, return_inverse=True)
    if len(ua) < 2 or len(ub) < 2:
        return 0.0
    n = len(a)
    joint = np.zeros((len(ua), len(ub)))
    np.add.at(joint, (ia, ib), 1.0)
    joint /= n
    pa, pb = joint.sum(1, keepdims=True), joint.sum(0, keepdims=True)
    with np.errstate(divide="ignore", invalid="ignore"):
        t = joint * np.log(joint / (pa @ pb))
    return float(np.nansum(t))


def entropy(col):
    v = [c for c in col if c in AA]
    if len(v) < 20:
        return np.nan
    _, cnt = np.unique(v, return_counts=True)
    p = cnt / cnt.sum()
    return float(-(p * np.log(p)).sum())


def swap_scores(MA, MB, pairs, species, pseudo=0.5):
    """Leave-two-out PMI compatibility matrix S(s, t) over the given pairs.

    Leave-two-out is the whole ballgame. If species s and t stayed in the
    frequency table, S(s, s) would be scored against a distribution that
    contains s's own residue combination, and the diagonal would win no matter
    what the biology is.

    Done naively this is n_species^2 x n_contacts rebuilds of a 20x20 table
    (~5 million for one interface), which is far too slow. But removing two
    observations only perturbs three quantities, so the whole leave-two-out
    grid closes in closed form. Writing r = rowsums, c = colsums, N = total:

        C'[a_s, b_t] = C[a_s, b_t] - [b_s == b_t] - [t != s][a_t == a_s]
        r'[a_s]      = r[a_s] - 1     - [t != s][a_t == a_s]
        c'[b_t]      = c[b_t] - [b_s == b_t] - [t != s]
        N'           = N - 1 - [t != s]

    which is three boolean n x n matrices and a gather - fully vectorised.
    """
    n = len(species)
    S = np.zeros((n, n))
    seen = np.zeros((n, n))
    used = 0
    idx = {c: k for k, c in enumerate(AA)}
    for i, j in pairs:
        a, b = MA[:, i], MB[:, j]
        ai = np.array([idx.get(x, -1) for x in a])
        bi = np.array([idx.get(x, -1) for x in b])
        ok = (ai >= 0) & (bi >= 0)
        if ok.sum() < 30 or len(set(ai[ok])) < 2 or len(set(bi[ok])) < 2:
            continue                      # invariant column: no information
        used += 1
        C = np.zeros((20, 20))
        np.add.at(C, (ai[ok], bi[ok]), 1.0)
        r, c, N = C.sum(1), C.sum(0), float(ok.sum())

        As = ai[:, None]                  # a_s, varies down rows
        Bt = bi[None, :]                  # b_t, varies across columns
        m_same_a = (ai[None, :] == ai[:, None])       # a_t == a_s
        m_same_b = (bi[:, None] == bi[None, :])       # b_s == b_t
        neq = ~np.eye(n, dtype=bool)

        num = C[np.clip(As, 0, 19), np.clip(Bt, 0, 19)] \
            - m_same_b - (neq & m_same_a)
        rr = r[np.clip(As, 0, 19)] - 1.0 - (neq & m_same_a)
        cc = c[np.clip(Bt, 0, 19)] - m_same_b - neq
        NN = N - 1.0 - neq

        num = num + pseudo
        rr = rr + 20 * pseudo
        cc = cc + 20 * pseudo
        NN = NN + 400 * pseudo

        valid = ok[:, None] & ok[None, :]
        with np.errstate(divide="ignore", invalid="ignore"):
            s_ij = np.log(num * NN / (rr * cc))
        S += np.where(valid, np.nan_to_num(s_ij), 0.0)
        seen += valid
    with np.errstate(divide="ignore", invalid="ignore"):
        S = np.where(seen > 0, S / np.maximum(seen, 1) * max(used, 1), np.nan)
    return S, used


def matched_noncontacts(MA, MB, contacts, rng, n_want):
    """Non-contact position pairs matched on both columns' entropies."""
    cset = set(contacts)
    ea = {i: entropy(MA[:, i]) for i in {i for i, _ in contacts}}
    eb = {j: entropy(MB[:, j]) for j in {j for _, j in contacts}}
    allA = [i for i in range(MA.shape[1]) if not np.isnan(entropy(MA[:, i]))]
    allB = [j for j in range(MB.shape[1]) if not np.isnan(entropy(MB[:, j]))]
    EA = {i: entropy(MA[:, i]) for i in allA}
    EB = {j: entropy(MB[:, j]) for j in allB}
    out = []
    for i, j in contacts:
        if np.isnan(ea.get(i, np.nan)) or np.isnan(eb.get(j, np.nan)):
            continue
        ca = [k for k in allA if abs(EA[k] - ea[i]) < 0.15]
        cb = [k for k in allB if abs(EB[k] - eb[j]) < 0.15]
        for _ in range(40):
            if not ca or not cb:
                break
            p = (rng.choice(ca), rng.choice(cb))
            if p not in cset:
                out.append(p)
                break
        if len(out) >= n_want:
            break
    return out


def native_vs_mismatch(S, min_species=20):
    """Diagonal (native) against off-diagonal (mismatched), per species.

    A species that is gapped at every informative contact scores against
    nothing and comes back all-NaN. Those rows are dropped rather than allowed
    to propagate: a single NaN in the diagonal would otherwise turn the whole
    interface's result into NaN, which is how the first run silently lost
    three of the four complex II interfaces.
    """
    S = np.asarray(S, float)
    good = ~np.all(np.isnan(S), axis=1) & ~np.all(np.isnan(S), axis=0)
    S = S[np.ix_(good, good)]
    n = S.shape[0]
    if n < min_species:
        return {}
    diag = np.diag(S).astype(float)
    off = np.array([np.nanmean(np.delete(S[i], i)) for i in range(n)])
    keep = np.isfinite(diag) & np.isfinite(off)
    if keep.sum() < min_species:
        return {}
    diag, off = diag[keep], off[keep]
    d = diag - off
    w = stats.wilcoxon(diag, off)
    return dict(n_species=int(keep.sum()), native=float(np.mean(diag)),
                mismatch=float(np.mean(off)), delta=float(np.mean(d)),
                frac_native_wins=float((d > 0).mean()),
                p=float(w.pvalue))


# --------------------------------------------------------------- driver

def analyse_complex(pdb_id, panel, species, rng, min_contacts=15):
    """Every inter-chain interface in one structure, scored and controlled."""
    seqs, chain_entity, coords = parse_cif(fetch_cif(pdb_id))
    ident = identify_chains(seqs, chain_entity, panel)
    if len(ident) < 2:
        return []
    # one chain per gene (complexes carry duplicate copies)
    by_gene = {}
    for ch, info in ident.items():
        if info["gene"] not in by_gene or info["ident"] > by_gene[info["gene"]]["ident"]:
            by_gene[info["gene"]] = dict(info, chain=ch)

    cols = {g: columns(panel, g, species) for g in by_gene}
    genes = sorted(by_gene)
    out = []
    for x in range(len(genes)):
        for y in range(x + 1, len(genes)):
            ga, gb = genes[x], genes[y]
            A, B = by_gene[ga], by_gene[gb]
            raw = chain_contacts(coords, A["chain"], B["chain"])
            pairs = [(A["map"][i], B["map"][j]) for i, j in raw
                     if i in A["map"] and j in B["map"]]
            pairs = sorted(set(pairs))
            if len(pairs) < min_contacts:
                continue
            MA, MB = cols[ga], cols[gb]
            S, used = swap_scores(MA, MB, pairs, species)
            if used < 5:
                out.append(dict(pdb=pdb_id, gene_a=ga, gene_b=gb,
                                n_contacts=len(pairs), n_informative=used,
                                note="interface too conserved to test"))
                continue
            res = native_vs_mismatch(S)
            if not res:
                out.append(dict(pdb=pdb_id, gene_a=ga, gene_b=gb,
                                n_contacts=len(pairs), n_informative=used,
                                note="too few species scoreable"))
                continue

            # A single random set of matched non-contacts is a noisy
            # control: on SUCLG1/SUCLG2 two different draws moved the excess
            # from +0.03 to +0.91, which is the difference between "no
            # coevolution" and "clear coevolution". Average several.
            nulls, deltas, fracs = [], [], []
            for k in range(N_NULL):
                nl = matched_noncontacts(MA, MB, pairs, rng, len(pairs))
                if not nl:
                    continue
                Sn, used_n = swap_scores(MA, MB, nl, species)
                rn = native_vs_mismatch(Sn) if used_n >= 5 else {}
                if rn:
                    nulls.append(nl)
                    deltas.append(rn["delta"])
                    fracs.append(rn["frac_native_wins"])
            null = nulls[-1] if nulls else []
            used_n = len(deltas)
            resn = dict(delta=float(np.mean(deltas)),
                        delta_sd=float(np.std(deltas)),
                        frac_native_wins=float(np.mean(fracs))) if deltas else {}

            # diagnostic: the non-contact control is only a control if it
            # really is matched on variability. A negative "excess" would
            # otherwise just mean the interface is more conserved.
            ec = np.nanmean([entropy(MA[:, i]) + entropy(MB[:, j])
                             for i, j in pairs])
            en = np.nanmean([entropy(MA[:, i]) + entropy(MB[:, j])
                             for i, j in null]) if null else np.nan
            mi_c = np.nanmean([_mi(MA[:, i], MB[:, j]) for i, j in pairs])
            mi_n = np.nanmean([_mi(MA[:, i], MB[:, j]) for i, j in null]) \
                if null else np.nan
            out.append(dict(
                pdb=pdb_id, gene_a=ga, gene_b=gb, n_contacts=len(pairs),
                n_informative=used, n_null=used_n,
                native=res["native"], mismatch=res["mismatch"],
                delta=res["delta"], frac_native_wins=res["frac_native_wins"],
                p=res["p"],
                null_delta=resn.get("delta", np.nan),
                null_sd=resn.get("delta_sd", np.nan),
                null_frac=resn.get("frac_native_wins", np.nan),
                n_null_draws=used_n,
                excess=res["delta"] - resn.get("delta", np.nan),
                mi_contact=mi_c, mi_null=mi_n,
                entropy_contact=ec, entropy_null=en, note=""))
    return out


def named_pair(panel, species, gene_a, gene_b, pdb_id, rng):
    """The 2x2 the proposal asks for, on named species, plus its control."""
    seqs, ce, coords = parse_cif(fetch_cif(pdb_id))
    ident = identify_chains(seqs, ce, panel)
    pick = {}
    for ch, info in ident.items():
        g = info["gene"]
        if g in (gene_a, gene_b) and (g not in pick
                                      or info["ident"] > pick[g]["ident"]):
            pick[g] = dict(info, chain=ch)
    if len(pick) < 2:
        return None
    A, B = pick[gene_a], pick[gene_b]
    raw = chain_contacts(coords, A["chain"], B["chain"])
    pairs = sorted({(A["map"][i], B["map"][j]) for i, j in raw
                    if i in A["map"] and j in B["map"]})
    MA, MB = columns(panel, gene_a, species), columns(panel, gene_b, species)
    S, used = swap_scores(MA, MB, pairs, species)
    null = matched_noncontacts(MA, MB, pairs, rng, len(pairs))
    Sn, _ = swap_scores(MA, MB, null, species)
    return dict(S=S, S_null=Sn, pairs=pairs, used=used, species=species)


def is_mt(gene):
    return str(gene).startswith("MT-")


def report():
    print("=" * 72)
    print("T24 - subunit-swap matrix: do interacting partners coevolve?")
    print("=" * 72)
    panel = t22_orthologs.build(verbose=False)
    log(f"panel: {len(panel)} genes")
    rng = np.random.default_rng(0)

    rows = []
    for pdb_id, label in COMPLEXES.items():
        sub = {g: d for g, d in panel.items()}
        sp = t22_orthologs.common_species(sub, 0.0)
        # a common species set per structure, so every interface in it is
        # scored on identical taxa
        try:
            res = analyse_complex(pdb_id, panel, sp, rng)
        except Exception as exc:                       # noqa: BLE001
            log(f"  {pdb_id} failed: {exc}")
            continue
        for r in res:
            r["complex"] = label
        rows += res
        log(f"  {pdb_id} {label}: {len(res)} interfaces")

    df = pd.DataFrame(rows)
    if df.empty:
        print("no interfaces recovered")
        return df
    df = df.sort_values("excess", ascending=False)

    print(f"\n{'complex':<34}{'A':<9}{'B':<9}{'cont':>6}{'inf':>5}"
          f"{'delta':>8}{'null':>8}{'excess':>9}{'MIc':>8}{'MIn':>8}")
    for _, r in df.iterrows():
        if r.get("note"):
            print(f"{r['complex']:<34}{r['gene_a']:<9}{r['gene_b']:<9}"
                  f"{r['n_contacts']:>6}{r['n_informative']:>5}"
                  f"   {r['note']}")
            continue
        print(f"{r['complex']:<34}{r['gene_a']:<9}{r['gene_b']:<9}"
              f"{r['n_contacts']:>6}{r['n_informative']:>5}"
              f"{r['delta']:>8.3f}{r['null_delta']:>8.3f}{r['excess']:>9.3f}"
              f"{r['mi_contact']:>8.4f}{r['mi_null']:>8.4f}")

    ok = df[df["excess"].notna()]
    if len(ok):
        w = stats.wilcoxon(ok["delta"], ok["null_delta"])
        print(f"\ninterfaces tested: {len(ok)}")
        print(f"  native-minus-mismatch at CONTACTS      median "
              f"{ok['delta'].median():.3f}")
        print(f"  same at matched NON-contacts (control) median "
              f"{ok['null_delta'].median():.3f}")
        print(f"  excess attributable to the interface   median "
              f"{ok['excess'].median():+.3f}   "
              f"(Wilcoxon p={w.pvalue:.3g}, n={len(ok)})")
        print(f"  interfaces where contacts beat the control: "
              f"{(ok['excess'] > 0).sum()}/{len(ok)}")
        print(f"  fraction of the raw signal that survives the control: "
              f"{ok['excess'].median() / ok['delta'].median():.0%}")

        # the sharp prediction: mtDNA evolves ~10x faster, so nuclear subunits
        # bolted onto an mtDNA-encoded core should be the ones forced to track
        mn = ok[[is_mt(a) != is_mt(b)
                 for a, b in zip(ok["gene_a"], ok["gene_b"])]]
        nn = ok[[not is_mt(a) and not is_mt(b)
                 for a, b in zip(ok["gene_a"], ok["gene_b"])]]
        if len(mn) >= 4 and len(nn) >= 4:
            u = stats.mannwhitneyu(mn["excess"], nn["excess"],
                                   alternative="greater")
            print(f"\n  mito-nuclear interfaces   n={len(mn):>2}  "
                  f"median excess {mn['excess'].median():+.3f}")
            print(f"  nuclear-nuclear interfaces n={len(nn):>2}  "
                  f"median excess {nn['excess'].median():+.3f}")
            print(f"  Mann-Whitney (mito-nuclear greater): p={u.pvalue:.4f}")

    df.to_csv(os.path.join(DERIVED, "t24_swap.tsv"), sep="\t", index=False)
    return df


if __name__ == "__main__":
    report()
