"""Stage 3 + Stage 4: the chemistry filter and the active-site geometry check.

This is the part of the pipeline that converts a *fold* statement into a
*chemistry* statement. A structural hit without this filter is not a candidate
enzyme, it is a shape.

Two ideas do the work:

1. Catalytic residues are mapped through the **structural** alignment, not a
   sequence alignment. At ~15% sequence identity a sequence alignment places
   them wrong; the structure superposition places them right.
2. Conservation is judged by **chemical class**, not identity. In the
   precedent (cGAS -> plant NTPs) Glu<->Asp substitution was accepted because
   both maintain an acidic active-site environment; anything non-acidic was
   fatal.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from structure import core_rmsd, kabsch_rmsd

# Chemical equivalence classes for "conservative" substitution.
CLASSES = {
    "acidic": set("DE"),
    "basic": set("KRH"),
    "amide": set("NQ"),
    "hydroxyl": set("STY"),
    "aromatic": set("FWY"),
    "hydrophobic": set("AVLIMFWCP"),
    "small": set("GAS"),
    "any": set("ACDEFGHIKLMNPQRSTVWY"),
}


def residue_map(qstart, tstart, qaln, taln):
    """Map query -> target residue index through a gapped alignment.

    Indices are 1-based positions in each chain's ordered residue list, which
    is the coordinate system Foldseek reports (qstart/tstart).
    Gap-aligned positions are simply absent from the mapping.
    """
    mapping = {}
    qi, ti = int(qstart), int(tstart)
    for qc, tc in zip(qaln, taln):
        if qc != "-" and tc != "-":
            mapping[qi] = ti
            qi += 1
            ti += 1
        elif qc != "-":
            qi += 1
        elif tc != "-":
            ti += 1
    return mapping


@dataclass
class CatalyticSite:
    """The query's catalytic residues, given as author residue numbers."""
    resnums: list
    klass: str = "acidic"          # required chemical class at every position
    label: str = ""

    def allowed(self):
        return CLASSES[self.klass]


def anchor_index(mapping, qi):
    """Target index for a query position, extrapolating across alignment gaps.

    A catalytic residue sitting in a gap column is common and is not evidence
    of anything; the nearest aligned position gives a usable anchor to search
    around.
    """
    if qi in mapping:
        return mapping[qi], 0
    if not mapping:
        return None, None
    qj = min(mapping, key=lambda k: abs(k - qi))
    return mapping[qj] + (qi - qj), abs(qi - qj)


def _assign(order, candidates):
    """Pick one target residue per catalytic position.

    Maximises the number of matched positions, then minimises the total
    registration offset, subject to: no target residue used twice, and target
    indices strictly increasing along the query. Catalytic sets are tiny
    (3-5 positions), so exhaustive search is both exact and instant.
    """
    best = {"n": -1, "cost": 0, "pick": {}}

    def recurse(i, used_idx, last_idx, pick, n, cost):
        if i == len(order):
            if n > best["n"] or (n == best["n"] and cost < best["cost"]):
                best.update(n=n, cost=cost, pick=dict(pick))
            return
        k = order[i]
        for ti, off, res in candidates[k]:
            if ti in used_idx or ti <= last_idx:
                continue
            pick[k] = (ti, off, res)
            recurse(i + 1, used_idx | {ti}, ti, pick, n + 1, cost + abs(off))
            del pick[k]
        # allow this position to go unmatched
        recurse(i + 1, used_idx, last_idx, pick, n, cost)

    recurse(0, frozenset(), -10**9, {}, 0, 0)
    return best["pick"]


@dataclass
class SiteCheck:
    passed: bool
    positions: list = field(default_factory=list)  # per-residue detail
    n_conserved: int = 0          # counts windowed matches
    n_conserved_exact: int = 0    # counts offset-0 matches only
    n_total: int = 0
    max_offset: int = 0
    local_rmsd: float = float("nan")      # all aligned pairs in the shell
    n_shell: int = 0
    core_rmsd: float = float("nan")       # outlier-trimmed -- judge by this
    n_core: int = 0
    core_frac: float = 0.0
    shell_plddt: float = None
    carboxylate_deltas: list = field(default_factory=list)
    note: str = ""

    def summary(self):
        found = "".join(p["target_aa"] or "-" for p in self.positions)
        return (f"{self.n_conserved}/{self.n_total} [{found}] "
                f"core={self.core_rmsd:.2f}A/{self.n_core} "
                f"shell={self.local_rmsd:.2f}A/{self.n_shell}")


def check_site(query, target, site, mapping, shell_radius=12.0, window=0):
    """Run the chemistry filter and the local geometry check for one hit.

    query, target : Structure
    site          : CatalyticSite defined on the query
    mapping       : query index -> target index, from residue_map()
    window        : allow the counterpart to sit +/- this many residues from
                    the structurally-aligned position. At low sequence identity
                    a structural alignment routinely misregisters a catalytic
                    motif by 1-2 residues; 0 means strict. Any non-zero offset
                    is recorded per position so the evidence stays visible.
    """
    allowed = site.allowed()
    positions, ok = [], True

    # Build per-position candidate lists, then choose one consistent
    # assignment: distinct target residues, in the same order as the query.
    # Without that constraint a windowed search lets two catalytic positions
    # claim the same residue, which looks like conservation but is an artefact.
    order = sorted(range(len(site.resnums)),
                   key=lambda k: query.resnum_to_index(site.resnums[k]) or 0)
    anchors, candidates = {}, {}
    for k in order:
        rn = site.resnums[k]
        qi = query.resnum_to_index(rn)
        anchor, gap = (anchor_index(mapping, qi) if qi else (None, None))
        anchors[k] = (anchor, gap)
        opts = []
        if anchor is not None:
            for off in sorted(range(-window, window + 1), key=abs):
                cand = target.by_index(anchor + off)
                if cand is not None and cand.one in allowed:
                    opts.append((anchor + off, off, cand))
        candidates[k] = opts
    chosen = _assign(order, candidates)

    cat_idx = []
    for k in range(len(site.resnums)):
        rn = site.resnums[k]
        qi = query.resnum_to_index(rn)
        qres = query.by_resnum(rn)
        anchor, gap = anchors[k]
        pick = chosen.get(k)
        hit_res, offset = (pick[2], pick[1]) if pick else (None, None)
        direct = target.by_index(anchor) if anchor is not None else None
        conserved = hit_res is not None
        if not conserved:
            ok = False
        if qi:
            cat_idx.append(qi)
        positions.append({
            "query_resnum": rn,
            "query_aa": qres.one if qres else None,
            "query_index": qi,
            "target_index": (anchor + offset) if hit_res is not None else anchor,
            "target_resnum": hit_res.resnum if hit_res else
                             (direct.resnum if direct else None),
            "target_aa": hit_res.one if hit_res else
                         (direct.one if direct else None),
            "conserved": conserved,
            "offset": offset,
            "gap_extrapolated": bool(gap),
            "plddt": (hit_res or direct).plddt if (hit_res or direct) else None,
        })

    chk = SiteCheck(
        passed=ok,
        positions=positions,
        n_conserved=sum(p["conserved"] for p in positions),
        n_conserved_exact=sum(1 for p in positions
                              if p["conserved"] and p["offset"] == 0),
        n_total=len(site.resnums),
        max_offset=max((abs(p["offset"]) for p in positions
                        if p["offset"] is not None), default=0),
    )

    # --- local superposition over the catalytic shell only ------------------
    # Global RMSD is expected to be poor for a divergent scaffold; the active
    # site is what must superimpose. See SKILL.md, Stage 4.
    cat_ca = [query.by_index(i).ca for i in cat_idx
              if query.by_index(i) is not None and query.by_index(i).ca is not None]
    if cat_ca and mapping:
        cat_ca = np.array(cat_ca)
        shell_q, shell_t = [], []
        for qi, ti in mapping.items():
            qres, tres = query.by_index(qi), target.by_index(ti)
            if qres is None or tres is None:
                continue
            if qres.ca is None or tres.ca is None:
                continue
            if np.min(np.linalg.norm(cat_ca - qres.ca, axis=1)) <= shell_radius:
                shell_q.append(qres.ca)
                shell_t.append(tres.ca)
        if len(shell_q) >= 3:
            chk.local_rmsd, chk.n_shell = kabsch_rmsd(shell_t, shell_q)
            chk.core_rmsd, chk.n_core, chk.core_frac = core_rmsd(shell_t, shell_q)
            idxs = [mapping[qi] for qi in mapping
                    if query.by_index(qi) is not None
                    and query.by_index(qi).ca is not None
                    and np.min(np.linalg.norm(cat_ca - query.by_index(qi).ca, axis=1)) <= shell_radius]
            chk.shell_plddt = target.plddt(idxs)

    # --- acidic-cluster geometry (needs side chains; local mode only) -------
    qcs = [query.by_resnum(rn) for rn in site.resnums]
    tcs = [target.by_index(p["target_index"]) if p["target_index"] else None
           for p in positions]
    for a in range(len(site.resnums)):
        for b in range(a + 1, len(site.resnums)):
            qa, qb = qcs[a], qcs[b]
            ta, tb = tcs[a], tcs[b]
            if not all([qa, qb, ta, tb]):
                continue
            qca, qcb = qa.carboxylate_c(), qb.carboxylate_c()
            tca, tcb = ta.carboxylate_c(), tb.carboxylate_c()
            if any(x is None for x in (qca, qcb, tca, tcb)):
                continue
            dq = float(np.linalg.norm(qca - qcb))
            dt = float(np.linalg.norm(tca - tcb))
            chk.carboxylate_deltas.append({
                "pair": (site.resnums[a], site.resnums[b]),
                "query_dist": round(dq, 2),
                "target_dist": round(dt, 2),
                "delta": round(dt - dq, 2),
            })
    return chk
