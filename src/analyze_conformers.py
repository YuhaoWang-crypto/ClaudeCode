"""
analyze_conformers.py
---------------------
Align a conformer pair, measure where the structure responds, and rank DISTAL
candidate regions (far from the active center) plus HINGE residues -- the places
to install a His-pair / His-triplet metal switch.

Method
------
1. Parse both states; match residues by (chain, resnum), keeping CA (and CB for
   side-chain direction). Mismatched/absent residues are dropped.
2. Rigid-core superposition (Kabsch + iterative outlier trimming): repeatedly
   superpose on the currently-inlier CA, recompute per-residue displacement, and
   drop residues that move more than `core_tol` * MAD above the core median. What
   remains is the rigid reference frame; everything else's displacement is then a
   genuine conformational response, not an alignment artifact. (Superposing on the
   whole protein would smear a domain motion across all residues.)
3. Per-residue Cα displacement `d_i` and distance-to-active-center `r_i`.
4. Distal candidate regions: contiguous runs of residues that are distal
   (r_i > distal_cutoff) AND responsive (d_i above the chosen percentile),
   merged and ranked by mean displacement.
5. Hinge residues: large local displacement gradient |d_{i+w} - d_{i-w}| while
   distal -- the loops that connect a moving domain to the rigid core. These are
   the highest-value clamp points for a metal switch.
"""
from __future__ import annotations
import numpy as np


# ----------------------------------------------------------------------------- parsing
def parse_chain(path: str, chain: str):
    """Return ordered dict-like: resnum -> {'resname','CA','CB'} for one chain."""
    import prody
    prody.confProDy(verbosity="none")
    ag = prody.parsePDB(path)
    sel = ag.select(f"protein and chain {chain}")
    if sel is None:                      # bundled/renamed chain fallback
        sel = ag.select("protein")
    res = {}
    for atom in sel:
        rn = int(atom.getResnum())
        name = atom.getName()
        if name not in ("CA", "CB"):
            continue
        entry = res.setdefault(rn, {"resname": atom.getResname(), "CA": None, "CB": None})
        entry[name] = np.array(atom.getCoords(), dtype=float)
    # keep only residues that have a CA
    return {k: v for k, v in sorted(res.items()) if v["CA"] is not None}


def match_residues(ref: dict, alt: dict):
    """Common residue numbers present (with CA) in both states."""
    common = sorted(set(ref) & set(alt))
    resnums = np.array(common, dtype=int)
    ref_ca = np.array([ref[r]["CA"] for r in common])
    alt_ca = np.array([alt[r]["CA"] for r in common])
    resnames = [ref[r]["resname"] for r in common]
    return resnums, resnames, ref_ca, alt_ca


# ----------------------------------------------------------------------------- superposition
def kabsch(mobile: np.ndarray, target: np.ndarray):
    """Rotation R and translation t minimizing ||mobile@R + t - target||."""
    mc, tc = mobile.mean(0), target.mean(0)
    H = (mobile - mc).T @ (target - tc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    R = (Vt.T @ D @ U.T).T            # apply as mobile@R
    t = tc - mc @ R
    return R, t


def iterative_core_superpose(mobile: np.ndarray, target: np.ndarray,
                             core_tol: float = 2.5, max_iter: int = 12):
    """Superpose `mobile` onto `target` on an auto-detected rigid core.

    Returns (mobile_aligned, displacement_per_residue, core_mask).
    """
    n = len(mobile)
    inliers = np.ones(n, dtype=bool)
    aligned = mobile.copy()
    for _ in range(max_iter):
        R, t = kabsch(mobile[inliers], target[inliers])
        aligned = mobile @ R + t
        d = np.linalg.norm(aligned - target, axis=1)
        med = np.median(d[inliers])
        mad = np.median(np.abs(d[inliers] - med)) + 1e-9
        new_inliers = d <= (med + core_tol * 1.4826 * mad)
        if new_inliers.sum() < max(8, int(0.15 * n)):   # never trim below 15%
            break
        if np.array_equal(new_inliers, inliers):
            break
        inliers = new_inliers
    return aligned, d, inliers


# ----------------------------------------------------------------------------- scoring
def smooth(x: np.ndarray, w: int = 2):
    if w <= 0:
        return x
    k = np.ones(2 * w + 1) / (2 * w + 1)
    return np.convolve(x, k, mode="same")


def analyze(ref_path, ref_chain, alt_path, alt_chain, active_site, cfg):
    """Run the full single-pair analysis. Returns a result dict."""
    ref = parse_chain(ref_path, ref_chain)
    alt = parse_chain(alt_path, alt_chain)
    resnums, resnames, ref_ca, alt_ca = match_residues(ref, alt)
    if len(resnums) < 10:
        raise RuntimeError(f"Only {len(resnums)} matched residues between the two states.")

    aligned, disp, core_mask = iterative_core_superpose(alt_ca, ref_ca)

    # active-center centroid (CA of active-site residues that exist in the match)
    asite_idx = [i for i, r in enumerate(resnums) if int(r) in set(active_site)]
    if asite_idx:
        center = ref_ca[asite_idx].mean(0)
    else:
        center = ref_ca.mean(0)
    r_to_active = np.linalg.norm(ref_ca - center, axis=1)

    disp_s = smooth(disp, w=cfg.get("smooth_w", 2))
    distal_cut = cfg["distal_cutoff_A"]
    disp_thr = np.percentile(disp_s, cfg["displacement_pctl"])

    distal = r_to_active > distal_cut
    responsive = disp_s >= disp_thr
    candidate = distal & responsive

    # contiguous candidate regions (in sequence)
    regions = _contiguous_regions(resnums, candidate, disp_s, r_to_active,
                                  min_len=cfg.get("min_region_len", 2))

    # hinge residues: local displacement gradient, distal, surface-ish
    w = cfg.get("hinge_w", 3)
    grad = np.zeros_like(disp_s)
    for i in range(len(disp_s)):
        lo, hi = max(0, i - w), min(len(disp_s) - 1, i + w)
        grad[i] = abs(disp_s[hi] - disp_s[lo]) / max(1, (hi - lo))
    hinge_score = grad * distal.astype(float)

    cb_lookup = {int(r): (ref[int(r)]["CB"] if ref[int(r)]["CB"] is not None else None)
                 for r in resnums}

    return {
        "resnums": resnums,
        "resnames": resnames,
        "ref_ca": ref_ca,
        "cb_lookup": cb_lookup,
        "alt_ca_aligned": aligned,
        "displacement": disp,
        "displacement_smooth": disp_s,
        "core_mask": core_mask,
        "r_to_active": r_to_active,
        "active_center": center,
        "candidate_mask": candidate,
        "regions": regions,
        "hinge_score": hinge_score,
        "disp_threshold": disp_thr,
        "displacement_pctl": cfg["displacement_pctl"],
        "distal_cutoff": distal_cut,
        "n_matched": len(resnums),
        "max_disp": float(disp.max()),
        "rmsd_core": float(np.sqrt((disp[core_mask] ** 2).mean())),
        "rmsd_all": float(np.sqrt((disp ** 2).mean())),
    }


def _contiguous_regions(resnums, mask, disp, rdist, min_len):
    regions = []
    i = 0
    n = len(resnums)
    while i < n:
        if mask[i]:
            j = i
            while j + 1 < n and mask[j + 1] and resnums[j + 1] == resnums[j] + 1:
                j += 1
            if (j - i + 1) >= min_len:
                regions.append({
                    "start": int(resnums[i]),
                    "end": int(resnums[j]),
                    "length": int(j - i + 1),
                    "mean_disp": float(disp[i:j + 1].mean()),
                    "max_disp": float(disp[i:j + 1].max()),
                    "min_dist_active": float(rdist[i:j + 1].min()),
                    "mean_dist_active": float(rdist[i:j + 1].mean()),
                })
            i = j + 1
        else:
            i += 1
    regions.sort(key=lambda d: d["mean_disp"], reverse=True)
    return regions


# ----------------------------------------------------------------------------- ground-truth check
def recovery_stats(result, known_allosteric, top_n_regions=20):
    """How many known allosteric residues fall inside the top candidate regions?"""
    if not known_allosteric:
        return None
    covered = set()
    for reg in result["regions"][:top_n_regions]:
        covered.update(range(reg["start"], reg["end"] + 1))
    hits = sorted(set(known_allosteric) & covered)
    return {
        "n_known": len(set(known_allosteric)),
        "n_recovered": len(hits),
        "recovered": hits,
        "fraction": len(hits) / max(1, len(set(known_allosteric))),
    }
