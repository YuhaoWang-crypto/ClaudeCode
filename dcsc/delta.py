"""
Brain-hole 2 — Δ-state signatures and Residual Pursuit ("recipe by vector
subtraction").

Library entry  = measured ΔS = S(stimulated) − S(control) in the 32-d program
                 space, plus its ternary signature (2 bits / program).
Query          = ΔS of (target sample − baseline sample).
Residual pursuit:
    r0 = ΔS
    loop: coarse candidates by ternary Hamming(sig(r), library)   [POPCNT-fast]
          pick the candidate with best cosine to r                 [exact rescoring]
          refit ALL chosen doses jointly by NNLS (orthogonal MP, so an early
          wrong dose is corrected later), r = ΔS − L_sel · dose
    stop when |r|/|ΔS| is small or stops improving.
A large final residual is reported, not hidden: it means a missing stimulus
or a non-additive interaction (synergy) that no linear recipe explains.
"""
import numpy as np
from scipy.optimize import nnls

from . import codec as C
from . import synth as S

STIMULI = {
    "IFNa":          {"IFNa": 1.2, "Antigen_pres": 0.6},
    "IFNb":          {"IFNa": 1.3, "Antigen_pres": 0.5, "Apoptosis": 0.3},
    "IFNg":          {"IFNg": 1.3, "Antigen_pres": 1.0, "IFNa": 0.4},
    "TNFa":          {"TNF_NFkB": 1.3, "Apoptosis": 0.5, "Inflammasome": 0.4},
    "IL1b":          {"TNF_NFkB": 0.9, "Inflammasome": 1.0, "IL6_JAK_STAT3": 0.5},
    "IL6":           {"IL6_JAK_STAT3": 1.2, "Complement": 0.5},
    "IL10":          {"IL10_STAT3": 1.2, "Antigen_pres": -0.8, "TNF_NFkB": -0.5},
    "IL4":           {"IL4_STAT6": 1.2, "Type2_cytokine": 0.8},
    "IL13":          {"IL4_STAT6": 1.0, "Type2_cytokine": 1.0},
    "TGFb1":         {"TGFb": 1.3, "EMT": 0.9, "Prolif_G2M": -0.6},
    "Hypoxia_1pct":  {"Hypoxia": 1.4, "Glycolysis": 1.0, "Oxphos": -0.7},
    "LPS":           {"TNF_NFkB": 1.2, "Inflammasome": 0.8, "IL6_JAK_STAT3": 0.7, "IFNa": 0.5},
    "Dexamethasone": {"TNF_NFkB": -0.9, "IL6_JAK_STAT3": -0.6, "Apoptosis": 0.4},
    "Rapamycin":     {"mTORC1": -1.2, "Autophagy": 0.9, "Prolif_S": -0.5},
    "Tunicamycin":   {"ER_stress_UPR": 1.5, "Apoptosis": 0.6},
    "Thapsigargin":  {"ER_stress_UPR": 1.3, "Heat_shock": 0.5, "NRF2_oxidative": 0.4},
    "Doxorubicin":   {"p53_DNAdamage": 1.4, "Prolif_G2M": -1.0, "Apoptosis": 0.6},
    "Nutlin3":       {"p53_DNAdamage": 1.2, "Prolif_S": -0.7},
    "CHIR99021":     {"Wnt": 1.3},
    "DAPT":          {"Notch": -1.2},
    "SAG":           {"Hedgehog": 1.2},
    "DHT":           {"Androgen": 1.3},
    "Estradiol":     {"Estrogen": 1.3, "Prolif_S": 0.4},
    "Simvastatin":   {"Cholesterol": 1.2, "mTORC1": -0.3},
    "Heat_42C":      {"Heat_shock": 1.5, "ER_stress_UPR": 0.4},
    "FICZ":          {"AhR": 1.3},
    "Sulforaphane":  {"NRF2_oxidative": 1.3},
    "Etoposide_sen": {"SASP": 1.2, "p53_DNAdamage": 0.6, "Prolif_G2M": -0.8},
    "C5a":           {"Complement": 1.0, "Inflammasome": 0.5},
    "Oligomycin":    {"Oxphos": -1.3, "Glycolysis": 0.8},
    "2DG":           {"Glycolysis": -1.2, "ER_stress_UPR": 0.4},
    "EGF":           {"Prolif_S": 0.9, "Prolif_G2M": 0.8, "mTORC1": 0.6},
    "Palbociclib":   {"Prolif_S": -1.2, "Prolif_G2M": -1.0},
    "IL2":           {"Prolif_S": 0.8, "IL6_JAK_STAT3": 0.3, "mTORC1": 0.5},
}
NAMES = list(STIMULI)
TAU_DELTA = 0.4


def effect_matrix():
    E = np.zeros((len(NAMES), S.P))
    for i, n in enumerate(NAMES):
        for prog, v in STIMULI[n].items():
            E[i, S.PROGRAMS.index(prog)] = v
    return E


def pseudobulk_delta(world, enc, leaf, a_target, rng, n=150, env_shift=None):
    """ΔS (40-d posterior mean) of n target cells vs n baseline cells."""
    z0 = np.zeros
    b = rng.integers(0, world.n_batches)
    e_t = np.zeros((n, S.NE)) if env_shift is None else np.tile(env_shift, (n, 1))
    Xc = world.sample_norm(np.full(n, leaf), z0((n, S.P)), z0((n, S.NM)), z0((n, S.NE)),
                           np.full(n, b), rng).mean(0)
    Xt = world.sample_norm(np.full(n, leaf), np.tile(a_target, (n, 1)), z0((n, S.NM)), e_t,
                           np.full(n, b), rng).mean(0)
    z = enc.encode(np.vstack([Xc, Xt]), leaf=[leaf, leaf], n_cells=n)["z"]
    naive = enc.naive_scores(np.vstack([Xc, Xt]), [leaf, leaf])
    return z[1] - z[0], naive[1] - naive[0]


def build_library(world, enc, rng, leaf=0, n=300):
    E = effect_matrix()
    L = np.array([pseudobulk_delta(world, enc, leaf, E[i], rng, n)[0][C.SL_MOD] for i in range(len(NAMES))])
    return L, signature(L)


def signature(D):
    """Ternary Δ-signature packed in the same 2-bitplane layout as the code."""
    t = np.sign(D) * (np.abs(D) > TAU_DELTA)
    z = np.zeros((len(D), S.K_LATENT))
    z[:, :S.P] = t
    code = np.zeros((len(D), C.CODE_BYTES), np.uint8)
    pos = np.zeros((len(D), C.MOD_CAP), bool)
    neg = np.zeros((len(D), C.MOD_CAP), bool)
    pos[:, :S.P], neg[:, :S.P] = t > 0, t < 0
    code[:, 16:32] = np.packbits(pos, axis=1, bitorder="little")
    code[:, 32:48] = np.packbits(neg, axis=1, bitorder="little")
    return code


def residual_pursuit(delta, L, Lsig, max_terms=3, prefilter=8, stop_frac=0.25,
                     min_gain=0.05, min_dose=0.2):
    from .index import mod_distance
    norm0 = np.linalg.norm(delta) + 1e-12
    sel, r, dose, hist = [], delta.copy(), np.zeros(0), []
    Ln = L / (np.linalg.norm(L, axis=1, keepdims=True) + 1e-12)
    for _ in range(max_terms):
        if Lsig is None:            # no coarse stage: rescore every entry
            cand = [c for c in range(len(L)) if c not in sel]
        else:
            ham = mod_distance(Lsig, signature(r[None])[0])
            cand = [c for c in np.argsort(ham, kind="stable") if c not in sel][:prefilter]
        best = max(cand, key=lambda c: Ln[c] @ r)
        if Ln[best] @ r <= 0:
            break
        trial = sel + [best]
        d, _ = nnls(L[trial].T, delta)
        frac = np.linalg.norm(delta - L[trial].T @ d) / norm0
        if hist and hist[-1] - frac < min_gain:
            break
        sel, dose, hist = trial, d, hist + [frac]
        r = delta - L[sel].T @ dose
        if frac < stop_frac:
            break
    keep = dose >= min_dose
    sel = [s for s, k in zip(sel, keep) if k]
    if sel:
        dose, _ = nnls(L[sel].T, delta)
        frac = np.linalg.norm(delta - L[sel].T @ dose) / norm0
    else:
        dose, frac = np.zeros(0), 1.0
    return dict(stimuli=sel, dose=dose, unexplained=frac, path=hist)


def equivalence(L, thr=0.9):
    Ln = L / np.linalg.norm(L, axis=1, keepdims=True)
    return (Ln @ Ln.T) > thr
