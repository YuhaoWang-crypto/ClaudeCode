"""
bioif.real.copula -- the cheap joint model, as a production edge.

§8.3 of INTEROP.md established that the p53 -> mutagenicity edge is
CONDITIONAL: it carries information only in compounds that are not already
flagged cytotoxic (OR 2.90 outside the cytotoxic stratum, OR 1.13 and p=0.77
inside it). The adapter built there handles that by switching between three
measured 2x2 tables, which works but throws away most of what is known: it
uses the cytotoxicity readout as a binary stratifier and nothing else.

A Gaussian copula over the Tox21 endpoint vector does the same job properly.
Fit 12 independent logistic marginals, then estimate the 66 free
correlations of a latent multivariate probit; P(p53 | fingerprint, MMP =
measured) then has a closed form. That gives, measured on held-out
scaffolds:

    ✅ p53 | fingerprint          AUROC 0.8081   AP 0.2940
    ✅ p53 | fingerprint + MMP    AUROC 0.8327   AP 0.3369
       gain from the measured MMP      +0.0246       +0.0430

on 1,512 test compounds with both labels observed, for **66 parameters**.
The joint-diffusion experiment in `diffusion.py` was run partly to see
whether a generative model could do better than this; it could not, and the
copula is why that comparison was worth making rather than assumed.

Three reasons this is the right shape for the contract, not just a better
number:

  * it is **closed form**, so the conditional is exact rather than sampled,
    and the edge costs microseconds instead of a GPU;
  * it conditions on **any** subset of measured endpoints, not only the one
    pair someone hard-coded. The stratified 2x2 answers "p53 given MMP"; the
    copula also answers "p53 given MMP and ATAD5", which is the same
    question the chain will ask next;
  * it degrades gracefully to the marginal when nothing is measured, which
    is the common case, so the edge has one implementation rather than two.

⚠️ It is still an association model over assay endpoints. A latent
correlation is not a mechanism, and `rho(p53, MMP) = 0.328` says these
readouts co-occur, not that cytotoxicity causes reporter activation.

NOTE ON DUPLICATION: the probit primitives below also exist in
`diffusion.py`, which is where they were first written and validated. This
module is their intended home; `diffusion.py` keeps its copy only while that
experiment is still being re-run, and should import from here once it
settles. The two implementations are identical and the test suite checks
they agree.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from scipy.special import ndtr, ndtri

from . import qsar, tox

SEED = 0
N_BITS = 2048
ENDPOINTS = tox.TOX21_ENDPOINTS
NE = len(ENDPOINTS)
P53_J = ENDPOINTS.index(tox.P53)
MMP_J = ENDPOINTS.index(tox.CYTOTOX_CONTROL)
ATAD5_J = ENDPOINTS.index(tox.DDR_CONTROL)

#: Fitting the correlations costs a few minutes of orthant integration, so
#: the fitted artefact is cached. It is a deterministic function of a
#: committed snapshot, so caching it hides nothing.
ARTEFACT = tox.SNAPSHOT / "copula_p53.npz"

# --------------------------------------------------------------------------
# Probit primitives
# --------------------------------------------------------------------------

_GL_X, _GL_W = np.polynomial.legendre.leggauss(48)
_GL_NODE = 0.5 * (_GL_X + 1.0)
_GL_WT = 0.5 * _GL_W


def bvn_sf(t1, t2, rho: float) -> np.ndarray:
    """
    P(Z1 > t1, Z2 > t2) for a standard bivariate normal with correlation rho.

    Integrated after the substitution v = Phi(u), which puts the quadrature
    on a fixed [0,1] grid and so vectorises over compounds, each of which has
    its own thresholds. Checked against scipy's multivariate_normal.cdf.
    """
    t1 = np.asarray(t1, float)
    t2 = np.asarray(t2, float)
    s = math.sqrt(max(1e-12, 1.0 - rho * rho))
    lo = ndtr(t1)[..., None]
    v = np.clip(lo + (1.0 - lo) * _GL_NODE, 1e-12, 1 - 1e-12)
    integ = ndtr(-(t2[..., None] - rho * ndtri(v)) / s)
    return ndtr(-t1) * (integ * _GL_WT).sum(-1)


def thresholds(P: np.ndarray) -> np.ndarray:
    """Latent cut points implied by per-compound marginal probabilities."""
    return ndtri(1.0 - np.clip(P, 1e-6, 1 - 1e-6))


def fit_corr(P: np.ndarray, Y: np.ndarray, M: np.ndarray, grid: int = 31,
             cap: int = 2500, seed: int = SEED) -> np.ndarray:
    """
    Pairwise composite-likelihood estimate of the latent correlation matrix,
    then projection onto the nearest PSD correlation matrix.

    Pairwise rather than full-likelihood: with per-compound thresholds the
    full 12-dimensional probit likelihood needs a 12-d orthant probability
    per compound, while the pairwise estimator is consistent and is all the
    conditionals actually use.
    """
    T = thresholds(P)
    rng = np.random.default_rng(seed)
    rhos = np.linspace(-0.95, 0.95, grid)
    R = np.eye(NE)
    for j in range(NE):
        for k in range(j + 1, NE):
            both = np.where(M[:, j] & M[:, k])[0]
            if len(both) > cap:
                both = rng.choice(both, cap, replace=False)
            if len(both) < 50:
                continue                      # too few pairs: leave at 0
            t1, t2 = T[both, j], T[both, k]
            y1, y2 = Y[both, j] == 1, Y[both, k] == 1
            q1, q2 = ndtr(-t1), ndtr(-t2)
            ll = np.empty(grid)
            for g, rho in enumerate(rhos):
                p11 = np.clip(bvn_sf(t1, t2, rho), 1e-12, 1.0)
                p10 = np.clip(q1 - p11, 1e-12, 1.0)
                p01 = np.clip(q2 - p11, 1e-12, 1.0)
                p00 = np.clip(1.0 - q1 - q2 + p11, 1e-12, 1.0)
                ll[g] = np.log(np.where(y1 & y2, p11,
                               np.where(y1 & ~y2, p10,
                               np.where(~y1 & y2, p01, p00)))).sum()
            g = int(np.argmax(ll))
            rho = rhos[g]
            if 0 < g < grid - 1:                          # parabolic refine
                a, b, c = ll[g - 1], ll[g], ll[g + 1]
                den = a - 2 * b + c
                if den < 0:
                    rho += 0.5 * (a - c) / den * (rhos[1] - rhos[0])
            R[j, k] = R[k, j] = float(np.clip(rho, -0.97, 0.97))
    w, V = np.linalg.eigh(R)
    R = V @ np.diag(np.clip(w, 1e-4, None)) @ V.T
    dg = np.sqrt(np.diag(R))
    return R / np.outer(dg, dg)


def conditional(P: np.ndarray, R: np.ndarray, j: int, k: int,
                y_k: np.ndarray) -> np.ndarray:
    """P(y_j = 1 | y_k = observed), closed form, one conditioning endpoint."""
    T = thresholds(P)
    tj, tk = T[:, j], T[:, k]
    rho = float(R[j, k])
    p11 = bvn_sf(tj, tk, rho)
    qj, qk = ndtr(-tj), ndtr(-tk)
    out = np.where(np.asarray(y_k) == 1,
                   p11 / np.clip(qk, 1e-9, None),
                   (qj - p11) / np.clip(1.0 - qk, 1e-9, None))
    return np.clip(out, 0.0, 1.0)


# --------------------------------------------------------------------------
# The fitted model
# --------------------------------------------------------------------------

@dataclass
class CopulaP53:
    """Marginal logistic coefficients plus the latent correlation matrix."""
    coef: np.ndarray            # (12, n_bits)
    intercept: np.ndarray       # (12,)
    R: np.ndarray               # (12, 12)
    q: dict                     # label-conditional conformal thresholds
    n_train: int = 0
    n_cal: int = 0
    note: str = ""
    #: fraction of prediction sets that are NOT a singleton, per regime.
    #: Same definition the plain QSAR adapter reports, so the registry can
    #: rank the two models of this hop against each other.
    abstain: dict = field(default_factory=dict)

    # -- marginals ---------------------------------------------------------
    def _marginals(self, X: np.ndarray) -> np.ndarray:
        z = X @ self.coef.T + self.intercept
        return 1.0 / (1.0 + np.exp(-z))

    def proba(self, smiles: list[str], given: dict[str, int] | None = None,
              j: int = P53_J) -> np.ndarray:
        """
        P(endpoint j active | structure, and any measured endpoints in `given`).

        `given` maps a Tox21 endpoint name to its measured 0/1 value. With
        nothing given this is the marginal; with one endpoint given it is the
        exact probit conditional. More than one is applied sequentially,
        which is an approximation and is flagged as such by
        `conditioning_is_exact`.
        """
        X, ok = qsar.featurize(list(smiles), n_bits=N_BITS)
        P = self._marginals(X)
        if not given:
            return P[:, j]
        out = P[:, j]
        for name, val in given.items():
            k = ENDPOINTS.index(name)
            Pk = P.copy()
            Pk[:, j] = out                 # carry the running conditional
            out = conditional(Pk, self.R, j, k,
                              np.full(len(X), int(val)))
        return out

    @staticmethod
    def conditioning_is_exact(given: dict | None) -> bool:
        """Exact for zero or one conditioning endpoint; sequential beyond."""
        return not given or len(given) <= 1

    # -- conformal ---------------------------------------------------------
    def predict_set(self, smiles: list[str],
                    given: dict[str, int] | None = None,
                    j: int = P53_J) -> list[frozenset]:
        key = "cond" if given else "marg"
        p1 = self.proba(smiles, given, j)
        out = []
        for p in p1:
            s = set()
            if p <= self.q[f"{key}_0"] + 1e-12:          # s(x,0) = 1 - P(0) = p
                s.add(0)
            if (1.0 - p) <= self.q[f"{key}_1"] + 1e-12:
                s.add(1)
            out.append(frozenset(s))
        return out

    def rho(self, a: str, b: str) -> float:
        return float(self.R[ENDPOINTS.index(a), ENDPOINTS.index(b)])

    def save(self, path: Path = ARTEFACT) -> Path:
        np.savez_compressed(path, coef=self.coef.astype(np.float32),
                            intercept=self.intercept.astype(np.float32),
                            R=self.R, n_train=self.n_train, n_cal=self.n_cal,
                            q_keys=np.array(list(self.q.keys())),
                            q_vals=np.array(list(self.q.values()), float),
                            a_keys=np.array(list(self.abstain.keys())),
                            a_vals=np.array(list(self.abstain.values()), float),
                            note=np.array(self.note))
        return path


def _conformal_q(scores: list[float], alpha: float) -> float:
    n = len(scores)
    if n == 0:
        return math.inf
    k = math.ceil((n + 1) * (1.0 - alpha))
    return math.inf if k > n else sorted(scores)[k - 1]


def fit(alpha: float = 0.1, seed: int = SEED) -> tuple[CopulaP53, dict]:
    """
    Fit on train scaffolds, estimate R from CROSS-FITTED train probabilities,
    calibrate conformal on the calibration fold, evaluate on test.

    Cross-fitting matters: in-sample marginals explain each label too well on
    their own, which drives the latent thresholds to extremes and attenuates
    the estimated correlation towards zero -- i.e. it would quietly destroy
    the very signal this model exists to use.
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import average_precision_score, roc_auc_score

    rows = tox.load_tox21()
    smiles = [r["smiles"] for r in rows]
    X, ok = qsar.featurize(smiles, n_bits=N_BITS)
    Y = np.zeros((len(rows), NE), dtype=np.int8)
    M = np.zeros((len(rows), NE), dtype=bool)
    for i, r in enumerate(rows):
        for jj, e in enumerate(ENDPOINTS):
            v = r[e]
            if v not in ("", "NA"):
                M[i, jj] = True
                Y[i, jj] = int(float(v))
    keep = np.where(ok)[0]
    X, Y, M = X[keep], Y[keep], M[keep]
    smiles = [smiles[i] for i in keep]
    scafs = np.array([qsar.scaffold(s) or f"__s{i}" for i, s in enumerate(smiles)],
                     dtype=object)
    tr, ca, te = qsar.scaffold_split(smiles, fracs=(0.5, 0.2, 0.3), seed=seed)

    coef = np.zeros((NE, N_BITS), dtype=np.float64)
    inter = np.zeros(NE)
    for jj in range(NE):
        rowsj = [i for i in tr if M[i, jj]]
        m = LogisticRegression(max_iter=4000, class_weight="balanced")
        m.fit(X[rowsj], Y[rowsj, jj])
        coef[jj], inter[jj] = m.coef_[0], m.intercept_[0]
    model = CopulaP53(coef, inter, np.eye(NE), {}, len(tr), len(ca))

    # cross-fitted train probabilities, scaffold-grouped
    uniq = sorted({scafs[i] for i in tr})
    rng = np.random.default_rng(seed)
    assign = dict(zip(uniq, rng.integers(0, 2, len(uniq))))
    P_oof = np.zeros((len(tr), NE))
    for f in (0, 1):
        fit_rows = [i for i in tr if assign[scafs[i]] != f]
        pred_pos = [a for a, i in enumerate(tr) if assign[scafs[i]] == f]
        if not pred_pos:
            continue
        c2 = np.zeros((NE, N_BITS)); i2 = np.zeros(NE)
        for jj in range(NE):
            rj = [i for i in fit_rows if M[i, jj]]
            mm = LogisticRegression(max_iter=4000, class_weight="balanced")
            mm.fit(X[rj], Y[rj, jj])
            c2[jj], i2[jj] = mm.coef_[0], mm.intercept_[0]
        z = X[[tr[a] for a in pred_pos]] @ c2.T + i2
        P_oof[pred_pos] = 1.0 / (1.0 + np.exp(-z))
    model.R = fit_corr(P_oof, Y[tr], M[tr], seed=seed)

    # conformal, label-conditional, both with and without a measured MMP
    P_ca = model._marginals(X[ca])
    cal_rows = [a for a, i in enumerate(ca) if M[i, P53_J]]
    y_ca = Y[[ca[a] for a in cal_rows], P53_J]
    p_marg = P_ca[cal_rows, P53_J]
    mmp_ok = [a for a in cal_rows if M[ca[a], MMP_J]]
    y_cond = Y[[ca[a] for a in mmp_ok], P53_J]
    p_cond = conditional(P_ca[mmp_ok], model.R, P53_J, MMP_J,
                         Y[[ca[a] for a in mmp_ok], MMP_J])
    q = {}
    for tag, pv, yv in (("marg", p_marg, y_ca), ("cond", p_cond, y_cond)):
        for cls in (0, 1):
            sc = [1.0 - (1.0 - p if cls == 0 else p)
                  for p, y in zip(pv, yv) if y == cls]
            q[f"{tag}_{cls}"] = _conformal_q(sc, alpha)
    model.q = q
    model.note = f"alpha={alpha} seed={seed} scaffold-split"

    # -- evaluation on test, scored once ----------------------------------
    P_te = model._marginals(X[te])
    both = [a for a, i in enumerate(te) if M[i, P53_J] and M[i, MMP_J]]
    y = Y[[te[a] for a in both], P53_J].astype(int)
    mmp = Y[[te[a] for a in both], MMP_J].astype(int)
    unc = P_te[both, P53_J]
    con = conditional(P_te[both], model.R, P53_J, MMP_J, mmp)
    ev = {
        "n_test_both": len(both),
        "prevalence": float(y.mean()),
        "auroc_marginal": float(roc_auc_score(y, unc)),
        "ap_marginal": float(average_precision_score(y, unc)),
        "auroc_conditional": float(roc_auc_score(y, con)),
        "ap_conditional": float(average_precision_score(y, con)),
        "rho_p53_mmp": model.rho(tox.P53, tox.CYTOTOX_CONTROL),
        "n_free_params": NE * (NE - 1) // 2,
        "n_train": len(tr), "n_cal": len(ca), "n_test": len(te),
    }
    ev["auroc_gain"] = ev["auroc_conditional"] - ev["auroc_marginal"]
    ev["ap_gain"] = ev["ap_conditional"] - ev["ap_marginal"]

    # Abstention rate, in BOTH regimes. This is what the registry ranks on,
    # and it has to be the same quantity the plain QSAR adapter reports
    # (fraction of prediction sets that are not a singleton) -- otherwise two
    # models of one hop would be compared on different scales, which is the
    # error this whole package exists to prevent.
    for tag, idx, given_vals in (("marg", both, None), ("cond", both, mmp)):
        sets = []
        for a, i in enumerate(idx):
            pp = con[a] if given_vals is not None else unc[a]
            s0 = pp <= model.q[f"{tag}_0"] + 1e-12
            s1 = (1.0 - pp) <= model.q[f"{tag}_1"] + 1e-12
            sets.append(int(s0) + int(s1))
        ev[f"abstain_{tag}"] = float(np.mean([z != 1 for z in sets]))
    model.abstain = {"marg": ev["abstain_marg"], "cond": ev["abstain_cond"]}
    return model, ev


def load(path: Path = ARTEFACT) -> CopulaP53:
    if not path.exists():
        raise FileNotFoundError(
            f"no fitted copula at {path}; run python3 -m bioif.real.copula")
    z = np.load(path, allow_pickle=False)
    q = dict(zip([str(k) for k in z["q_keys"]], z["q_vals"].tolist()))
    ab = dict(zip([str(k) for k in z["a_keys"]], z["a_vals"].tolist())) \
        if "a_keys" in z else {}
    return CopulaP53(z["coef"].astype(np.float64), z["intercept"].astype(np.float64),
                     z["R"], q, int(z["n_train"]), int(z["n_cal"]),
                     str(z["note"]), ab)


def report(ev: dict | None = None) -> str:
    if ev is None:
        _, ev = fit()
    L = ["The cheap joint model: p53 conditioned on a measured cytotoxicity call",
         f"  scaffold split train/cal/test = {ev['n_train']}/{ev['n_cal']}/"
         f"{ev['n_test']}",
         f"  scored on {ev['n_test_both']} test compounds with BOTH p53 and "
         f"MMP observed (prevalence {ev['prevalence']:.3f})",
         f"  fitted rho(p53, MMP) = {ev['rho_p53_mmp']:.3f}   "
         f"free parameters: {ev['n_free_params']}",
         "",
         f"  {'predictor':<34}{'AUROC':>8}{'AP':>8}",
         f"  {'p53 | fingerprint':<34}{ev['auroc_marginal']:>8.4f}"
         f"{ev['ap_marginal']:>8.4f}",
         f"  {'p53 | fingerprint + measured MMP':<34}"
         f"{ev['auroc_conditional']:>8.4f}{ev['ap_conditional']:>8.4f}",
         f"  {'gain from the measurement':<34}{ev['auroc_gain']:>+8.4f}"
         f"{ev['ap_gain']:>+8.4f}",
         "",
         f"  conformal abstention rate: {ev['abstain_marg']:.3f} marginal, "
         f"{ev['abstain_cond']:.3f} conditioned",
         "",
         "Reading:",
         f"  {ev['n_free_params']} correlation parameters turn the cytotoxicity readout from a",
         "  binary stratifier into a conditioning variable, and buy "
         f"{ev['auroc_gain']:+.3f} AUROC",
         f"  / {ev['ap_gain']:+.3f} AP. The joint-diffusion model in diffusion.py was",
         "  given the same information through inpainting and gained +0.000.",
         "",
         "  ⚠️ Still an association over assay endpoints: rho says these readouts",
         "  co-occur, not that cytotoxicity causes reporter activation."]
    return "\n".join(L)


if __name__ == "__main__":
    m, ev = fit()
    m.save()
    print(report(ev))
    print(f"\nfitted artefact written to {ARTEFACT}")
