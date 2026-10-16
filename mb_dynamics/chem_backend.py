"""SWV current for conformation x multi-species redox chemistry (round 2: stepwise MB chemistry).

Chemistries (species along a chain, oxidized -> reduced):
  C0  O --2e--> R                       concerted apparent 2e- Butler-Volmer (identical to round 1)
  C1  O --e--> I --e--> R               EE; protonation fast, folded into apparent potentials
  C2  O --e--> I --H+--> IH --e--> R    E-C-E; semiquinone protonation with finite rate k_p, fixed K_p
Electron-transfer edges in conformational state i run at k * w_i (w_i = <exp(-beta (z - z_min))> in i);
the chemical (protonation) edge is not distance-dependent.  O, I, IH, R share pi and L.

Every edge obeys detailed balance and the landscape is shared, so the joint chain is reversible with
respect to P(E) = pi (x) q(E).  With u = p - P(E) and v = P^{-1/2} u, a constant-potential half-pulse
obeys dv/dt = -G^T G v, where G stacks one row per edge:
  conformational edge (i, j), species s:  sqrt(F_ij) (e_is / sqrt(pi_i) - e_js / sqrt(pi_j))
  reaction edge a -> b in conformation i: sqrt(k_ab(i)) e_ia - sqrt(k_ba(i)) e_ib
(no q in either row).  The spectrum comes from the SVD of G (accurate when exchange is many orders faster
than reaction).  On a potential change p is continuous, so u jumps by P_old - P_new.  Current per
electron-transfer edge = k_ab u_ia - k_ba u_ib (equilibrium cancels exactly); I = (Q/2) sum n_e (...).
q is handled in log space, so no species population under- or overflows.
"""
import numpy as np

import swv_backend as sb

FN1 = sb.F_CONST / (sb.R_GAS * sb.T_DEFAULT)      # F/RT for one electron


class Chemistry:
    """name in {C0, C1, C2}; params dict: see build_edges."""

    def __init__(self, name, K_p=10.0):
        self.name = name
        self.K_p = K_p
        self.species = {"C0": ["O", "R"], "C1": ["O", "I", "R"], "C2": ["O", "I", "IH", "R"]}[name]
        self.S = len(self.species)

    def edges(self, p, E):
        """Per-edge (a, b, is_et, n_e, log k_ab, log k_ba, log(q_b/q_a)) on the potential array E (V)."""
        a = p["alpha"]
        out = []
        if self.name == "C0":
            eta = E - p["E0"]
            lk = np.log(p["k"])
            out.append((0, 1, True, 2, lk - a * 2 * FN1 * eta, lk + (1 - a) * 2 * FN1 * eta, -2 * FN1 * eta))
            return out
        if self.name == "C1":
            E1, E2 = p["Em"] + p["dE"] / 2, p["Em"] - p["dE"] / 2
            chain = [(0, 1, p["k1"], E1), (1, 2, p["k2"], E2)]
        else:
            lnK = np.log(self.K_p)
            s = 2 * p["Em"] - lnK / FN1           # E1 + E2'
            E1, E2 = (s + p["dE"]) / 2, (s - p["dE"]) / 2
            chain = [(0, 1, p["k1"], E1), (2, 3, p["k2"], E2)]
            lkp = np.log(p["kp"])
            out.append((1, 2, False, 0, lkp + 0 * E, lkp - lnK + 0 * E, lnK + 0 * E))
        for ia, ib, k, Ee in chain:
            eta = E - Ee
            lk = np.log(k)
            out.append((ia, ib, True, 1, lk - a * FN1 * eta, lk + (1 - a) * FN1 * eta, -FN1 * eta))
        return out

    def logq(self, p, E):
        """log equilibrium fractions (K, S) along the chain."""
        rel = np.zeros((len(E), self.S))
        for a, b, _, _, _, _, lr in sorted(self.edges(p, E), key=lambda e: e[0]):
            rel[:, b] = rel[:, a] + lr
        m = rel.max(1, keepdims=True)
        return rel - (m + np.log(np.exp(rel - m).sum(1, keepdims=True)))


def bank(kin, chem, p, mv):
    E = mv / 1000.0
    K, M, S = len(E), kin.M, chem.S
    edges = chem.edges(p, E)
    lq = chem.logq(p, E)
    lpi = np.log(kin.pi)
    w = kin.w
    Bc = np.kron(kin.B, np.eye(S)) if kin.B.shape[0] else np.zeros((0, M * S))
    nr = M * len(edges)
    G = np.zeros((K, Bc.shape[0] + nr, M * S))
    G[:, :Bc.shape[0], :] = Bc[None]
    c = np.zeros((K, M * S))
    r = Bc.shape[0]
    for (a, b, is_et, n_e, lkab, lkba, _) in edges:
        for i in range(M):
            lw = np.log(w[i]) if is_et else 0.0
            G[:, r, i * S + a] = np.exp(0.5 * (lkab + lw))
            G[:, r, i * S + b] = -np.exp(0.5 * (lkba + lw))
            r += 1
            if is_et:
                c[:, i * S + a] += 0.5 * n_e * np.exp(lkab + lw + 0.5 * (lpi[i] + lq[:, a]))
                c[:, i * S + b] -= 0.5 * n_e * np.exp(lkba + lw + 0.5 * (lpi[i] + lq[:, b]))
    _, sig, Wt = np.linalg.svd(G, full_matrices=False)
    lP = (lpi[None, :, None] + lq[:, None, :]).reshape(K, M * S)     # log(pi_i q_s), index i*S + s
    return dict(rate=sig ** 2, W=np.transpose(Wt, (0, 2, 1)), c=c, lP=lP)


def simulate(kin, chem, p, wf_list, Q=1.0):
    mv_all = np.unique(np.concatenate([np.r_[w.mv_fwd, w.mv_bwd, int(round(w.e_start * 1000))] for w in wf_list]))
    bk = bank(kin, chem, p, mv_all)
    idx = {int(m): i for i, m in enumerate(mv_all)}
    out = []
    for wf in wf_list:
        h = wf.h
        decay = np.exp(-bk["rate"] * h)
        j_prev = idx[int(round(wf.e_start * 1000))]
        v = np.zeros(kin.M * chem.S)
        nk = len(wf.e_stair)
        i_f, i_b = np.empty(nk), np.empty(nk)
        for k in range(nk):
            for half, mvk in ((0, wf.mv_fwd[k]), (1, wf.mv_bwd[k])):
                j = idx[int(mvk)]
                if j != j_prev:
                    lo, ln = bk["lP"][j_prev], bk["lP"][j]
                    v = v * np.exp(0.5 * (lo - ln)) + np.exp(lo - 0.5 * ln) - np.exp(0.5 * ln)
                    j_prev = j
                W = bk["W"][j]
                v = W @ (decay[j] * (W.T @ v))
                cur = np.dot(bk["c"][j], v)
                (i_f if half == 0 else i_b)[k] = cur
        out.append(dict(E=wf.e_stair, fwd=Q * i_f, bwd=Q * i_b, net=Q * (i_f - i_b)))
    return out


def full_generator(kin, chem, p, E):
    """Independent dense column generator of the joint chain at one potential (for tests)."""
    M, S = kin.M, chem.S
    sq = np.sqrt(kin.pi)
    L = (kin.S * sq[:, None]) / sq[None, :]
    Gm = np.kron(L, np.eye(S))
    for (a, b, is_et, n_e, lkab, lkba, _) in chem.edges(p, np.array([E])):
        for i in range(M):
            f = kin.w[i] if is_et else 1.0
            kab, kba = np.exp(lkab[0]) * f, np.exp(lkba[0]) * f
            ia, ib = i * S + a, i * S + b
            Gm[ib, ia] += kab
            Gm[ia, ia] -= kab
            Gm[ia, ib] += kba
            Gm[ib, ib] -= kba
    return Gm


def current_from_p(kin, chem, p, E, pvec):
    I = 0.0
    S = chem.S
    for (a, b, is_et, n_e, lkab, lkba, _) in chem.edges(p, np.array([E])):
        if not is_et:
            continue
        for i in range(kin.M):
            I += 0.5 * n_e * kin.w[i] * (np.exp(lkab[0]) * pvec[i * S + a] - np.exp(lkba[0]) * pvec[i * S + b])
    return I


def equilibrium(kin, chem, p, E):
    lq = chem.logq(p, np.array([E]))[0]
    return (kin.pi[:, None] * np.exp(lq)[None, :]).ravel()
