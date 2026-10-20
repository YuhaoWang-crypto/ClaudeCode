"""
Two-stage retrieval: compact codes for the coarse pass, full vectors for the
exact pass.

    stage 1  score every library item using ONLY its compact code
             -> shortlist of M candidates
    stage 2  rescore the M candidates on their full-precision vectors
             (cosine for k-NN; residual pursuit for "cocktail" queries)

Stage-1 code formats (all built from the same K-dim program basis W):

    "ternary"  2 bits / program   {-1,0,+1} at |z| > tau_j      (K/4 bytes)
    "int8"     8 bits / program   per-dimension scaled int8     (K bytes)
    "pq"       faiss product quantizer on the full vector       (m bytes)

and two ways to compare a query with a code:

    symmetric   quantize the query too, compare code to code (Hamming)
    asymmetric  keep the query at full precision, compare it to the
                decoded code (ADC). Costs nothing extra at query time and
                is the default.

For pursuit, the shortlist is refreshed on the current residual at every
step, so a weak second component that the first one masked can still enter.

Only stage 2 touches full vectors: they can live on disk / object storage and
be fetched for M items, while the codes stay in RAM.
"""
import time

import numpy as np

from . import delta as DL


def _pack_ternary(T):
    """T in {-1,0,1}^(n,K) -> uint8 bit-planes (pos | neg), K padded to 64."""
    K = T.shape[1]
    Kp = ((K + 63) // 64) * 64
    pos = np.zeros((len(T), Kp), bool)
    neg = np.zeros((len(T), Kp), bool)
    pos[:, :K], neg[:, :K] = T > 0, T < 0
    return np.hstack([np.packbits(pos, 1, bitorder="little"), np.packbits(neg, 1, bitorder="little")])


class TwoStageIndex:
    def __init__(self, full, K=32, code="int8", basis=None, tau_q=0.67, pq_m=8, fit_rows=None):
        """full: (N, G) float32 vectors kept for stage 2 (library Δ or cell profiles).
        basis: (G, K) orthonormal program loadings; default = top-K right
        singular vectors of `full` (or of full[fit_rows])."""
        self.full = np.ascontiguousarray(full, dtype=np.float32)
        self.norm = np.linalg.norm(self.full, axis=1) + 1e-12
        self.code = code
        if code == "pq":
            import faiss
            self.pad = (-full.shape[1]) % pq_m           # PQ needs d % m == 0
            d = full.shape[1] + self.pad
            self.pq = faiss.IndexPQ(d, pq_m, 8, faiss.METRIC_INNER_PRODUCT)
            Xn = np.hstack([self.full / self.norm[:, None], np.zeros((len(self.full), self.pad), np.float32)])
            self.pq.train(Xn if fit_rows is None else Xn[fit_rows])
            # ADC scores q·decode(code) are exactly what faiss computes; decoding
            # once keeps evaluation simple (timing uses the native faiss search)
            self.Xhat = self.pq.sa_decode(self.pq.sa_encode(Xn))
            self.pq.add(Xn)
            self.bytes = pq_m
            return
        if basis is None:
            Xf = self.full if fit_rows is None else self.full[fit_rows]
            _, _, vt = np.linalg.svd(Xf, full_matrices=False)
            basis = vt[:K].T
        self.W = basis.astype(np.float32)
        Z = self.full @ self.W
        if code == "ternary":
            # per-program threshold: tau_q-quantile of |z| (so ~1-tau_q of bits are non-zero)
            self.tau = np.quantile(np.abs(Z), tau_q, axis=0)
            self.T = (np.sign(Z) * (np.abs(Z) > self.tau)).astype(np.int8)
            self.bits = _pack_ternary(self.T)
            self.Tnorm = np.sqrt((self.T.astype(np.float32) ** 2).sum(1)) + 1e-6
            self.bytes = K // 4                                # 2 bits per program
        elif code == "int8":
            self.scale = np.abs(Z).max(0) / 127 + 1e-12
            self.Q = np.clip(np.round(Z / self.scale), -127, 127).astype(np.int8)
            self.Zhat = self.Q.astype(np.float32) * self.scale
            self.Zhat_norm = np.linalg.norm(self.Zhat, axis=1) + 1e-12
            self.bytes = K
        else:
            raise ValueError(code)

    # ---- stage 1 ----------------------------------------------------------
    def coarse_scores(self, q, symmetric=False):
        """Higher = better. Uses codes only (+ the K-dim projection of q)."""
        if self.code == "pq":
            qn = np.concatenate([q / (np.linalg.norm(q) + 1e-12), np.zeros(self.pad)]).astype(np.float32)
            return self.Xhat @ qn
        z = q @ self.W
        if self.code == "ternary":
            if symmetric:
                tq = (np.sign(z) * (np.abs(z) > self.tau)).astype(np.int8)
                qb = _pack_ternary(tq[None])[0]
                return -np.bitwise_count(self.bits ^ qb).sum(1).astype(np.float32)
            return (self.T @ z) / self.Tnorm                   # asymmetric
        return (self.Zhat @ z) / self.Zhat_norm                # int8 ADC (cosine)

    def shortlist(self, q, M, exclude=(), symmetric=False):
        s = self.coarse_scores(q, symmetric)
        if len(exclude):
            s = s.copy()
            s[list(exclude)] = -np.inf
        M = min(M, len(s))
        idx = np.argpartition(-s, M - 1)[:M]
        return idx[np.argsort(-s[idx])]

    # ---- stage 2 ----------------------------------------------------------
    def search(self, q, k=10, M=100, symmetric=False):
        """k-NN by cosine on full vectors, among a stage-1 shortlist of M."""
        cand = self.shortlist(q, M, symmetric=symmetric)
        cos = (self.full[cand] @ q) / (self.norm[cand] * (np.linalg.norm(q) + 1e-12))
        o = np.argsort(-cos)[:k]
        return cand[o], cos[o]

    def code_only(self, q, k=10, symmetric=False):
        return self.shortlist(q, k, symmetric=symmetric)

    def pursuit(self, delta, M=16, refresh=True, symmetric=False, **kw):
        """Residual pursuit whose candidates come from stage 1.
        refresh=True: new shortlist from the residual at every step."""
        static = None if refresh else self.shortlist(delta, M, symmetric=symmetric)

        def cand(r, sel):
            if static is not None:
                return [int(c) for c in static]
            return [int(c) for c in self.shortlist(r, M + len(sel), exclude=sel, symmetric=symmetric)]
        return DL.residual_pursuit(delta, self.full, None, candidates=cand, **kw)


def exact_knn(full, q, k=10):
    n = np.linalg.norm(full, axis=1) + 1e-12
    cos = (full @ q) / (n * (np.linalg.norm(q) + 1e-12))
    idx = np.argpartition(-cos, k - 1)[:k]
    return idx[np.argsort(-cos[idx])]


def timeit(fn, reps=20):
    fn()
    t = time.perf_counter()
    for _ in range(reps):
        fn()
    return (time.perf_counter() - t) / reps
