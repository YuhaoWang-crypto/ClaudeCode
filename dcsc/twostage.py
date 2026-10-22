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


class HybridIndex:
    """PQ coarse pass -> program-code rescoring. No full vectors at query time.

    Per item in RAM: PQ code (pq_m bytes) + int8 program code (K bytes).

        coarse_space  "raw"      PQ over the z-scored gene vector (default)
                      "program"  PQ over the K float program coordinates
        coarse        "pq"       exhaustive ADC scan over all PQ codes
                      "ivfpq"    inverted file: only nprobe of nlist cells
                                 are scanned -> sublinear in library size
        rescore       "int8"     cosine on the int8 program code (default)
                      "float"    cosine on float program coordinates
                      "none"     keep the coarse order (PQ only)

    `full` is used only at build time (basis, codes); `add()` streams more
    items in chunks without keeping them.
    """

    def __init__(self, full, K=32, pq_m=32, coarse="pq", coarse_space="raw", nlist=1024,
                 fit_rows=None, basis=None, add_init=True):
        import faiss
        full = np.asarray(full, dtype=np.float32)
        fit = full if fit_rows is None else full[fit_rows]
        if basis is None:
            _, _, vt = np.linalg.svd(fit, full_matrices=False)
            basis = vt[:K].T
        self.W = basis.astype(np.float32)
        self.K, self.pq_m, self.coarse_space = K, pq_m, coarse_space
        Zfit = fit @ self.W
        self.scale = np.abs(Zfit).max(0) / 127 + 1e-12
        cdata = self._coarse_vec(fit)
        d = cdata.shape[1]
        if coarse == "pq":
            self.index = faiss.IndexPQ(d, pq_m, 8, faiss.METRIC_INNER_PRODUCT)
        else:
            q = faiss.IndexFlatIP(d)
            self.index = faiss.IndexIVFPQ(q, d, nlist, pq_m, 8, faiss.METRIC_INNER_PRODUCT)
            self._quantizer = q          # keep a reference (faiss does not own it)
        self.index.train(cdata)
        self.Q = np.zeros((0, K), np.int8)
        self.Zf = np.zeros((0, K), np.float32)   # float coords, only for rescore="float"
        if add_init:          # False: `full` is only a training sample; stream items with add()
            self.add(full)

    def _coarse_vec(self, X):
        V = X @ self.W if self.coarse_space == "program" else X
        V = V / (np.linalg.norm(V, axis=1, keepdims=True) + 1e-12)
        pad = (-V.shape[1]) % self.pq_m
        return np.ascontiguousarray(np.hstack([V, np.zeros((len(V), pad), np.float32)]), dtype=np.float32)

    def add(self, X, chunk=50_000):
        """Stream items in; codes are appended in chunks (linear-time build)."""
        if not hasattr(self, "_chunks"):
            self._chunks = []
        for s in range(0, len(X), chunk):
            x = np.asarray(X[s:s + chunk], dtype=np.float32)
            z = x @ self.W
            self._chunks.append((np.clip(np.round(z / self.scale), -127, 127).astype(np.int8), z))
            self.index.add(self._coarse_vec(x))
        self._dirty = True

    def _finalize(self):
        if not getattr(self, "_dirty", False):
            return
        if self._chunks:
            self.Q = np.vstack([self.Q] + [c[0] for c in self._chunks])
            self.Zf = np.vstack([self.Zf] + [c[1] for c in self._chunks])
            self._chunks = []
        self.Zhat = self.Q.astype(np.float32) * self.scale
        self.Zhat /= np.linalg.norm(self.Zhat, axis=1, keepdims=True) + 1e-12
        self.Zf_n = self.Zf / (np.linalg.norm(self.Zf, axis=1, keepdims=True) + 1e-12)
        self._dirty = False

    def __len__(self):
        self._finalize()
        return len(self.Q)

    @property
    def bytes_per_item(self):
        return self.pq_m + self.K

    def search(self, q, k=10, M=1000, nprobe=16, rescore="int8"):
        if hasattr(self.index, "nprobe"):
            self.index.nprobe = nprobe
        self._finalize()
        _, I = self.index.search(self._coarse_vec(q[None]), M)
        cand = I[0][I[0] >= 0]
        if rescore == "none":
            return cand[:k]
        z = q @ self.W
        z /= np.linalg.norm(z) + 1e-12
        S = self.Zhat if rescore == "int8" else self.Zf_n
        s = S[cand] @ z
        return cand[np.argsort(-s)[:k]]

    def program_scan(self, q, k=10):
        """Reference: exhaustive cosine over all int8 program codes (O(N))."""
        self._finalize()
        z = q @ self.W
        s = self.Zhat @ (z / (np.linalg.norm(z) + 1e-12))
        idx = np.argpartition(-s, k - 1)[:k]
        return idx[np.argsort(-s[idx])]
