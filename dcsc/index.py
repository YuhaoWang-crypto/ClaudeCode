"""
Search over 64-byte DCSC codes.

Three engines, same answers:
  * numpy  — np.bitwise_count on the four 64-bit MOD words (portable);
  * C      — a 10-line POPCNT kernel compiled at runtime (what "hundreds of
             millions of comparisons per second per core" actually needs);
  * faiss  — IndexBinaryFlat / IndexBinaryHNSW on the 256-bit MOD block.

Query modes (brain-hole 1):
  "mod"      match programs only -> cross-tissue functional prototypes
  "same_id"  fix the identity path, then match programs -> same cell type,
             different state
  "full"     identity mismatch penalty + programs + Mem/Env L1
"""
import ctypes
import os
import subprocess
import tempfile

import numpy as np

from . import codec as C

MOD_WORDS = slice(2, 6)     # bytes 16..47


def mod_distance(codes, q, mask=None):
    """Ternary Hamming: 0 vs ±1 costs 1, +1 vs -1 costs 2."""
    x = codes.view(np.uint64)[:, MOD_WORDS] ^ q.view(np.uint64)[MOD_WORDS]
    if mask is not None:
        x = x & mask.view(np.uint64)
    return np.bitwise_count(x).sum(1, dtype=np.int32)


def _i8(codes, blk):
    a, b = C.BLOCKS[blk]
    return codes[:, a:b].view(np.int8).astype(np.int16)


def distances(codes, q, mode="mod", mask=None, depth=2, w_slow=1 / 40):
    d = mod_distance(codes, q, mask)
    if mode == "same_id":
        d = np.where(np.all(codes[:, :depth] == q[:depth], 1), d, 10_000)
    elif mode == "full":
        id_pen = 4 * (codes[:, 0] != q[0]) + 2 * (codes[:, 1] != q[1])
        slow = (np.abs(_i8(codes, "MEM") - _i8(q[None], "MEM")).sum(1)
                + np.abs(_i8(codes, "ENV") - _i8(q[None], "ENV")).sum(1))
        d = d + id_pen + w_slow * slow
    return d


def search(codes, q, k=10, **kw):
    d = distances(codes, q, **kw)
    k = min(k, len(d))
    idx = np.argpartition(d, k - 1)[:k]
    idx = idx[np.argsort(d[idx], kind="stable")]
    return idx, d[idx]


# ---------------------------------------------------------------------------
_C_SRC = r"""
#include <stdint.h>
void hamming_mod(const uint64_t* codes, long n, const uint64_t* q,
                 const uint64_t* m, uint8_t* out) {
    uint64_t q0=q[2], q1=q[3], q2=q[4], q3=q[5];
    for (long i = 0; i < n; i++) {
        const uint64_t* c = codes + 8*i;
        out[i] = __builtin_popcountll((c[2]^q0)&m[0]) + __builtin_popcountll((c[3]^q1)&m[1])
               + __builtin_popcountll((c[4]^q2)&m[2]) + __builtin_popcountll((c[5]^q3)&m[3]);
    }
}
"""
_lib = None


def c_kernel():
    """Compile the POPCNT kernel once; returns None if no compiler."""
    global _lib
    if _lib is not None:
        return _lib
    try:
        d = tempfile.mkdtemp(prefix="dcsc_")
        src, so = os.path.join(d, "k.c"), os.path.join(d, "k.so")
        with open(src, "w") as f:
            f.write(_C_SRC)
        subprocess.run(["gcc", "-O3", "-march=native", "-shared", "-fPIC", src, "-o", so],
                       check=True, capture_output=True)
        lib = ctypes.CDLL(so)
        lib.hamming_mod.restype = None
        lib.hamming_mod.argtypes = [ctypes.c_void_p, ctypes.c_long, ctypes.c_void_p,
                                    ctypes.c_void_p, ctypes.c_void_p]
        _lib = lib
    except Exception:
        _lib = False
    return _lib


def mod_distance_c(codes, q, mask=None):
    lib = c_kernel()
    if not lib:
        return mod_distance(codes, q, mask).astype(np.uint8)
    codes = np.ascontiguousarray(codes)
    m = np.full(4, np.uint64(0xFFFFFFFFFFFFFFFF)) if mask is None else np.ascontiguousarray(mask).view(np.uint64)
    qq = np.ascontiguousarray(q).view(np.uint64)
    out = np.empty(len(codes), np.uint8)
    lib.hamming_mod(codes.ctypes.data, ctypes.c_long(len(codes)), qq.ctypes.data,
                    m.ctypes.data, out.ctypes.data)
    return out


# ---------------------------------------------------------------------------
def faiss_mod_index(codes, kind="flat", M=32):
    import faiss
    a, b = C.BLOCKS["MOD"]
    nbits = (b - a) * 8
    idx = faiss.IndexBinaryFlat(nbits) if kind == "flat" else faiss.IndexBinaryHNSW(nbits, M)
    idx.add(np.ascontiguousarray(codes[:, a:b]))
    return idx


def mod_block(codes):
    a, b = C.BLOCKS["MOD"]
    return np.ascontiguousarray(codes[:, a:b])
