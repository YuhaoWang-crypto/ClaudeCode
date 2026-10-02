"""
Locate the BootLoops toolkit (github.com/BootLoops-ai/bootloops) and expose
its BALLER ball-arithmetic package (fail-closed printing + the dual-path
halt-not-average tripwire).

Search order: $BOOTLOOPS_DIR, then ../bootloops beside this repository (the
sibling layout BootLoops' own INSTALL.md assumes). Install with
`bash scripts/setup_bootloops.sh`. Returns None when it is absent, so callers
can skip by name instead of crashing.
"""
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
_cached = None


def bootloops_dir():
    cands = [os.environ.get("BOOTLOOPS_DIR"),
             os.path.join(os.path.dirname(_REPO), "bootloops")]
    for c in cands:
        if c and os.path.isdir(os.path.join(c, "tools", "baller")):
            return os.path.abspath(c)
    return None


def baller():
    """Import BALLER (and run its sha-pin verify, fail-closed). None if absent."""
    global _cached
    if _cached is not None:
        return _cached
    root = bootloops_dir()
    if root is None:
        return None
    p = os.path.join(root, "tools", "baller")
    if p not in sys.path:
        sys.path.insert(0, p)
    import baller as _b
    v = _b.verify()
    if v.get("vendor_bad"):
        raise RuntimeError(f"BootLoops baller pin check failed: {v['vendor_bad']}")
    _cached = _b
    return _b
