"""
CGMacros (PhysioNet, open access, CC BY-NC-SA 4.0)
  Das A, Kerr D, Glantz N, Bevier W, Santiago R, Gutierrez-Osuna R,
  Mortazavi BJ. "CGMacros: a pilot scientific dataset for personalized
  nutrition and diet monitoring." Sci Data 12:1557 (2025).
  doi:10.1038/s41597-025-05851-7 ; https://physionet.org/content/cgmacros/1.0.0/

45 adults (15 healthy, 16 prediabetes, 14 T2D), ~10 days of 1-min CGM (Libre +
Dexcom) with logged meals = ~30 natural, macronutrient-annotated
PERTURBATIONS per person -> a perturbation-relaxation dataset.

The archive is a single 627 MB zip (mostly food photos). We read only the
central directory and the CSV members via HTTP Range requests (~10 MB total),
caching them under figures/_data_cache/cgmacros/.
"""
import io
import os
import struct
import subprocess
import zipfile
import zlib

from netjac import DATA_CACHE

URL = ("https://physionet.org/files/cgmacros/1.0.0/"
       "CGMacros_dateshifted365.zip")
DIR = os.path.join(DATA_CACHE, "cgmacros")


def _range(start, end, tries=5):
    import time
    for k in range(tries):
        r = subprocess.run(["curl", "-sf", "--max-time", "120", "-r",
                            f"{start}-{end}", URL], capture_output=True)
        if r.returncode == 0 and len(r.stdout) == end - start + 1:
            return r.stdout
        time.sleep(2 ** (k + 1))
    raise RuntimeError(f"range {start}-{end} failed (curl exit {r.returncode})")


def _size():
    out = subprocess.run(["curl", "-sfI", "--max-time", "60", URL],
                         check=True, capture_output=True, text=True).stdout
    for line in out.splitlines():
        if line.lower().startswith("content-length"):
            return int(line.split(":")[1])
    raise RuntimeError("no content-length")


class _RemoteZip(io.RawIOBase):
    """Seekable file-like over HTTP ranges, enough for zipfile to parse the
    central directory (read-ahead blocks of 1 MB)."""

    def __init__(self):
        self.size = _size()
        self.pos = 0
        self.cache = {}

    def seekable(self):
        return True

    def readable(self):
        return True

    def seek(self, off, whence=0):
        self.pos = {0: off, 1: self.pos + off, 2: self.size + off}[whence]
        return self.pos

    def tell(self):
        return self.pos

    def readinto(self, b):
        n = min(len(b), self.size - self.pos)
        if n <= 0:
            return 0
        blk = 1 << 20
        out = bytearray()
        p = self.pos
        while len(out) < n:
            k = p // blk
            if k not in self.cache:
                self.cache[k] = _range(k * blk, min((k + 1) * blk,
                                                    self.size) - 1)
            chunk = self.cache[k][p - k * blk:]
            take = min(len(chunk), n - len(out))
            out += chunk[:take]
            p += take
        b[:n] = out
        self.pos += n
        return n


def _extract_member(info):
    # local header is 30 bytes + name + extra; read header to get lengths
    h = _range(info.header_offset, info.header_offset + 29)
    nlen, xlen = struct.unpack("<HH", h[26:30])
    start = info.header_offset + 30 + nlen + xlen
    raw = _range(start, start + info.compress_size - 1)
    if info.compress_type == zipfile.ZIP_DEFLATED:
        return zlib.decompress(raw, -15)
    if info.compress_type == zipfile.ZIP_STORED:
        return raw
    raise RuntimeError(f"unsupported compression {info.compress_type}")


def fetch(force=False):
    """Download bio.csv + every CGMacros-0XX.csv. Returns the cache dir."""
    os.makedirs(DIR, exist_ok=True)
    marker = os.path.join(DIR, ".complete")
    if os.path.exists(marker) and not force:
        return DIR
    z = zipfile.ZipFile(_RemoteZip())
    want = [i for i in z.infolist() if i.filename.lower().endswith(".csv")]
    for info in want:
        dest = os.path.join(DIR, os.path.basename(info.filename))
        if os.path.exists(dest) and os.path.getsize(dest) > 0:
            continue
        with open(dest, "wb") as fh:
            fh.write(_extract_member(info))
    open(marker, "w").write("\n".join(i.filename for i in want))
    return DIR


def load():
    """Return (bio DataFrame, {participant_id: CGM DataFrame})."""
    import glob
    import pandas as pd
    d = fetch()
    bio = pd.read_csv(os.path.join(d, "bio.csv"))
    traces = {}
    for p in sorted(glob.glob(os.path.join(d, "CGMacros-*.csv"))):
        pid = int(os.path.basename(p).split("-")[1].split(".")[0])
        traces[pid] = pd.read_csv(p)
    return bio, traces
