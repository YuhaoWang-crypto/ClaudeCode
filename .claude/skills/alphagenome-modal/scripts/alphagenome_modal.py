"""
AlphaGenome (genomicsxai PyTorch port) variant-effect scanning on Modal.

Model code : github.com/genomicsxai/alphagenome-pytorch  (PyTorch port of DeepMind
             AlphaGenome, 450M params; API: AlphaGenome.from_pretrained / .predict)
Weights    : Modal Volume "alphagenome-weights" ->
             gtca__alphagenome_pytorch/model_all_folds.safetensors (+ config.json)
             (gtca/alphagenome_pytorch on HF; converted from official JAX checkpoint;
              includes per-head track_means -> proper experimental-space unscaling)
App        : "alphagenome-pytorch"   GPU: A100-80GB, model loaded ONCE per container.
Window     : 131,072 bp, variant centered. Human organism_index = 0. One-hot (B,S,4).
hg38 seq   : UCSC REST (api.genome.ucsc.edu), block-fetched + cached per container.
Tissue map : port ships data/track_metadata_human.parquet -> TrackMetadataCatalog,
             so raw track index -> assay/biosample/ontology IS resolvable here.

FIDELITY CAVEAT: community PyTorch port of AlphaGenome (weights converted from the
official JAX checkpoint), NOT DeepMind's official served model. Deltas are
directional/relative signals; treat magnitudes as relative, not calibrated official
AlphaGenome variant scores.
"""

import modal

WINDOW = 131_072
WEIGHTS_DIR = "/weights/gtca__alphagenome_pytorch"
SAFETENSORS = f"{WEIGHTS_DIR}/model_all_folds.safetensors"

HEADS = ("rna_seq", "atac", "dnase", "cage", "splice_sites", "splice_site_usage")
GENOME_HEADS = ("rna_seq", "atac", "dnase", "cage")
GENOME_HALF = 512
SPLICE_HALF = 64

image = (
    modal.Image.debian_slim(python_version="3.12")
    .apt_install("git")
    .pip_install(
        "torch==2.5.1",
        "numpy<2",
        "pandas>=2.0",
        "pyarrow>=15",
        "safetensors>=0.4.3",
        "requests",
    )
    .pip_install(
        "alphagenome-pytorch @ git+https://github.com/genomicsxai/alphagenome-pytorch.git",
        extra_options="--no-deps",
    )
)

app = modal.App("alphagenome-pytorch")
weights_vol = modal.Volume.from_name("alphagenome-weights")


@app.cls(image=image, volumes={"/weights": weights_vol}, gpu="A100-80GB",
         timeout=3000, scaledown_window=60)
class AlphaGenomeModel:
    @modal.enter()
    def load(self):
        import warnings, torch
        from alphagenome_pytorch import AlphaGenome

        self.torch = torch
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            model = AlphaGenome.from_pretrained(SAFETENSORS, device="cuda")
            self.load_warnings = [str(x.message) for x in w]
        model.eval()
        self.model = model

        # builtin track metadata catalog (human) for tissue mapping
        self.catalog = None
        self.catalog_err = None
        try:
            from alphagenome_pytorch.named_outputs import TrackMetadataCatalog
            self.catalog = TrackMetadataCatalog.load_builtin(0)
            self.catalog_outputs = list(self.catalog.outputs(0))
        except Exception as e:
            self.catalog_err = repr(e)
            self.catalog_outputs = []

        self._cache = {}  # per-chrom list of (start0, end0, seq)
        self.load_report = {
            "model": "genomicsxai/alphagenome-pytorch",
            "weights": "gtca model_all_folds.safetensors",
            "load_warnings": self.load_warnings,
            "catalog_loaded": self.catalog is not None,
            "catalog_outputs": self.catalog_outputs,
            "catalog_err": self.catalog_err,
        }

    # ---- hg38 sequence, block-fetched from UCSC, cached ------------------
    def _fetch_ucsc(self, chrom, start0, end0):
        import requests
        url = ("https://api.genome.ucsc.edu/getData/sequence"
               f"?genome=hg38;chrom={chrom};start={start0};end={end0}")
        last = None
        for _ in range(4):
            try:
                r = requests.get(url, timeout=90)
                j = r.json()
                if r.status_code == 200 and "dna" in j:
                    return j["dna"].upper()
                last = j
            except Exception as e:
                last = repr(e)
        raise RuntimeError(f"UCSC fetch failed {chrom}:{start0}-{end0}: {last}")

    def get_window(self, chrom, pos, W=WINDOW):
        var0 = pos - 1
        need_s = var0 - W // 2
        need_e = need_s + W
        blocks = self._cache.setdefault(chrom, [])
        for (bs, be, seq) in blocks:
            if bs <= need_s and need_e <= be:
                return seq[need_s - bs: need_e - bs], W // 2, "cache"
        bs = max(0, var0 - (3 * W) // 2)
        be = bs + 3 * W
        seq = self._fetch_ucsc(chrom, bs, be)
        blocks.append((bs, be, seq))
        return seq[need_s - bs: need_e - bs], W // 2, "fetch"

    @staticmethod
    def _onehot(seq):
        import numpy as np
        m = {"A": 0, "C": 1, "G": 2, "T": 3}
        idx = np.fromiter((m.get(b, -1) for b in seq), dtype=np.int64, count=len(seq))
        oh = np.zeros((len(seq), 4), dtype=np.float32)
        valid = idx >= 0
        oh[np.arange(len(seq))[valid], idx[valid]] = 1.0
        return oh

    def _predict(self, onehot_np):
        torch = self.torch
        x = torch.from_numpy(onehot_np).unsqueeze(0).to("cuda")
        with torch.no_grad():
            out = self.model.predict(x, 0, heads=HEADS, channels_last=True)
        return out

    def _local_1bp(self, out, head, idx, half):
        """Return CPU float tensor (P, T) local window at 1bp for `head`."""
        v = out[head]
        if head == "splice_sites":
            t = v["probs"] if isinstance(v, dict) else v      # (B,S,5)
        elif head == "splice_site_usage":
            t = v["predictions"] if isinstance(v, dict) else v  # (B,S,T) sigmoid
        else:  # genome tracks: dict{res:(B,S,T)}
            t = v[1] if isinstance(v, dict) else v
        t = t[0]                                              # (S, T)
        S = t.shape[0]
        c = idx if idx < S else S // 2
        lo, hi = max(0, c - half), min(S, c + half + 1)
        return t[lo:hi].float().cpu(), (c - lo)

    def _tissue(self, head, tidx):
        if self.catalog is None:
            return None
        # try model head name and a couple of aliases against catalog outputs
        aliases = {
            "splice_sites": ["splice_sites", "splice_site", "splice_sites_classification"],
            "splice_site_usage": ["splice_site_usage", "splice_sites_usage", "splice_usage"],
        }.get(head, [head])
        for name in aliases:
            try:
                tracks = self.catalog.get_tracks(name, organism=0)
            except Exception:
                tracks = ()
            if tracks and tidx < len(tracks):
                tk = tracks[tidx]
                fields = {}
                for f in ("track_name", "biosample_name", "ontology_curie",
                          "assay", "biosample_type", "strand"):
                    try:
                        val = tk.get(f) if hasattr(tk, "get") else None
                    except Exception:
                        val = None
                    if val not in (None, "", "nan"):
                        fields[f] = val
                if fields:
                    return {"output": name, "n_tracks": len(tracks), **fields}
        return {"output": None, "note": "no metadata for this head (likely splice-class or padding track)"}

    def score_one(self, chrom, pos, ref, alt):
        torch = self.torch
        win, idx, src = self.get_window(chrom, pos)
        obs_ref = win[idx]
        ref_ok = (obs_ref == ref.upper())
        alt_win = win[:idx] + alt.upper() + win[idx + 1:]

        ref_out = self._predict(self._onehot(win))
        # slice ref locals, then free full ref output
        ref_loc = {}
        for h in HEADS:
            half = SPLICE_HALF if h.startswith("splice") else GENOME_HALF
            ref_loc[h] = self._local_1bp(ref_out, h, idx, half)
        del ref_out
        torch.cuda.empty_cache()

        alt_out = self._predict(self._onehot(alt_win))
        res = {}
        for h in HEADS:
            half = SPLICE_HALF if h.startswith("splice") else GENOME_HALF
            a, cc = self._local_1bp(alt_out, h, idx, half)
            r, _ = ref_loc[h]
            d = (a - r).abs()                     # (P, T)
            per_track_max = d.max(dim=0).values   # (T,)
            top = int(per_track_max.argmax().item())
            res[h] = {
                "mean_abs": round(float(d.mean().item()), 6),
                "at_var": round(float(d[cc].mean().item()), 6),
                "max_abs": round(float(d.max().item()), 6),
                "top_track_idx": top,
                "top_track_delta": round(float(per_track_max[top].item()), 6),
                "tissue": self._tissue(h, top),
            }
        del alt_out
        torch.cuda.empty_cache()
        return {"obs_ref": obs_ref, "ref_match": ref_ok,
                "seq_source": src, "window_var_index": idx, "heads": res}

    @modal.method()
    def batch_variant_effect(self, variants: list):
        """Load model ONCE (in .load); loop variants; reuse cached hg38 blocks."""
        import time
        results = []
        t0 = time.time()
        for v in variants:
            meta = {k: v.get(k) for k in
                    ("gene", "rsid", "chrom", "pos", "ref", "alt", "klass", "status")}
            try:
                out = self.score_one(v["chrom"], int(v["pos"]), v["ref"], v["alt"])
                out.update(meta)
            except Exception as e:
                out = {**meta, "error": repr(e)}
            results.append(out)
        self.load_report["elapsed_sec"] = round(time.time() - t0, 1)
        return {"load_report": self.load_report, "results": results}

    @modal.method()
    def variant_effect_at(self, chrom, pos, ref, alt):
        return self.score_one(chrom, int(pos), ref, alt)


# ---------------------------------------------------------------------------
# Curated hg38 variant table (INHBE / ACVR1C / GPR75) + driver
# ---------------------------------------------------------------------------
VARIANTS = [
    # gene, rsid, chrom, 1-based hg38 pos, ref, alt, class, status
    # ---- INHBE  chr12 (+strand), tx ENST00000266646, 2 exons ----
    ["INHBE", "rs375342858",  "chr12", 57455835, "G", "T", "splice_donor",        "documented|ANCHOR"],
    ["INHBE", "rs1420168629", "chr12", 57455837, "G", "A", "splice_donor_region", "documented"],
    ["INHBE", "rs1870821812", "chr12", 57456092, "A", "C", "splice_acceptor",     "documented"],
    ["INHBE", "rs150777893",  "chr12", 57456093, "G", "C", "splice_acceptor",     "documented"],
    ["INHBE", "rs866560678",  "chr12", 57455339, "C", "T", "5_prime_UTR",         "documented"],
    ["INHBE", "rs1408776511", "chr12", 57455322, "G", "A", "5_prime_UTR",         "documented"],
    # ---- ACVR1C  chr2 (-strand), tx ENST00000243349, 9 exons ----
    ["ACVR1C", "rs55920843",  "chr2", 157556189, "T", "C", "missense_CODING",     "documented|ANCHOR"],
    ["ACVR1C", "rs887404649", "chr2", 157587421, "G", "A", "splice_region",       "documented"],
    ["ACVR1C", "rs1205245099","chr2", 157587429, "A", "G", "splice_polypyrimidine","documented"],
    ["ACVR1C", "rs1688099936","chr2", 157556336, "C", "T", "splice_region",       "documented"],
    ["ACVR1C", "rs934393720", "chr2", 157628868, "C", "T", "5_prime_UTR",         "documented"],
    ["ACVR1C", "rs895955903", "chr2", 157628877, "G", "A", "5_prime_UTR",         "documented"],
    # ---- GPR75  chr2 (-strand), tx ENST00000394705, 2 exons ----
    ["GPR75", "PROMO_DEMO",   "chr2", 53860050, "C", "G", "promoter_upstream",    "demonstrative|ANCHOR"],
    ["GPR75", "rs1244093517", "chr2", 53859823, "C", "T", "splice_donor_5th_base","documented"],
    ["GPR75", "rs1678332197", "chr2", 53859822, "T", "C", "splice_donor_region",  "documented"],
    ["GPR75", "rs1321192968", "chr2", 53854865, "C", "A", "splice_region",        "documented"],
    ["GPR75", "rs753524658",  "chr2", 53859835, "G", "A", "5_prime_UTR",          "documented"],
    ["GPR75", "rs1041290706", "chr2", 53859834, "G", "T", "5_prime_UTR",          "documented"],
]


def _as_dicts():
    keys = ["gene", "rsid", "chrom", "pos", "ref", "alt", "klass", "status"]
    return [dict(zip(keys, v)) for v in VARIANTS]


@app.local_entrypoint()
def main():
    import json
    model = AlphaGenomeModel()
    out = model.batch_variant_effect.remote(_as_dicts())
    print("LOAD_REPORT=" + json.dumps(out["load_report"]))
    print("RESULTS=" + json.dumps(out["results"]))
