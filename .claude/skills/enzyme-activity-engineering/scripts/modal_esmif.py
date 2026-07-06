import modal, json

app = modal.App("is621-esmif")
# Pin torch 2.2.1+cu121 so PyG wheels resolve cleanly
img = (
    modal.Image.debian_slim(python_version="3.11")
    .pip_install("torch==2.2.1", "torchvision==0.17.1",
                 index_url="https://download.pytorch.org/whl/cu121")
    .pip_install("torch-scatter", "torch-sparse", "torch-cluster", "torch-geometric",
                 find_links="https://data.pyg.org/whl/torch-2.2.1+cu121.html")
    .pip_install("fair-esm==2.0.0", "biotite==0.39.0", "numpy<2")
    .env({"TORCH_HOME": "/models/torch"})
)
vol = modal.Volume.from_name("is621-models", create_if_missing=True)
AAS = "ACDEFGHIKLMNPQRSTVWY"


@app.function(image=img, gpu="A10G", timeout=2400, volumes={"/models": vol})
def score(pdb_text: str):
    import torch, esm, numpy as np
    import esm.inverse_folding as esmif
    open("/tmp/a.pdb", "w").write(pdb_text)
    model, alphabet = esm.pretrained.esm_if1_gvp4_t16_142M_UR50()
    model = model.eval().cuda()
    coords, seq = esmif.util.load_coords("/tmp/a.pdb", "A")  # N,CA,C coords + modeled seq
    device = next(model.parameters()).device
    batch_converter = esmif.util.CoordBatchConverter(alphabet)
    batch = [(coords, None, seq)]
    coords_t, confidence, strs, tokens, padding_mask = batch_converter(batch, device=device)
    prev = tokens[:, :-1]
    with torch.no_grad():
        logits, _ = model.forward(coords_t, padding_mask, confidence, prev)  # [B, vocab, L]
    logp = torch.log_softmax(logits, dim=1)[0].cpu().numpy()  # [vocab, L]
    tok = {aa: alphabet.get_idx(aa) for aa in AAS}
    # scores aligned to modeled seq (which may start at res 4 of native)
    out = []
    for i, wt in enumerate(seq):
        col = logp[:, i]
        wt_ll = float(col[tok[wt]]) if wt in tok else None
        muts = {m: float(col[tok[m]] - col[tok[wt]]) for m in AAS if m != wt and wt in tok}
        out.append({"idx": i, "wt": wt, "wt_ll": wt_ll, "muts": muts})
    vol.commit()
    return {"model": "esm_if1_gvp4_t16_142M_UR50", "modeled_seq": seq, "per_pos": out}


@app.local_entrypoint()
def main():
    out = score.remote(open("data/8WT9_chainA.pdb").read())
    with open("data/esmif_scores.json", "w") as f:
        json.dump(out, f)
    print("WROTE data/esmif_scores.json  modeled_len", len(out["modeled_seq"]))
    print("modeled seq head:", out["modeled_seq"][:50])
