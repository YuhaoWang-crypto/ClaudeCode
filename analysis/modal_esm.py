import modal, json

app = modal.App("is621-esm")
vol = modal.Volume.from_name("is621-models", create_if_missing=True)
img = (
    modal.Image.debian_slim(python_version="3.11")
    .pip_install("torch", "fair-esm==2.0.0", "numpy")
    .env({"TORCH_HOME": "/models/torch"})
)

SEQ = ("MEHELHYIGIDTAKEKLDVDVLRPDGRHRTKKFANTTKGHDELVSWLKGHKIDHAHICIEATGTYMEP"
       "VAECLYDAGYIVSVINPALGKAFAQSEGLRNKTDTVDARMLAEFCRQKRPAAWEAPHPLERALRALVVRH"
       "QALTDMHTQELNRTETAREVQRPSIDAHLLWLEAELKRLEKQIKDLTDDDPDMKHRRKLLESIPGIGEKT"
       "SAVLLAYIGLKDRFAHARQFAAFAGLTPRRYESGSSVRGASRMSKAGHVSLRRALYMPAMVATSKTEWGR"
       "AFRDRLAANGKKGKVILGAMMRKLAQVAYGVLKSGVPFDASRHNPVAA")
AAS = "ACDEFGHIKLMNPQRSTVWY"


@app.function(image=img, gpu="A10G", timeout=3600, volumes={"/models": vol})
def score(model_name: str):
    import torch, esm, numpy as np
    torch.set_grad_enabled(False)
    loader = getattr(esm.pretrained, model_name)
    model, alphabet = loader()
    model = model.eval().cuda()
    bc = alphabet.get_batch_converter()
    L = len(SEQ)
    mask_idx = alphabet.mask_idx
    aa_tok = {aa: alphabet.get_idx(aa) for aa in AAS}

    # masked-marginal: mask each position, read log-probs there
    scores = {}  # pos(1-based)-> {mut: delta_llr}
    wt_ll = np.zeros(L)
    B = 12
    positions = list(range(L))
    for s in range(0, L, B):
        batch_pos = positions[s:s + B]
        data = []
        for p in batch_pos:
            mseq = list(SEQ)
            mseq[p] = "<mask>"  # placeholder; we edit tokens after tokenizing
            data.append((f"p{p}", SEQ))
        _, _, toks = bc(data)
        toks = toks.cuda()
        for j, p in enumerate(batch_pos):
            toks[j, p + 1] = mask_idx  # +1 for BOS
        logits = model(toks)["logits"]
        logp = torch.log_softmax(logits, dim=-1)
        for j, p in enumerate(batch_pos):
            wt = SEQ[p]
            row = logp[j, p + 1]
            wt_ll[p] = row[aa_tok[wt]].item()
            d = {}
            for mut in AAS:
                if mut == wt:
                    continue
                d[mut] = float(row[aa_tok[mut]].item() - row[aa_tok[wt]].item())
            scores[p + 1] = {"wt": wt, "muts": d}
    vol.commit()
    return {"model": model_name, "length": L, "scores": scores}


@app.local_entrypoint()
def main():
    out = {}
    for mn in ["esm1v_t33_650M_UR90S_1", "esm2_t33_650M_UR50D"]:
        print("scoring", mn, "...")
        out[mn] = score.remote(mn)
        print("  done", mn, "positions:", len(out[mn]["scores"]))
    with open("data/esm_scores.json", "w") as f:
        json.dump(out, f)
    print("WROTE data/esm_scores.json")
