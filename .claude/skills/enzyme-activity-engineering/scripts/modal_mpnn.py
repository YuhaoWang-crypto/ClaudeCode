import modal, json

app = modal.App("is621-mpnn")
img = (
    modal.Image.debian_slim(python_version="3.11")
    .apt_install("git")
    .pip_install("torch", "numpy")
    .run_commands("git clone --depth 1 https://github.com/dauparas/ProteinMPNN.git /ProteinMPNN")
)
MPNN_ALPHABET = "ACDEFGHIKLMNPQRSTVWYX"


@app.function(image=img, gpu="T4", timeout=1200)
def score(pdb_text: str, weights: str = "v_48_020.pt"):
    import sys, torch, numpy as np
    sys.path.insert(0, "/ProteinMPNN")
    from protein_mpnn_utils import ProteinMPNN, parse_PDB, tied_featurize, StructureDatasetPDB
    dev = torch.device("cuda")
    open("/tmp/a.pdb", "w").write(pdb_text)
    pdb = parse_PDB("/tmp/a.pdb")
    ds = StructureDatasetPDB(pdb, truncate=None, max_length=20000)
    ckpt = torch.load(f"/ProteinMPNN/vanilla_model_weights/{weights}", map_location=dev)
    model = ProteinMPNN(num_letters=21, node_features=128, edge_features=128,
                        hidden_dim=128, num_encoder_layers=3, num_decoder_layers=3,
                        k_neighbors=ckpt["num_edges"])
    model.to(dev); model.load_state_dict(ckpt["model_state_dict"]); model.eval()
    chains = [it[-1:] for it in list(pdb[0]) if it[:9] == "seq_chain"]
    cid = {pdb[0]["name"]: (chains, [])}
    batch = [ds[0]]
    feats = tied_featurize(batch, dev, cid, None, None, None, None, None, ca_only=False)
    X, S, mask, lengths, chain_M, chain_enc = feats[0], feats[1], feats[2], feats[3], feats[4], feats[5]
    residue_idx = feats[12]  # order: ...,chain_M_pos(10),omit_AA_mask(11),residue_idx(12),...
    with torch.no_grad():
        logp = model.unconditional_probs(X, mask, residue_idx, chain_enc)  # [B,L,21]
    logp = logp[0].cpu().numpy()
    S = S[0].cpu().numpy()
    seq = "".join(MPNN_ALPHABET[i] for i in S)
    return {"weights": weights, "seq": seq, "log_p": logp.tolist(), "alphabet": MPNN_ALPHABET}


@app.local_entrypoint()
def main():
    out = score.remote(open("data/8WT9_chainA.pdb").read())
    with open("data/mpnn_scores.json", "w") as f:
        json.dump(out, f)
    print("WROTE data/mpnn_scores.json  seqlen", len(out["seq"]))
    print("seq head:", out["seq"][:50])
