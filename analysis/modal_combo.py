import modal, json, itertools

app = modal.App("is621-combo")
vol = modal.Volume.from_name("is621-models", create_if_missing=True)
img = (modal.Image.debian_slim(python_version="3.11")
       .pip_install("torch", "fair-esm==2.0.0", "numpy")
       .env({"TORCH_HOME": "/models/torch"}))

SEQ = ("MEHELHYIGIDTAKEKLDVDVLRPDGRHRTKKFANTTKGHDELVSWLKGHKIDHAHICIEATGTYMEP"
       "VAECLYDAGYIVSVINPALGKAFAQSEGLRNKTDTVDARMLAEFCRQKRPAAWEAPHPLERALRALVVRH"
       "QALTDMHTQELNRTETAREVQRPSIDAHLLWLEAELKRLEKQIKDLTDDDPDMKHRRKLLESIPGIGEKT"
       "SAVLLAYIGLKDRFAHARQFAAFAGLTPRRYESGSSVRGASRMSKAGHVSLRRALYMPAMVATSKTEWGR"
       "AFRDRLAANGKKGKVILGAMMRKLAQVAYGVLKSGVPFDASRHNPVAA")
SITES = {168: "K", 193: "R", 231: "Y", 248: "K", 257: "A"}  # native pos -> mutant aa


@app.function(image=img, gpu="A10G", timeout=2400, volumes={"/models": vol})
def score(model_name: str, combos: list):
    import torch, esm, numpy as np
    torch.set_grad_enabled(False)
    model, alphabet = getattr(esm.pretrained, model_name)()
    model = model.eval().cuda()
    bc = alphabet.get_batch_converter()
    mask_idx = alphabet.mask_idx
    tok = {aa: alphabet.get_idx(aa) for aa in "ACDEFGHIKLMNPQRSTVWY"}

    def apply_muts(positions):
        s = list(SEQ)
        for p in positions:
            s[p - 1] = SITES[p]
        return "".join(s)

    def logp_at(seq_str, pos_list_mask):
        """mask each position in pos_list_mask (1-based) of seq_str, one seq per mask; return logsoftmax rows."""
        data = [(f"m{k}", seq_str) for k in range(len(pos_list_mask))]
        _, _, toks = bc(data)
        toks = toks.cuda()
        for k, p in enumerate(pos_list_mask):
            toks[k, p] = mask_idx  # +1 BOS offset handled: p is 1-based residue -> token index p
        lg = model(toks)["logits"]
        return torch.log_softmax(lg, dim=-1)

    out = {}
    for combo in combos:
        combo = list(combo)
        mut_seq = apply_muts(combo)
        add_terms, epi_terms = [], []
        # WT-context masked-marginals (additive baseline): mask pos i in WT seq
        lp_wt = logp_at(SEQ, combo)   # token index = residue number (1-based) because BOS shifts by +1 -> pos i residue at token i?
        # NOTE: residue r (1-based) -> token index r (BOS at 0). mask index passed as r.
        for k, p in enumerate(combo):
            wt = SEQ[p - 1]; mt = SITES[p]
            row = lp_wt[k, p]
            add_terms.append(float(row[tok[mt]] - row[tok[wt]]))
        # MUT-context masked-marginals: mask pos i in the full mutant seq (others carry mutations)
        lp_mut = logp_at(mut_seq, combo)
        for k, p in enumerate(combo):
            wt = SEQ[p - 1]; mt = SITES[p]
            row = lp_mut[k, p]
            epi_terms.append(float(row[tok[mt]] - row[tok[wt]]))
        key = "+".join(f"{SEQ[p-1]}{p}{SITES[p]}" for p in combo)
        out[key] = {"combo": combo, "additive": sum(add_terms),
                    "epistatic": sum(epi_terms),
                    "epistasis": sum(epi_terms) - sum(add_terms),
                    "add_terms": add_terms, "epi_terms": epi_terms}
    vol.commit()
    return {"model": model_name, "scores": out}


@app.local_entrypoint()
def main():
    pos = list(SITES.keys())
    combos = []
    for r in range(1, 6):
        combos += [list(c) for c in itertools.combinations(pos, r)]
    print(f"scoring {len(combos)} variants (incl. singles)")
    res = {}
    for mn in ["esm1v_t33_650M_UR90S_1", "esm2_t33_650M_UR50D"]:
        print("  ", mn)
        res[mn] = score.remote(mn, combos)
    json.dump(res, open("data/combo_scores.json", "w"))
    print("WROTE data/combo_scores.json")
