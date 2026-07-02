"""
EGFR Multi-Model Pipeline on Modal Cloud GPU
Runs: ESM-2 650M embedding + ESMFold structure + GenMol drug generation
Saves results to /tmp/egfr_results.json
"""
import modal, json

app = modal.App("egfr-pipeline")

image = (
    modal.Image.debian_slim(python_version="3.11")
    .pip_install("torch", "transformers", "sentencepiece", "requests", "numpy")
)

# EGFR kinase domain fragment (72 aa) for structure prediction
EGFR_PEPTIDE = "KVLGSGAFGTVYKGLWIPEGEKVKIPVAIKELREATSPKANKEILDEAYVMASVDNPHVCRLLGICL"

# Known EGFR-related sequences for embedding comparison
SEQUENCES = {
    "EGFR_kinase":   EGFR_PEPTIDE,
    "EGFR_short":    "KVLGSGAFGTVYKGLWIPEGEKVKIPVAIKELR",
    "HER2_fragment": "MELAALCRWGLLLALLPPGAASTQVCTGTDMKLRLPASPETHLDMLRHLYQGCQVVQGNLELTYLPTNASLSFLNDQFTRGVNKSQQPRSKNLGLHQKSSQTQEGPEGQLEYKFPETLPYCAHSSVSAGRGTSSSGTPQHQDPHRRTSGPTQFCQESGGTNRPFPSGPWDPNYPNSTGRLDNALRQFLPPQDLCLGLLHYMDL",
    "KRAS_fragment":  "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSY",
    "p53_fragment":   "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQ",
}

@app.function(gpu="T4", image=image, timeout=600)
def run_esm2_embeddings():
    import torch
    from transformers import AutoTokenizer, AutoModel
    from torch.nn.functional import cosine_similarity

    print("=== ESM-2 650M Embeddings ===")
    tokenizer = AutoTokenizer.from_pretrained("facebook/esm2_t33_650M_UR50D")
    model = AutoModel.from_pretrained("facebook/esm2_t33_650M_UR50D").cuda().eval()

    names = list(SEQUENCES.keys())
    seqs = list(SEQUENCES.values())

    inputs = tokenizer(seqs, return_tensors="pt", padding=True, truncation=True, max_length=512)
    inputs = {k: v.cuda() for k, v in inputs.items()}

    with torch.no_grad():
        out = model(**inputs)

    mask = inputs["attention_mask"].unsqueeze(-1).float()
    emb = (out.last_hidden_state * mask).sum(1) / mask.sum(1)

    # Pairwise cosine similarity matrix
    n = len(names)
    sim_matrix = {}
    for i in range(n):
        for j in range(n):
            sim = cosine_similarity(emb[i:i+1], emb[j:j+1]).item()
            sim_matrix[f"{names[i]}|{names[j]}"] = round(sim, 4)

    print(f"GPU: {torch.cuda.get_device_name(0)}")
    print(f"Embedding dim: {emb.shape[1]}")
    for name in names:
        print(f"  {name}: {len(SEQUENCES[name])} aa")

    return {
        "model": "ESM-2 650M",
        "gpu": torch.cuda.get_device_name(0),
        "embedding_dim": emb.shape[1],
        "sequences": {n: {"length": len(s), "name": n} for n, s in SEQUENCES.items()},
        "similarity_matrix": sim_matrix,
        "sequence_names": names,
    }


@app.function(image=image, timeout=120, secrets=[modal.Secret.from_name("nvidia-nim-key")])
def run_esmfold():
    import os, requests

    api_key = os.environ["NGC_API_KEY"]
    url = "https://health.api.nvidia.com/v1/biology/meta/esmfold"
    headers = {"Authorization": f"Bearer {api_key}", "Content-Type": "application/json"}

    print(f"=== ESMFold: predicting structure for {len(EGFR_PEPTIDE)} aa EGFR fragment ===")
    r = requests.post(url, headers=headers,
                      json={"sequence": EGFR_PEPTIDE}, timeout=90)
    r.raise_for_status()
    data = r.json()
    pdb = data.get("pdbs", [""])[0]

    # Extract b-factors (pLDDT confidence scores) from PDB
    plddt_scores = []
    for line in pdb.split("\n"):
        if line.startswith("ATOM"):
            try:
                plddt_scores.append(float(line[60:66]))
            except:
                pass

    mean_plddt = sum(plddt_scores) / len(plddt_scores) if plddt_scores else 0
    high_conf = sum(1 for s in plddt_scores if s >= 70) / len(plddt_scores) * 100 if plddt_scores else 0

    print(f"PDB length: {len(pdb)} chars")
    print(f"ATOM records: {pdb.count(chr(10)+'ATOM')}")
    print(f"Mean pLDDT: {mean_plddt:.1f}")
    print(f"High-confidence residues (≥70): {high_conf:.1f}%")

    return {
        "model": "ESMFold (NIM)",
        "sequence": EGFR_PEPTIDE,
        "sequence_length": len(EGFR_PEPTIDE),
        "pdb_length": len(pdb),
        "atom_count": pdb.count("\nATOM"),
        "mean_plddt": round(mean_plddt, 2),
        "high_confidence_pct": round(high_conf, 1),
        "pdb_preview": pdb[:500],
    }


@app.function(image=image, timeout=180, secrets=[modal.Secret.from_name("nvidia-nim-key")])
def run_genmol():
    import os, requests

    api_key = os.environ["NGC_API_KEY"]
    url = "https://health.api.nvidia.com/v1/biology/nvidia/genmol/generate"
    headers = {"Authorization": f"Bearer {api_key}", "Content-Type": "application/json"}

    print("=== GenMol: de novo drug-like molecule generation ===")
    payload = {
        "smiles": "[*{25-35}]",   # de novo, drug-sized
        "num_molecules": 30,
        "temperature": "1.0",
        "noise": "1.0",
        "step_size": 1,
        "scoring": "QED",
        "unique": True,
    }
    r = requests.post(url, headers=headers, json=payload, timeout=150)
    r.raise_for_status()
    data = r.json()
    mols = sorted(data.get("molecules", []), key=lambda m: m["score"], reverse=True)

    print(f"Generated: {len(mols)} molecules")
    for i, m in enumerate(mols[:5], 1):
        print(f"  #{i} QED={m['score']:.4f}  {m['smiles']}")

    return {
        "model": "GenMol (NIM)",
        "total_generated": len(mols),
        "molecules": [{"rank": i+1, "smiles": m["smiles"], "qed": round(m["score"], 4)}
                      for i, m in enumerate(mols)],
        "top5": mols[:5],
        "mean_qed": round(sum(m["score"] for m in mols) / len(mols), 4) if mols else 0,
    }


@app.local_entrypoint()
def main():
    import json, pathlib

    print("Starting EGFR multi-model pipeline on Modal...")

    # Run all 3 in parallel
    esm2_fut = run_esm2_embeddings.spawn()
    fold_fut  = run_esmfold.spawn()
    mol_fut   = run_genmol.spawn()

    results = {
        "esm2":    esm2_fut.get(),
        "esmfold": fold_fut.get(),
        "genmol":  mol_fut.get(),
    }

    out = pathlib.Path("egfr_project/data/modal_results.json")
    out.write_text(json.dumps(results, indent=2))
    print(f"\n✅ Saved to {out}")
    print(f"  ESM-2: {results['esm2']['embedding_dim']}D embeddings for {len(results['esm2']['sequences'])} proteins")
    print(f"  ESMFold: pLDDT={results['esmfold']['mean_plddt']}, conf={results['esmfold']['high_confidence_pct']}%")
    print(f"  GenMol: {results['genmol']['total_generated']} molecules, mean QED={results['genmol']['mean_qed']}")
