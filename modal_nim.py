"""NVIDIA NIM API calls from Modal (open internet, no proxy)"""
import modal

app = modal.App("nvidia-nim-api")

image = modal.Image.debian_slim(python_version="3.11").pip_install("requests")

@app.function(image=image, timeout=120, secrets=[modal.Secret.from_name("nvidia-nim-key")])
def test_esm2_nim():
    import os, requests

    api_key = os.environ["NGC_API_KEY"]
    headers = {"Authorization": f"Bearer {api_key}", "Content-Type": "application/json"}
    payload = {
        "sequence": "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEK",
        "pooling": "mean",
    }
    # Try multiple known endpoint variants
    candidates = [
        "https://health.api.nvidia.com/v1/biology/meta/esm2-650m/embeddings",
        "https://integrate.api.nvidia.com/v1/biology/meta/esm2-650m/embeddings",
        "https://health.api.nvidia.com/v1/biology/nvidia/esm2-650m/embeddings",
        "https://health.api.nvidia.com/v1/biology/meta/esm2/embeddings",
    ]
    for url in candidates:
        print(f"Trying: {url}")
        r = requests.post(url, headers=headers, json=payload, timeout=30)
        print(f"  Status: {r.status_code}")
        if r.status_code == 200:
            data = r.json()
            emb = data.get("embedding", [])
            print(f"  Embedding dim: {len(emb)}")
            print(f"  First 5: {[round(v, 4) for v in emb[:5]]}")
            return {"dim": len(emb), "url": url}
        elif r.status_code != 404:
            print(f"  Body: {r.text[:200]}")
    return {"dim": 0, "error": "no working endpoint found"}


@app.function(image=image, timeout=120, secrets=[modal.Secret.from_name("nvidia-nim-key")])
def test_esmfold_nim():
    import os, requests

    api_key = os.environ["NGC_API_KEY"]
    url = "https://health.api.nvidia.com/v1/biology/meta/esmfold"
    headers = {"Authorization": f"Bearer {api_key}", "Content-Type": "application/json"}
    # 短序列避免超时
    payload = {"sequence": "MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNG"}
    print("Calling NIM ESMFold API...")
    r = requests.post(url, headers=headers, json=payload, timeout=60)
    print(f"Status: {r.status_code}")
    r.raise_for_status()
    data = r.json()
    pdb = data.get("pdbs", [""])[0]
    atom_count = pdb.count("\nATOM")
    print(f"PDB output: {len(pdb)} chars, {atom_count} ATOM records")
    print(pdb[:300])
    return {"pdb_len": len(pdb), "atom_count": atom_count}


@app.function(image=image, timeout=120, secrets=[modal.Secret.from_name("nvidia-nim-key")])
def test_genmol_nim():
    import os, requests

    api_key = os.environ["NGC_API_KEY"]
    url = "https://health.api.nvidia.com/v1/biology/nvidia/genmol/generate"
    headers = {"Authorization": f"Bearer {api_key}", "Content-Type": "application/json"}
    payload = {
        "smiles": "[*{20-30}]",
        "num_molecules": 10,
        "temperature": "1.0",
        "noise": "1.0",
        "step_size": 1,
        "scoring": "QED",
        "unique": True,
    }
    print("Calling NIM GenMol API (de novo generation)...")
    r = requests.post(url, headers=headers, json=payload, timeout=120)
    print(f"Status: {r.status_code}")
    r.raise_for_status()
    data = r.json()
    mols = data.get("molecules", [])
    mols_sorted = sorted(mols, key=lambda m: m["score"], reverse=True)
    print(f"Generated {len(mols_sorted)} molecules")
    for i, m in enumerate(mols_sorted[:5], 1):
        print(f"  #{i:2d} QED={m['score']:.4f}  {m['smiles']}")
    return mols_sorted


@app.local_entrypoint()
def main():
    print("="*60)
    print("NVIDIA NIM APIs via Modal (open internet)")
    print("="*60)

    print("\n[1/3] ESM-2 650M Embedding (probing endpoints)...")
    r1 = test_esm2_nim.remote()
    print(f"Done. dim={r1.get('dim', 'N/A')}, url={r1.get('url', r1.get('error'))}")

    print("\n[2/3] ESMFold Structure Prediction...")
    r2 = test_esmfold_nim.remote()
    print(f"Done. PDB={r2['pdb_len']} chars, {r2['atom_count']} ATOM records")

    print("\n[3/3] GenMol De Novo Molecule Generation...")
    r3 = test_genmol_nim.remote()
    print(f"Done. {len(r3)} molecules generated")

    print("\n✅ NIM API calls via Modal completed!")
