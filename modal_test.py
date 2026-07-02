import modal

app = modal.App("bionemo-test")

@app.function(gpu="T4", timeout=120)
def check_gpu():
    import subprocess, sys
    result = subprocess.run(["nvidia-smi"], capture_output=True, text=True)
    print(result.stdout)

    import torch
    print(f"PyTorch: {torch.__version__}")
    print(f"CUDA available: {torch.cuda.is_available()}")
    if torch.cuda.is_available():
        print(f"GPU: {torch.cuda.get_device_name(0)}")
        print(f"VRAM: {torch.cuda.get_device_properties(0).total_memory / 1e9:.1f} GB")
    return torch.cuda.get_device_name(0) if torch.cuda.is_available() else "no gpu"

@app.local_entrypoint()
def main():
    result = check_gpu.remote()
    print(f"\nGPU confirmed: {result}")
