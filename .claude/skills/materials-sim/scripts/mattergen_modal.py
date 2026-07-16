"""
MatterGen on Modal (GPU).

在 Modal 云端 GPU 上安装 MatterGen 并运行无条件晶体生成 demo，
把生成的 CIF 结果拉回本地。

用法:
    # 凭据已通过环境变量 MODAL_TOKEN_ID / MODAL_TOKEN_SECRET 提供
    modal run mattergen_modal.py                      # 默认: mattergen_base, batch 16
    modal run mattergen_modal.py --model-name dft_mag_density \
        --batch-size 16 --extra "--properties_to_condition_on={'dft_mag_density':0.15} --diffusion_guidance_factor=2.0"
"""
import modal

# ---- 镜像: python3.10 + 克隆 mattergen + uv 安装(cu118, 适配 GPU) ----
image = (
    modal.Image.debian_slim(python_version="3.10")
    .apt_install("git", "build-essential")
    .run_commands(
        "git clone --depth 1 https://github.com/microsoft/mattergen.git /mattergen",
        "pip install uv",
        # editable 安装; pyproject 已固定 torch==2.2.1+cu118 与 pyg cu118 wheels
        "cd /mattergen && uv pip install --system -e .",
    )
    # HF 下载缓存目录
    .env({"HF_HOME": "/root/.cache/huggingface"})
)

app = modal.App("mattergen-demo", image=image)

# 持久化生成结果的卷(可选)
vol = modal.Volume.from_name("mattergen-results", create_if_missing=True)


@app.function(gpu="T4", timeout=3600, volumes={"/out": vol})
def generate(
    model_name: str = "mattergen_base",
    batch_size: int = 16,
    num_batches: int = 1,
    extra: str = "",
) -> dict:
    """在 GPU 上运行 mattergen-generate, 返回 CIF zip 的字节内容。"""
    import subprocess
    import pathlib

    import torch

    print("torch:", torch.__version__, "cuda available:", torch.cuda.is_available())
    if torch.cuda.is_available():
        print("GPU:", torch.cuda.get_device_name(0))

    results_path = f"/out/{model_name}"
    cmd = [
        "mattergen-generate",
        results_path,
        f"--pretrained-name={model_name}",
        f"--batch_size={batch_size}",
        f"--num_batches={num_batches}",
    ]
    if extra:
        cmd += extra.split()

    print("运行:", " ".join(cmd))
    subprocess.run(cmd, cwd="/mattergen", check=True)

    vol.commit()

    out = {}
    for name in ("generated_crystals_cif.zip", "generated_crystals.extxyz"):
        p = pathlib.Path(results_path) / name
        if p.exists():
            out[name] = p.read_bytes()
            print(f"生成文件 {name}: {p.stat().st_size} bytes")
    return out


@app.local_entrypoint()
def main(
    model_name: str = "mattergen_base",
    batch_size: int = 16,
    num_batches: int = 1,
    extra: str = "",
):
    import pathlib

    out = generate.remote(
        model_name=model_name,
        batch_size=batch_size,
        num_batches=num_batches,
        extra=extra,
    )
    dest = pathlib.Path("modal_results") / model_name
    dest.mkdir(parents=True, exist_ok=True)
    for name, data in out.items():
        (dest / name).write_bytes(data)
        print(f"已保存到本地: {dest / name} ({len(data)} bytes)")
    if not out:
        print("警告: 未返回任何结果文件")
