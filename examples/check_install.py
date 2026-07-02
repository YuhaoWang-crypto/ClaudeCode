"""Sanity check: verify the stack is installed and report GPU/token status.

Runs WITHOUT downloading any gated model weights, so it works before you have
set up HuggingFace access. Run:  python examples/check_install.py
"""

import importlib
import os

try:  # pick up HF_TOKEN from a local .env if present
    from dotenv import load_dotenv

    load_dotenv()
except Exception:  # noqa: BLE001
    pass


def check(name: str, dist: str | None = None) -> str:
    try:
        mod = importlib.import_module(name)
        ver = getattr(mod, "__version__", None)
        if ver is None:
            from importlib.metadata import version as _v

            try:
                ver = _v(dist or name)
            except Exception:  # noqa: BLE001
                ver = "?"
        return f"  ok   {name:16s} {ver}"
    except Exception as exc:  # noqa: BLE001
        return f"  MISS {name:16s} ({exc})"


def main() -> None:
    print("Dependency check:")
    print(check("torch"))
    print(check("ase"))
    print(check("fairchem", dist="fairchem-core"))
    print(check("huggingface_hub"))
    print(check("pymatgen") + "   (optional; ASE already reads CIF)")

    print("\nCompute:")
    try:
        import torch

        cuda = torch.cuda.is_available()
        print(f"  CUDA available : {cuda}")
        if cuda:
            print(f"  GPU            : {torch.cuda.get_device_name(0)}")
        else:
            print("  Running on CPU — fine for small tests, slow for big cells.")
    except Exception as exc:  # noqa: BLE001
        print(f"  torch not usable: {exc}")

    tok = os.environ.get("HF_TOKEN") or os.environ.get("HUGGING_FACE_HUB_TOKEN")
    print("\nHuggingFace access:")
    print("  HF token       :", "set" if tok else "NOT set")
    print("  Needed to download gated UMA / OMAT24 weights.")
    print("  1) create account, 2) accept license at")
    print("     https://huggingface.co/facebook/UMA")
    print("  3) `huggingface-cli login` or `export HF_TOKEN=...`")


if __name__ == "__main__":
    main()
