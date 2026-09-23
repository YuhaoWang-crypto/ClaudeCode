"""Create (or update) a Hugging Face Docker Space that serves Laya.

Usage:
    pip install -U huggingface_hub
    export HF_TOKEN=hf_xxx            # write-scoped token
    python deploy.py                  # Space name defaults to "laya-api"
    python deploy.py --name my-laya --private

Prints the Space URL and the generated LAYA_API_KEY (stored as a Space secret).
Re-running keeps the existing key unless --new-key is passed.
"""
import argparse
import os
import secrets
from pathlib import Path

from huggingface_hub import HfApi

HERE = Path(__file__).resolve().parent


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--name", default="laya-api")
    p.add_argument("--private", action="store_true")
    p.add_argument("--new-key", action="store_true", help="rotate LAYA_API_KEY")
    p.add_argument("--hardware", default=None,
                   help="e.g. cpu-upgrade, t4-small (paid); default free cpu-basic")
    args = p.parse_args()

    token = os.environ.get("HF_TOKEN")
    if not token:
        raise SystemExit("Set HF_TOKEN to a write-scoped Hugging Face token.")

    api = HfApi(token=token)
    user = api.whoami()["name"]
    repo_id = f"{user}/{args.name}"

    exists = api.repo_exists(repo_id, repo_type="space")
    api.create_repo(repo_id, repo_type="space", space_sdk="docker",
                    private=args.private, exist_ok=True)

    key = None
    if not exists or args.new_key:
        key = secrets.token_hex(32)
        api.add_space_secret(repo_id, "LAYA_API_KEY", key)

    api.upload_folder(repo_id=repo_id, repo_type="space", folder_path=str(HERE),
                      allow_patterns=["Dockerfile", "README.md"],
                      commit_message="Deploy Laya API")
    if args.hardware:
        api.request_space_hardware(repo_id, args.hardware)

    host = repo_id.replace("/", "-").replace("_", "-").lower()
    print(f"Space:    https://huggingface.co/spaces/{repo_id}")
    print(f"Endpoint: https://{host}.hf.space/v1/systemone")
    if key:
        print(f"LAYA_API_KEY={key}   (save this; it is only shown once)")
    else:
        print("LAYA_API_KEY unchanged (pass --new-key to rotate).")


if __name__ == "__main__":
    main()
