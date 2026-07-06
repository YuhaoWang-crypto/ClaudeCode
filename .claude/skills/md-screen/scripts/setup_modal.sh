#!/usr/bin/env bash
# One-time setup for running mdscreen on Modal from this sandboxed environment.
# Installs the Modal client + proxy support and verifies authentication.
set -euo pipefail

echo "==> Installing Modal client and proxy support"
pip install --quiet modal pyyaml 'python-socks' 'modal[api-proxy-support]'

echo "==> Checking credentials"
: "${MODAL_TOKEN_ID:?MODAL_TOKEN_ID not set}"
: "${MODAL_TOKEN_SECRET:?MODAL_TOKEN_SECRET not set}"

echo "==> Verifying auth with a trivial remote call"
python3 - <<'PY'
import modal
app = modal.App("mdscreen-auth-check")

@app.function()
def ping():
    return "pong"

with modal.enable_output():
    with app.run():
        assert ping.remote() == "pong"
print("AUTH OK — Modal is reachable and functions run remotely.")
PY

echo "==> Done. You can now: modal run modal_app.py::protein --pdb 1ubq --ns 1"
