"""Minimal pure-PyTorch stand-in for the one torch_cluster function gRNAde's
featurizer uses (knn_graph). Exact for a single (non-batched) point cloud."""
import torch

def knn_graph(x, k, batch=None, loop=False, flow="source_to_target", **kw):
    assert batch is None or int(batch.max()) == 0, "shim supports a single graph only"
    n = x.size(0)
    d = torch.cdist(x, x)
    if not loop:
        d = d + torch.eye(n, device=x.device) * float("inf")
    kk = min(k, n - (0 if loop else 1))
    idx = d.topk(kk, largest=False, dim=-1).indices          # [n, kk] neighbours of each node
    col = torch.arange(n, device=x.device).repeat_interleave(kk)
    row = idx.reshape(-1)
    return torch.stack([row, col], dim=0) if flow == "source_to_target" \
        else torch.stack([col, row], dim=0)
