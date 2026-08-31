"""Minimal pure-PyTorch stand-in for torch_scatter.scatter_add (the only function
gRNAde's src/layers.py imports). Exact for arbitrary dim / broadcasting index."""
import torch

def scatter_add(src, index, dim=-1, out=None, dim_size=None):
    dim = dim % src.dim()
    if dim_size is None:
        dim_size = int(index.max()) + 1 if index.numel() > 0 else 0
    idx = index
    if idx.dim() != src.dim():
        shape = [1] * src.dim(); shape[dim] = -1
        idx = idx.view(shape).expand_as(src)
    if out is None:
        size = list(src.shape); size[dim] = dim_size
        out = torch.zeros(size, dtype=src.dtype, device=src.device)
    return out.scatter_add_(dim, idx, src)
