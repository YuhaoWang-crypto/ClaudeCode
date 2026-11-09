"""Dual-head residual network f_theta(s) = (p, v), as in AlphaGo Zero.

Paper:  conv(256,3x3)-BN-ReLU stem, 19 or 39 residual blocks,
        policy head: conv(2,1x1)-BN-ReLU-FC(n*n+1)
        value  head: conv(1,1x1)-BN-ReLU-FC(256)-ReLU-FC(1)-tanh
Here:   identical shape, just narrower / shallower (CPU budget).
"""
import torch
import torch.nn as nn
import torch.nn.functional as F


class ResBlock(nn.Module):
    def __init__(self, ch):
        super().__init__()
        self.c1 = nn.Conv2d(ch, ch, 3, padding=1, bias=False)
        self.b1 = nn.BatchNorm2d(ch)
        self.c2 = nn.Conv2d(ch, ch, 3, padding=1, bias=False)
        self.b2 = nn.BatchNorm2d(ch)

    def forward(self, x):
        y = F.relu(self.b1(self.c1(x)))
        y = self.b2(self.c2(y))
        return F.relu(x + y)          # skip connection


class AlphaZeroNet(nn.Module):
    def __init__(self, n=7, in_planes=17, filters=64, blocks=4, value_hidden=64):
        super().__init__()
        self.n = n
        self.cfg = dict(n=n, in_planes=in_planes, filters=filters,
                        blocks=blocks, value_hidden=value_hidden)
        self.stem = nn.Conv2d(in_planes, filters, 3, padding=1, bias=False)
        self.stem_bn = nn.BatchNorm2d(filters)
        self.tower = nn.Sequential(*[ResBlock(filters) for _ in range(blocks)])
        # policy head
        self.p_conv = nn.Conv2d(filters, 2, 1, bias=False)
        self.p_bn = nn.BatchNorm2d(2)
        self.p_fc = nn.Linear(2 * n * n, n * n + 1)
        # value head
        self.v_conv = nn.Conv2d(filters, 1, 1, bias=False)
        self.v_bn = nn.BatchNorm2d(1)
        self.v_fc1 = nn.Linear(n * n, value_hidden)
        self.v_fc2 = nn.Linear(value_hidden, 1)

    def forward(self, x):
        x = F.relu(self.stem_bn(self.stem(x)))
        x = self.tower(x)
        p = F.relu(self.p_bn(self.p_conv(x))).flatten(1)
        p = self.p_fc(p)                                   # logits
        v = F.relu(self.v_bn(self.v_conv(x))).flatten(1)
        v = torch.tanh(self.v_fc2(F.relu(self.v_fc1(v)))).squeeze(1)
        return p, v


def build(cfg):
    return AlphaZeroNet(n=cfg.board_size, in_planes=2 * cfg.history + 1,
                        filters=cfg.filters, blocks=cfg.blocks,
                        value_hidden=cfg.value_hidden)


def load(path, map_location="cpu"):
    ck = torch.load(path, map_location=map_location, weights_only=False)
    net = AlphaZeroNet(**ck["net_cfg"])
    net.load_state_dict(ck["state_dict"])
    net.eval()
    return net, ck
