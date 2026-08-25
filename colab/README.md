# Running the heavy jobs on Colab

## What I can and cannot do

I **cannot drive a Colab runtime**. There is no mechanism by which I hold an
interactive session on your notebook, and the link you sent returns a Google
sign-in page to me, so I cannot read its contents either.

There are three ways to connect the two, and they differ mostly in how much
trust they ask of you:

| pattern | how it works | what you hand over | verdict |
|---|---|---|---|
| **exchange through this repo** | I write a script here, you run three lines in Colab, the outputs come back as small files you paste or push | nothing | **recommended** |
| tunnel the runtime | `cloudflared`/`ngrok` + sshd in the notebook, I connect in | a live shell on your Google account's VM | possible, but it hands me arbitrary code execution against your account and Colab's terms discourage it — I would rather you didn't |
| paste | you run things and paste output | nothing | fine for small checks |

The first pattern needs nothing from you but a copy-paste, and everything stays
auditable. The scripts below are written for it.

## Why Colab at all — the two jobs this container genuinely cannot do

**1. Validate the full submission** (`run_vcc_prep.py`).

`vcc prep` holds the whole matrix as an in-memory CSR. Measured here:

| slice | cells | stored values | peak RSS | result |
|---|---|---|---|---|
| partial context | 40,000 | 238 M | 6.9 GB | passes |
| one context | 120,000 | 715 M | 14.0 GB | **SIGKILL** |
| full submission | 360,000 | 2,078 M | est. 40–60 GB | not attemptable |

This box has 15 GB. A Colab **High-RAM** runtime has ~51 GB and an A100
runtime ~83 GB, either of which clears it. The submission is inside the
scorer's own 4.75e9 cap, so this is a local-validator limit, not a spec
problem — but running it is the last unticked check.

**2. Test SE-600M as a context encoder** (`run_se_embed.py`).

The one use of Arc's SE-600M I have not been able to test. Its protein
embeddings helped slightly (+0.0012, and only for the 26 panel targets with no
measurement at all); its *cell* encoder is a different question and needs a GPU
plus the 2.7 GB checkpoint.

Be warned what the measurements already predict: the cross-context axis is the
one a context encoder would act on, and it is close to the measurement floor —
two runs of the same K562 screen agree at r = 0.319 against this model's 0.328.
The `transferability-prior-eval` skill also measured an *oracle* context feature
scoring **worse** than no context feature at all. So this is worth testing and
worth not expecting much from.

The script does the GPU part only — embed control cells, save a small `.npz` —
and leaves the analysis here, including the track-resolution gate that decides
whether the encoder can even tell these cell lines apart.

## Running them

Open a new Colab notebook, pick the runtime each script names, and paste:

```python
!git clone --depth 1 -b claude/virtual-cell-model-prediction-qu31ry \
    https://github.com/yuhaowang-crypto/claudecode.git /content/repo
%cd /content/repo
!python colab/run_vcc_prep.py        # or colab/run_se_embed.py
```

Each script prints a short block at the end marked `=== PASTE BACK ===`. That
block is all I need; paste it into the chat.

## About your VCC token

`run_vcc_prep.py` needs the official controls bundle, which needs your API key.
Put it in Colab's **Secrets** panel (the key icon in the left sidebar) as
`VCC_TOKEN` — the script reads it from there. Do not paste it into the chat, and
do not paste it into a notebook cell that gets saved.

The token you posted earlier in this conversation should be treated as
compromised and rotated: `vcc auth rotate` in your own terminal.
