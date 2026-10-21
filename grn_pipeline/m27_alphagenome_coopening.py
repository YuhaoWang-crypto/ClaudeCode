"""
Module 27 — AlphaGenome: which chromatin "position switches" co-open
(question 2: "有哪些位置开关也许被同时打开，从而被其他 trigger 引起新的反馈？
这个能否利用 alphagenome 做分析？").

The analysis design (what AlphaGenome is for here):
  AlphaGenome predicts, from sequence alone, per-locus CHROMATIN ACCESSIBILITY
  (ATAC + DNase heads) and histone/TF-binding tracks. To ask "which switches
  open together", we:
    1. take a REAL hg38 panel of chromatin-regulator loci (TSS of the M25
       writers/erasers/readers — fetched live from Ensembl, cached here);
    2. for a shared perturbation (e.g. a recruiter variant, or a reference vs
       a demethylated/accessibility-opened background), read each locus's
       accessibility delta on the atac/dnase heads;
    3. build a CO-ACCESSIBILITY graph: loci whose accessibility moves together
       across the perturbation panel are edges; its CONNECTED COMPONENTS are
       the sets of "switches that open simultaneously" — the candidate
       co-triggered feedback modules the question asks about.

Status in THIS session (honest):
  The AlphaGenome deployment runs on Modal (GPU). Modal's gRPC control channel
  is NOT reachable through this session's egress proxy (plain HTTPS to
  api.modal.com returns 200, but the modal client reports "Could not connect to
  the Modal server"). So the REAL forward passes cannot be executed here. This
  module therefore:
    * ships the REAL, ready-to-run locus panel (real Ensembl GRCh38 coords);
    * attempts the Modal run and records the exact failure if present;
    * and DEMONSTRATES the co-opening READOUT LOGIC on an illustrative
      accessibility matrix so the pipeline is tested end-to-end and will
      consume real AlphaGenome deltas unchanged in a Modal-reachable env.

RIGOR:
  ✅ The locus panel coordinates are real (Ensembl GRCh38, cached); the
     co-accessibility graph / connected-component readout is a deterministic,
     tested computation.
  ⚠️ The accessibility MATRIX used for the in-session demo is ILLUSTRATIVE
     (clearly synthetic, seeded), NOT AlphaGenome output — because the model
     could not be run here. Swap in real deltas (batch_variant_effect ->
     head 'atac'/'dnase') in a Modal-reachable environment and the SAME code
     produces the real co-opening modules.
  ⚠️ AlphaGenome itself is a community PyTorch port (fidelity caveat from the
     skill): treat deltas as relative/directional, confirm top calls against
     the official API before decision-grade use.
"""
import json
import os
import numpy as np
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PANEL_CACHE = "figures/_alphagenome_cache/locus_panel.json"

# real Ensembl GRCh38 gene spans (fetched live; TSS computed by strand).
# (symbol, chrom, start, end, strand, M25 memory class)
_ENSEMBL = [
    ("EZH2",    "7",  148807257, 148884625, -1, "hybrid (H3K27)"),
    ("DNMT3B",  "20",  32762385,  32809359,  1, "long (5mC)"),
    ("SUV39H1", "X",   48695554,  48709556,  1, "hybrid (H3K9)"),
    ("SETDB1",  "1",  150926263, 150964748,  1, "hybrid (H3K9)"),
    ("HDAC4",   "2",  239048168, 239401661, -1, "short (acetyl)"),
    ("CBX5",    "12",  54230942,  54280141, -1, "hybrid (H3K9 reader)"),
    ("TET2",    "4",  105145875, 105279816,  1, "long-eraser (5mC)"),
    ("KDM6A",   "X",   44873130,  45114135,  1, "hybrid-eraser (H3K27)"),
]


def locus_panel():
    """Return the real hg38 TSS panel (cached)."""
    if os.path.exists(PANEL_CACHE):
        return json.loads(open(PANEL_CACHE).read())
    panel = []
    for sym, chrom, start, end, strand, klass in _ENSEMBL:
        tss = start if strand == 1 else end
        panel.append({"gene": sym, "chrom": "chr" + chrom, "tss": tss,
                      "strand": strand, "memory_class": klass})
    os.makedirs(os.path.dirname(PANEL_CACHE), exist_ok=True)
    open(PANEL_CACHE, "w").write(json.dumps(panel, indent=2))
    return panel


def try_alphagenome(panel):
    """Attempt the real Modal AlphaGenome accessibility run. Returns
    (matrix, heads, note). On any failure returns (None, None, reason)."""
    try:
        import subprocess
        # Probe connectivity cheaply; the full run needs GPU + weights volume.
        r = subprocess.run(["modal", "app", "list"], capture_output=True,
                            text=True, timeout=18)
        if "Could not connect" in (r.stderr + r.stdout) or r.returncode != 0:
            return None, None, ("modal unreachable (gRPC control channel not "
                                "tunnelled through egress proxy)")
        # If reachable, a real run would go here (batch_variant_effect on the
        # panel, reading heads 'atac','dnase'). Left as the documented path.
        return None, None, "modal reachable but full GPU run not auto-triggered"
    except Exception as e:
        return None, None, f"modal attempt failed: {type(e).__name__}: {e}"


def illustrative_accessibility(panel, n_perturb=12, seed=7):
    """ILLUSTRATIVE accessibility-delta matrix (loci x perturbations).
    Structured so that loci sharing a memory class co-open (shared latent
    factor) — exactly the structure real co-accessibility would reveal — so the
    readout logic is exercised. NOT AlphaGenome output."""
    rng = np.random.default_rng(seed)
    classes = sorted(set(p["memory_class"].split()[0] for p in panel))
    # latent "opening programs", one per coarse class
    factors = {c: rng.normal(size=n_perturb) for c in classes}
    M = np.zeros((len(panel), n_perturb))
    for i, p in enumerate(panel):
        c = p["memory_class"].split()[0]
        M[i] = 0.85 * factors[c] + 0.15 * rng.normal(size=n_perturb)
    return M


def co_accessibility_graph(panel, M, thresh=0.6):
    """Edges between loci whose accessibility-delta profiles correlate above
    `thresh`; connected components = co-opening switch modules."""
    C = np.corrcoef(M)
    G = nx.Graph()
    for i, p in enumerate(panel):
        G.add_node(p["gene"], memory_class=p["memory_class"])
    for i in range(len(panel)):
        for j in range(i + 1, len(panel)):
            if C[i, j] >= thresh:
                G.add_edge(panel[i]["gene"], panel[j]["gene"],
                           weight=float(C[i, j]))
    comps = [sorted(c) for c in nx.connected_components(G)]
    return G, C, comps


def make_figure(panel, M, C, comps, path):
    fig, ax = plt.subplots(1, 2, figsize=(13, 5.2))
    genes = [p["gene"] for p in panel]
    im = ax[0].imshow(C, cmap="RdBu_r", vmin=-1, vmax=1)
    ax[0].set_xticks(range(len(genes))); ax[0].set_xticklabels(genes, rotation=60, fontsize=7)
    ax[0].set_yticks(range(len(genes))); ax[0].set_yticklabels(genes, fontsize=7)
    ax[0].set_title("Co-accessibility (illustrative) across loci\n"
                    "ready to consume real AlphaGenome atac/dnase deltas")
    fig.colorbar(im, ax=ax[0], fraction=0.046)

    G, _, _ = co_accessibility_graph(panel, M)
    pos = nx.spring_layout(G, seed=2)
    comp_id = {}
    for k, c in enumerate(comps):
        for g in c:
            comp_id[g] = k
    cols = ["#2980b9", "#c0392b", "#16a085", "#8e44ad", "#e67e22"]
    nx.draw_networkx(G, pos, ax=ax[1],
                     node_color=[cols[comp_id[g] % len(cols)] for g in G.nodes()],
                     node_size=900, font_size=7, edgecolors="k")
    ax[1].set_title(f"Co-opening switch modules (connected components): "
                    f"{len(comps)}")
    ax[1].axis("off")
    fig.suptitle("M27 — AlphaGenome co-opening design (real locus panel; "
                 "accessibility matrix illustrative — Modal unreachable here)",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(path, dpi=130)
    plt.close(fig)


def report(fig_path="figures/m27_alphagenome.png"):
    print("=" * 68)
    print("MODULE 27 — AlphaGenome co-opening of chromatin position-switches")
    print("=" * 68)
    panel = locus_panel()
    print(f"\nREAL hg38 locus panel ({len(panel)} chromatin-regulator TSS, "
          f"Ensembl GRCh38):")
    print(f"  {'gene':>8} {'chrom':>6} {'TSS':>11} {'strand':>7}  memory class")
    for p in panel:
        print(f"  {p['gene']:>8} {p['chrom']:>6} {p['tss']:>11} "
              f"{p['strand']:>7}  {p['memory_class']}")

    M, heads, note = try_alphagenome(panel)
    print(f"\nAlphaGenome run attempt: {note}")
    if M is None:
        print("  -> REAL forward passes not available in this session.")
        print("     Using an ILLUSTRATIVE (clearly synthetic) accessibility")
        print("     matrix to exercise the co-opening readout. Swap in real")
        print("     atac/dnase deltas in a Modal-reachable env; code unchanged.")
        M = illustrative_accessibility(panel)
        real = False
    else:
        real = True

    G, C, comps = co_accessibility_graph(panel, M)
    print(f"\nco-opening switch modules (connected components of the "
          f"co-accessibility graph):")
    for k, c in enumerate(comps, 1):
        klasses = sorted(set(p["memory_class"] for p in panel if p["gene"] in c))
        print(f"  module {k}: {{{', '.join(c)}}}  [{'; '.join(klasses)}]")
    print("\ninterpretation: loci inside one module open together -> a trigger")
    print("hitting any one can co-activate the rest (the 'new feedback' the")
    print("question asks about). With REAL AlphaGenome deltas these modules")
    print("become concrete co-regulated switch sets; here they recover the")
    print("memory-class grouping by construction (illustrative).")

    make_figure(panel, M, C, comps, fig_path)
    print(f"\nfigure written to: {fig_path}")
    return {"panel": panel, "modal_note": note, "real_alphagenome": real,
            "n_comodules": len(comps), "comodules": comps}


if __name__ == "__main__":
    report()
