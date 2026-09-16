"""Step 1 -- define the search space and lock the validation set.

The search space is the SURFY in-silico surfaceome (Bausch-Fluck et al.,
PNAS 2018).  Each surface protein carries its own topology string, so the
extracellular domain length and the number of transmembrane passes are
parsed per gene rather than assumed.

The pre-registered validation set (clinically validated LUAD antigens) and
the negative controls are written to disk *before* any expression data is
touched, so the recall numbers in step 7 cannot be tuned after the fact.
"""

from __future__ import annotations

import json
import re
import urllib.request
import warnings

import pandas as pd

from . import config as C

warnings.filterwarnings("ignore", category=UserWarning)

EXPERIMENTAL_SOURCES = {"pos. trainingset", "GPI (UniProt)"}
_SEG = re.compile(r"(SP|NC|TM|CY):(\d+)-(\d+)")


def _download(url: str, path):
    if path.exists():
        return path
    print(f"  downloading {url}")
    urllib.request.urlretrieve(url, path)
    return path


def parse_topology(topology: str) -> dict:
    """Turn 'SP:1-24;NC:25-305;TM:306-330;CY:331-362' into measurements."""
    out = {
        "n_tm": 0,
        "ecd_length": 0,
        "max_ecd_segment": 0,
        "cyto_length": 0,
        "has_signal_peptide": False,
        "topology_class": "unknown",
    }
    if not isinstance(topology, str) or not topology:
        return out
    segs = [(k, int(a), int(b)) for k, a, b in _SEG.findall(topology)]
    if not segs:
        return out
    nc = [b - a + 1 for k, a, b in segs if k == "NC"]
    cy = [b - a + 1 for k, a, b in segs if k == "CY"]
    out["n_tm"] = sum(1 for k, _, _ in segs if k == "TM")
    out["ecd_length"] = sum(nc)
    out["max_ecd_segment"] = max(nc) if nc else 0
    out["cyto_length"] = sum(cy)
    out["has_signal_peptide"] = any(k == "SP" for k, _, _ in segs)
    order = [k for k, _, _ in segs]
    if out["n_tm"] == 0:
        out["topology_class"] = "no TM (GPI/peripheral)"
    elif out["n_tm"] == 1:
        first_non_sp = next((k for k in order if k != "SP"), "")
        out["topology_class"] = "type I" if first_non_sp == "NC" else "type II"
    else:
        out["topology_class"] = f"multi-pass ({out['n_tm']} TM)"
    return out


def run() -> pd.DataFrame:
    print("[S1] search space: SURFY in-silico surfaceome")
    _download(C.SURFY_URL, C.SURFY_FILE)
    raw = pd.read_excel(C.SURFY_FILE, sheet_name="SurfaceomeMasterTable", skiprows=1)
    surface = raw[raw["Surfaceome Label"] == "surface"].copy()
    print(f"  SURFY master table: {len(raw)} proteins, {len(surface)} labelled surface")

    n_labelled = len(surface)
    surface["gene"] = surface["UniProt gene"].astype(str).str.split().str[0]
    surface = surface[surface["gene"].notna() & (surface["gene"] != "nan")]
    n_with_symbol = len(surface)
    surface["is_experimental"] = surface["Surfaceome Label Source"].isin(EXPERIMENTAL_SOURCES)
    topo = surface["topology"].apply(parse_topology).apply(pd.Series)
    surface = pd.concat([surface.reset_index(drop=True), topo.reset_index(drop=True)], axis=1)

    # One isoform per gene: keep the entry with the longest extracellular domain,
    # which is the one an antibody would be raised against.
    surface = surface.sort_values(["gene", "ecd_length", "length"], ascending=[True, False, False])
    first = surface.groupby("gene", as_index=False).first()
    exp_any = surface.groupby("gene")["is_experimental"].any().rename("surface_experimental")
    genes = first.merge(exp_any, on="gene")

    genes["accessible_sites"] = (
        genes["peps with accessible noncyt. nxst"].fillna(0)
        + genes["peps with accessible noncyt. Trp"].fillna(0)
        + genes["peps with accessible noncyt. Tyr"].fillna(0)
    )
    genes["cspa_confirmed"] = genes["CSPA category"].astype(str).str.contains("high confidence")
    genes["surface_evidence"] = genes["surface_experimental"].map(
        {True: "experimental", False: "machine_learning"}
    )

    out = genes[[
        "gene", "UniProt accession", "UniProt description", "Ensembl gene",
        "surface_evidence", "cspa_confirmed", "topology", "topology_class",
        "n_tm", "ecd_length", "max_ecd_segment", "cyto_length",
        "has_signal_peptide", "accessible_sites", "Membranome Almen main-class",
        "UniProt subcellular", "MachineLearning score",
    ]].rename(columns={
        "UniProt accession": "uniprot",
        "UniProt description": "protein_name",
        "Ensembl gene": "ensembl_gene",
        "Membranome Almen main-class": "almen_class",
        "UniProt subcellular": "uniprot_subcellular",
        "MachineLearning score": "ml_score",
    })
    out = out.sort_values("gene").reset_index(drop=True)
    out.to_csv(C.RESULTS / "s1_surfaceome.csv", index=False)

    n_exp = int((out.surface_evidence == "experimental").sum())
    n_ml = int((out.surface_evidence == "machine_learning").sum())
    print(f"  search space: {len(out)} surface genes "
          f"({n_exp} experimental evidence, {n_ml} machine-learning prediction)")

    # ---- lock the validation set before anything is ranked -----------------
    registry = {
        "locked_before_ranking": True,
        "positives": {
            g: {"rationale": r, "in_search_space": bool((out.gene == g).any())}
            for g, r in C.VALIDATION_POSITIVES.items()
        },
        "negative_controls": {
            g: {"rationale": r, "in_search_space": bool((out.gene == g).any())}
            for g, r in C.NEGATIVE_CONTROLS.items()
        },
        "recall_k": list(C.RECALL_K),
        "negative_control_fail_rank": C.NEGATIVE_CONTROL_FAIL_RANK,
    }
    (C.RESULTS / "s1_validation_registry.json").write_text(json.dumps(registry, indent=2))
    in_space = sum(v["in_search_space"] for v in registry["positives"].values())
    neg_in = sum(v["in_search_space"] for v in registry["negative_controls"].values())
    print(f"  validation set locked: {in_space}/{len(C.VALIDATION_POSITIVES)} positives "
          f"and {neg_in}/{len(C.NEGATIVE_CONTROLS)} negative controls inside the search space")

    C.write_provenance("s1", {
        "source": C.SURFY_URL,
        "proteins_in_master_table": int(len(raw)),
        "surface_labelled_proteins": int(n_labelled),
        "surface_labelled_with_gene_symbol": int(n_with_symbol),
        "unique_surface_genes": int(len(out)),
        "experimental_evidence_genes": n_exp,
        "machine_learning_genes": n_ml,
        "cspa_high_confidence_genes": int(out.cspa_confirmed.sum()),
        "positives_in_search_space": in_space,
        "negative_controls_in_search_space": neg_in,
    })
    return out


if __name__ == "__main__":
    run()
