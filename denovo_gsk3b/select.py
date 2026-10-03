"""Chemotype-aware top-N selection from the filter-cascade survivors.

Why not simply take the N highest-scoring survivors: the GA concentrates
probability mass on whichever chemotype the surrogate oracle scores best, so the
top of the ranking is dozens of decorations of a single core. Bemis-Murcko
scaffold uniqueness does not fix this -- swapping a pendant ring changes the
Murcko scaffold while leaving the recognition motif identical.

Selection therefore clusters the survivors by fingerprint similarity (Butina,
Tanimoto 0.4 within-cluster) and takes the best-scoring representative of each
of the N highest-scoring clusters. The result is N distinct chemotypes rather
than N analogues, at the cost of some score at the margin.
"""

from __future__ import annotations

from rdkit import Chem, DataStructs
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.ML.Cluster import Butina

from .filters import fingerprint, properties

BUTINA_CUTOFF = 0.6     # distance; 0.6 => members within Tanimoto 0.4 of each other


def murcko(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ""
    try:
        return MurckoScaffold.MurckoScaffoldSmiles(mol=mol)
    except Exception:
        return ""


def cluster(records: list[dict], cutoff: float = BUTINA_CUTOFF):
    """Butina-cluster records by Morgan/Tanimoto. Returns (clusters, records)."""
    paired = [(r, fingerprint(r["smiles"])) for r in records]
    paired = [(r, f) for r, f in paired if f is not None]
    recs = [r for r, _ in paired]
    fps = [f for _, f in paired]
    if not recs:
        return [], []
    dists: list[float] = []
    for i in range(1, len(fps)):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
        dists.extend(1.0 - s for s in sims)
    clusters = Butina.ClusterData(dists, len(fps), cutoff, isDistData=True)
    return clusters, recs


def select_representatives(records: list[dict], n: int = 10,
                           cutoff: float = BUTINA_CUTOFF,
                           id_prefix: str = "GSK3B-DN") -> list[dict]:
    """Best-scoring member of each of the n highest-scoring chemotype clusters."""
    clusters, recs = cluster(records, cutoff=cutoff)
    if not clusters:
        return []
    ranked = sorted(clusters,
                    key=lambda c: max(recs[i]["score"] for i in c),
                    reverse=True)
    chosen: list[dict] = []
    for rank, members in enumerate(ranked[:n], start=1):
        best = max(members, key=lambda i: recs[i]["score"])
        rec = dict(recs[best])
        rec["cluster_rank"] = rank
        rec["cluster_size"] = len(members)
        rec["scaffold"] = murcko(rec["smiles"])
        rec.update(properties(rec["smiles"]))
        chosen.append(rec)
    for i, rec in enumerate(chosen, start=1):
        rec["design_id"] = f"{id_prefix}-{i:02d}"
    return chosen
