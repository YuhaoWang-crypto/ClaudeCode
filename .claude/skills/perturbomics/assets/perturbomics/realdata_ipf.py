"""REAL-DATA example: reverse an IPF disease signature with drugs + CRISPR KOs.

Downloads three genuine perturbation libraries from Enrichr (no login), builds
a consensus idiopathic-pulmonary-fibrosis signature, then ranks single-agent
reversers (drug + CRISPR in one table) and the best drug+CRISPR combinations.

    python3 -m perturbomics.realdata_ipf [cache_dir]

Data (cached in cache_dir, ~45 MB total, CC-BY / open):
    Disease_Perturbations_from_GEO_up/_down     (CREEDS disease DE signatures)
    LINCS_L1000_Chem_Pert_up/_down              (L1000 drug perturbations)
    LINCS_L1000_CRISPR_KO_Consensus_Sigs        (LINCS CRISPR KO signatures)

RIGOR: every score is a deterministic WTCS / coverage computation. "reverses IPF
in cells" and "this pair is synergistic" are HYPOTHESES to validate in a fibrosis
model — the ranked table is a prioritised starting point, not a result.
"""

from __future__ import annotations

import random
import sys

from . import enrichr, rank_reversers, best_combinations, Signature

DISEASE = "idiopathic pulmonary fibrosis"
N_DRUG_SAMPLE = 5000     # subsample the 33k L1000 drug profiles for speed
SEED = 0


def _drug_name(term: str) -> str:
    # "CPC001 HA1E 24H-cilnidipine-10.0" -> "cilnidipine"
    tail = term.split("-", 1)[1] if "-" in term else term
    return tail.rsplit("-", 1)[0]


def main(cache_dir: str = "enrichr_cache", disease: str = DISEASE,
         label: str | None = None) -> dict:
    label = label or "".join(w[:1].upper() for w in disease.split()) or "DIS"
    random.seed(SEED)
    print(f"Loading real perturbation libraries into {cache_dir}/ ...")

    # 1. consensus disease signature over its GEO replicates
    dis_up = enrichr.load_library("Disease_Perturbations_from_GEO_up", cache_dir)
    dis_dn = enrichr.load_library("Disease_Perturbations_from_GEO_down", cache_dir)
    ipf = enrichr.consensus_signature(
        dis_up, dis_dn, match=lambda t: disease.lower() in t.lower(),
        name=f"{label}_consensus", min_recurrence=2)
    print(f"[disease] {ipf}  (consensus of {ipf.meta['n_replicates']} replicates)")

    # 2. real L1000 drug signatures (subsampled)
    cu = enrichr.load_library("LINCS_L1000_Chem_Pert_up", cache_dir)
    cd = enrichr.load_library("LINCS_L1000_Chem_Pert_down", cache_dir)
    drugs_all = enrichr.paired_signatures(cu, cd, "drug", name_fn=_drug_name)
    drugs = random.sample(drugs_all, min(N_DRUG_SAMPLE, len(drugs_all)))
    print(f"[drugs]  {len(drugs_all)} L1000 profiles -> sampled {len(drugs)}")

    # 3. real CRISPR KO signatures
    ck = enrichr.load_library("LINCS_L1000_CRISPR_KO_Consensus_Sigs", cache_dir)
    crispr = enrichr.crispr_ko_signatures(ck)
    print(f"[crispr] {len(crispr)} CRISPR KO signatures")

    library = drugs + crispr

    # 4. single-agent reversers (most negative WTCS = strongest corrector)
    rev = rank_reversers(ipf, library, k=100)
    print(f"\n=== Top single reversers of {label} (drug + CRISPR, one table) ===")
    print(rev.head(10).to_string(index=False,
          formatters={"connectivity": "{:+.3f}".format}))
    print("\n--- top CRISPR-KO reversers ---")
    print(rev[rev.modality == "crispr_ko"].head(6).to_string(index=False,
          formatters={"connectivity": "{:+.3f}".format}))

    # 5. best drug + CRISPR combinations (binarize disease -> coverage in [0,1])
    ipf_bin = Signature.from_ranked(ipf.top(150, "up"), ipf.top(150, "down"),
                                    name="IPF_consensus", modality="disease")
    drug_names = (rev[(rev.modality == "drug") & (rev.connectivity < 0)]
                  .drop_duplicates("name").head(30)["name"])
    ko_names = rev[(rev.modality == "crispr_ko")
                   & (rev.connectivity < 0)].head(15)["name"]
    keep, pool, seen = set(drug_names) | set(ko_names), [], set()
    for s in library:
        if s.name in keep and s.name not in seen:
            pool.append(s); seen.add(s.name)
    combos = best_combinations(ipf_bin, pool, k_genes=300,
                               require_cross_modality=True, top_n=8)
    print(f"\n=== Best drug + CRISPR-KO combinations (coverage of {label} genes) ===")
    print(combos[["pert_a", "pert_b", "coverage", "complementarity",
                  "combo_score"]].to_string(index=False,
                  float_format="{:.3f}".format))

    print("\nRIGOR: scores are deterministic; therapeutic/synergy claims are "
          "HYPOTHESES to validate experimentally.")
    return {"disease": ipf, "reversers": rev, "combinations": combos,
            "n_drugs": len(drugs_all), "n_crispr": len(crispr)}


if __name__ == "__main__":
    # usage: python3 -m perturbomics.realdata_ipf [cache_dir] [disease substring]
    cd = sys.argv[1] if len(sys.argv) > 1 else "enrichr_cache"
    dz = sys.argv[2] if len(sys.argv) > 2 else DISEASE
    main(cd, dz)
