#!/usr/bin/env python3
"""
pathway-biomarker-triangulation
--------------------------------
Cross-validate a pathway's core biomarkers with TWO orthogonal, evidence-grounded
methods, then hand off to a THIRD (dynamical) method:

  Method 1  ATTRIBUTION  = measurability x (assoc_w*disease_assoc + cent_w*centrality)
                           -> greedy minimal cover >= coverage_thr (minimal essential unit)
  Method 2  RWR/DIFFUSION = personalized PageRank (restart r) from the seeds on the
                           STRING network (pure topology; no disease/measurability)
  Method 3  DYNAMICAL     = the user's `network-biomarker` skill (grn_pipeline):
                           irreducible core (automorphism/fibration) + CRNT delta +
                           DNB / critical-slowing early-warning biomarker. Documented
                           hand-off (see SKILL.md); not executed here.

Convergence across Method 1 (evidence) and Method 2 (topology) = the most defensible
core; Method 3 adds switch-capacity + a dynamical tipping-point observable.

Public live deps: STRING API, Open Targets Platform GraphQL. Local: networkx, matplotlib.
Everything is best-effort with graceful fallbacks (a target missing from Open Targets
gets assoc=0 with a warning -- exactly the INHBE case that must be rescued by other
evidence, see SKILL.md).

Usage:
  python triangulate.py --seeds GPR75,INHBE,ACVR1C,ACVR2B,MSTN,GDF15,GFRAL \
      --disease EFO_0001073 --out out/ [--topk 13] [--string-score 700]
"""
import argparse, json, sys, time, urllib.request, urllib.parse, pathlib
import networkx as nx

STRING = "https://string-db.org/api/tsv/network"
OT = "https://api.platform.opentargets.org/api/v4/graphql"

# ---- measurability priors: analyte class -> weight (edit per project) ----
MEAS_DEFAULT = 0.5
MEAS = {  # secreted/circulating ligands & genotype = 1.0; receptors/occupancy 0.7;
    # phospho-node (needs biopsy) 0.8; intracellular-only 0.3
}
SECRETED = {"INHBE","INHBA","INHBB","MSTN","GDF8","GDF11","GDF15","BMP8A","BMP8B",
            "FST","FSTL3","LEP","ADIPOQ","NODAL","INHA","INHBC"}
RECEPTOR = {"ACVR1C","ACVR2B","ACVR2A","ACVR1B","ACVR1","TGFBR1","GFRAL","RET",
            "GPR75","BMPR2","TGFBR2"}
PHOSPHO  = {"SMAD2","SMAD3","SMAD4","SMAD7"}

def meas_weight(g):
    if g in MEAS: return MEAS[g]
    if g in SECRETED: return 1.0
    if g in PHOSPHO:  return 0.8
    if g in RECEPTOR: return 0.7
    return 0.3  # intracellular-only default

def http_get(url, tries=2, timeout=20):
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as r: return r.read().decode()
        except Exception:
            if i==tries-1: raise
            time.sleep(1.5*(i+1))

def http_post_json(url, payload, tries=2, timeout=20):
    data=json.dumps(payload).encode()
    for i in range(tries):
        try:
            req=urllib.request.Request(url, data=data, headers={"Content-Type":"application/json"})
            with urllib.request.urlopen(req, timeout=timeout) as r: return json.loads(r.read().decode())
        except Exception:
            if i==tries-1: raise
            time.sleep(1.5*(i+1))

# ---------- STRING network (seeds + first shell) ----------
def string_network(seeds, species=9606, score=700, add_nodes=30):
    ids="%0d".join(seeds)
    url=f"{STRING}?identifiers={ids}&species={species}&required_score={score}&add_nodes={add_nodes}"
    txt=http_get(url); G=nx.Graph()
    for ln in txt.strip().splitlines()[1:]:
        f=ln.split("\t")
        if len(f)>=6:
            a,b,sc=f[2],f[3],float(f[5]); G.add_edge(a,b,weight=sc)
    return G

# ---------- Open Targets disease association (disease-centric, one call) ----------
def ot_disease_targets(disease_id, size=3000):
    """Return {approvedSymbol: association_score} for a disease (EFO/MONDO id).
    One GraphQL call; OT rows are score-desc so the top `size` covers the relevant targets.
    Missing targets (e.g. INHBE for obesity) simply aren't in the map -> assoc 0 -> rescue
    via other evidence (see SKILL.md)."""
    q={"query":"query($d:String!,$n:Int!){disease(efoId:$d){id name "
             "associatedTargets(page:{index:0,size:$n}){count rows{score target{approvedSymbol}}}}}",
       "variables":{"d":disease_id,"n":size}}
    d=http_post_json(OT,q,timeout=40)
    dis=(d.get("data") or {}).get("disease")
    if not dis:
        raise RuntimeError(f"Open Targets: no disease node for '{disease_id}' "
                           f"(use a valid MONDO/EFO id, e.g. obesity = MONDO_0011122)")
    at=dis["associatedTargets"]
    return {r["target"]["approvedSymbol"]: float(r["score"]) for r in at["rows"]}, dis["name"], at["count"]

# ---------- Method 1: attribution + greedy minimal cover ----------
def attribution(nodes, G, assoc, assoc_w=0.65, cent_w=0.35):
    deg=dict(G.degree());
    bet=nx.betweenness_centrality(G, weight="weight") if G.number_of_edges() else {}
    dmax=max(deg.values()) if deg else 1
    out={}
    for g in nodes:
        c = 0.6*bet.get(g,0.0) + 0.4*(deg.get(g,0)/dmax if dmax else 0.0)
        a = assoc.get(g,0.0)
        out[g]=dict(assoc=a, centrality=c, meas=meas_weight(g),
                    attribution=meas_weight(g)*(assoc_w*a + cent_w*c))
    return out

def greedy_cover(attr, thr=0.80):
    tot=sum(v["attribution"] for v in attr.values()) or 1.0
    ordered=sorted(attr, key=lambda g:-attr[g]["attribution"])
    core=[]; cum=0.0
    for g in ordered:
        core.append(g); cum+=attr[g]["attribution"]/tot
        if cum>=thr: break
    return core, cum, ordered

# ---------- Method 2: RWR (personalized PageRank) ----------
def rwr(G, seeds, r=0.3):
    seeds=[s for s in seeds if s in G]
    if not seeds: return {}
    pers={n:(1.0/len(seeds) if n in seeds else 0.0) for n in G}
    return nx.pagerank(G, alpha=1-r, personalization=pers, weight="weight")

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--seeds", required=True)
    ap.add_argument("--disease", default="EFO_0001073", help="EFO/MONDO id (default obesity)")
    ap.add_argument("--out", default="out")
    ap.add_argument("--topk", type=int, default=0, help="RWR top-k (0 = match minimal-set size)")
    ap.add_argument("--string-score", type=int, default=700)
    ap.add_argument("--coverage", type=float, default=0.80)
    ap.add_argument("--assoc-w", type=float, default=0.65)
    ap.add_argument("--rwr-r", type=float, default=0.30)
    a=ap.parse_args()
    seeds=[s.strip() for s in a.seeds.split(",") if s.strip()]
    outdir=pathlib.Path(a.out); outdir.mkdir(parents=True, exist_ok=True)

    print(f"[1/5] STRING network for {len(seeds)} seeds ...", flush=True)
    G=string_network(seeds, score=a.string_score)
    print(f"      nodes={G.number_of_nodes()} edges={G.number_of_edges()}")
    pool=sorted(set(G.nodes())|set(seeds))

    print(f"[2/5] Open Targets association to {a.disease} (one disease-centric query) ...", flush=True)
    dmap,dname,dcount=ot_disease_targets(a.disease)
    print(f"      disease='{dname}' associated_targets={dcount} (mapped {len(dmap)})", flush=True)
    assoc={g: dmap.get(g,0.0) for g in pool}
    for g in seeds:
        if assoc.get(g,0.0)==0.0:
            print(f"      WARN {g}: OT assoc=0 (not in obesity associations) -> rescue via other evidence "
                  f"(e.g. AlphaGenome splice / pLOF)", flush=True)

    print("[3/5] Method 1: attribution + greedy minimal cover ...", flush=True)
    attr=attribution(pool,G,assoc,assoc_w=a.assoc_w,cent_w=1-a.assoc_w)
    core,cov,ordered=greedy_cover(attr,thr=a.coverage)
    print(f"      minimal core: {len(core)}/{len(pool)} nodes, coverage={cov:.3f}")

    print("[4/5] Method 2: RWR diffusion ...", flush=True)
    rw=rwr(G,seeds,r=a.rwr_r)
    k=a.topk or len(core)
    rwr_top=[g for g,_ in sorted(rw.items(), key=lambda x:-x[1])[:k]]

    shared=[g for g in core if g in rwr_top]
    print(f"      RWR top-{k}; overlap with Method1 core = {len(shared)} "
          f"(Jaccard={len(shared)/len(set(core)|set(rwr_top)):.2f})")

    # write TSV
    tsv=outdir/"triangulation.tsv"
    with open(tsv,"w") as f:
        f.write("gene\tassoc\tcentrality\tmeas\tattribution\trwr\tin_m1_core\tin_rwr_top\n")
        for g in ordered:
            f.write(f"{g}\t{attr[g]['assoc']:.3f}\t{attr[g]['centrality']:.3f}\t{attr[g]['meas']:.2f}"
                    f"\t{attr[g]['attribution']:.3f}\t{rw.get(g,0):.4f}"
                    f"\t{int(g in core)}\t{int(g in rwr_top)}\n")

    # figure
    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        gs=ordered[:min(20,len(ordered))]; ys=range(len(gs))
        cols=["#0b6b5e" if g in core else "#c0c7d0" for g in gs]
        fig,ax=plt.subplots(1,2,figsize=(11,6))
        ax[0].barh(list(ys),[attr[g]["attribution"] for g in gs],color=cols)
        ax[0].set_yticks(list(ys)); ax[0].set_yticklabels(gs,fontsize=8); ax[0].invert_yaxis()
        ax[0].set_title(f"Method 1 attribution (green=minimal core {len(core)}/{len(pool)})")
        # cumulative
        tot=sum(attr[g]["attribution"] for g in ordered) or 1
        cum=[]; s=0
        for g in ordered: s+=attr[g]["attribution"]/tot; cum.append(s)
        ax[1].plot(range(1,len(cum)+1),cum,marker="o",ms=3); ax[1].axhline(a.coverage,ls="--",c="r")
        ax[1].set_title(f"Cumulative coverage; cutoff {a.coverage:.2f} -> {len(core)} nodes")
        ax[1].set_xlabel("nodes (sorted)"); ax[1].set_ylabel("cum. attribution")
        fig.suptitle("Biomarker triangulation: Method 1 (evidence) + Method 2 (RWR)")
        fig.tight_layout(); fig.savefig(outdir/"triangulation.png",dpi=140)
    except Exception as e:
        print("      figure skipped:",e)

    summary=dict(seeds=seeds, disease=a.disease,
                 pool=len(pool), m1_core=core, m1_coverage=round(cov,3),
                 rwr_top=rwr_top, shared_core=shared,
                 jaccard=round(len(shared)/len(set(core)|set(rwr_top)),3),
                 note="Method 3 (dynamical: irreducible core + CRNT delta + DNB) -> run the "
                      "network-biomarker skill (grn_pipeline) on this pathway; see SKILL.md")
    (outdir/"triangulation.json").write_text(json.dumps(summary,indent=2))
    print("[5/5] wrote:", tsv, outdir/'triangulation.png', outdir/'triangulation.json')
    print("\nSHARED ROBUST CORE (Method1 ∩ Method2):", ", ".join(shared))
    print("Method3 hand-off:", summary["note"])

if __name__=="__main__":
    main()
