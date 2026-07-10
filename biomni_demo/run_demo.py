"""
Biomni (snap-stanford/biomni) demo script.

Demonstrates two things that work WITHOUT any LLM API key:
  1. Initializing the A1 agent (skipping the ~11 GB S3 data-lake download).
  2. Calling standalone biomedical tool functions that hit public web APIs
     directly (PubMed, arXiv, AlphaFold, RCSB PDB).

The full autonomous agent loop (`agent.go(...)`) additionally needs an LLM key,
e.g. `export ANTHROPIC_API_KEY=sk-ant-...`. See README.md.

Run:  python run_demo.py
"""

import textwrap


def hr(title):
    print("\n" + "=" * 70)
    print(title)
    print("=" * 70)


def demo_agent_init():
    hr("1. Initialize the A1 agent (no data-lake download, no retriever)")
    from biomni.agent import A1

    # expected_data_lake_files=[] -> skip the ~11 GB S3 download.
    # use_tool_retriever=False   -> skip the embedding-model download.
    agent = A1(
        path="/tmp/biomni_data",
        llm="claude-sonnet-4-5",
        use_tool_retriever=False,
        expected_data_lake_files=[],
    )
    total = sum(len(v) for v in agent.module2api.values())
    print(f"\n-> Agent ready. {len(agent.module2api)} tool modules, {total} tool APIs registered.")
    return agent


def demo_pubmed():
    hr("2. Literature: query_pubmed  (live NCBI E-utilities, no key)")
    from biomni.tool.literature import query_pubmed

    out = query_pubmed("BRCA1 breast cancer targeted therapy", max_papers=2)
    print(textwrap.shorten(out.replace("\n", " "), width=700))


def demo_arxiv():
    hr("3. Literature: query_arxiv  (live arXiv API, no key)")
    from biomni.tool.literature import query_arxiv

    out = query_arxiv("protein language model drug discovery", max_papers=2)
    print(textwrap.shorten(out.replace("\n", " "), width=700))


def demo_alphafold():
    hr("4. Database: query_alphafold  (live EBI AlphaFold API, no key)")
    from biomni.tool.database import query_alphafold

    r = query_alphafold("P38398", endpoint="prediction")  # P38398 = BRCA1
    entry = r["result"][0]
    for k in ("modelEntityId", "modelCreatedDate", "globalMetricValue", "latestVersion"):
        print(f"  {k}: {entry.get(k)}")


def demo_pdb():
    hr("5. Database: query_pdb_identifiers  (live RCSB PDB API, no key)")
    from biomni.tool.database import query_pdb_identifiers

    r = query_pdb_identifiers(["1TUP"], return_type="entry")  # p53-DNA complex
    d = r["detailed_results"][0]["data"]
    cit = d["citation"][0]
    print(f"  Structure 1TUP: {cit['journal_abbrev']} {cit['journal_volume']}:{cit['page_first']} "
          f"(PubMed {cit['pdbx_database_id_PubMed']})")


if __name__ == "__main__":
    demo_agent_init()
    demo_pubmed()
    demo_arxiv()
    demo_alphafold()
    demo_pdb()
    hr("All no-key demos completed.")
