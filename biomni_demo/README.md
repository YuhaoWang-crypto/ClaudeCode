# Biomni (snap-stanford/biomni) — install & demo

Notes from installing [`snap-stanford/biomni`](https://github.com/snap-stanford/biomni)
(the Stanford SNAP general-purpose biomedical AI agent) and running working demos
in this environment.

- **Package:** `biomni` 0.0.8 (PyPI)
- **Python:** 3.11
- **Result:** Installed ✅, A1 agent initializes with **218 tool APIs across 21
  domains** ✅, standalone biomedical tools run against live public APIs ✅.

## Install

The published wheel under-declares its dependencies (metadata lists only
`pydantic`, `langchain`, `python-dotenv`), so the runtime deps have to be added
explicitly:

```bash
python3 -m venv venv && source venv/bin/activate
pip install biomni \
    langgraph langchain-anthropic langchain-openai langchain-community \
    pandas numpy scikit-learn tqdm requests beautifulsoup4 \
    biopython pubchempy PyPDF2 langchain-text-splitters \
    googlesearch-python pymed arxiv
```

## Run the demo

```bash
python run_demo.py
```

### What it does (no LLM key required)

1. **Initializes the `A1` agent** while skipping the heavy downloads:
   - `expected_data_lake_files=[]` skips the **~11 GB** S3 data-lake download.
   - `use_tool_retriever=False` skips the embedding-model download.
2. Calls standalone tool functions that hit **public web APIs directly**:
   - `query_pubmed` — NCBI E-utilities (BRCA1 breast-cancer literature)
   - `query_arxiv` — arXiv API (protein-LM drug-discovery papers)
   - `query_alphafold` — EBI AlphaFold API (BRCA1 / UniProt P38398 structure)
   - `query_pdb_identifiers` — RCSB PDB API (1TUP, the p53–DNA complex)

Sample output:

```
1. Initialize the A1 agent ...
-> Agent ready. 21 tool modules, 218 tool APIs registered.
4. Database: query_alphafold  (live EBI AlphaFold API, no key)
  modelEntityId: AF-P38398-F1
  globalMetricValue: 41.59
5. Database: query_pdb_identifiers  (live RCSB PDB API, no key)
  Structure 1TUP: Science 265:346 (PubMed 8023157)
```

## Running the full autonomous agent (needs an LLM key)

The 218-tool autonomous reasoning loop (`agent.go(...)`) and the natural-language
`query_*` tools (e.g. `query_uniprot(prompt=...)`, which use an LLM to translate a
prompt into an API call) require an LLM API key. In this sandbox Claude Code
authenticates via OAuth, so no plain `ANTHROPIC_API_KEY` is exposed — those paths
report an auth error. To run them, supply a key:

```bash
export ANTHROPIC_API_KEY=sk-ant-...      # or OPENAI_API_KEY / GEMINI_API_KEY
```

```python
from biomni.agent import A1
agent = A1(path="./data", llm="claude-sonnet-4-5")   # first real run downloads the data lake
agent.go("Plan a CRISPR screen to study T-cell exhaustion and suggest 32 genes to perturb.")
```

Default LLM is `claude-sonnet-4-5` (configurable via `BiomniConfig` / the `A1(llm=..., source=...)`
constructor; supported sources: OpenAI, AzureOpenAI, Anthropic, Ollama, Gemini, Bedrock, Groq, Custom).

## Tool domains registered by A1

biochemistry, bioengineering, bioimaging, biophysics, cancer_biology, cell_biology,
database (40 tools), genetics, genomics (19), glycoengineering, immunology,
literature, microbiology, molecular_biology (18), pathology, pharmacology (23),
physiology, support_tools, synthetic_biology, systems_biology, lab_automation.
