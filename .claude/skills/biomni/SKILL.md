---
name: biomni
description: >-
  Install and run snap-stanford/biomni — Stanford SNAP's general-purpose
  biomedical AI agent (the A1 agent + 218 tool APIs across 21 domains:
  genomics, pharmacology, molecular biology, literature, structural biology,
  clinical trials, etc.). Use when the user wants to set up biomni, initialize
  the A1 agent without the ~11 GB data-lake download, call a specific biomni
  tool (PubMed / arXiv / AlphaFold / PDB / ChEMBL / UniProt / Ensembl / KEGG
  lookups), run the autonomous agent loop (agent.go), or debug biomni install /
  dependency / API-key errors. Distinguishes tools that work with NO LLM key
  from those that need one.
---

# Biomni — install & run the Stanford biomedical agent

[`snap-stanford/biomni`](https://github.com/snap-stanford/biomni) is a
general-purpose biomedical AI agent. Its `A1` agent wires **218 tool APIs
across 21 domains** into an LLM reasoning loop (LangGraph). Many of those tools
are also plain Python functions you can call directly.

Verified working: **biomni 0.0.8**, Python 3.11.

## The two things to get right first

1. **The PyPI wheel under-declares its dependencies** (metadata lists only
   `pydantic`, `langchain`, `python-dotenv`). A bare `pip install biomni` will
   `ModuleNotFoundError` on first real use. Install the full set below.
2. **`A1()` downloads a ~11 GB S3 data lake on first init** unless you stop it.
   Pass `expected_data_lake_files=[]` to skip it for demos / tool-only use.

## Install

```bash
python3 -m venv venv && source venv/bin/activate
pip install biomni \
    langgraph langchain-anthropic langchain-openai langchain-community \
    pandas numpy scikit-learn tqdm requests beautifulsoup4 \
    biopython pubchempy PyPDF2 langchain-text-splitters \
    googlesearch-python pymed arxiv
```

Add more only when an import error names them (biomni lazy-imports per tool).

## Initialize the agent WITHOUT the big download

```python
from biomni.agent import A1
agent = A1(
    path="./biomni_data",
    llm="claude-sonnet-4-5",
    use_tool_retriever=False,      # skips embedding-model download
    expected_data_lake_files=[],   # skips the ~11 GB S3 data-lake download
)
# agent.module2api -> dict of 21 modules; ~218 tool APIs total
```

## Key rule: which tools need an LLM API key?

This is the #1 source of confusion. See `reference/tools.md` for the full map.

- **NO key needed** — tools that take *structured* args and hit a public API
  directly: `query_pubmed`, `query_arxiv`, `query_alphafold(uniprot_id, ...)`,
  `query_pdb_identifiers([...])`, `blast_sequence`, plus the MCP-style database
  endpoints. Great for demos and deterministic lookups.
- **Key REQUIRED** — (a) the autonomous loop `agent.go("...")`, and (b) the
  natural-language `query_*(prompt="...")` tools (e.g.
  `query_uniprot(prompt=...)`, `query_pubchem(prompt=...)`) that use an LLM to
  translate the prompt into an API call. Without a key these raise
  *"Could not resolve authentication method."*

```bash
export ANTHROPIC_API_KEY=sk-ant-...   # or OPENAI_API_KEY / GEMINI_API_KEY
```

> In a Claude Code remote sandbox there is usually **no plain `ANTHROPIC_API_KEY`**
> (auth is OAuth via file descriptor, not reusable by subprocesses). Demo the
> no-key tools; tell the user to export a key for the full loop.

## Run the bundled demo (no key)

```bash
python assets/run_demo.py
```

Exercises: A1 init (skipping download) → `query_pubmed` → `query_arxiv` →
`query_alphafold` (BRCA1 / P38398) → `query_pdb_identifiers` (1TUP, p53–DNA).
Expected tail:

```
-> Agent ready. 21 tool modules, 218 tool APIs registered.
  Structure 1TUP: Science 265:346 (PubMed 8023157)
```

## Run the full autonomous agent (needs a key)

```python
from biomni.agent import A1
agent = A1(path="./data", llm="claude-sonnet-4-5")   # first REAL run downloads the data lake
agent.go("Plan a CRISPR screen to study T-cell exhaustion; propose 32 genes to perturb.")
```

Default LLM `claude-sonnet-4-5`. Sources supported by `A1(source=...)`:
OpenAI, AzureOpenAI, Anthropic, Ollama, Gemini, Bedrock, Groq, Custom.
Config also via env: `BIOMNI_LLM`, `BIOMNI_PATH`, `BIOMNI_TIMEOUT_SECONDS`, etc.

## Troubleshooting

| Symptom | Cause / fix |
|---|---|
| `ModuleNotFoundError: Bio / PyPDF2 / googlesearch / pymed` | wheel under-declares deps — install the full set above |
| Init hangs / fills disk | the ~11 GB data-lake download — pass `expected_data_lake_files=[]` |
| `Could not resolve authentication method` | no LLM key — export one, or use a no-key tool instead |
| slow first `A1()` with retriever | embedding-model download — `use_tool_retriever=False` for tool-only use |
