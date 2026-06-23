# Small-Molecule Annotation Pipeline

A minimum-viable cheminformatics pipeline that annotates a list of small molecules
through eight sequential steps, from structure curation to ADMET risk.

## Pipeline steps

| # | Module | Sources |
|---|--------|---------|
| 1 | Structure standardization | RDKit (MolStandardize, PAINS/Brenk/NIH alerts) |
| 2 | Experimental targets | ChEMBL · BindingDB · PubChem BioAssay · GtoPdb |
| 3 | Known drug mechanisms | DrugCentral · Repurposing Hub |
| 4 | Predicted targets | SwissTargetPrediction · DrugCLIP (optional) |
| 5 | Pathways & disease context | Reactome · Open Targets |
| 6 | Chemical probe quality | Chemical Probes Portal |
| 7 | Functional signatures | LINCS / iLINCS · CMap (optional) |
| 8 | ADMET & toxicity | ADMETlab 2.0 · CompTox (EPA) · ProTox-3.0 |

## Quick start

```bash
# 1. Install dependencies
pip install -r small_molecule_pipeline/requirements.txt

# 2. Run with the bundled example molecules
python small_molecule_pipeline/pipeline.py \
  --input  small_molecule_pipeline/data/input/example_molecules.smi \
  --output-dir small_molecule_pipeline/data/output

# 3. Run only specific steps
python small_molecule_pipeline/pipeline.py \
  --input molecules.smi --steps 1,2,3

# 4. Use a YAML config
cp small_molecule_pipeline/config_example.yaml config.yaml
# edit config.yaml …
python small_molecule_pipeline/pipeline.py --config config.yaml
```

## Input format

One SMILES per line (lines starting with `#` are ignored).
Also accepts `.csv` / `.tsv` with SMILES in the first column.

```
CC(=O)Oc1ccccc1C(=O)O
Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1
```

## Output

Each step writes two files to the output directory:

```
data/output/
├── 01_rdkit_standardization.json / .csv
├── 02_experimental_targets.json / .csv
├── 03_drug_mechanisms.json / .csv
├── 04_predicted_targets.json / .csv
├── 05_pathways_diseases.json / .csv
├── 06_chemical_probes.json / .csv
├── 07_lincs_cmap.json / .csv
├── 08_admet_toxicity.json / .csv
└── pipeline_summary.json
```

## Optional features

| Feature | How to enable |
|---------|--------------|
| DrugCLIP target prediction | `export DRUGCLIP_MODEL_PATH=/path/to/ckpt` |
| CMap / Clue.io queries | `export CLUE_API_KEY=your_key` |
| CompTox authenticated access | `export COMPTOX_API_KEY=your_key` |
| Local Repurposing Hub TSV | `export REPURPOSING_HUB_PATH=/path/to/file.tsv` |
| Run specific steps only | `export PIPELINE_STEPS=1,2,3` |

## Structure

```
small_molecule_pipeline/
├── pipeline.py          # Orchestrator / CLI
├── config.py            # Configuration dataclass
├── requirements.txt
├── config_example.yaml
├── modules/
│   ├── step1_rdkit.py
│   ├── step2_experimental_targets.py
│   ├── step3_drug_mechanisms.py
│   ├── step4_prediction.py
│   ├── step5_pathways.py
│   ├── step6_probes.py
│   ├── step7_lincs.py
│   └── step8_admet.py
├── utils/
│   ├── api_client.py    # Rate-limited HTTP with retry
│   └── io_utils.py      # File I/O helpers
└── data/
    ├── input/example_molecules.smi
    └── output/          # Generated at runtime
```
