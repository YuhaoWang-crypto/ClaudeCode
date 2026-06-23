"""Small-molecule annotation pipeline — main orchestrator.

Usage
-----
  python pipeline.py --input molecules.smi --output-dir results/
  python pipeline.py --input molecules.smi --steps 1,2,3   # run subset
  python pipeline.py --input molecules.smi --config cfg.yaml
"""

import argparse
import json
import logging
import sys
import time
from pathlib import Path

# Ensure package root is on sys.path when run directly
sys.path.insert(0, str(Path(__file__).parent))

from config import PipelineConfig
from utils.io_utils import read_smiles_list, save_results, load_step_results
from modules import (
    step1_rdkit,
    step2_experimental_targets,
    step3_drug_mechanisms,
    step4_prediction,
    step5_pathways,
    step6_probes,
    step7_lincs,
    step8_admet,
)

STEP_NAMES = {
    1: "01_rdkit_standardization",
    2: "02_experimental_targets",
    3: "03_drug_mechanisms",
    4: "04_predicted_targets",
    5: "05_pathways_diseases",
    6: "06_chemical_probes",
    7: "07_lincs_cmap",
    8: "08_admet_toxicity",
}


def setup_logging(level: str) -> None:
    logging.basicConfig(
        format="%(asctime)s  %(levelname)-8s  %(name)s — %(message)s",
        datefmt="%H:%M:%S",
        level=getattr(logging, level.upper(), logging.INFO),
    )


def run_pipeline(cfg: PipelineConfig) -> None:
    logger = logging.getLogger("pipeline")
    logger.info("=" * 60)
    logger.info("Small-Molecule Annotation Pipeline")
    logger.info("Input  : %s", cfg.input_path)
    logger.info("Output : %s", cfg.output_dir)
    logger.info("Steps  : %s", cfg.steps)
    logger.info("=" * 60)

    smiles_list = read_smiles_list(cfg.input_path)
    if not smiles_list:
        logger.error("No molecules loaded — aborting.")
        sys.exit(1)

    t0 = time.time()

    # ── Step 1 ──────────────────────────────────────────────────────────────
    s1_results = []
    if 1 in cfg.steps:
        logger.info("\n▶ Step 1: RDKit standardization & structural alerts")
        s1_results = step1_rdkit.run(smiles_list)
        save_results(s1_results, cfg.output_dir, STEP_NAMES[1])
    else:
        s1_results = load_step_results(cfg.output_dir, STEP_NAMES[1])

    valid_mols = [r for r in s1_results if r.get("passed")]
    logger.info(
        "  → %d / %d molecules passed Step 1", len(valid_mols), len(s1_results)
    )

    # ── Step 2 ──────────────────────────────────────────────────────────────
    s2_results = []
    if 2 in cfg.steps:
        logger.info("\n▶ Step 2: Experimental targets (ChEMBL / BindingDB / PubChem / GtoPdb)")
        s2_results = step2_experimental_targets.run(s1_results)
        save_results(s2_results, cfg.output_dir, STEP_NAMES[2])
    else:
        s2_results = load_step_results(cfg.output_dir, STEP_NAMES[2])

    # ── Step 3 ──────────────────────────────────────────────────────────────
    s3_results = []
    if 3 in cfg.steps:
        logger.info("\n▶ Step 3: Known drug mechanisms (DrugCentral / Repurposing Hub)")
        s3_results = step3_drug_mechanisms.run(s1_results)
        save_results(s3_results, cfg.output_dir, STEP_NAMES[3])
    else:
        s3_results = load_step_results(cfg.output_dir, STEP_NAMES[3])

    # ── Step 4 ──────────────────────────────────────────────────────────────
    s4_results = []
    if 4 in cfg.steps:
        logger.info("\n▶ Step 4: Predicted targets (SwissTargetPrediction / DrugCLIP)")
        s4_results = step4_prediction.run(s1_results)
        save_results(s4_results, cfg.output_dir, STEP_NAMES[4])
    else:
        s4_results = load_step_results(cfg.output_dir, STEP_NAMES[4])

    # ── Step 5 ──────────────────────────────────────────────────────────────
    s5_results = []
    if 5 in cfg.steps:
        logger.info("\n▶ Step 5: Pathways & disease context (Reactome / Open Targets)")
        s5_results = step5_pathways.run(s1_results, s2_results, s3_results)
        save_results(s5_results, cfg.output_dir, STEP_NAMES[5])
    else:
        s5_results = load_step_results(cfg.output_dir, STEP_NAMES[5])

    # ── Step 6 ──────────────────────────────────────────────────────────────
    s6_results = []
    if 6 in cfg.steps:
        logger.info("\n▶ Step 6: Chemical Probes Portal quality check")
        s6_results = step6_probes.run(s1_results)
        save_results(s6_results, cfg.output_dir, STEP_NAMES[6])
    else:
        s6_results = load_step_results(cfg.output_dir, STEP_NAMES[6])

    # ── Step 7 ──────────────────────────────────────────────────────────────
    s7_results = []
    if 7 in cfg.steps:
        logger.info("\n▶ Step 7: LINCS / CMap transcriptional signatures")
        s7_results = step7_lincs.run(s1_results)
        save_results(s7_results, cfg.output_dir, STEP_NAMES[7])
    else:
        s7_results = load_step_results(cfg.output_dir, STEP_NAMES[7])

    # ── Step 8 ──────────────────────────────────────────────────────────────
    s8_results = []
    if 8 in cfg.steps:
        logger.info("\n▶ Step 8: ADMET / toxicity / interference risk")
        s8_results = step8_admet.run(s1_results)
        save_results(s8_results, cfg.output_dir, STEP_NAMES[8])
    else:
        s8_results = load_step_results(cfg.output_dir, STEP_NAMES[8])

    # ── Summary report ───────────────────────────────────────────────────────
    elapsed = time.time() - t0
    _write_summary(
        cfg, s1_results, s2_results, s3_results, s4_results,
        s5_results, s6_results, s7_results, s8_results, elapsed,
    )
    logger.info("\n✓ Pipeline complete in %.1f s", elapsed)


def _write_summary(cfg, s1, s2, s3, s4, s5, s6, s7, s8, elapsed):
    logger = logging.getLogger("pipeline.summary")

    summary = {
        "input_file": cfg.input_path,
        "n_input": len(s1),
        "n_valid": sum(1 for r in s1 if r.get("passed")),
        "n_duplicates": sum(1 for r in s1 if r.get("error") == "duplicate"),
        "n_invalid_smiles": sum(1 for r in s1 if r.get("error") == "invalid_smiles"),
        "n_with_alerts": sum(
            1 for r in s1 if r.get("passed") and r.get("n_alerts", 0) > 0
        ),
        "n_with_experimental_targets": sum(
            1 for r in s2 if r.get("n_targets", 0) > 0
        ),
        "n_known_drugs": sum(1 for r in s3 if r.get("is_known_drug")),
        "n_with_predictions": sum(
            1 for r in s4 if r.get("predicted_targets")
        ),
        "n_probes": sum(1 for r in s6 if r.get("is_probe")),
        "n_with_lincs": sum(1 for r in s7 if r.get("has_lincs_data")),
        "risk_distribution": {
            "high": sum(1 for r in s8 if r.get("risk_level") == "high"),
            "medium": sum(1 for r in s8 if r.get("risk_level") == "medium"),
            "low": sum(1 for r in s8 if r.get("risk_level") == "low"),
        },
        "elapsed_seconds": round(elapsed, 1),
    }

    out = Path(cfg.output_dir) / "pipeline_summary.json"
    with open(out, "w") as f:
        json.dump(summary, f, indent=2)

    logger.info("\n─── SUMMARY ───────────────────────────────────────")
    for k, v in summary.items():
        if k != "elapsed_seconds":
            logger.info("  %-35s %s", k, v)
    logger.info("  %-35s %.1f s", "elapsed", elapsed)
    logger.info("─────────────────────────────────────────────────")
    logger.info("Full summary written to %s", out)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Small-molecule annotation pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--input", "-i", help="Input SMILES file (.smi, .csv, .tsv)")
    p.add_argument("--output-dir", "-o", help="Output directory (default: data/output)")
    p.add_argument(
        "--steps",
        help="Comma-separated step numbers to run (e.g. 1,2,3). Default: all",
    )
    p.add_argument("--config", help="YAML config file")
    p.add_argument("--log-level", default="INFO", help="Logging level")
    return p.parse_args()


def main() -> None:
    args = parse_args()

    if args.config:
        cfg = PipelineConfig.from_yaml(args.config)
    else:
        cfg = PipelineConfig.from_env()

    if args.input:
        cfg.input_path = args.input
    if args.output_dir:
        cfg.output_dir = args.output_dir
    if args.steps:
        cfg.steps = [int(s.strip()) for s in args.steps.split(",")]
    if args.log_level:
        cfg.log_level = args.log_level

    setup_logging(cfg.log_level)
    run_pipeline(cfg)


if __name__ == "__main__":
    main()
