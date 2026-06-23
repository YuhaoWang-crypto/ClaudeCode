"""Pipeline configuration — override via environment variables or config.yaml."""

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional


@dataclass
class PipelineConfig:
    # I/O
    input_path: str = "data/input/molecules.smi"
    output_dir: str = "data/output"

    # Steps to run (all enabled by default)
    steps: List[int] = field(default_factory=lambda: [1, 2, 3, 4, 5, 6, 7, 8])

    # Step 2 — experimental targets
    chembl_max_activities: int = 100

    # Step 4 — prediction
    stp_organism: str = "Homo sapiens"
    drugclip_model_path: str = ""  # empty = disabled
    drugclip_target_db: str = ""

    # Step 7 — LINCS/CMap
    clue_api_key: str = ""  # empty = CMap disabled

    # Step 8 — ADMET
    comptox_api_key: str = ""  # empty = unauthenticated (rate-limited)

    # Misc
    log_level: str = "INFO"
    max_workers: int = 1  # sequential by default (public API rate limits)

    @classmethod
    def from_env(cls) -> "PipelineConfig":
        cfg = cls()
        cfg.input_path = os.environ.get("INPUT_PATH", cfg.input_path)
        cfg.output_dir = os.environ.get("OUTPUT_DIR", cfg.output_dir)
        cfg.drugclip_model_path = os.environ.get("DRUGCLIP_MODEL_PATH", "")
        cfg.drugclip_target_db = os.environ.get("DRUGCLIP_TARGET_DB", "")
        cfg.clue_api_key = os.environ.get("CLUE_API_KEY", "")
        cfg.comptox_api_key = os.environ.get("COMPTOX_API_KEY", "")
        cfg.log_level = os.environ.get("LOG_LEVEL", cfg.log_level)
        steps_env = os.environ.get("PIPELINE_STEPS", "")
        if steps_env:
            cfg.steps = [int(s.strip()) for s in steps_env.split(",") if s.strip()]
        return cfg

    @classmethod
    def from_yaml(cls, path: str) -> "PipelineConfig":
        import yaml  # optional dependency

        with open(path) as f:
            raw = yaml.safe_load(f)
        cfg = cls()
        for k, v in raw.items():
            if hasattr(cfg, k):
                setattr(cfg, k, v)
        # Env vars take precedence over YAML
        env = cls.from_env()
        for k in ["drugclip_model_path", "drugclip_target_db", "clue_api_key", "comptox_api_key"]:
            if getattr(env, k):
                setattr(cfg, k, getattr(env, k))
        return cfg
