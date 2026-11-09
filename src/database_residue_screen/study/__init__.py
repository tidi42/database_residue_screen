"""Config-driven study orchestration with immutable attempts and per-stage provenance (see docs/study-engine.md)."""
from .common import REDACTED, StudyError
from .config import StudyConfig, load_config
from .engine import (PIPELINE_REFERENCES, Pipeline, Stage, StageContext, UpstreamAttempt, available_pipelines, load_pipeline,
                     register_pipeline, register_pipeline_reference)
from .identity import CodeIdentity, CodeIndex, collect_code_identity
from .runner import run_study
from .status import format_status, study_status, verify_stage_code

__all__ = [
    "PIPELINE_REFERENCES", "REDACTED", "CodeIdentity", "CodeIndex", "Pipeline", "Stage", "StageContext", "StudyConfig", "StudyError",
    "UpstreamAttempt", "available_pipelines", "collect_code_identity", "format_status", "load_config", "load_pipeline",
    "register_pipeline", "register_pipeline_reference", "run_study", "study_status", "verify_stage_code",
]
