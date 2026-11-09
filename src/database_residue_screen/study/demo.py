"""Built-in demonstration pipeline: four small stages on a text file of integers (no BLAST, network or R).

This module only registers stages; their code lives in :mod:`demo_stages`, so editing the registration (or the engine)
does not change any stage's code closure.

    from database_residue_screen.study.demo import write_demo_study
    write_demo_study("runs/demo-input")
    # residue-screen run-study runs/demo-input/study.toml --output runs/demo --working
"""
from pathlib import Path

from . import demo_stages
from .engine import Pipeline

PIPELINE = Pipeline("demo", "Ingest integers, scale them, describe them and write a short report")
PIPELINE.section("numbers", {
    "type": "object",
    "required": ["source"],
    "properties": {"source": {"type": "string", "format": "path", "description": "Text file with one integer per line"}},
})
PIPELINE.section("scaling", {"type": "object", "properties": {"factor": {"type": "integer", "minimum": 1, "default": 1}}})
PIPELINE.section("report", {"type": "object", "properties": {"heading": {"type": "string", "default": "Demonstration report"}}})
PIPELINE.stage(config_sections=("numbers",), input_files=lambda config: [config["numbers"]["source"]])(demo_stages.ingest)
PIPELINE.stage(depends_on=("ingest",), config_sections=("scaling",))(demo_stages.scale)
PIPELINE.stage(depends_on=("ingest",))(demo_stages.describe)
PIPELINE.stage(depends_on=("scale", "describe"), config_sections=("report",))(demo_stages.report)


def write_demo_study(directory: Path, *, study_id: str = "demo", numbers=(3, 1, 4, 1, 5, 9, 2, 6), factor: int = 2) -> Path:
    """Write ``study.toml`` and ``numbers.txt`` for the demonstration pipeline into a new directory."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=False)
    (directory / "numbers.txt").write_text("".join(f"{number}\n" for number in numbers), encoding="utf-8")
    config = directory / "study.toml"
    config.write_text(f"""# Demonstration study of the study engine.
[study]
id = "{study_id}"
title = "Demonstration study"
pipeline = "demo"

[numbers]
source = "numbers.txt"   # resolved against this file's directory

[scaling]
factor = {factor}

[report]
heading = "Demonstration report"
""", encoding="utf-8")
    return config
