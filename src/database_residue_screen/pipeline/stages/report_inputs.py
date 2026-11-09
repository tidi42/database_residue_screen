"""Output layout of the acquisition and analysis stages of the `screen` pipeline, as read by the reporting stages.

Every entry names the directory (relative to the stage's attempt output) that holds a marker file; ``directory``
returns it after checking that the marker exists, so a changed layout fails with a precise message.
"""
from pathlib import Path

from ...contracts import ContractError, read_json

LAYOUT = {
    "panel": (".", "panel.json"),
    "controls": ("controls", "control_qc.json"),
    "design": ("search_panel", "panel.json"),
    "search": ("search", "complete.json"),
    "metadata": ("metadata", "metadata.json"),
    "taxonomy": ("analysis", "analysis.json"),
    "insdc": ("insdc", "insdc.json"),
    "positional": (".", "positional_sensitivity.tsv"),
    "record_cases": ("record_cases", "record_cases.tsv"),
    "census": ("census", "census.json"),
    "exact_scan": ("scan", "scan.json"),
    "scan_validation": ("scan_validation", "scan_concordance_by_query.tsv"),
    "cap_review": (".", "cap_review.tsv"),
    "target_limit": ("comparison", "target_limit.tsv"),
    "benchmark": (".", "benchmark.json"),
    "inference": ("precision", "precision.json"),
    "validation": (".", "validation.json"),
    "metrics": (".", "metrics.json"),
    "reports": (".", "reports.json"),
}


def directory(stage: str, output: Path) -> Path:
    """The directory of `stage`'s attempt output that holds its marker file."""
    relative, marker = LAYOUT[stage]
    found = Path(output) / relative
    if not (found / marker).is_file():
        raise ContractError(f"Stage {stage}: {marker} missing in {found}")
    return found


def file(stage: str, output: Path) -> Path:
    """The marker file of `stage` (used where a report reads a single table or record)."""
    return directory(stage, output) / LAYOUT[stage][1]


def precision_specification(output: Path) -> Path:
    path = Path(output) / "specification" / "precision_specification.json"
    if not path.is_file():
        raise ContractError(f"Stage inference: precision_specification.json missing in {path.parent}")
    return path


def panel_definition(output: Path) -> Path:
    """The panel definition named in panel.json, checked against its recorded SHA-256."""
    from ...reporting.common import sha256
    record = read_json(file("panel", output))
    definition = Path(record["definition"])
    if not definition.is_file() or sha256(definition) != record["definition_sha256"]:
        raise ContractError("The panel definition named in panel.json is missing or changed")
    return definition


def benchmark_directories(output: Path) -> dict[str, Path]:
    """Benchmark runs by database-size label (runs/<label>, each holding benchmark.json)."""
    sizes = read_json(file("benchmark", output))["sizes"]
    result = {}
    for label in sizes:
        run = Path(output) / "runs" / label
        if not (run / "benchmark.json").is_file():
            raise ContractError(f"Stage benchmark: benchmark.json missing for size {label}")
        result[label] = run
    if not result:
        raise ContractError("Stage benchmark: no database size was run")
    return result
