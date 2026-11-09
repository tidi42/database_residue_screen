"""Shared helpers of the reporting layer: deterministic tables and JSON, verified inputs and output manifests.

Values are written at full precision (``repr`` of floats); rounding happens only at display. Manifests list input
and output SHA-256 digests but no paths or times, so identical inputs give byte-identical output directories.
"""
import csv
import hashlib
import json
import math
from pathlib import Path

from ..contracts import ContractError, read_json

NE = "not_estimated"
GROUPS = ("ONT", "Illumina", "BGI/DNBSEQ", "PacBio")
GLOBAL = "Global union"
ENDPOINTS = ("strict", "tolerant")
MEASURES = ("detection", "breadth")
EXCLUDED_STATES = ("censored_untestable", "chance_limited", "partially_score_censored")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 22), b""):
            digest.update(block)
    return digest.hexdigest()


def text(value, boolean=("yes", "no")) -> str:
    """Full-precision text: floats by repr, booleans as the given pair, None as not_estimated."""
    if value is None:
        return NE
    if isinstance(value, bool):
        return boolean[0] if value else boolean[1]
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ContractError(f"Non-finite value in a report table: {value}")
        return repr(value)
    return str(value)


def ratio(numerator, denominator):
    """numerator / denominator as a float, or None (not estimable) for an empty denominator."""
    return numerator / denominator if denominator else None


def rate(numerator, denominator) -> str:
    return text(ratio(numerator, denominator))


def mean(values):
    values = list(values)
    return sum(values) / len(values) if values else None


def write_tsv(path: Path, rows: list[dict], fields=None, *, boolean=("yes", "no")) -> None:
    """Write a new tab-separated table; every row must have exactly the declared columns."""
    rows = list(rows)
    if fields is None:
        if not rows:
            raise ContractError(f"{Path(path).name}: an empty table needs declared columns")
        fields = list(rows[0])
    fields = list(fields)
    for row in rows:
        if list(row) != fields and set(row) != set(fields):
            raise ContractError(f"{Path(path).name}: inconsistent columns {sorted(set(row) ^ set(fields))}")
    with Path(path).open("x", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(fields)
        for row in rows:
            writer.writerow([text(row[name], boolean) for name in fields])


def read_tsv(path: Path) -> list[dict]:
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_json(path: Path, value) -> None:
    Path(path).write_text(json.dumps(value, indent=1, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")


def verified_files(directory: Path, record: dict, key: str = "files") -> dict:
    """Check the SHA-256 of every file listed under ``record[key]`` (paths relative to ``directory``)."""
    for name, checksum in record.get(key, {}).items():
        path = Path(directory, name)
        if Path(name).is_absolute() or ".." in Path(name).parts or not path.is_file() or sha256(path) != checksum:
            raise ContractError(f"Changed or missing input artifact: {path}")
    return record


def load_analysis(directory: Path) -> dict:
    """analysis.json of a completed analysis, after checking every listed table."""
    return verified_files(directory, read_json(Path(directory) / "analysis.json"))


def new_directory(path: Path) -> Path:
    """Create an output directory; an existing one must be empty (a stage's prepared output directory)."""
    path = Path(path)
    if path.exists() and (not path.is_dir() or any(path.iterdir())):
        raise ContractError(f"Output directory already exists and is not empty: {path}")
    path.mkdir(parents=True, exist_ok=True)
    return path


def group_key(group: str):
    return (GROUPS.index(group), group) if group in GROUPS else (len(GROUPS) + (group == GLOBAL), group)


def ordered_groups(groups) -> list[str]:
    return sorted(set(groups), key=group_key)


def write_manifest(destination: Path, kind: str, inputs: dict, statement: str, **extra) -> dict:
    """manifest.json: kind, input and output digests and declared facts; deterministic (no paths or times)."""
    destination = Path(destination)
    outputs = {str(path.relative_to(destination)): sha256(path) for path in sorted(destination.rglob("*"))
               if path.is_file() and path.name != "manifest.json"}
    manifest = dict(kind=kind, statement=statement, rounding="Full precision in all tables; rounding only at display",
                    inputs={name: sha256(path) for name, path in sorted(inputs.items()) if path is not None},
                    outputs=outputs, **extra)
    write_json(destination / "manifest.json", manifest)
    return manifest


def verify_manifest(directory: Path, kind: str | None = None) -> dict:
    """Read a reporting manifest and check every output it lists."""
    manifest = read_json(Path(directory) / "manifest.json")
    if kind is not None and manifest.get("kind") != kind:
        raise ContractError(f"{directory} is not a {kind} output (kind {manifest.get('kind')!r})")
    return verified_files(directory, manifest, "outputs")
