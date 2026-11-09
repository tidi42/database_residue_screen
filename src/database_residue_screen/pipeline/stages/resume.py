"""Continuation of long stages after a failed or interrupted attempt with an identical request.

The engine starts a new attempt after an unfinished one with an identical request. Long stages copy the completed,
self-verifying units of such earlier attempts (search batches, scanned volumes) into the new attempt and continue
there; the copies are recorded and the earlier attempts are never modified.
"""
import json
import re
import shutil
from pathlib import Path

ATTEMPT = re.compile(r"^attempt_([1-9][0-9]*)$")


def _json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def earlier_attempts(context) -> list[Path]:
    """Output directories of the earlier unfinished (failed or interrupted) attempts of this stage whose request
    digest equals this attempt's, latest first."""
    digest = _json(context.directory / "request.json")["request_digest"]
    found = []
    for sibling in context.directory.parent.iterdir():
        match = ATTEMPT.match(sibling.name)
        if match is None or sibling == context.directory or not sibling.is_dir():
            continue
        if (sibling / "stage.json").is_file() and _json(sibling / "stage.json").get("status") == "completed":
            continue
        if (sibling / "request.json").is_file() and _json(sibling / "request.json").get("request_digest") == digest:
            found.append((int(match.group(1)), sibling / "output"))
    return [path for _, path in sorted(found, reverse=True) if path.is_dir()]


def copy_file(source: Path, target: Path) -> None:
    """Copy content only (earlier attempts may be read-only)."""
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, target)


def copy_directory(source: Path, target: Path) -> int:
    count = 0
    for path in sorted(source.rglob("*")):
        if path.is_file() and not path.is_symlink():
            copy_file(path, target / path.relative_to(source))
            count += 1
    return count


def sources_with(context, relative: str, marker: str) -> list[Path]:
    """`relative` directories of the earlier unfinished attempts that contain `marker` and the same marker content."""
    candidates = [earlier / relative for earlier in earlier_attempts(context) if (earlier / relative / marker).is_file()]
    if not candidates:
        return []
    reference = (candidates[0] / marker).read_bytes()
    return [path for path in candidates if (path / marker).read_bytes() == reference]
