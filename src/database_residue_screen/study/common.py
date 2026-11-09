"""Shared helpers of the study engine: errors, canonical JSON, hashing, atomic writes and UTC timestamps."""
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path

from ..contracts import ContractError

REDACTED = "<redacted>"


class StudyError(ContractError):
    """A study configuration, run directory or provenance contract was violated."""


def canonical_json(value) -> bytes:
    """Sorted-key, whitespace-free UTF-8 JSON; NaN and infinity are rejected."""
    try:
        return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False).encode("utf-8")
    except (TypeError, ValueError) as error:
        raise StudyError(f"Value cannot be represented as canonical JSON: {error}") from error


def json_sha256(value) -> str:
    return hashlib.sha256(canonical_json(value)).hexdigest()


def file_sha256(path: Path) -> str:
    with open(path, "rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def listing_sha256(entries) -> str:
    """Digest of a sha256sum-style listing: one '<sha256>  <relative path>' line per file, sorted by path."""
    text = "".join(f"{digest}  {path}\n" for path, digest in sorted(entries))
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def utc_now() -> datetime:
    return datetime.now(timezone.utc)


def utc_text(moment: datetime) -> str:
    return moment.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%S.%fZ")


def parse_utc(text: str) -> datetime:
    return datetime.strptime(text, "%Y-%m-%dT%H:%M:%S.%fZ").replace(tzinfo=timezone.utc)


def read_json(path: Path) -> dict:
    def unique(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise StudyError(f"{path}: duplicate JSON key {key!r}")
            result[key] = value
        return result

    try:
        value = json.loads(Path(path).read_text(encoding="utf-8"), object_pairs_hook=unique)
    except json.JSONDecodeError as error:
        raise StudyError(f"{path}: invalid JSON: {error}") from error
    if not isinstance(value, dict):
        raise StudyError(f"{path}: expected a JSON object")
    return value


def write_json_atomic(path: Path, value) -> None:
    """Write JSON through a synced temporary file and an atomic rename."""
    text = json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        with open(temporary, "w", encoding="utf-8") as handle:
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()
    if os.name == "posix":
        descriptor = os.open(path.parent, os.O_RDONLY)
        try:
            os.fsync(descriptor)
        finally:
            os.close(descriptor)
