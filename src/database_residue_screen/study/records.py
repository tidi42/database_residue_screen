"""Run-directory records: attempts, output listings, verification, write protection and request differences."""
import json
import os
import re
import stat
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping

from .common import StudyError, file_sha256, json_sha256, listing_sha256, parse_utc, read_json, utc_now
from .engine import OUTPUT

STAGE_RECORD, REQUEST_RECORD, STAGE_LOG = "stage.json", "request.json", "stage.log"
RUN_RECORD, RUN_LOG, RUN_LOCK = "run.json", "run.log", "run.lock"
ATTEMPT_NAME = re.compile(r"^attempt_([1-9][0-9]*)$")
ATTEMPT_REFERENCE = re.compile(r"^[0-9]{2,}_[a-z][a-z0-9_]*/attempt_[1-9][0-9]*$")
STAGE_DIRECTORY = re.compile(r"^([0-9]{2,})_([a-z][a-z0-9_]*)$")
_IDENTIFIER = re.compile(r"^[A-Za-z_][A-Za-z0-9_-]*$")


@dataclass(frozen=True, eq=False)
class Attempt:
    number: int
    directory: Path
    record: dict | None
    started: dict | None

    @property
    def status(self) -> str:
        return self.record["status"] if self.record else "incomplete"

    @property
    def request(self) -> dict | None:
        source = self.record or self.started
        return source.get("request") if source else None

    @property
    def request_digest(self) -> str | None:
        source = self.record or self.started
        return source.get("request_digest") if source else None

    @property
    def code(self) -> dict | None:
        source = self.record or self.started
        return source.get("code") if source else None

    @property
    def mode(self) -> str | None:
        source = self.record or self.started
        return source.get("mode") if source else None

    def checked_request(self) -> dict | None:
        """The recorded request, after checking that it still matches its recorded digest."""
        request = self.request
        if request is not None and json_sha256(request) != self.request_digest:
            raise StudyError(f"{self.directory}: recorded request does not match its digest; the attempt record was modified")
        return request


def scan_attempts(stage_directory: Path) -> list[Attempt]:
    attempts = []
    if stage_directory.is_dir():
        for child in stage_directory.iterdir():
            match = ATTEMPT_NAME.match(child.name)
            if match and child.is_dir() and not child.is_symlink():
                record = read_json(child / STAGE_RECORD) if (child / STAGE_RECORD).is_file() else None
                started = read_json(child / REQUEST_RECORD) if (child / REQUEST_RECORD).is_file() else None
                attempts.append(Attempt(int(match.group(1)), child, record, started))
    return sorted(attempts, key=lambda attempt: attempt.number)


def current_attempt(attempts: list[Attempt]) -> Attempt | None:
    """The attempt in force: the latest completed one."""
    return max((attempt for attempt in attempts if attempt.status == "completed"), key=lambda attempt: attempt.number, default=None)


def latest_requested(attempts: list[Attempt]) -> Attempt | None:
    """The latest attempt that recorded a request (any status)."""
    return max((attempt for attempt in attempts if attempt.request is not None), key=lambda attempt: attempt.number, default=None)


def describe_attempt(attempt: Attempt, running: bool = False) -> dict:
    record = attempt.record or {}
    started = (attempt.record or attempt.started or {}).get("started_utc")
    status, duration = attempt.status, record.get("duration_seconds")
    if status == "incomplete" and running:
        status = "running"
        duration = round((utc_now() - parse_utc(started)).total_seconds(), 3) if started else None
    return dict(attempt=attempt.number, status=status, mode=attempt.mode, request_digest=attempt.request_digest,
                code_closure_sha256=((attempt.request or {}).get("code") or {}).get("closure_sha256"),
                stage_json_sha256=file_sha256(attempt.directory / STAGE_RECORD) if attempt.record else None,
                started_utc=started, ended_utc=record.get("ended_utc"), duration_seconds=duration)


def summarise_attempts(directory: str, attempts: list[Attempt], running: bool = False) -> dict:
    """Run-record entry of one stage: status of the latest attempt, the current attempt and every attempt."""
    current = current_attempt(attempts)
    listed = [describe_attempt(attempt, running and index == len(attempts) - 1) for index, attempt in enumerate(attempts)]
    return dict(directory=directory, status=listed[-1]["status"] if listed else "not_started",
                current_attempt=current.number if current else None, latest_attempt=attempts[-1].number if attempts else None,
                attempts=listed)


def _entries(output: Path) -> dict[str, Path]:
    found = {}
    for root, directories, files in os.walk(output):
        names = list(files) + [name for name in directories if os.path.islink(os.path.join(root, name))]
        directories[:] = sorted(name for name in directories if not os.path.islink(os.path.join(root, name)))
        for name in names:
            path = Path(root, name)
            found[path.relative_to(output).as_posix()] = path
    return found


def list_outputs(output: Path, strict: bool = True) -> list[dict]:
    """Every output file with relative path, size and SHA-256. Strict listings reject links and special files."""
    entries = []
    for relative, path in sorted(_entries(output).items()):
        info = os.lstat(path)
        if stat.S_ISREG(info.st_mode):
            try:
                entries.append(dict(path=relative, size=info.st_size, sha256=file_sha256(path)))
                continue
            except OSError as error:
                if strict:
                    raise
                entries.append(dict(path=relative, size=info.st_size, error=str(error)))
                continue
        kind = "symbolic link" if stat.S_ISLNK(info.st_mode) else "special file"
        if strict:
            raise StudyError(f"Stage outputs must be regular files; {relative} is a {kind}")
        entries.append(dict(path=relative, size=info.st_size, kind=kind))
    return entries


def outputs_digest(outputs: list[dict]) -> str:
    return listing_sha256((entry["path"], entry.get("sha256", "-")) for entry in outputs)


def verify_attempt(directory: Path, record: Mapping, full: bool) -> list[str]:
    """Differences between a completed attempt on disk and its stage.json (sizes always, content if ``full``)."""
    expected = {entry["path"]: entry for entry in record["outputs"]}
    actual = _entries(directory / OUTPUT)
    problems = [f"missing output {name}" for name in sorted(set(expected) - set(actual))]
    problems += [f"unexpected output {name}" for name in sorted(set(actual) - set(expected))]
    checks = [(name, actual[name], expected[name]) for name in sorted(set(expected) & set(actual))]
    for key in ("log", "request_json"):
        entry = record.get(key)
        if entry:
            checks.append((entry["path"], directory / entry["path"], entry))
    for name, path, entry in checks:
        try:
            info = os.lstat(path)
        except FileNotFoundError:
            problems.append(f"missing {name}")
            continue
        if not stat.S_ISREG(info.st_mode) or info.st_size != entry["size"]:
            problems.append(f"size of {name} changed")
        elif full and file_sha256(path) != entry["sha256"]:
            problems.append(f"content of {name} changed")
    return problems


def make_read_only(directory: Path) -> None:
    """Remove write permission from a finished attempt (files first, then directories, bottom-up)."""
    for root, directories, files in os.walk(directory, topdown=False):
        for name in files:
            path = os.path.join(root, name)
            if not os.path.islink(path):
                os.chmod(path, stat.S_IMODE(os.lstat(path).st_mode) & ~0o222)
        if not os.path.islink(root):
            os.chmod(root, stat.S_IMODE(os.lstat(root).st_mode) & ~0o222)


def input_digest(path: Path) -> dict:
    """SHA-256 and size of an input file, or a digest over the sorted listing of an input directory."""
    if path.is_file():
        return dict(size=path.stat().st_size, sha256=file_sha256(path))
    if path.is_dir():
        entries, size = [], 0
        for root, directories, files in os.walk(path):
            directories.sort()
            for name in files:
                file = Path(root, name)
                entries.append((file.relative_to(path).as_posix(), file_sha256(file)))
                size += file.stat().st_size
        return dict(size=size, files=len(entries), sha256=listing_sha256(entries))
    raise StudyError(f"Input does not exist: {path}")


def _join(prefix: str, key) -> str:
    key = str(key)
    if _IDENTIFIER.match(key):
        return f"{prefix}.{key}" if prefix else key
    return f"{prefix}[{json.dumps(key)}]"


def changed_fields(old, new, prefix: str = "") -> list[str]:
    """Dotted paths of the leaves that differ between two JSON values."""
    if isinstance(old, dict) and isinstance(new, dict):
        result = []
        for key in sorted(set(old) | set(new), key=str):
            path = _join(prefix, key)
            if key not in old:
                result.append(f"{path} (added)")
            elif key not in new:
                result.append(f"{path} (removed)")
            else:
                result.extend(changed_fields(old[key], new[key], path))
        return result
    return [] if old == new else [prefix or "(value)"]


def changed_keys(old: Mapping | None, new: Mapping | None) -> list[str]:
    """Keys whose values differ between two flat mappings (e.g. closure listings)."""
    old, new = old or {}, new or {}
    return sorted(key for key in set(old) | set(new) if old.get(key) != new.get(key))
