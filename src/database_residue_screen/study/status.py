"""Status of a study run directory, with optional re-verification of outputs and of stage code."""
import os
from pathlib import Path
from typing import Mapping

from .common import StudyError, file_sha256, listing_sha256, read_json
from .identity import PACKAGE_DIRECTORY, CodeIndex, environment_identity, git_state, tool_identity
from .logs import lock_held
from .records import (RUN_LOCK, RUN_RECORD, STAGE_RECORD, changed_fields, changed_keys, current_attempt, scan_attempts,
                      summarise_attempts, verify_attempt)

PASSED = ("ok", "-")


def study_status(run_directory: Path, *, verify: bool = False, verify_code: bool = False, package_dir: Path | None = None,
                 environ: Mapping[str, str] | None = None) -> dict:
    """Per-stage status from the attempt records on disk. ``verify`` re-hashes the outputs of every current attempt;
    ``verify_code`` checks each current attempt against the current package source and environment."""
    directory = Path(run_directory).resolve()
    path = directory / RUN_RECORD
    if not path.is_file():
        raise StudyError(f"{directory} is not a study run directory (no {RUN_RECORD})")
    document = read_json(path)
    locked = lock_held(directory / RUN_LOCK)
    rows = []
    for stage in document["pipeline"]["stages"]:
        attempts = scan_attempts(directory / stage["directory"])
        running = locked and document["stages"].get(stage["name"], {}).get("status") == "running"
        row = dict(stage=stage["name"], **summarise_attempts(stage["directory"], attempts, running))
        current = current_attempt(attempts)
        if verify and current is not None:
            row["problems"] = verify_attempt(current.directory, current.record, full=True)
            row["verified"] = not row["problems"]
        rows.append(row)
    report = dict(run_directory=str(directory), study_id=document["study_id"], title=document["title"],
                  pipeline=document["pipeline"]["name"], mode=document["mode"], config_sha256=document["config_sha256"],
                  updated_utc=document["updated_utc"], running=locked, stages=rows)
    if verify_code:
        report["code"] = verify_stage_code(directory, document, package_dir=package_dir, environ=environ)
    return report


def verify_stage_code(directory: Path, document: dict, *, package_dir: Path | None = None,
                      environ: Mapping[str, str] | None = None) -> dict:
    """For every stage: was its current attempt produced by the stage code now in the package (same import closure),
    with the same tools, dependencies, interpreter and platform, and from the current attempts of its dependencies?"""
    environ = os.environ if environ is None else environ
    package_dir = Path(package_dir or PACKAGE_DIRECTORY).resolve()
    index, environment, git = CodeIndex(package_dir), environment_identity(package_dir), git_state(package_dir, environ)
    currents, references = {}, {}
    for stage in document["pipeline"]["stages"]:
        currents[stage["name"]] = current = current_attempt(scan_attempts(directory / stage["directory"]))
        if current is not None:
            references[stage["name"]] = (f"{stage['directory']}/attempt_{current.number}", file_sha256(current.directory / STAGE_RECORD))
    rows = []
    for stage in document["pipeline"]["stages"]:
        current = currents[stage["name"]]
        if current is None:
            rows.append(dict(stage=stage["name"], attempt=None, mode=None, code="NOT RUN", tools="-", dependencies="-", python="-",
                             platform="-", upstream="-", ok=False, details=["no completed attempt"]))
            continue
        request, recorded = current.record["request"], current.record.get("code") or {}
        bound, details = request.get("code"), []
        if bound is None or "module" not in recorded:
            rows.append(dict(stage=stage["name"], attempt=current.number, mode=current.mode, code="UNKNOWN", tools="-", dependencies="-",
                             python="-", platform="-", upstream="-", ok=False, details=["attempt record has no stage code block"]))
            continue
        try:
            closure = index.closure(recorded["module"])
        except StudyError as error:
            code = "UNKNOWN"
            details.append(str(error))
        else:
            code = "ok" if listing_sha256(closure.items()) == bound["closure_sha256"] else "CHANGED"
            if code == "CHANGED":
                details.append("changed code: " + ", ".join(changed_keys(recorded["closure"], closure)))
        tools = []
        for name, tool in sorted(bound["tools"].items()):
            try:
                now = tool_identity(name, tool, environ)
            except StudyError as error:
                tools.append(f"{name} ({error})")
                continue
            if now != tool:
                tools.append(f"{name} ({', '.join(changed_fields(tool, now))})")
        if tools:
            details.append("changed tools: " + "; ".join(tools))
        dependencies = changed_keys(bound["dependencies"], environment["dependencies"])
        if dependencies:
            details.append("changed dependencies: " + ", ".join(f"{name} {bound['dependencies'].get(name)} -> "
                                                                f"{environment['dependencies'].get(name)}" for name in dependencies))
        for key in ("python", "platform"):
            if bound[key] != environment[key]:
                details.append(f"{key}: {' '.join(map(str, bound[key].values()))} -> {' '.join(map(str, environment[key].values()))}")
        stale = [name for name, reference in request["upstream"].items()
                 if references.get(name) != (reference["attempt"], reference["stage_json_sha256"])]
        if stale:
            details.append("not computed from the current attempt of: " + ", ".join(stale))
        row = dict(stage=stage["name"], attempt=current.number, mode=current.mode, code=code,
                   tools="-" if not bound["tools"] else ("CHANGED" if tools else "ok"),
                   dependencies="CHANGED" if dependencies else "ok",
                   python="ok" if bound["python"] == environment["python"] else "CHANGED",
                   platform="ok" if bound["platform"] == environment["platform"] else "CHANGED",
                   upstream="-" if not request["upstream"] else ("STALE" if stale else "ok"), details=details)
        row["ok"] = all(row[key] in PASSED for key in ("code", "tools", "dependencies", "python", "platform", "upstream"))
        rows.append(row)
    return dict(package_directory=str(package_dir), git=dict(commit=git.get("commit"), clean=git.get("clean")) if git["available"]
                else dict(commit=None, clean=None), ok=all(row["ok"] for row in rows), stages=rows)


def format_duration(seconds: float | None) -> str:
    if seconds is None:
        return "-"
    if seconds < 60:
        return f"{seconds:.1f}s"
    minutes, second = divmod(int(round(seconds)), 60)
    if minutes < 60:
        return f"{minutes}m{second:02d}s"
    hours, minute = divmod(minutes, 60)
    if hours < 24:
        return f"{hours}h{minute:02d}m"
    days, hour = divmod(hours, 24)
    return f"{days}d{hour:02d}h"


def _table(rows: list[list[str]]) -> list[str]:
    widths = [max(len(row[index]) for row in rows) for index in range(len(rows[0]))]
    return ["  ".join(cell.ljust(width) for cell, width in zip(row, widths)).rstrip() for row in rows]


def format_status(report: dict) -> str:
    """Status table: STATUS and DURATION describe the latest attempt, CURRENT is the attempt in force."""
    verified = any("verified" in row for row in report["stages"])
    rows = [["STAGE", "STATUS", "CURRENT", "ATTEMPTS", "DURATION", "REQUEST"] + (["OUTPUTS"] if verified else [])]
    for row in report["stages"]:
        latest = row["attempts"][-1] if row["attempts"] else {}
        rows.append([row["stage"], row["status"], str(row["current_attempt"] or "-"), str(len(row["attempts"])),
                     format_duration(latest.get("duration_seconds")), (latest.get("request_digest") or "-")[:12]]
                    + ([("ok" if row["verified"] else "MODIFIED") if "verified" in row else "-"] if verified else []))
    lines = [f"Study {report['study_id']}: {report['title']}",
             f"Run {report['run_directory']} | pipeline {report['pipeline']} | last invocation in {report['mode']} mode | config "
             f"{report['config_sha256'][:12]} | updated {report['updated_utc']}" + (" | invocation running" if report["running"] else ""),
             *_table(rows)]
    lines += [f"{row['stage']}: {problem}" for row in report["stages"] for problem in row.get("problems", [])]
    if "code" in report:
        code = report["code"]
        git = f"git {code['git']['commit'][:12]}{'' if code['git']['clean'] else ' (dirty tree)'}" if code["git"]["commit"] else "no git"
        lines += ["", f"Code verification against {code['package_directory']} ({git}): {'all identical' if code['ok'] else 'MISMATCH'}"]
        table = [["STAGE", "ATTEMPT", "MODE", "CODE", "TOOLS", "DEPENDENCIES", "PYTHON", "PLATFORM", "UPSTREAM"]]
        table += [[row["stage"], str(row["attempt"] or "-"), row["mode"] or "-", row["code"], row["tools"], row["dependencies"],
                   row["python"], row["platform"], row["upstream"]] for row in code["stages"]]
        lines += _table(table) + [f"{row['stage']}: {detail}" for row in code["stages"] for detail in row["details"]]
    return "\n".join(lines)
