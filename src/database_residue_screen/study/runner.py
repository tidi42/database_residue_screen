"""Execution of a study: per-stage requests, rerun decisions, immutable attempts and the run record.

Run directory layout::

    RUN_DIR/run.json                            run record, rewritten atomically after every stage
    RUN_DIR/run.log                             run-level log (UTC)
    RUN_DIR/run.lock                            advisory lock held by the running invocation
    RUN_DIR/NN_stage/attempt_K/request.json     request, written when the attempt starts
    RUN_DIR/NN_stage/attempt_K/stage.log        attempt log (UTC)
    RUN_DIR/NN_stage/attempt_K/output/          the only directory the stage writes
    RUN_DIR/NN_stage/attempt_K/stage.json       final record; the attempt is then made read-only

A stage's current attempt is its latest completed attempt.
"""
import copy
import json
import logging
import os
import time
import traceback
from dataclasses import replace
from pathlib import Path
from typing import Iterable, Mapping

from .common import StudyError, canonical_json, file_sha256, json_sha256, read_json, utc_now, utc_text, write_json_atomic
from .config import StudyConfig, read_config_file, redactor, resolve_config
from .engine import OUTPUT, Pipeline, Stage, StageContext, UpstreamAttempt, load_pipeline
from .identity import CodeIdentity, collect_code_identity, git_ignored
from .logs import PACKAGE_LOGGER, RunLock, attempt_logging, run_logging
from .records import (ATTEMPT_REFERENCE, REQUEST_RECORD, RUN_LOCK, RUN_LOG, RUN_RECORD, STAGE_DIRECTORY, STAGE_LOG, STAGE_RECORD, Attempt,
                      changed_fields, changed_keys, current_attempt, input_digest, latest_requested, list_outputs, make_read_only,
                      outputs_digest, scan_attempts, summarise_attempts, verify_attempt)

SCHEMA_VERSION = "2.0.0"
CHANGED_ADVICE = "Use --rerun-changed to run a new attempt (earlier attempts are kept), or use a new run directory."


def run_study(config_path: Path, output: Path, *, until: str | None = None, only: Iterable[str] = (), mode: str = "production",
              rerun_changed: bool = False, environ: Mapping[str, str] | None = None, package_dir: Path | None = None,
              console=None) -> dict:
    """Run or resume the selected stages of the configured pipeline in run directory ``output``.

    A stage whose current attempt has an identical request is skipped; a stage whose latest attempt failed or is
    incomplete with an identical request gets a new attempt. A changed request stops the invocation before any stage
    executes, naming the changed fields, unless ``rerun_changed`` is set: then the stage gets a new attempt and its
    selected descendants rerun because their upstream digests change."""
    environ = os.environ if environ is None else environ
    raw = read_config_file(config_path)
    pipeline = load_pipeline(raw.pipeline)
    stages = pipeline.ordered()
    config = resolve_config(raw, pipeline.sections(), environ=environ, secrets=pipeline.secrets)
    for stage in stages:
        missing = [tool for tool in stage.tools if tool not in config.values["tools"]]
        if missing:
            raise StudyError(f"Stage {stage.name!r} uses tool(s) {', '.join(missing)}; configure them in [tools]")
    only = (only,) if isinstance(only, str) else tuple(only)
    execute, verify = pipeline.select(until, only)
    identity = collect_code_identity(mode, config.values["tools"], package_dir=package_dir, pipeline=pipeline, environ=environ)
    directory = Path(output).resolve()
    _check_location(directory, identity, environ)
    directory.mkdir(parents=True, exist_ok=True)
    selection = dict(until=until, only=list(only), rerun_changed=bool(rerun_changed))
    with RunLock(directory / RUN_LOCK), run_logging(directory / RUN_LOG, config.secret_values, console) as logger:
        try:
            return _Run(directory, pipeline, stages, config, identity, logger, bool(rerun_changed)).invoke(execute, verify, selection)
        except BaseException as error:
            logger.error("invocation stopped: %s", str(error) or type(error).__name__)
            raise


def _check_location(directory: Path, identity: CodeIdentity, environ: Mapping[str, str]) -> None:
    if directory.exists() and not directory.is_dir():
        raise StudyError(f"Run directory path exists and is not a directory: {directory}")
    if directory.is_dir() and not (directory / RUN_RECORD).exists():
        foreign = sorted(child.name for child in directory.iterdir() if child.name not in {RUN_LOCK, RUN_LOG})
        if foreign:
            raise StudyError(f"{directory} is not empty and has no {RUN_RECORD} ({', '.join(foreign[:5])}); use a new or empty run directory")
    root = identity.git_root
    if root is not None and directory.is_relative_to(root) and not git_ignored(root, directory, environ):
        raise StudyError(f"Run directory {directory} is inside the git working tree {root} but not ignored by git; its outputs would "
                         "make the tree dirty. Use a git-ignored location (such as runs/) or a directory outside the repository.")


def _summary(result, redact):
    if result is None:
        return None
    if not isinstance(result, Mapping):
        raise StudyError(f"A stage returns None or a JSON-serialisable mapping, not {type(result).__name__}")
    return _redacted(json.loads(canonical_json(dict(result))), redact)


def _redacted(value, redact):
    if isinstance(value, str):
        return redact(value)
    if isinstance(value, list):
        return [_redacted(item, redact) for item in value]
    if isinstance(value, dict):
        return {redact(key): _redacted(item, redact) for key, item in value.items()}
    return value


def _current_digest(path: str) -> dict | None:
    try:
        return input_digest(Path(path))
    except (StudyError, OSError):
        return None


def _file_changed(path: str, sha256: str) -> bool:
    try:
        return file_sha256(Path(path)) != sha256
    except OSError:
        return True


class _Run:
    def __init__(self, directory: Path, pipeline: Pipeline, stages: list[Stage], config: StudyConfig, identity: CodeIdentity,
                 logger: logging.Logger, rerun_changed: bool):
        self.directory, self.pipeline, self.stages, self.config, self.identity, self.logger = (
            directory, pipeline, stages, config, identity, logger)
        self.rerun_changed = rerun_changed
        self.redact = redactor(config.secret_values)
        self.directories: dict[str, str] = {}
        self.orphans: list[str] = []
        self.code: dict[str, tuple[dict, dict]] = {}
        self.document = None
        self.current = None
        self.failed = None
        self.verified = set()

    def attempts(self, stage: Stage) -> list[Attempt]:
        return scan_attempts(self.directory / self.directories[stage.name])

    def invoke(self, execute: set[str], verify: set[str], selection: dict) -> dict:
        previous = self.load_record()
        self.assign_directories(previous)
        self.logger.info("invocation in %s: study %s, pipeline %s, %s mode, config %s, git %s, selection %s", self.directory,
                         self.config.study_id, self.pipeline.name, self.identity.mode, self.config.sha256[:12],
                         (self.identity.recorded["git"].get("commit") or "unavailable")[:12], json.dumps(selection))
        if self.orphans:
            self.logger.warning("directories of stages no longer in pipeline %s are kept untouched: %s", self.pipeline.name,
                                ", ".join(self.orphans))
        self.code = {stage.name: self.identity.stage_code(stage.run, stage.tools) for stage in self.stages}
        scope = execute | verify
        resolved, requests = self.preflight(scope, verify)
        incomplete = [stage.name for stage in self.stages if stage.name in verify and stage.name not in resolved]
        if incomplete:
            raise StudyError(f"--only {', '.join(selection['only'])} needs completed dependencies; not complete: {', '.join(incomplete)}. "
                             "Run them first (for example with --until) or include them in --only.")
        self.begin(previous, selection)
        actions = {}
        try:
            for stage in self.stages:
                if stage.name not in scope or stage.name in verify:
                    actions[stage.name] = "verified" if stage.name in verify else "not_selected"
                    continue
                if stage.name not in resolved:
                    request = requests.get(stage.name) or self.request(stage, resolved)
                    decision, reference, fields = self.decide(stage, self.attempts(stage), request)
                    if decision != "skip":
                        if decision == "changed" and not self.rerun_changed:
                            raise StudyError(self.changed_message(stage, reference, fields))
                        if decision == "changed":
                            self.logger.info("stage %s: request changed from %s/attempt_%d (%s); running a new attempt", stage.name,
                                             self.directories[stage.name], reference.number, "; ".join(fields))
                        elif decision == "retry":
                            self.logger.info("stage %s: retrying the request of %s/attempt_%d (%s)", stage.name,
                                             self.directories[stage.name], reference.number, reference.status)
                        for dependency in stage.depends_on:
                            self.verify(resolved[dependency], full=True)
                        resolved[stage.name] = self.execute(stage, request, resolved)
                        actions[stage.name] = "executed"
                        continue
                    resolved[stage.name] = self.upstream_attempt(stage, reference)
                actions[stage.name] = "skipped"
                self.logger.info("stage %s: skipped; current attempt %s has an identical request", stage.name, resolved[stage.name].relative)
        except BaseException as error:
            self.end("failed", actions, error)
            raise
        self.end("completed", actions)
        self.logger.info("invocation completed: executed %s; skipped %s",
                         ", ".join(name for name, action in actions.items() if action == "executed") or "none",
                         ", ".join(name for name, action in actions.items() if action == "skipped") or "none")
        return dict(run_directory=str(self.directory), study_id=self.config.study_id, mode=self.identity.mode,
                    config_sha256=self.config.sha256, rerun_changed=self.rerun_changed,
                    stages=[dict(stage=name, action=actions[name], status=entry["status"], current_attempt=entry["current_attempt"],
                                 request_digest=next((item["request_digest"] for item in entry["attempts"]
                                                      if item["attempt"] == entry["current_attempt"]), None))
                            for name, entry in self.document["stages"].items()])

    def load_record(self) -> dict | None:
        path = self.directory / RUN_RECORD
        if not path.exists():
            return None
        document = read_json(path)
        if document.get("schema_version") != SCHEMA_VERSION:
            raise StudyError(f"{path}: unsupported run record version {document.get('schema_version')!r} (expected {SCHEMA_VERSION})")
        if document["study_id"] != self.config.study_id:
            raise StudyError(f"Run directory {self.directory} belongs to study {document['study_id']!r}, not {self.config.study_id!r}")
        if document["pipeline"]["name"] != self.pipeline.name:
            raise StudyError(f"Run directory {self.directory} belongs to pipeline {document['pipeline']['name']!r}, not {self.pipeline.name!r}")
        return document

    def assign_directories(self, previous: dict | None) -> None:
        """Keep every stage's directory once assigned; stages new to the pipeline get the next free number."""
        recorded = {entry["name"]: entry["directory"] for entry in (previous or {}).get("pipeline", {}).get("stages", [])}
        existing: dict[str, list[str]] = {}
        for child in sorted(self.directory.iterdir()):
            if child.is_dir():
                match = STAGE_DIRECTORY.match(child.name)
                if match is None:
                    raise StudyError(f"Unexpected directory {child} in the run directory")
                existing.setdefault(match.group(2), []).append(child.name)
        for name, directories in existing.items():
            if len(directories) > 1 or (name in recorded and directories[0] != recorded[name]):
                raise StudyError(f"Run directory layout disagrees with {RUN_RECORD} for stage {name!r}: {', '.join(directories)}")
        assigned = {**{name: directories[0] for name, directories in existing.items()}, **recorded}
        following = max((int(name.split("_", 1)[0]) for name in assigned.values()), default=0)
        for position, stage in enumerate(self.stages, 1):
            if stage.name not in assigned:
                following += 1
                assigned[stage.name] = f"{position if not existing and not recorded else following:02d}_{stage.name}"
            self.directories[stage.name] = assigned[stage.name]
        self.orphans = sorted(names[0] for name, names in existing.items() if name not in self.directories)

    def preflight(self, scope: set[str], verify: set[str]) -> tuple[dict, dict]:
        """Decide every stage whose request is computable now; nothing executes before all checks pass. Without
        --rerun-changed every stage with attempts is checked, selected or not."""
        resolved, requests = {}, {}
        for stage in self.stages:
            attempts = self.attempts(stage)
            if stage.name not in scope and (self.rerun_changed or not attempts):
                continue
            self.check_references(stage, attempts)
            pending = [name for name in stage.depends_on if name not in resolved]
            if pending:
                if not self.rerun_changed and latest_requested(attempts) is not None:
                    raise StudyError(f"Stage {stage.name!r}: upstream stage(s) {', '.join(pending)} have no current attempt matching their "
                                     f"request, so this stage's request will change (changed fields: "
                                     f"{', '.join(f'upstream.{name}' for name in pending)}). {CHANGED_ADVICE}")
                continue
            requests[stage.name] = request = self.request(stage, resolved)
            decision, reference, fields = self.decide(stage, attempts, request)
            if decision == "skip":
                resolved[stage.name] = self.upstream_attempt(stage, reference)
                self.verify(resolved[stage.name], full=False)
            elif decision == "changed" and (not self.rerun_changed or stage.name in verify):
                raise StudyError(self.changed_message(stage, reference, fields, stage.name in verify))
        return resolved, requests

    def decide(self, stage: Stage, attempts: list[Attempt], request: dict) -> tuple[str, Attempt | None, list[str]]:
        """skip (current attempt identical), retry (latest failed or incomplete attempt identical), new or changed."""
        digest = json_sha256(request)
        for attempt in attempts:
            attempt.checked_request()
        current, latest = current_attempt(attempts), latest_requested(attempts)
        if current is not None and current.request_digest == digest and (self.identity.mode == "working" or current.mode == "production"):
            return "skip", current, []
        if latest is not None and latest.status != "completed" and latest.request_digest == digest:
            return "retry", latest, []
        reference = current or latest
        if reference is None:
            return "new", None, []
        if reference.request_digest == digest:
            return "changed", reference, ["mode (production mode does not reuse a working-mode attempt)"]
        fields = changed_fields(reference.request, request)
        if "code.closure_sha256" in fields:
            modules = changed_keys((reference.code or {}).get("closure"), self.code[stage.name][1]["closure"])
            if modules:
                fields[fields.index("code.closure_sha256")] = f"code.closure_sha256 (changed modules: {', '.join(modules)})"
        return "changed", reference, fields

    def changed_message(self, stage: Stage, reference: Attempt, fields: list[str], required: bool = False) -> str:
        which = "current" if reference.status == "completed" else "latest"
        subject = f"Stage {stage.name!r}" + (" (required by --only)" if required else "")
        advice = "Select it as well (or use --until) together with --rerun-changed to rerun it." if required else CHANGED_ADVICE
        return (f"{subject}: the request differs from its {which} attempt {self.directories[stage.name]}/attempt_{reference.number} "
                f"({reference.status}); changed fields: {', '.join(fields)}. {advice}")

    def stage_config(self, stage: Stage) -> dict:
        return {name: self.config.section(name) for name in stage.config_sections}

    def request(self, stage: Stage, resolved: Mapping[str, UpstreamAttempt]) -> dict:
        upstream = {name: dict(attempt=resolved[name].relative, stage_json_sha256=resolved[name].stage_json_sha256)
                    for name in stage.depends_on}
        return dict(stage=stage.name, config={name: self.config.recorded_section(name) for name in stage.config_sections},
                    upstream=upstream, inputs=self.inputs(stage), code=self.code[stage.name][0])

    def inputs(self, stage: Stage) -> dict:
        if stage.input_files is None:
            return {}
        try:
            items = list(stage.input_files(self.stage_config(stage)))
        except StudyError:
            raise
        except Exception as error:
            raise StudyError(f"Stage {stage.name!r}: input_files failed: {error}") from error
        result = {}
        for item in items:
            path = Path(self.config.directory, item).resolve()
            try:
                result[str(path)] = input_digest(path)
            except StudyError as error:
                raise StudyError(f"Stage {stage.name!r}: {error}") from error
        return dict(sorted(result.items()))

    def check_references(self, stage: Stage, attempts: list[Attempt]) -> None:
        """Every recorded attempt must still find the exact upstream attempts it was computed from."""
        for attempt in attempts:
            for name, reference in ((attempt.checked_request() or {}).get("upstream") or {}).items():
                target = str(reference.get("attempt", ""))
                path = self.directory / target / STAGE_RECORD
                current = file_sha256(path) if ATTEMPT_REFERENCE.match(target) and path.is_file() else None
                if current != reference.get("stage_json_sha256"):
                    state = "is missing" if current is None else "has changed"
                    raise StudyError(f"Stage {stage.name!r}: attempt_{attempt.number} used upstream {target}, whose stage.json {state}; "
                                     f"changed fields: upstream.{name}.stage_json_sha256. The run directory was modified; attempts "
                                     "are immutable. Use a new run directory.")

    def upstream_attempt(self, stage: Stage, attempt: Attempt) -> UpstreamAttempt:
        return UpstreamAttempt(stage.name, attempt.number, attempt.directory, f"{self.directories[stage.name]}/attempt_{attempt.number}",
                               file_sha256(attempt.directory / STAGE_RECORD), copy.deepcopy(attempt.record))

    def verify(self, attempt: UpstreamAttempt, full: bool) -> None:
        if full and attempt.directory in self.verified:
            return
        problems = verify_attempt(attempt.directory, attempt.record, full)
        if problems:
            raise StudyError(f"Completed attempt {attempt.relative} was modified after completion ({'; '.join(problems[:10])}). "
                             "Completed attempts are immutable; use a new run directory.")
        if full:
            self.verified.add(attempt.directory)

    def execute(self, stage: Stage, request: dict, resolved: Mapping[str, UpstreamAttempt]) -> UpstreamAttempt:
        bound, code = self.code[stage.name]
        fresh = self.identity.stage_code(stage.run, stage.tools, fresh=True)[1]
        if fresh["closure"] != code["closure"]:
            raise StudyError(f"Stage {stage.name!r}: its code changed on disk during this invocation "
                             f"({', '.join(changed_keys(code['closure'], fresh['closure']))}); start a new invocation")
        number = max((attempt.number for attempt in self.attempts(stage)), default=0) + 1
        relative = f"{self.directories[stage.name]}/attempt_{number}"
        directory = self.directory / relative
        directory.parent.mkdir(exist_ok=True)
        directory.mkdir()
        (directory / OUTPUT).mkdir()
        digest, started, clock = json_sha256(request), utc_now(), time.monotonic()
        base = dict(schema_version=SCHEMA_VERSION, study_id=self.config.study_id, pipeline=self.pipeline.name, stage=stage.name,
                    attempt=number, directory=relative, mode=self.identity.mode, started_utc=utc_text(started))
        write_json_atomic(directory / REQUEST_RECORD, dict(base, pid=os.getpid(), request=request, request_digest=digest, code=code,
                                                          code_identity=self.identity.recorded))
        self.current = stage.name
        self.save()
        error = summary = None
        tools = {name: bound["tools"][name]["path"] for name in stage.tools}
        with attempt_logging(directory / STAGE_LOG, self.config.secret_values):
            self.logger.info("stage %s: attempt %d started in %s (request %s, code %s)", stage.name, number, relative, digest[:12],
                             bound["closure_sha256"][:12])
            log = logging.getLogger(f"{PACKAGE_LOGGER}.study.stage.{stage.name}")
            context = StageContext(name=stage.name, study_id=self.config.study_id, title=self.config.title, attempt=number,
                                   directory=directory, config=self.stage_config(stage),
                                   upstream={name: replace(resolved[name], record=copy.deepcopy(resolved[name].record))
                                             for name in stage.depends_on},
                                   tools=tools, mode=self.identity.mode, run_directory=self.directory, log=log)
            try:
                summary = _summary(stage.run(context), self.redact)
                changed = [path for path, entry in request["inputs"].items() if _current_digest(path) != entry]
                if changed:
                    raise StudyError("input changed while the stage was running: " + ", ".join(changed))
                after = self.identity.stage_code(stage.run, stage.tools, fresh=True)[1]
                if after["closure"] != code["closure"]:
                    raise StudyError("stage code changed while the stage was running: " + ", ".join(changed_keys(code["closure"], after["closure"])))
                changed = [name for name in stage.tools if _file_changed(bound["tools"][name]["resolved_path"], bound["tools"][name]["sha256"])]
                if changed:
                    raise StudyError("tool changed while the stage was running: " + ", ".join(changed))
                outputs = list_outputs(directory / OUTPUT)
            except BaseException as caught:
                error = caught
                log.error("stage failed", exc_info=caught)
                outputs = list_outputs(directory / OUTPUT, strict=False)
            ended, duration = utc_now(), round(time.monotonic() - clock, 6)
            self.logger.info("stage %s: attempt %d %s after %.3f s", stage.name, number, "failed" if error else "completed", duration)
        log_path, request_path = directory / STAGE_LOG, directory / REQUEST_RECORD
        record = dict(base, status="failed" if error else "completed", ended_utc=utc_text(ended), duration_seconds=duration,
                      request=request, request_digest=digest, code=code, code_identity=self.identity.recorded, outputs=outputs,
                      outputs_sha256=outputs_digest(outputs),
                      log=dict(path=STAGE_LOG, size=log_path.stat().st_size, sha256=file_sha256(log_path)),
                      request_json=dict(path=REQUEST_RECORD, size=request_path.stat().st_size, sha256=file_sha256(request_path)),
                      summary=summary, error=None if error is None else self.error_record(error), host=self.identity.host)
        write_json_atomic(directory / STAGE_RECORD, record)
        make_read_only(directory)
        self.current = None
        self.save()
        if error is not None:
            self.failed = stage.name
            if not isinstance(error, Exception):
                raise error
            raise StudyError(f"Stage {stage.name!r} failed in attempt {number}: {record['error']['message']}. Log: {log_path}") from error
        self.verified.add(directory)
        return UpstreamAttempt(stage.name, number, directory, relative, file_sha256(directory / STAGE_RECORD), record)

    def error_record(self, error: BaseException) -> dict:
        return dict(type=type(error).__name__, message=self.redact(str(error)) or type(error).__name__,
                    traceback=self.redact("".join(traceback.format_exception(error))))

    def begin(self, previous: dict | None, selection: dict) -> None:
        now, recorded = utc_text(utc_now()), self.identity.recorded
        invocation = dict(started_utc=now, ended_utc=None, status="running", pid=os.getpid(), mode=self.identity.mode,
                          selection=selection, config_sha256=self.config.sha256, git_commit=recorded["git"].get("commit"),
                          git_clean=recorded["git"].get("clean"), package_modules_sha256=recorded["package"]["modules_sha256"],
                          executed=[], skipped=[], failed=None, error=None)
        stages = [dict(name=stage.name, directory=self.directories[stage.name], depends_on=list(stage.depends_on),
                       config_sections=list(stage.config_sections), tools=list(stage.tools), description=stage.description,
                       module=self.code[stage.name][1]["module"]) for stage in self.stages]
        self.document = dict(
            schema_version=SCHEMA_VERSION, study_id=self.config.study_id, title=self.config.title,
            created_utc=previous["created_utc"] if previous else now, updated_utc=now, mode=self.identity.mode,
            pipeline=dict(name=self.pipeline.name, module=self.pipeline.module, description=self.pipeline.description, stages=stages,
                          orphan_directories=self.orphans),
            config_file=str(self.config.path), config_file_sha256=self.config.file_sha256, config_sha256=self.config.sha256,
            config=self.config.recorded, environment=self.config.environment,
            code_identity=dict(recorded, **self.identity.environment, tools=self.identity.tools), code_files=self.identity.files,
            host=self.identity.host, stages={}, invocations=(previous["invocations"] if previous else []) + [invocation])
        self.save()

    def save(self) -> None:
        self.document["updated_utc"] = utc_text(utc_now())
        self.document["stages"] = {stage.name: summarise_attempts(self.directories[stage.name], self.attempts(stage),
                                                                  stage.name == self.current) for stage in self.stages}
        write_json_atomic(self.directory / RUN_RECORD, self.document)

    def end(self, status: str, actions: dict, error: BaseException | None = None) -> None:
        if self.document is None:
            return
        self.document["invocations"][-1].update(
            ended_utc=utc_text(utc_now()), status=status, failed=self.failed,
            executed=[name for name, action in actions.items() if action == "executed"],
            skipped=[name for name, action in actions.items() if action == "skipped"],
            error=None if error is None else (self.redact(str(error)) or type(error).__name__))
        self.current = None
        self.save()
