"""Stage model and pipeline registry of the study engine.

A pipeline is a named set of stages plus the JSON schemas of the configuration sections they read::

    PIPELINE = Pipeline("paper3")
    PIPELINE.section("search", {"type": "object", "required": ["queries"],
                                "properties": {"queries": {"type": "string", "format": "path"},
                                               "threads": {"type": "integer", "minimum": 1, "default": 8}}})

    @PIPELINE.stage(depends_on=("panel",), config_sections=("search",), tools=("blastn",),
                    input_files=lambda config: [config["search"]["queries"]])
    def search(context):
        ...  # read context.config, context.upstream["panel"].output, context.tool("blastn"); write below context.output

A study configuration selects the pipeline with ``[study] pipeline = "<name>"`` (a registered name or a
``module:attribute`` reference); see :data:`PIPELINE_REFERENCES` and :func:`register_pipeline`.
"""
import copy
import importlib
import logging
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType
from typing import Callable, Iterable, Mapping

from jsonschema import Draft202012Validator
from jsonschema.exceptions import SchemaError

from .common import StudyError
from .config import BUILTIN_SECTIONS, ENVIRONMENT_NAME, SECTION_NAME

STAGE_NAME = re.compile(r"^[a-z][a-z0-9_]{0,63}$")
PIPELINE_NAME = re.compile(r"^[A-Za-z][A-Za-z0-9_.-]{0,63}$")
TOOL_NAME = re.compile(r"^[A-Za-z][A-Za-z0-9_.+-]{0,63}$")
OUTPUT = "output"

PIPELINE_REFERENCES: dict[str, str] = {"demo": "database_residue_screen.study.demo:PIPELINE",
                                       "screen": "database_residue_screen.pipeline.screen:PIPELINE"}
"""Pipeline names usable in ``[study] pipeline``, mapped to the ``module:attribute`` holding the :class:`Pipeline`.
Modules are imported only when a configuration names them."""
_REGISTERED: dict[str, "Pipeline"] = {}


def _names(value, what: str, pattern) -> tuple[str, ...]:
    names = (value,) if isinstance(value, str) else tuple(value)
    for name in names:
        if not isinstance(name, str) or not pattern.match(name):
            raise StudyError(f"Invalid {what}: {name!r}")
    if len(set(names)) != len(names):
        raise StudyError(f"Duplicate {what} in {names}")
    return names


@dataclass(frozen=True)
class Stage:
    """One step of a study. ``run(context)`` writes only below ``context.output`` and may return a small
    JSON-serialisable mapping recorded as the attempt summary. ``input_files(config)`` receives the declared
    configuration sections and returns the external files or directories (typically paths named in those sections)
    whose SHA-256 is bound into the request."""
    name: str
    run: Callable[["StageContext"], Mapping | None]
    depends_on: tuple[str, ...] = ()
    config_sections: tuple[str, ...] = ()
    input_files: Callable[[Mapping], Iterable] | None = None
    tools: tuple[str, ...] = ()
    description: str = ""

    def __post_init__(self):
        if not isinstance(self.name, str) or not STAGE_NAME.match(self.name):
            raise StudyError(f"Invalid stage name {self.name!r}: use lowercase letters, digits and underscores")
        if not callable(self.run) or (self.input_files is not None and not callable(self.input_files)):
            raise StudyError(f"Stage {self.name!r}: run and input_files must be callables")
        object.__setattr__(self, "depends_on", _names(self.depends_on, "dependency", STAGE_NAME))
        object.__setattr__(self, "config_sections", _names(self.config_sections, "section name", SECTION_NAME))
        object.__setattr__(self, "tools", _names(self.tools, "tool name", TOOL_NAME))
        if self.name in self.depends_on:
            raise StudyError(f"Stage {self.name!r} depends on itself")


class Pipeline:
    """Registry of stages and configuration sections; orders stages topologically (ties: registration order)."""

    def __init__(self, name: str, description: str = "", *, secrets: Iterable[str] = (), module: str | None = None):
        if not isinstance(name, str) or not PIPELINE_NAME.match(name):
            raise StudyError(f"Invalid pipeline name {name!r}")
        self.name = name
        self.description = description
        self.secrets = _names(tuple(secrets), "secret environment variable", ENVIRONMENT_NAME)
        self.module = module or sys._getframe(1).f_globals.get("__name__", "")
        self._stages: dict[str, Stage] = {}
        self._sections: dict[str, dict] = {}

    def __repr__(self) -> str:
        return f"Pipeline({self.name!r}, stages={list(self._stages)})"

    def section(self, name: str, schema: dict) -> dict:
        """Declare the JSON schema of configuration section ``[name]``. Mark path values with ``"format": "path"``;
        object schemas reject unknown keys unless they declare ``additionalProperties``; ``default`` values are filled
        in and recorded."""
        if not isinstance(name, str) or not SECTION_NAME.match(name) or name in BUILTIN_SECTIONS:
            raise StudyError(f"Invalid or reserved section name {name!r}")
        if name in self._sections:
            raise StudyError(f"Section [{name}] is already declared in pipeline {self.name!r}")
        try:
            Draft202012Validator.check_schema(schema)
        except SchemaError as error:
            raise StudyError(f"Invalid schema for section [{name}]: {error.message}") from error
        self._sections[name] = copy.deepcopy(schema)
        return schema

    def register(self, stage: Stage) -> Stage:
        if not isinstance(stage, Stage):
            raise StudyError(f"Expected a Stage, got {type(stage).__name__}")
        if stage.name in self._stages:
            raise StudyError(f"Stage {stage.name!r} is already registered in pipeline {self.name!r}")
        self._stages[stage.name] = stage
        return stage

    def stage(self, name: str | None = None, *, depends_on=(), config_sections=(), input_files=None, tools=(),
              description: str | None = None):
        """Decorator registering ``function(context)`` as a stage named after the function by default."""
        def decorate(function):
            summary = description if description is not None else (function.__doc__ or "").strip().split("\n")[0]
            self.register(Stage(name or function.__name__, function, depends_on, config_sections, input_files, tools, summary))
            return function
        return decorate

    def sections(self) -> dict[str, dict]:
        """Schemas of the declared sections read by at least one stage, in declaration order."""
        used = set()
        for stage in self._stages.values():
            for section in stage.config_sections:
                if section not in BUILTIN_SECTIONS and section not in self._sections:
                    raise StudyError(f"Stage {stage.name!r} reads undeclared section [{section}] of pipeline {self.name!r}")
                used.add(section)
        return {name: copy.deepcopy(schema) for name, schema in self._sections.items() if name in used}

    def ordered(self) -> list[Stage]:
        for stage in self._stages.values():
            unknown = [name for name in stage.depends_on if name not in self._stages]
            if unknown:
                raise StudyError(f"Stage {stage.name!r} depends on unknown stage(s): {', '.join(unknown)}")
        self.sections()
        order, done, remaining = [], set(), dict(self._stages)
        while remaining:
            ready = next((stage for stage in remaining.values() if set(stage.depends_on) <= done), None)
            if ready is None:
                raise StudyError(f"Pipeline {self.name!r} has a dependency cycle among: {', '.join(remaining)}")
            order.append(ready)
            done.add(ready.name)
            del remaining[ready.name]
        if not order:
            raise StudyError(f"Pipeline {self.name!r} has no stages")
        return order

    def ancestors(self, name: str) -> set[str]:
        result, pending = set(), list(self._stages[name].depends_on)
        while pending:
            current = pending.pop()
            if current not in result:
                result.add(current)
                pending.extend(self._stages[current].depends_on)
        return result

    def select(self, until: str | None = None, only: Iterable[str] = ()) -> tuple[set[str], set[str]]:
        """Return (stages to execute, completed dependencies to verify). ``until`` selects a stage and everything it
        depends on; ``only`` selects exactly the named stages, whose other dependencies must already be complete."""
        only = (only,) if isinstance(only, str) else tuple(only)
        if until and only:
            raise StudyError("--until and --only are mutually exclusive")
        unknown = [name for name in ([until] if until else []) + list(only) if name not in self._stages]
        if unknown:
            raise StudyError(f"Unknown stage(s) {', '.join(unknown)}; pipeline {self.name!r} has: {', '.join(self._stages)}")
        if until:
            return self.ancestors(until) | {until}, set()
        if only:
            execute = set(only)
            return execute, set().union(*(self.ancestors(name) for name in only)) - execute
        return set(self._stages), set()


def register_pipeline(pipeline: Pipeline, *, replace: bool = False) -> Pipeline:
    """Make an in-process pipeline available to ``[study] pipeline = "<name>"``."""
    if not isinstance(pipeline, Pipeline):
        raise StudyError(f"Expected a Pipeline, got {type(pipeline).__name__}")
    if not replace and _REGISTERED.get(pipeline.name, pipeline) is not pipeline:
        raise StudyError(f"Pipeline {pipeline.name!r} is already registered")
    _REGISTERED[pipeline.name] = pipeline
    return pipeline


def register_pipeline_reference(name: str, reference: str, *, replace: bool = False) -> None:
    """Map a pipeline name to the ``module:attribute`` defining it; the module is imported on first use."""
    module, _, attribute = reference.partition(":")
    if not PIPELINE_NAME.match(name) or not module or not attribute.isidentifier():
        raise StudyError(f"Invalid pipeline reference {name!r} -> {reference!r}; expected 'package.module:ATTRIBUTE'")
    if not replace and PIPELINE_REFERENCES.get(name, reference) != reference:
        raise StudyError(f"Pipeline name {name!r} already refers to {PIPELINE_REFERENCES[name]}")
    PIPELINE_REFERENCES[name] = reference


def available_pipelines() -> list[str]:
    return sorted(set(_REGISTERED) | set(PIPELINE_REFERENCES))


def load_pipeline(name: str) -> Pipeline:
    """Resolve a registered name, a :data:`PIPELINE_REFERENCES` entry or a ``module:attribute`` reference."""
    if name in _REGISTERED:
        return _REGISTERED[name]
    reference = PIPELINE_REFERENCES.get(name, name if ":" in name else None)
    if reference is None:
        raise StudyError(f"Unknown pipeline {name!r}; available: {', '.join(available_pipelines())} "
                         "(or name a 'package.module:ATTRIBUTE')")
    module_name, _, attribute = reference.partition(":")
    try:
        module = importlib.import_module(module_name)
    except ImportError as error:
        raise StudyError(f"Cannot import pipeline module {module_name!r}: {error}") from error
    pipeline = getattr(module, attribute, None)
    if not isinstance(pipeline, Pipeline):
        raise StudyError(f"{reference} is not a Pipeline")
    if ":" not in name and pipeline.name != name:
        raise StudyError(f"{reference} defines pipeline {pipeline.name!r}, not {name!r}")
    return pipeline


@dataclass(frozen=True, eq=False)
class UpstreamAttempt:
    """The exact completed attempt of a dependency that a downstream stage reads."""
    stage: str
    attempt: int
    directory: Path
    relative: str
    stage_json_sha256: str
    record: Mapping

    @property
    def output(self) -> Path:
        return self.directory / OUTPUT

    def path(self, relative: str) -> Path:
        """An output file of this attempt; it must be listed in its stage.json."""
        if relative not in {entry["path"] for entry in self.record["outputs"]}:
            raise StudyError(f"{self.relative} has no recorded output {relative!r}")
        return self.output / relative


class StageContext:
    """Everything a running stage may use. ``config`` holds only the sections the stage declared (secrets
    substituted); ``upstream`` holds only its direct dependencies."""

    def __init__(self, *, name: str, study_id: str, title: str, attempt: int, directory: Path, config: dict,
                 upstream: Mapping[str, UpstreamAttempt], tools: Mapping[str, str], mode: str, run_directory: Path,
                 log: logging.Logger):
        self.name, self.study_id, self.title, self.attempt = name, study_id, title, attempt
        self.directory, self.output, self.run_directory = directory, directory / OUTPUT, run_directory
        self.config, self.upstream, self.mode, self.log = config, MappingProxyType(dict(upstream)), mode, log
        self._tools = dict(tools)

    def tool(self, name: str) -> str:
        """Absolute path of a declared external tool, exactly as identified and hashed for this run."""
        if name not in self._tools:
            raise StudyError(f"Stage {self.name!r} did not declare tool {name!r}")
        return self._tools[name]
