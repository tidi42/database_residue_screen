"""Study configuration: TOML, environment substitution, secret redaction, path resolution, validation and hashing.

One TOML file configures a study. ``[study]`` names the study and its pipeline, the optional ``[tools]`` table names
external programs whose identity is recorded, and every other table is a section declared by the pipeline and read
by at least one of its stages. String values may reference environment variables as ``${NAME}`` (``$$`` is a literal
dollar sign). Schema properties marked ``"format": "path"`` are resolved against the directory of the file.
"""
import copy
import hashlib
import math
import os
import re
import tomllib
from dataclasses import dataclass
from datetime import date, datetime, time
from pathlib import Path
from typing import Iterable, Mapping

from jsonschema import Draft202012Validator
from jsonschema.exceptions import SchemaError

from .common import REDACTED, StudyError, json_sha256

SECTION_NAME = re.compile(r"^[a-z][a-z0-9_]{0,63}$")
ENVIRONMENT_NAME = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")
BUILTIN_SECTIONS = ("study", "tools")
SECRET_WORDS = frozenset({"KEY", "APIKEY", "TOKEN", "SECRET", "PASSWORD", "PASSWD", "PASSPHRASE", "CREDENTIAL", "CREDENTIALS"})
NCBI_TOOLS = frozenset({
    "blastn", "blastp", "blastx", "tblastn", "tblastx", "psiblast", "rpsblast", "rpstblastn", "deltablast", "blastdbcmd",
    "blastdbcheck", "blastdb_aliastool", "blast_formatter", "makeblastdb", "makembindex", "makeprofiledb", "dustmasker",
    "segmasker", "windowmasker", "convert2blastmask", "blastn_vdb", "tblastn_vdb", "blast_formatter_vdb", "blast_vdb_cmd",
})
_REFERENCE = re.compile(r"\$(?:(\$)|\{([^{}]*)\}|(\{))")

STUDY_SCHEMA = {
    "type": "object",
    "additionalProperties": False,
    "required": ["id", "title", "pipeline"],
    "properties": {
        "id": {"type": "string", "pattern": "^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$"},
        "title": {"type": "string", "minLength": 1},
        "pipeline": {"type": "string", "minLength": 1},
        "secrets": {"type": "array", "uniqueItems": True, "items": {"type": "string", "pattern": ENVIRONMENT_NAME.pattern}},
    },
}
TOOLS_SCHEMA = {
    "type": "object",
    "propertyNames": {"pattern": "^[A-Za-z][A-Za-z0-9_.+-]{0,63}$"},
    "additionalProperties": {"oneOf": [
        {"type": "string", "minLength": 1},
        {"type": "object", "additionalProperties": False, "required": ["command"],
         "properties": {"command": {"type": "string", "minLength": 1}, "version_args": {"type": "array", "items": {"type": "string"}}}},
    ]},
}


def is_secret(name: str, declared: Iterable[str] = ()) -> bool:
    """Declared secrets and variables named like credentials (e.g. NCBI_API_KEY, GITHUB_TOKEN) are never recorded."""
    return name in set(declared) or any(word in SECRET_WORDS for word in name.upper().split("_"))


def default_version_args(command: str) -> list[str]:
    return ["-version"] if Path(command).name in NCBI_TOOLS else ["--version"]


@dataclass(frozen=True)
class ConfigFile:
    path: Path
    data: dict
    sha256: str

    @property
    def pipeline(self) -> str:
        study = self.data.get("study")
        name = study.get("pipeline") if isinstance(study, dict) else None
        if not isinstance(name, str) or not name:
            raise StudyError(f'{self.path}: [study] must name the pipeline, for example pipeline = "demo"')
        if "$" in name:
            raise StudyError(f"{self.path}: [study] pipeline must be a literal name")
        return name


@dataclass(frozen=True, eq=False)
class StudyConfig:
    """A validated configuration. ``values`` holds substituted secrets; ``recorded`` is the redacted form that is hashed
    and written to provenance records."""
    path: Path
    file_sha256: str
    values: dict
    recorded: dict
    sha256: str
    environment: dict
    secret_values: tuple

    @property
    def directory(self) -> Path:
        return self.path.parent

    @property
    def study_id(self) -> str:
        return self.values["study"]["id"]

    @property
    def title(self) -> str:
        return self.values["study"]["title"]

    @property
    def pipeline(self) -> str:
        return self.values["study"]["pipeline"]

    @property
    def tools(self) -> dict:
        return copy.deepcopy(self.values["tools"])

    def section(self, name: str) -> dict:
        return copy.deepcopy(self.values[name])

    def recorded_section(self, name: str) -> dict:
        return copy.deepcopy(self.recorded[name])

    def redact(self, text: str) -> str:
        return redactor(self.secret_values)(text)


def redactor(secret_values: Iterable[str]):
    values = sorted({value for value in secret_values if value}, key=len, reverse=True)

    def redact(text: str) -> str:
        for value in values:
            text = text.replace(value, REDACTED)
        return text
    return redact


def read_config_file(path: Path) -> ConfigFile:
    path = Path(path).resolve()
    try:
        content = path.read_bytes()
    except OSError as error:
        raise StudyError(f"Cannot read study configuration {path}: {error}") from error
    try:
        data = tomllib.loads(content.decode("utf-8"))
    except (UnicodeDecodeError, tomllib.TOMLDecodeError) as error:
        raise StudyError(f"{path}: invalid TOML: {error}") from error
    return ConfigFile(path, data, hashlib.sha256(content).hexdigest())


def config_schema(sections: Mapping[str, dict]) -> dict:
    """JSON Schema of a configuration: [study], [tools] and the declared sections, all closed to unknown keys."""
    properties = {"study": copy.deepcopy(STUDY_SCHEMA), "tools": copy.deepcopy(TOOLS_SCHEMA)}
    required = ["study"]
    for name, schema in sections.items():
        if name in BUILTIN_SECTIONS or not SECTION_NAME.match(name):
            raise StudyError(f"Invalid or reserved section name: {name!r}")
        section = copy.deepcopy(schema)
        section.setdefault("type", "object")
        if section["type"] != "object":
            raise StudyError(f"Section [{name}] must be described by an object schema")
        properties[name] = _closed(section)
        if section.get("required"):
            required.append(name)
    schema = {"type": "object", "additionalProperties": False, "required": required, "properties": properties}
    try:
        Draft202012Validator.check_schema(schema)
    except SchemaError as error:
        raise StudyError(f"Invalid section schema: {error.message}") from error
    return schema


def resolve_config(raw: ConfigFile, sections: Mapping[str, dict], *, environ: Mapping[str, str] | None = None,
                   secrets: Iterable[str] = ()) -> StudyConfig:
    environ = os.environ if environ is None else environ
    declared = _declared_secrets(raw) | set(secrets)
    used: dict[str, str] = {}
    actual, recorded = _substitute(raw.data, "", environ, declared, used)
    secret_values = tuple(environ[name] for name in used if is_secret(name, declared))
    schema = config_schema(sections)
    for tree in (actual, recorded):
        for name, section_schema in schema["properties"].items():
            if name not in tree and name not in schema["required"]:
                tree[name] = {}
            _apply_defaults(tree.get(name), section_schema)
    _validate(actual, schema, redactor(secret_values))
    for tree in (actual, recorded):
        tree["tools"] = _normalise_tools(tree["tools"], raw.path.parent)
        for name, section_schema in schema["properties"].items():
            if name not in BUILTIN_SECTIONS:
                tree[name] = _resolve_paths(tree[name], section_schema, raw.path.parent, name)
    return StudyConfig(raw.path, raw.sha256, actual, recorded, json_sha256(recorded), dict(sorted(used.items())), secret_values)


def load_config(path: Path, sections: Mapping[str, dict], *, environ: Mapping[str, str] | None = None,
                secrets: Iterable[str] = ()) -> StudyConfig:
    return resolve_config(read_config_file(path), sections, environ=environ, secrets=secrets)


def _declared_secrets(raw: ConfigFile) -> set[str]:
    study = raw.data.get("study")
    names = study.get("secrets", []) if isinstance(study, dict) else []
    if not isinstance(names, list) or not all(isinstance(name, str) and ENVIRONMENT_NAME.match(name) for name in names):
        raise StudyError("[study] secrets must be a list of environment variable names")
    return set(names)


def _child(where: str, key) -> str:
    return f"{where}.{key}" if where else str(key)


def _substitute_text(text: str, where: str, environ, declared, used) -> tuple[str, str]:
    actual, recorded, position = [], [], 0
    for match in _REFERENCE.finditer(text):
        actual.append(text[position:match.start()])
        recorded.append(text[position:match.start()])
        position = match.end()
        if match.group(1):
            actual.append("$")
            recorded.append("$")
            continue
        if match.group(3) is not None:
            raise StudyError(f"{where}: unterminated '${{' in {text!r}")
        name = match.group(2)
        if not ENVIRONMENT_NAME.match(name):
            raise StudyError(f"{where}: invalid environment variable reference ${{{name}}}")
        value = environ.get(name)
        if not value:
            raise StudyError(f"{where}: environment variable {name} is {'empty' if value == '' else 'not set'}")
        shown = REDACTED if is_secret(name, declared) else value
        used[name] = shown
        actual.append(value)
        recorded.append(shown)
    actual.append(text[position:])
    recorded.append(text[position:])
    return "".join(actual), "".join(recorded)


def _substitute(value, where: str, environ, declared, used):
    if isinstance(value, str):
        return _substitute_text(value, where, environ, declared, used)
    if isinstance(value, bool | int):
        return value, value
    if isinstance(value, float):
        if not math.isfinite(value):
            raise StudyError(f"{where}: non-finite number {value}")
        return value, value
    if isinstance(value, dict):
        pairs = {key: _substitute(item, _child(where, key), environ, declared, used) for key, item in value.items()}
        return {key: pair[0] for key, pair in pairs.items()}, {key: pair[1] for key, pair in pairs.items()}
    if isinstance(value, list):
        pairs = [_substitute(item, f"{where}[{index}]", environ, declared, used) for index, item in enumerate(value)]
        return [pair[0] for pair in pairs], [pair[1] for pair in pairs]
    if isinstance(value, datetime | date | time):
        raise StudyError(f"{where}: TOML dates and times are not supported; write them as quoted strings")
    raise StudyError(f"{where}: unsupported value of type {type(value).__name__}")


def _closed(schema):
    """Reject unknown keys in every object schema that does not say otherwise."""
    if not isinstance(schema, dict):
        return schema
    if (schema.get("type") == "object" or "properties" in schema) and not any(
            key in schema for key in ("additionalProperties", "patternProperties", "unevaluatedProperties")):
        schema["additionalProperties"] = False
    for item in schema.get("properties", {}).values():
        _closed(item)
    for key in ("additionalProperties", "items"):
        if isinstance(schema.get(key), dict):
            _closed(schema[key])
    return schema


def _apply_defaults(value, schema) -> None:
    if not isinstance(schema, dict):
        return
    if isinstance(value, dict):
        properties = schema.get("properties", {})
        for key, item in properties.items():
            if key not in value and isinstance(item, dict) and "default" in item:
                value[key] = copy.deepcopy(item["default"])
            if key in value:
                _apply_defaults(value[key], item)
        extra = schema.get("additionalProperties")
        if isinstance(extra, dict):
            for key in value:
                if key not in properties:
                    _apply_defaults(value[key], extra)
    elif isinstance(value, list) and isinstance(schema.get("items"), dict):
        for item in value:
            _apply_defaults(item, schema["items"])


def _resolve_paths(value, schema, directory: Path, where: str):
    if not isinstance(schema, dict):
        return value
    if schema.get("format") == "path" and isinstance(value, str):
        if not value:
            raise StudyError(f"{where}: empty path")
        return str(Path(directory, value).resolve())
    if isinstance(value, dict):
        properties, extra = schema.get("properties", {}), schema.get("additionalProperties")
        return {key: _resolve_paths(item, properties.get(key, extra), directory, _child(where, key)) for key, item in value.items()}
    if isinstance(value, list) and isinstance(schema.get("items"), dict):
        return [_resolve_paths(item, schema["items"], directory, f"{where}[{index}]") for index, item in enumerate(value)]
    return value


def _normalise_tools(tools: dict, directory: Path) -> dict:
    result = {}
    for name, specification in sorted(tools.items()):
        if isinstance(specification, str):
            specification = {"command": specification}
        command = specification["command"]
        if os.sep in command or (os.altsep and os.altsep in command):
            command = str(Path(directory, command).resolve())
        result[name] = {"command": command, "version_args": list(specification.get("version_args", default_version_args(command)))}
    return result


def _location(path) -> str:
    return "".join(f"[{part}]" if isinstance(part, int) else ("." if index else "") + str(part) for index, part in enumerate(path))


def _describe(error) -> str:
    where = _location(error.absolute_path)
    if error.validator == "additionalProperties" and isinstance(error.instance, dict) and isinstance(error.schema, dict):
        known = error.schema.get("properties", {})
        patterns = error.schema.get("patternProperties", {})
        unknown = sorted(key for key in error.instance if key not in known and not any(re.search(pattern, key) for pattern in patterns))
        allowed = ", ".join(sorted(known)) or "none"
        if not where:
            sections = ", ".join(f"[{key}]" for key in unknown)
            return f"unknown section(s) {sections} (allowed: {allowed})"
        return f"unknown key(s) in [{where}]: {', '.join(unknown)} (allowed: {allowed})"
    return f"{where or 'configuration'}: {error.message}"


def _validate(instance: dict, schema: dict, redact) -> None:
    errors = sorted(Draft202012Validator(schema).iter_errors(instance),
                    key=lambda error: ([str(part) for part in error.absolute_path], error.message))
    if errors:
        raise StudyError("Invalid study configuration: " + "; ".join(redact(_describe(error)) for error in errors))
