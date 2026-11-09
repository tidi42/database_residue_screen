import csv
import json
import math
from importlib.resources import files
from pathlib import Path

from jsonschema import Draft202012Validator


class ContractError(ValueError):
    pass


def table_schema(strings: tuple, integers: tuple = (), numbers: tuple = ()) -> dict:
    properties = {name: {"type": "string"} for name in strings}
    properties.update({name: {"type": "integer"} for name in integers})
    properties.update({name: {"type": "number"} for name in numbers})
    return dict(type="object", additionalProperties=False, required=list(properties), properties=properties)


def schema(name: str) -> dict:
    resource = files("database_residue_screen").joinpath("schemas", name + ".json")
    value = json.loads(resource.read_text(encoding="utf-8"))
    Draft202012Validator.check_schema(value)
    return value


def validate_record(record: dict, definition: dict) -> None:
    errors = sorted(Draft202012Validator(definition).iter_errors(record), key=lambda error: str(error.path))
    if errors:
        raise ContractError("; ".join(f"{'.'.join(map(str, error.path)) or 'record'}: {error.message}" for error in errors))


def read_json(path: Path, definition: dict | None = None) -> dict:
    def unique_object(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise ContractError(f"Duplicate JSON key: {key}")
            result[key] = value
        return result

    def reject_constant(value):
        raise ContractError(f"Non-finite JSON number: {value}")

    value = json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=unique_object, parse_constant=reject_constant)
    if not isinstance(value, dict):
        raise ContractError("Expected a JSON object")
    if definition is not None:
        validate_record(value, definition)
    return value


def read_table(path: Path, definition: dict) -> list[dict]:
    columns = list(definition["properties"])
    result = []
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or len(reader.fieldnames) != len(set(reader.fieldnames)) or set(reader.fieldnames) != set(columns):
            raise ContractError(f"{path.name}: expected exactly these columns: {', '.join(columns)}")
        for line, row in enumerate(reader, 2):
            if None in row or any(value is None for value in row.values()):
                raise ContractError(f"{path.name}:{line}: wrong field count")
            for name, rule in definition["properties"].items():
                try:
                    if rule.get("type") == "integer":
                        row[name] = int(row[name])
                    elif rule.get("type") == "number":
                        row[name] = float(row[name])
                        if not math.isfinite(row[name]):
                            raise ValueError("non-finite value")
                except ValueError as error:
                    raise ContractError(f"{path.name}:{line}:{name}: {error}") from error
            try:
                validate_record(row, definition)
            except ContractError as error:
                raise ContractError(f"{path.name}:{line}: {error}") from error
            result.append(row)
    return result


def write_table(path: Path, rows: list[dict], definition: dict) -> None:
    for row in rows:
        validate_record(row, definition)
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(definition["properties"]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
