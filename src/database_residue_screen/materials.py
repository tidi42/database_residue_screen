"""Run-time study materials derived from a built panel and its pinned sources.

The study configuration stores rules, never sequences: the panel subset (by canonical parent ID or physical query ID),
census flank elements (a segment or slice of a named panel query), barcode catalogues (records of pinned sources read
by the panel builder's parsers) and benchmark templates (named panel queries). Every derivation is checked
(segments of a family must agree, stripped suffixes must match, catalogue coverage checks must hold) and fails
instead of guessing.
"""
import hashlib
import re

from .constructs import Library
from .contracts import ContractError, table_schema, validate_record
from .sequences import canonical_id, reverse_complement

QUERY_NAME = {"type": "string", "pattern": "^[A-Za-z0-9_.:+-]+$"}
SUBSET_SCHEMA = {"type": "object", "additionalProperties": False, "required": ["parents", "query_ids"],
                 "properties": {"parents": {"type": "integer", "minimum": 0},
                                "query_ids": {"type": "array", "uniqueItems": True, "items": QUERY_NAME}}}
ELEMENT_RULE_SCHEMA = {
    "type": "object", "additionalProperties": False, "required": ["name", "query"],
    "properties": {"name": QUERY_NAME, "query": QUERY_NAME, "segment": {"type": "string", "pattern": "^[a-z_]+$"},
                   "start": {"type": "integer", "minimum": 1}, "end": {"type": "integer", "minimum": 1},
                   "strip_suffix": {"type": "object", "additionalProperties": False, "required": ["query", "reverse_complement"],
                                    "properties": {"query": QUERY_NAME, "reverse_complement": {"type": "boolean"}}},
                   "reverse_complement": {"type": "boolean"}, "same_in_family": {"type": "boolean"}, "description": {"type": "string"}}}
CATALOGUE_PARSERS = ("dorado_barcodes", "ont_native_table", "ont_barcode_table", "fasta")
CATALOGUE_RULE_SCHEMA = {
    "type": "object", "additionalProperties": False, "required": ["set", "source", "parser", "label"],
    "properties": {"set": QUERY_NAME, "source": {"type": "string", "minLength": 1}, "parser": {"enum": list(CATALOGUE_PARSERS)},
                   "label": QUERY_NAME, "length": {"type": "integer", "minimum": 1}, "name_pattern": {"type": "string", "minLength": 1}}}
CATALOGUE_CHECK_SCHEMA = {"oneOf": [
    {"type": "object", "additionalProperties": False, "required": ["set", "family", "segment"],
     "properties": {"set": QUERY_NAME, "family": {"type": "string", "minLength": 1}, "segment": {"type": "string", "pattern": "^[a-z_]+$"}}},
    {"type": "object", "additionalProperties": False, "required": ["set", "source", "parser"],
     "properties": {"set": QUERY_NAME, "source": {"type": "string", "minLength": 1}, "parser": {"enum": list(CATALOGUE_PARSERS)},
                    "length": {"type": "integer", "minimum": 1}, "name_pattern": {"type": "string", "minLength": 1}}}]}
TEMPLATE_RULE_SCHEMA = {
    "type": "object", "additionalProperties": False, "required": ["name", "query", "repeat_from"],
    "properties": {"name": {"type": "string", "pattern": "^[A-Za-z0-9_]+$"}, "query": QUERY_NAME,
                   "repeat_from": {"type": "integer", "minimum": 1}, "description": {"type": "string"}}}
SEGMENT_SCHEMA = table_schema(("query_id", "kind", "element", "verification"), ("segment_index", "start", "end"))
CATALOGUE_SCHEMA = table_schema(("set", "name", "sequence"))


def sha256_text(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def parse_segments(text: str) -> list[dict]:
    """`front:1-14;barcode:15-38=NB01` -> [{kind, start, end, element}] (one-based, inclusive)."""
    segments = []
    for item in text.split(";"):
        match = re.fullmatch(r"([a-z_]+):(\d+)-(\d+)(?:=(.+))?", item)
        if match is None:
            raise ContractError(f"Invalid segment {item!r}")
        segments.append(dict(kind=match.group(1), start=int(match.group(2)), end=int(match.group(3)), element=match.group(4) or ""))
    return segments


def select_subset(sequences: dict[str, str], parents: int = 0, query_ids=()) -> tuple[list[str], dict]:
    """Physical query IDs of the study panel and a record of the selection.

    `parents` > 0 keeps the canonical parents with the smallest parent IDs (canonical sequence SHA-256) and every
    physical query of them; `query_ids` keeps exactly those physical queries; neither keeps the whole panel."""
    query_ids = list(query_ids)
    if parents < 0 or (parents and query_ids):
        raise ContractError("Select either a number of canonical parents or explicit query IDs, not both")
    by_parent = {}
    for name, sequence in sequences.items():
        by_parent.setdefault(canonical_id(sequence), []).append(name)
    if query_ids:
        unknown = sorted(set(query_ids) - set(sequences))
        if unknown:
            raise ContractError("Subset names queries that are not in the panel: " + ", ".join(unknown))
        selected, mode = sorted(query_ids), "query_ids"
    elif parents:
        if parents > len(by_parent):
            raise ContractError(f"Subset asks for {parents} canonical parents; the panel has {len(by_parent)}")
        chosen = sorted(by_parent)[:parents]
        selected, mode = sorted(name for parent in chosen for name in by_parent[parent]), "parents"
    else:
        selected, mode = sorted(sequences), "complete_panel"
    order = {name: index for index, name in enumerate(sequences)}
    selected.sort(key=order.__getitem__)
    chosen_parents = sorted({canonical_id(sequences[name]) for name in selected})
    record = dict(mode=mode, requested_parents=parents, requested_query_ids=sorted(query_ids), panel_queries=len(sequences),
                  panel_parents=len(by_parent), selected_queries=len(selected), selected_parents=len(chosen_parents),
                  query_ids=selected, parent_ids=chosen_parents,
                  rule=("complete panel" if mode == "complete_panel" else "canonical parents with the smallest parent IDs (SHA-256 of the "
                        "canonical sequence) and all their physical queries" if mode == "parents" else "explicit physical query IDs"))
    return selected, record


def segment_rows(entries: dict[str, dict], query_ids) -> list[dict]:
    rows = []
    for name in query_ids:
        for index, segment in enumerate(parse_segments(entries[name]["segments"]), 1):
            rows.append(dict(query_id=name, segment_index=index, kind=segment["kind"], start=segment["start"], end=segment["end"],
                             element=segment["element"], verification=entries[name]["verification"]))
    return rows


def _piece(rule: dict, name: str, entries: dict, sequences: dict) -> str:
    sequence = sequences[name]
    if "segment" in rule:
        found = [segment for segment in parse_segments(entries[name]["segments"]) if segment["kind"] == rule["segment"]]
        if len(found) != 1:
            raise ContractError(f"Element {rule['name']}: query {name} has {len(found)} '{rule['segment']}' segments")
        sequence = sequence[found[0]["start"] - 1:found[0]["end"]]
    if "start" in rule or "end" in rule:
        start, end = rule.get("start", 1), rule.get("end", len(sequence))
        if not 1 <= start <= end <= len(sequence):
            raise ContractError(f"Element {rule['name']}: slice {start}-{end} is outside its {len(sequence)}-nt source")
        sequence = sequence[start - 1:end]
    if "strip_suffix" in rule:
        suffix = sequences[rule["strip_suffix"]["query"]]
        suffix = reverse_complement(suffix) if rule["strip_suffix"]["reverse_complement"] else suffix
        if len(sequence) <= len(suffix) or not sequence.endswith(suffix):
            raise ContractError(f"Element {rule['name']}: query {name} does not end with the declared suffix query")
        sequence = sequence[:-len(suffix)]
    return reverse_complement(sequence) if rule.get("reverse_complement") else sequence


def derive_elements(rules: list[dict], entries: dict[str, dict], sequences: dict[str, str]) -> list[dict]:
    """Census flank elements from the built panel; `same_in_family` requires every query of the family to agree."""
    elements, seen = [], set()
    for rule in rules:
        validate_record(rule, ELEMENT_RULE_SCHEMA)
        if rule["name"] in seen:
            raise ContractError(f"Duplicate flank element {rule['name']}")
        seen.add(rule["name"])
        names = [rule["query"]] + ([rule["strip_suffix"]["query"]] if "strip_suffix" in rule else [])
        missing = [name for name in names if name not in sequences]
        if missing:
            raise ContractError(f"Element {rule['name']} names queries that are not in the panel: {', '.join(missing)}")
        sequence = _piece(rule, rule["query"], entries, sequences)
        family = entries[rule["query"]]["family"]
        checked = [rule["query"]]
        if rule.get("same_in_family"):
            checked = [name for name, entry in entries.items() if entry["family"] == family]
            differing = [name for name in checked if _piece(rule, name, entries, sequences) != sequence]
            if differing:
                raise ContractError(f"Element {rule['name']} differs between queries of family {family}: {', '.join(differing[:5])}")
        elements.append(dict(name=rule["name"], sequence=sequence, length=len(sequence), sha256=sha256_text(sequence),
                             source_query=rule["query"], family=family, queries_checked=len(checked),
                             rule={key: value for key, value in rule.items() if key != "description"}))
    return elements


def _records(library: Library, rule: dict) -> dict[str, str]:
    kind = rule["parser"]
    if kind == "dorado_barcodes":
        records = library.parsed("dorado", rule["source"]).barcodes
    elif kind == "ont_native_table":
        records = library.parsed("ont", rule["source"]).native
    elif kind == "ont_barcode_table":
        records = library.parsed("ont", rule["source"]).rapid_table
    else:
        records = library.parsed("fasta", rule["source"])
    selected = {name: sequence.upper() for name, sequence in records.items()
                if ("length" not in rule or len(sequence) == rule["length"])
                and ("name_pattern" not in rule or re.fullmatch(rule["name_pattern"], name))}
    if not selected:
        raise ContractError(f"Catalogue rule {rule['parser']} on {rule['source']} selects no record")
    if any(set(sequence) - set("ACGT") for sequence in selected.values()):
        raise ContractError(f"Catalogue records of {rule['source']} must be unambiguous ACGT")
    return selected


def derive_catalogue(rules: list[dict], checks: list[dict], library: Library, entries: dict[str, dict],
                     sequences: dict[str, str]) -> tuple[list[dict], dict]:
    """Barcode catalogue rows (set, name, sequence) from pinned sources, and the coverage checks that must hold."""
    rows, seen = [], set()
    for rule in rules:
        validate_record(rule, CATALOGUE_RULE_SCHEMA)
        for name, sequence in sorted(_records(library, rule).items()):
            key = rule["set"], f"{rule['label']}:{name}", sequence
            if key not in seen:
                seen.add(key)
                rows.append(dict(set=key[0], name=key[1], sequence=sequence))
    sets = {}
    for row in rows:
        sets.setdefault(row["set"], set()).add(row["sequence"])
    results = []
    for check in checks:
        validate_record(check, CATALOGUE_CHECK_SCHEMA)
        catalogue = sets.get(check["set"], set())
        if "family" in check:
            members = [name for name, entry in entries.items() if entry["family"] == check["family"]]
            if not members:
                raise ContractError(f"Catalogue check names an unknown panel family: {check['family']}")
            required = {_piece(dict(name="check", segment=check["segment"]), name, entries, sequences) for name in members}
            origin = f"panel family {check['family']} segment {check['segment']}"
        else:
            required = set(_records(library, check).values())
            origin = f"{check['parser']} records of {check['source']}"
        absent = sorted(sequence for sequence in required if sequence not in catalogue and reverse_complement(sequence) not in catalogue)
        if absent:
            raise ContractError(f"Catalogue set {check['set']} lacks {len(absent)} of {len(required)} sequences of the {origin}")
        results.append(dict(check, required_sequences=len(required), status="passed"))
    summary = dict(sets={name: dict(entries=sum(row["set"] == name for row in rows), distinct_sequences=len(values))
                         for name, values in sorted(sets.items())},
                   checks=results, sources=sorted({rule["source"] for rule in rules}),
                   matching="Census matches catalogue entries in either orientation")
    return rows, summary


def derive_templates(rules: list[dict], sequences: dict[str, str]) -> list[dict]:
    """Benchmark templates: a named panel query extended by repeating its tail from `repeat_from` (one-based)."""
    templates, seen = [], set()
    for rule in rules:
        validate_record(rule, TEMPLATE_RULE_SCHEMA)
        if rule["name"] in seen or rule["query"] not in sequences:
            raise ContractError(f"Benchmark template {rule['name']}: duplicate name or query {rule['query']} not in the panel")
        seen.add(rule["name"])
        sequence = sequences[rule["query"]]
        if rule["repeat_from"] > len(sequence):
            raise ContractError(f"Benchmark template {rule['name']}: repeat_from exceeds the query length")
        templates.append(dict(name=rule["name"], query=rule["query"], sequence=sequence, repeat_from=rule["repeat_from"],
                              length=len(sequence), sha256=sha256_text(sequence)))
    return templates
