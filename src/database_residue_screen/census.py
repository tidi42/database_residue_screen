"""Flank census of retained strict reference hits (descriptive).

For declared query families every retained strict HSP window is retrieved from the local database with fixed flanks,
checked to contain the physical query exactly at the HSP coordinates, oriented to the family's declared orientation
and classified by ordered flank rules (secondary: at most one mismatch per flank element). The census specification
is data: per family its member queries (conditions or physical query IDs), the classification orientation, ordered
classes and the flank lengths below which an unclassified window is `insufficient_flank`. Class rules name flank
elements and barcode catalogues that are derived at run time from the built panel and its pinned sources (see
materials.py); neither a specification nor this module contains a sequence. Unresolved identifiers and failed
retrievals are listed, never dropped. Counts are windows and accessions; the census assigns context, not biological
origin, and adds no search.
"""
import hashlib
import json
import os
import subprocess
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

from .analysis import ANNOTATED_EVENT_SCHEMA
from .contracts import ContractError, read_json, read_table, table_schema, validate_record, write_table
from .design import digest
from .panels import LOGICAL_SCHEMA
from .sequences import read_fasta, reverse_complement

RESERVED_CLASSES = ("insufficient_flank", "none", "integrity_failed", "retrieval_failed", "unresolved_identifier")
NAME = {"type": "string", "pattern": "^[A-Za-z0-9_.:+-]+$"}
TEST_SCHEMA = {"oneOf": [
    {"type": "object", "additionalProperties": False, "required": ["side", "element"],
     "properties": {"side": {"enum": ["upstream", "downstream"]}, "element": NAME, "offset": {"type": "integer", "minimum": 0}}},
    {"type": "object", "additionalProperties": False, "required": ["side", "catalogue", "length"],
     "properties": {"side": {"enum": ["upstream", "downstream", "either"]}, "catalogue": NAME,
                    "length": {"type": "integer", "minimum": 1}, "offset": {"type": "integer", "minimum": 0}}},
    {"type": "object", "additionalProperties": False, "required": ["taxon"], "properties": {"taxon": {"type": "integer", "minimum": 1}}},
]}
CLASS_SCHEMA = {"type": "object", "additionalProperties": False, "required": ["name", "require"],
                "properties": {"name": NAME, "description": {"type": "string"},
                               "require": {"type": "array", "minItems": 1, "items": TEST_SCHEMA}}}
FAMILY_SCHEMA = {"type": "object", "additionalProperties": False, "required": ["name", "orientation", "required_flank", "classes"],
                 "properties": {"name": NAME, "conditions": {"type": "array", "uniqueItems": True, "items": {"type": "string", "minLength": 1}},
                                "query_ids": {"type": "array", "uniqueItems": True, "items": NAME},
                                "orientation": {"enum": ["query", "reverse_complement"]}, "note": {"type": "string"},
                                "required_flank": {"type": "object", "additionalProperties": False, "required": ["upstream", "downstream"],
                                                   "properties": {"upstream": {"type": "integer", "minimum": 0},
                                                                  "downstream": {"type": "integer", "minimum": 0}}},
                                "classes": {"type": "array", "minItems": 1, "items": CLASS_SCHEMA}}}
SPECIFICATION_SCHEMA = {"type": "object", "additionalProperties": False, "required": ["flank", "families"],
                        "properties": {"flank": {"type": "integer", "minimum": 1},
                                       "families": {"type": "array", "minItems": 1, "items": FAMILY_SCHEMA}}}
WINDOW_SCHEMA = table_schema(("family", "physical_id", "accession_version", "strand", "status", "census_class", "census_class_hamming1",
                              "upstream", "downstream", "taxids", "testability"),
                             ("lower", "upper", "record_length", "query_length", "upstream_length", "downstream_length"))
SUMMARY_SCHEMA = table_schema(("family", "rule", "census_class"), ("windows", "accessions", "family_windows", "family_accessions"))
CATALOGUE_SCHEMA = table_schema(("set", "name", "sequence"))


def validate_specification(specification: dict) -> dict:
    """Schema plus the cross-field rules of a census specification; returns it unchanged."""
    validate_record(specification, SPECIFICATION_SCHEMA)
    names = [family["name"] for family in specification["families"]]
    if len(set(names)) != len(names):
        raise ContractError("Census family names must be unique")
    for family in specification["families"]:
        if not family.get("conditions") and not family.get("query_ids"):
            raise ContractError(f"Census family {family['name']} needs member conditions or query IDs")
        classes = [rule["name"] for rule in family["classes"]]
        if len(set(classes)) != len(classes) or set(classes) & set(RESERVED_CLASSES):
            raise ContractError(f"Census family {family['name']}: class names must be unique and not one of {', '.join(RESERVED_CLASSES)}")
        if max(family["required_flank"].values()) > specification["flank"]:
            raise ContractError(f"Census family {family['name']}: required flank exceeds the retrieved flank")
    return specification


def specification_references(specification: dict) -> dict[str, set]:
    """Names of the flank elements, catalogue sets and taxa that the class rules use."""
    found = dict(elements=set(), catalogues=set(), taxa=set())
    for family in specification["families"]:
        for rule in family["classes"]:
            for test in rule["require"]:
                if "element" in test:
                    found["elements"].add(test["element"])
                elif "catalogue" in test:
                    found["catalogues"].add(test["catalogue"])
                else:
                    found["taxa"].add(test["taxon"])
    return found


def load_elements(path: Path) -> dict[str, str]:
    """Flank elements written by the panel stage: {"elements": [{"name", "sequence", ...}]}."""
    elements = {}
    for item in read_json(path).get("elements", []):
        sequence = item.get("sequence", "")
        if not sequence or set(sequence) - set("ACGT") or item.get("name") in elements:
            raise ContractError(f"{path.name}: invalid or duplicate flank element {item.get('name')!r}")
        elements[item["name"]] = sequence
    return elements


def matches(observed: str, expected: str, mismatches: int) -> bool:
    return len(observed) == len(expected) and sum(a != b for a, b in zip(observed, expected)) <= mismatches


def in_catalogue(observed: str, catalogue: set[str], mismatches: int) -> bool:
    if mismatches == 0:
        return observed in catalogue
    return any(matches(observed, item, mismatches) for item in catalogue if len(item) == len(observed))


def flank_piece(up: str, down: str, side: str, length: int, offset: int = 0) -> str | None:
    """`length` bases of a flank, `offset` bases away from the query; None if the flank is too short."""
    if side == "upstream":
        return up[len(up) - offset - length:len(up) - offset] if len(up) >= offset + length else None
    return down[offset:offset + length] if len(down) >= offset + length else None


def holds(test: dict, up: str, down: str, elements: dict, catalogues: dict, lineage, mismatches: int) -> bool:
    if "taxon" in test:
        return bool(lineage(test["taxon"]))
    if "element" in test:
        expected = elements[test["element"]]
        observed = flank_piece(up, down, test["side"], len(expected), test.get("offset", 0))
        return observed is not None and matches(observed, expected, mismatches)
    catalogue = catalogues.get(test["catalogue"], set())
    for side in (("upstream", "downstream") if test["side"] == "either" else (test["side"],)):
        observed = flank_piece(up, down, side, test["length"], test.get("offset", 0))
        if observed is not None and (in_catalogue(observed, catalogue, mismatches)
                                     or in_catalogue(reverse_complement(observed), catalogue, mismatches)):
            return True
    return False


def classify(family: dict, up: str, down: str, elements: dict, catalogues: dict, lineage=lambda taxon: False,
             mismatches: int = 0) -> str:
    """First class whose every test holds; otherwise `insufficient_flank` when a flank is shorter than required, else `none`.
    Catalogue tests accept either orientation; `mismatches` applies to every element and catalogue comparison."""
    for rule in family["classes"]:
        if all(holds(test, up, down, elements, catalogues, lineage, mismatches) for test in rule["require"]):
            return rule["name"]
    required = family["required_flank"]
    return "insufficient_flank" if len(up) < required["upstream"] or len(down) < required["downstream"] else "none"


def in_lineage(taxids: str, taxonomy, taxon: int) -> bool:
    if taxonomy is None:
        return False
    for raw in taxids.replace(",", ";").split(";"):
        if raw.strip().isdigit():
            resolved, _ = taxonomy.resolve(int(raw))
            if resolved is not None and taxon in taxonomy.lineage(resolved):
                return True
    return False


def window_requests(analysis: Path, panel: Path, specification: dict) -> list[dict]:
    report = read_json(analysis / "analysis.json")
    for name, checksum in report["files"].items():
        if digest(analysis / name) != checksum:
            raise ContractError("Analysis artifact changed: " + name)
    definition = read_json(panel / "panel.json")
    for name, checksum in definition["files"].items():
        if digest(panel / name) != checksum:
            raise ContractError("Panel artifact changed: " + name)
    physical = read_fasta(panel / "physical_queries.fasta")
    logical = read_table(analysis / "logical_queries.tsv", LOGICAL_SCHEMA)
    members = defaultdict(list)
    for row in logical:
        if row["role"] != "real":
            continue
        for family in specification["families"]:
            if row["condition"] in family.get("conditions", ()) or row["logical_id"] in family.get("query_ids", ()):
                members[row["sequence_id"]].append((family["name"], row))
    policy = report["analysis_policy"]
    testability = {row["query_length"]: row["strict_testability"] for row in policy.get("detectability", [])}
    requests, seen = [], set()
    for event in read_table(analysis / "real_events.tsv", ANNOTATED_EVENT_SCHEMA):
        if event["strict"] != "yes":
            continue
        for family, row in members.get(event["query_id"], []):
            lower, upper = min(event["sstart"], event["send"]), max(event["sstart"], event["send"])
            canonical_plus = event["sstart"] <= event["send"]
            physical_plus = canonical_plus == (row["orientation"] == "+")
            key = family, row["logical_id"], event["accession_version"], lower, upper
            if key in seen:
                continue
            seen.add(key)
            sequence = physical[row["logical_id"]]
            requests.append(dict(family=family, physical_id=row["logical_id"], query=sequence, accession_version=event["accession_version"],
                                 lower=lower, upper=upper, record_length=event["slen"], strand="+" if physical_plus else "-",
                                 taxids=event.get("taxids", ""), testability=testability.get(len(sequence), "not_calibrated")))
    return requests


def retrieve(requests: list[dict], prefix: Path, blastdbcmd: str, scratch: Path, flank: int = 50, batch: int = 2000) -> None:
    resolvable = [item for item in requests if not item["accession_version"].startswith("unresolved:")]
    for item in requests:
        if item["accession_version"].startswith("unresolved:"):
            item.update(status="unresolved_identifier", window="")
    for offset in range(0, len(resolvable), batch):
        chunk = resolvable[offset:offset + batch]
        entries = scratch / f"entries_{offset:07d}.txt"
        entries.write_text("".join(f"{item['accession_version']} {max(1, item['lower'] - flank)}-"
                                   f"{min(item['record_length'], item['upper'] + flank)} plus\n" for item in chunk), encoding="ascii")
        process = subprocess.run([blastdbcmd, "-db", str(prefix), "-entry_batch", str(entries), "-target_only", "-outfmt", "%a %s"],
                                 capture_output=True, text=True, env=dict(os.environ, LC_ALL="C", BLAST_USAGE_REPORT="false"))
        lines = [line for line in process.stdout.splitlines() if line.strip()]
        if process.returncode or len(lines) != len(chunk):
            for item in chunk:
                item.setdefault("status", "retrieval_failed")
                item.setdefault("window", "")
            raise ContractError(f"blastdbcmd batch retrieval returned {len(lines)} of {len(chunk)} windows: {process.stderr.strip()[:300]}")
        for item, line in zip(chunk, lines):
            accession, _, window = line.partition(" ")
            if accession.split(".")[0] != item["accession_version"].split(".")[0]:
                raise ContractError("Batch retrieval order or identity mismatch for " + item["accession_version"])
            item.update(window=window.upper(), status="retrieved")


def orient_and_classify(item: dict, family: dict, elements: dict, catalogues: dict, taxonomy, flank: int = 50) -> dict:
    if item["status"] != "retrieved":
        return dict(item, census_class=item["status"], census_class_hamming1=item["status"], upstream="", downstream="")
    start = max(1, item["lower"] - flank)
    left = item["window"][:item["lower"] - start]
    core = item["window"][item["lower"] - start:item["upper"] - start + 1]
    right = item["window"][item["upper"] - start + 1:]
    if item["strand"] == "-":
        left, core, right = reverse_complement(right), reverse_complement(core), reverse_complement(left)
    if core != item["query"]:
        return dict(item, status="integrity_failed", census_class="integrity_failed", census_class_hamming1="integrity_failed",
                    upstream=left, downstream=right)
    if family["orientation"] == "reverse_complement":
        left, right = reverse_complement(right), reverse_complement(left)

    def lineage(taxon):
        return in_lineage(item["taxids"], taxonomy, taxon)

    return dict(item, status="verified", upstream=left, downstream=right,
                census_class=classify(family, left, right, elements, catalogues, lineage, 0),
                census_class_hamming1=classify(family, left, right, elements, catalogues, lineage, 1))


def run_census(analysis: Path, panel: Path, prefix: Path, blastdbcmd: str, specification: dict, elements: dict[str, str],
               catalogue: Path, destination: Path, taxonomy_directory: Path | None = None) -> dict:
    """Census of the strict reference windows of `analysis` for the families of `specification`; `elements` maps the
    flank element names of the class rules to sequences and `catalogue` (TSV set/name/sequence) holds the barcode sets."""
    validate_specification(specification)
    catalogues = defaultdict(set)
    for row in read_table(catalogue, CATALOGUE_SCHEMA):
        catalogues[row["set"]].add(row["sequence"].upper())
    used = specification_references(specification)
    missing = sorted(used["elements"] - set(elements)) + sorted(used["catalogues"] - set(catalogues))
    if missing:
        raise ContractError("Census rules name flank elements or catalogue sets that were not supplied: " + ", ".join(missing))
    if any(set(sequence) - set("ACGT") or not sequence for sequence in elements.values()):
        raise ContractError("Flank elements must be nonempty ACGT sequences")
    if used["taxa"] and taxonomy_directory is None:
        raise ContractError("Census rules with taxon tests need a pinned taxonomy snapshot")
    taxonomy = None
    if taxonomy_directory is not None:
        from .taxonomy import load_taxonomy
        taxonomy = load_taxonomy(taxonomy_directory)[0]
    families = {family["name"]: family for family in specification["families"]}
    flank = specification["flank"]
    requests = window_requests(analysis, panel, specification)
    destination.mkdir(parents=True, exist_ok=False)
    scratch = destination / "retrieval"
    scratch.mkdir()
    retrieve(requests, prefix, blastdbcmd, scratch, flank)
    rows = [orient_and_classify(item, families[item["family"]], elements, catalogues, taxonomy, flank) for item in requests]
    windows = [dict(family=r["family"], physical_id=r["physical_id"], accession_version=r["accession_version"], strand=r["strand"],
                    status=r["status"], census_class=r["census_class"], census_class_hamming1=r["census_class_hamming1"],
                    upstream=r["upstream"], downstream=r["downstream"], taxids=r["taxids"], testability=r["testability"],
                    lower=r["lower"], upper=r["upper"], record_length=r["record_length"], query_length=len(r["query"]),
                    upstream_length=len(r["upstream"]), downstream_length=len(r["downstream"])) for r in rows]
    windows.sort(key=lambda r: (r["family"], r["physical_id"], r["accession_version"], r["lower"]))
    write_table(destination / "census_windows.tsv", windows, WINDOW_SCHEMA)
    summary = []
    for family in sorted({row["family"] for row in windows}):
        selected = [row for row in windows if row["family"] == family]
        family_accessions = len({row["accession_version"] for row in selected})
        for rule, field in (("exact", "census_class"), ("hamming1", "census_class_hamming1")):
            counts, accessions = Counter(), defaultdict(set)
            for row in selected:
                counts[row[field]] += 1
                accessions[row[field]].add(row["accession_version"])
            summary.extend(dict(family=family, rule=rule, census_class=name, windows=count, accessions=len(accessions[name]),
                                family_windows=len(selected), family_accessions=family_accessions) for name, count in sorted(counts.items()))
    write_table(destination / "census_summary.tsv", summary, SUMMARY_SCHEMA)
    report = dict(kind="flank_census_v1", created_utc=datetime.now(timezone.utc).isoformat(), analysis=str(analysis.resolve()),
                  analysis_sha256=digest(analysis / "analysis.json"), panel_sha256=digest(panel / "panel.json"),
                  database_prefix=str(prefix), blastdbcmd=blastdbcmd, blastdbcmd_sha256=digest(Path(blastdbcmd)),
                  catalogue_sha256=digest(catalogue), taxonomy=str(taxonomy_directory) if taxonomy_directory else None,
                  flank=flank, specification=specification,
                  elements={name: dict(length=len(sequence), sha256=hashlib.sha256(sequence.encode("ascii")).hexdigest())
                            for name, sequence in sorted(elements.items()) if name in used["elements"]},
                  catalogue_sets={name: len(values) for name, values in sorted(catalogues.items())},
                  windows=len(windows), status_counts=dict(Counter(row["status"] for row in windows)),
                  integrity_complete=all(row["status"] in {"verified", "unresolved_identifier"} for row in windows),
                  code_sha256=digest(Path(__file__)),
                  files={name: digest(destination / name) for name in ("census_windows.tsv", "census_summary.tsv")},
                  interpretation="Descriptive context of retained strict hits; classes are exact flank rules, not proof of origin")
    (destination / "census.json").write_text(json.dumps(report, indent=1) + "\n", encoding="utf-8")
    return report
