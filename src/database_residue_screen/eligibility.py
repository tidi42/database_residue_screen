import csv
from pathlib import Path

from .contracts import ContractError, read_table, schema
from .sequences import accession_version


ANNOTATION_FIELDS = ("molecule", "origin", "targeted", "organelle", "identity_status", "blast_taxids", "record_taxids",
                     "metadata_source", "metadata_date", "context")
RECORD_PROPERTIES = {
    "molecule": {"enum": ["DNA", "RNA", "other", "unknown", "conflicting"]},
    "origin": {"enum": ["synthetic", "vector", "artificial_other", "engineered", "unknown", "conflicting"]},
    "targeted": {"enum": ["yes", "no", "unknown", "conflicting"]},
    "organelle": {"enum": ["yes", "no", "unknown", "conflicting"]},
    "context": {"enum": ["cellular", "viral", "unknown", "conflicting"]},
}


def record_schema() -> dict:
    definition = schema("subjects")
    definition["properties"].update(RECORD_PROPERTIES)
    definition["properties"]["topology"]["enum"].append("conflicting")
    definition["properties"]["record_type"]["enum"].extend(["other", "conflicting"])
    definition["required"] = list(definition["properties"])
    return definition


def read_records(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as handle:
        columns = next(csv.reader(handle, delimiter="\t"), [])
    definition = record_schema() if "molecule" in columns else schema("subjects")
    return read_table(path, definition)


def annotation_defaults(row: dict) -> dict:
    return dict(molecule="unknown", origin="unknown", targeted="unknown", organelle="unknown", context="unknown",
                identity_status="unresolved_id" if row["accession_version"].startswith("unresolved:") else "from_search",
                blast_taxids=row.get("blast_taxids", "" if row.get("metadata_status") == "exact_version_esummary" else row["taxids"]),
                record_taxids="", metadata_source="", metadata_date="")


def annotate_record(row: dict, subject: dict | None, problem: str = "") -> dict:
    result = dict(row, **annotation_defaults(row))
    result.update(record_type="unknown", topology="unknown", metadata_status="not_retrieved")
    if problem in {"length_conflict", "version_mismatch", "conflicting"}:
        result.update(identity_status=problem, metadata_status=problem)
        return result
    if subject is None:
        return result
    if subject["accession_version"] != row["accession_version"]:
        result.update(identity_status="version_mismatch", metadata_status="version_mismatch")
        return result
    if subject["length"] != row["slen"]:
        result.update(identity_status="length_conflict", metadata_status="length_conflict")
        return result
    result.update(record_type=subject["record_type"], topology=subject["topology"], identity_status="exact_version",
                  metadata_status="exact_version_record", record_taxids=subject["taxids"],
                  metadata_source=subject["metadata_source"], metadata_date=subject["metadata_date"])
    result.update({name: subject.get(name, "unknown") for name in RECORD_PROPERTIES})
    return result


def record_eligibility(row: dict, population: str) -> tuple[str, str]:
    if population not in {"all_records", "linear_genomic"}:
        raise ContractError("Unknown record population: " + population)
    try:
        if accession_version(row["accession_version"]) != row["accession_version"]:
            return "unresolved", "identifier_not_normalised"
    except ContractError:
        return "unresolved", "unresolved_identifier"
    identity = row.get("identity_status", "from_search")
    if identity in {"unresolved_id", "version_mismatch", "length_conflict", "conflicting"}:
        return "unresolved", identity
    if population == "all_records":
        return "included", "identified_deposited_record"
    for field in ("molecule", "record_type", "topology", "origin", "targeted"):
        if row.get(field) == "conflicting":
            return "unresolved", field + "_conflicting"
    if row.get("record_type") in {"transcript", "amplicon", "synthetic", "vector"}:
        return "excluded", "record_type_" + row["record_type"]
    if row.get("origin") in {"synthetic", "vector", "artificial_other", "engineered"}:
        return "excluded", "origin_" + row["origin"]
    if row.get("targeted") == "yes":
        return "excluded", "targeted_record"
    if row.get("molecule") == "RNA":
        return "excluded", "rna_record"
    if row.get("topology") == "circular":
        return "excluded", "circular_record"
    for field, expected in (("molecule", "DNA"), ("record_type", "genomic"), ("topology", "linear")):
        if row.get(field) != expected:
            return "unresolved", field + "_not_established"
    return "included", "linear_genomic_dna_known_exclusions_removed"


def event_eligibility(events: list[dict], populations: list[str]) -> dict:
    decisions, lengths = {}, {}
    length_conflicts = set()
    for row in events:
        accession = row["accession_version"]
        if "slen" in row:
            if accession in lengths and lengths[accession] != row["slen"]:
                length_conflicts.add(accession)
            lengths[accession] = row["slen"]
        for population in populations:
            key = accession, population
            decision = record_eligibility(row, population)
            decisions[key] = ("unresolved", "conflicting_record_evidence") if key in decisions and decisions[key] != decision else decision
    for accession in length_conflicts:
        for population in populations:
            decisions[accession, population] = "unresolved", "length_conflict"
    return decisions