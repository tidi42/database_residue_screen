"""Result bundle: one analysable search from historical and new completed executions.

The bundle combines the original executions behind a historical fixed-threshold
search (classified a priori by a hash-bound calibration migration) with new
searches that carry their runtime calibration, for example a corrected construct
panel and a bridge set of re-executed historical queries. Nothing is re-run and
no historical file is changed: every bundle batch is a verified projection of
completed source query blocks.

Sources must be scientifically equivalent: the same accepted database, inventory,
BLAST binary, parser and search settings. Environment variables may differ only
in telemetry (BLAST_USAGE_REPORT). Per-length calibration rows must be identical
wherever historical and new lengths overlap, and every query searched by more
than one source must have identical alignments (the bridge check), otherwise the
assembly stops.
"""
from pathlib import Path

from .contracts import ContractError, read_json, read_table
from .design import digest
from .search import STATUS_SCHEMA
from .subset_search import BUNDLE_CONTRACT, panel_queries, source_records, write_projection
from .testability_migration import completed_search, validate_migration

EQUIVALENCE_FIELDS = ("mode", "database_acceptance", "database_prefix", "database_evidence_sha256", "database_metadata_conflicts",
                      "database_identity_verified", "executable", "executable_version", "executable_sha256", "settings", "parser_sha256")
TELEMETRY_ONLY_ENVIRONMENT = {"BLAST_USAGE_REPORT"}
CALIBRATION_HEADER = ("method", "threshold_mode", "maximum_evalue", "informativeness_limit", "testable_chance_limit", "chance_model")


def science_environment(request: dict) -> dict:
    return {key: value for key, value in request["environment"].items() if key not in TELEMETRY_ONLY_ENVIRONMENT}


def require_equivalent_fields(left: dict, right: dict) -> None:
    differing = [name for name in EQUIVALENCE_FIELDS if left.get(name) != right.get(name)]
    if differing or science_environment(left) != science_environment(right) or left.get("study_binding") != right.get("study_binding"):
        raise ContractError("Bundle sources differ in scientific contract: " + ", ".join(differing or ["environment/study binding"]))


def require_equivalent_science(bundle: dict, original: dict) -> None:
    """Batch-level check used when validating a bundle projection against its original execution."""
    require_equivalent_fields(bundle, original)
    provenance = bundle["calibration_provenance"]
    rows = {row["query_length"]: row for row in bundle["length_calibration"]["lengths"]}
    if original.get("length_calibration"):
        calibration = original["length_calibration"]
        if (any(calibration[key] != bundle["length_calibration"][key] for key in CALIBRATION_HEADER)
                or original.get("calibration_code_sha256") != bundle["calibration_code_sha256"]
                or any(rows.get(row["query_length"]) != row for row in calibration["lengths"])):
            raise ContractError("Bundle calibration disagrees with a runtime-calibrated source")
    elif object_sha(original) not in provenance["migration"]["covered_source_request_sha256"]:
        raise ContractError("Uncalibrated source is not covered by the bundle's calibration migration")


def object_sha(request: dict) -> str:
    from .search import object_digest
    return object_digest(request)


def assemble_bundle(panel: Path, historical: Path, migration: Path, new_searches: list[Path], destination: Path,
                    protocol: Path | None = None) -> dict:
    queries = panel_queries(panel)
    historical_request, historical_complete = completed_search(historical)
    migrated = validate_migration(migration, historical, historical_request, historical_complete)
    if historical_request.get("assembly_contract"):
        originals = []
        for entry in historical_request["source_searches"]:
            directory = Path(entry["directory"])
            if digest(directory / "request.json") != entry["request_sha256"] or read_json(directory / "request.json") != entry["request"]:
                raise ContractError("Historical original execution changed since its assembly")
            originals.append(directory)
    else:
        originals = [historical]
    calibrated = []
    for search in new_searches:
        request, _ = completed_search(search)
        calibration = request.get("length_calibration")
        if not calibration or calibration["threshold_mode"] != "fixed" or calibration["maximum_evalue"] is not None:
            raise ContractError("New bundle sources need a fixed-threshold runtime calibration: " + str(search))
        if any(calibration[key] != migrated["length_calibration"][key] for key in CALIBRATION_HEADER):
            raise ContractError("New source calibration header differs from the historical migration")
        if request.get("calibration_code_sha256") != migrated["calibration_code_sha256"]:
            raise ContractError("New source calibration code differs from the migration")
        calibrated.append((search, request))
    lengths = {len(sequence) for sequence in queries.values()}
    merged, origin = {}, {}
    for row in migrated["length_calibration"]["lengths"]:
        merged[row["query_length"]], origin[row["query_length"]] = row, ["historical_migration"]
    for search, request in calibrated:
        for row in request["length_calibration"]["lengths"]:
            length = row["query_length"]
            if length in merged and merged[length] != row:
                raise ContractError(f"Calibration rows disagree for length {length}")
            merged.setdefault(length, row)
            origin.setdefault(length, []).append(str(search.resolve()))
    if not lengths <= merged.keys():
        raise ContractError(f"Bundle calibration misses lengths {sorted(lengths - merged.keys())}")
    covered_requests = [object_sha(read_json(directory / "request.json")) for directory in originals]
    template = dict(read_json(originals[0] / "request.json"))
    for name in ("reuse_source", "assembly_contract", "source_searches", "assembler_sha256", "execution_scope"):
        template.pop(name, None)
    template.update(length_calibration=dict(calibrated[0][1]["length_calibration"] if calibrated else migrated["length_calibration"],
                                            lengths=[merged[length] for length in sorted(lengths)]),
                    calibration_code_sha256=migrated["calibration_code_sha256"],
                    calibration_provenance=dict(
                        migration=dict(path=str(migration.resolve()), sha256=digest(migration), source_search=str(historical.resolve()),
                                       covered_source_request_sha256=covered_requests),
                        runtime_sources=[dict(directory=str(search.resolve()), request_sha256=digest(search / "request.json"))
                                         for search, _ in calibrated],
                        per_length={str(length): origin[length] for length in sorted(lengths)}))

    def check(left, right):
        require_equivalent_fields(left, right)

    for search, request in calibrated:
        check(template, request)
    requests, entries, covered = source_records(originals + [search for search, _ in calibrated], queries, check)
    if covered != queries.keys():
        raise ContractError(f"Bundle panel has {len(queries.keys() - covered)} queries without completed evidence")
    overlap = sum(len(entry["query_ids"]) for entry in entries)
    request = dict(template, queries_sha256=digest(panel / "search_queries.fasta"), query_count=len(queries),
                   assembly_contract=BUNDLE_CONTRACT, source_searches=requests, assembler_sha256=digest(Path(__file__)),
                   science_equivalence=dict(fields=list(EQUIVALENCE_FIELDS), telemetry_only_environment=sorted(TELEMETRY_ONLY_ENVIRONMENT),
                                            overlapping_queries="identical alignment multisets required; earlier source owns the query"),
                   protocol_sha256=digest(protocol) if protocol else None,
                   execution_scope="Query-level bundle of completed historical and new executions; source requests keep batch sizes, "
                                   "threads and times; no BLAST was re-run for the bundle")
    if overlap != len(queries):
        raise ContractError("Bundle query ownership is not one-to-one")
    report = write_projection(destination, request, requests, entries, queries,
                              "Bundle of verified completed executions; retained results, not exhaustive retrieval")
    for batch in report["batches"]:
        directory = destination / batch["directory"]
        source = read_json(directory / "source_batch.json")
        original = Path(requests[source["source_index"]]["directory"]) / source["batch_directory"]
        projected = {row["query_id"]: row for row in read_table(directory / "status.tsv", STATUS_SCHEMA)}
        recorded = {row["query_id"]: row for row in read_table(original / "status.tsv", STATUS_SCHEMA)}
        if any(dict(recorded[name], batch_id=row["batch_id"]) != row for name, row in projected.items()):
            raise ContractError("Current parser disagrees with the status recorded by the original execution: " + str(original))
    return report
