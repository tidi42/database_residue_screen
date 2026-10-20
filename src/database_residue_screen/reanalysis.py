import csv
from collections import Counter, defaultdict
from pathlib import Path

from .contracts import ContractError, read_json, read_table, schema, table_schema, write_table
from .design import digest
from .eligibility import ANNOTATION_FIELDS, annotation_defaults, record_eligibility
from .example import save_json
from .legacy import REAL_BLAST_FILES
from .matching import parse_legacy_hsp, strict_hsp, tolerant_hsp
from .panels import LOGICAL_SCHEMA
from .outcomes import CONTRAST_SCHEMA, ELIGIBILITY_SCHEMA, MEMBERSHIP_SCHEMA, OUTCOME_SCHEMA, build_outcomes, project_real_events
from .policy import applicability, load_policy, primary_event
from .search import STATUS_SCHEMA, iter_hsps, object_digest, validate_completed_batch
from .sequences import accession_version, canonical, canonical_id, read_fasta


EVENT_SCHEMA = table_schema(("source_file", "query_id", "parent_id", "group", "condition", "raw_subject", "accession_version",
                             "canonical_strand", "strict", "tolerant", "record_type", "topology", "metadata_status", "taxids", "title"),
                            ("source_row", "qlen", "qstart", "qend", "sstart", "send", "slen", "length", "distance_to_end"), ("pident",))
QUALITY_SCHEMA = table_schema(("group",), ("identity_bin", "coverage_bin", "hsp_rows"))
QUERY_SCHEMA = table_schema(("query_id", "parent_id", "group", "condition", "construct_class", "provenance_status", "completion"),
                            ("length", "retained_hsp_rows", "legacy_perfect_rows", "strict_rows", "tolerant_rows", "strict_accessions", "tolerant_accessions"))
GROUP_SCHEMA = table_schema(("group",), ("physical_queries", "condition_canonical_units", "global_canonical_union", "legacy_positive_ids",
                                        "strict_positive_ids", "tolerant_positive_ids", "strict_rows", "tolerant_rows", "strict_accessions", "strict_orientation_deduplicated_events"))
DENSITY_SCHEMA = table_schema(("group", "accession_version"), ("strict_raw_hsps", "strict_orientation_deduplicated_events"))
POSITION_SCHEMA = table_schema(("group", "population", "fold_n2"), ("window", "minimum_length", "events", "terminal_events"), ("expected_n2",))
OVERLAP_SCHEMA = table_schema(("groups",), ("accessions",))
ANNOTATED_EVENT_SCHEMA = dict(EVENT_SCHEMA, properties=dict(EVENT_SCHEMA["properties"], **{name: {"type": "string"} for name in ANNOTATION_FIELDS}),
                              required=list(EVENT_SCHEMA["properties"]) + list(ANNOTATION_FIELDS))


def read_events(directory: Path, report: dict) -> list[dict]:
    return read_table(directory / "events.tsv", ANNOTATED_EVENT_SCHEMA if report.get("event_contract") == "annotated_v1" else EVENT_SCHEMA)


def event(row: dict, metadata: dict, sequence: str, source: str, source_row: int) -> dict:
    query_direction = 1 if row["qend"] >= row["qstart"] else -1
    subject_direction = 1 if row["send"] >= row["sstart"] else -1
    canonical_direction = 1 if sequence == canonical(sequence) else -1
    try:
        accession = accession_version(row["sseqid"])
    except ContractError:
        accession = "unresolved:" + row["sseqid"]
    return dict(source_file=source, source_row=source_row, query_id=row["qseqid"], parent_id=canonical_id(sequence),
                group=metadata.get("group", ""), condition=metadata.get("condition", ""), raw_subject=row["sseqid"], accession_version=accession,
                canonical_strand="+" if query_direction * subject_direction * canonical_direction == 1 else "-",
                strict="yes" if strict_hsp(row) else "no", tolerant="yes" if tolerant_hsp(row) else "no",
                record_type="unknown", topology="unknown", metadata_status="not_retrieved", taxids=row["staxids"], title=row["stitle"],
                **{key: row[key] for key in ("qlen", "qstart", "qend", "sstart", "send", "slen", "length", "pident")},
                distance_to_end=min(row["sstart"] - 1, row["send"] - 1, row["slen"] - max(row["sstart"], row["send"])))


def orientation_key(row: dict) -> tuple:
    return row["parent_id"], row["accession_version"], min(row["sstart"], row["send"]), max(row["sstart"], row["send"]), row["canonical_strand"]


def position_summary(events: list[dict], group: str, window: int, minimum_length: int, primary: bool = False) -> dict:
    selected = [row for row in events if row["group"] == group and row["strict"] == "yes" and row["slen"] >= minimum_length
                and (not primary or (record_eligibility(row, "linear_genomic")[0] == "included" if "identity_status" in row else
                                     row["topology"] == "linear" and row["record_type"] == "genomic"))]
    selected = list({orientation_key(row): row for row in selected}.values())
    expected = sum(min(1.0, 2 * (window + 1) / (row["slen"] - row["length"] + 1)) for row in selected)
    observed = sum(row["distance_to_end"] <= window for row in selected)
    return dict(group=group, population="eligible_linear_genomic" if primary else "all_record_ends_geometry_only",
                window=window, minimum_length=minimum_length, events=len(selected), terminal_events=observed,
                expected_n2=expected, fold_n2=str(observed / expected) if expected else "not_estimated")


def derive_event_tables(events: list[dict], groups: list[str], destination: Path, prefix: str = "", policy: dict | None = None) -> dict:
    # These tables summarise strict support. Keep the complete retained inventory
    # in events.tsv/real_events.tsv, but apply the same a-priori gate as outcomes
    # before building primary event, accession and positional summaries.
    if policy and "detectability" in policy:
        events = [row for row in events if applicability(row["qlen"], "strict", policy) == "operational_retained_search"]
    if policy:
        events = [row for row in events if primary_event(row, policy)]
    density, unique = Counter(), defaultdict(set)
    memberships = defaultdict(set)
    for row in events:
        if row["strict"] != "yes":
            continue
        key = row["group"], row["accession_version"]
        density[key] += 1
        unique[key].add(orientation_key(row))
        if not row["accession_version"].startswith("unresolved:"):
            memberships[row["accession_version"]].add(row["group"])
    write_table(destination / (prefix + "density.tsv"), [dict(group=key[0], accession_version=key[1], strict_raw_hsps=count,
                                                   strict_orientation_deduplicated_events=len(unique[key])) for key, count in sorted(density.items())], DENSITY_SCHEMA)
    windows = policy["position"]["windows"] if policy else (100, 500, 1000)
    minimum_lengths = policy["position"]["minimum_lengths"] if policy else (1000, 10000, 50000)
    positions = [position_summary(events, group, window, minimum, primary) for group in groups
                 for window in windows for minimum in minimum_lengths for primary in (False, True)]
    write_table(destination / (prefix + "positions.tsv"), positions, POSITION_SCHEMA)
    overlaps = Counter(tuple(sorted(members)) for members in memberships.values())
    overlap_rows = []
    for mask in range(1, 2 ** len(groups)):
        selected = tuple(sorted(group for index, group in enumerate(groups) if mask & (1 << index)))
        overlap_rows.append(dict(groups=";".join(selected), accessions=overlaps[selected]))
    write_table(destination / (prefix + "overlap.tsv"), overlap_rows, OVERLAP_SCHEMA)
    return dict(strict_distinct_accessions=len(memberships), unknown_metadata_events=sum(row["metadata_status"] == "not_retrieved" for row in events),
                biological_event_collapse=False, orientation_deduplication="canonical_query_accession_footprint_strand",
                positional_inference="No inference on unknown topology/type; geometric diagnostics are separate")


def reanalyse_legacy(root: Path, query_fasta: Path, query_metadata: Path, destination: Path) -> dict:
    sequences = read_fasta(query_fasta)
    metadata = {row["query_id"]: row for row in read_table(query_metadata, schema("queries"))}
    sources = [root / "blast" / name for name in REAL_BLAST_FILES]
    source_hashes = {str(path.relative_to(root)): digest(path) for path in sources}
    query_rows = {identifier: dict(query_id=identifier, parent_id=canonical_id(sequence), group=metadata[identifier]["group"],
                                  condition=metadata[identifier]["condition"], construct_class=metadata[identifier]["construct_class"],
                                  provenance_status=metadata[identifier]["provenance_status"], completion="unknown_not_zero", length=len(sequence),
                                  retained_hsp_rows=0, legacy_perfect_rows=0, strict_rows=0, tolerant_rows=0,
                                  strict_accessions=0, tolerant_accessions=0) for identifier, sequence in sequences.items()
                  if metadata[identifier]["role"] != "reference"}
    quality, strict_accessions, tolerant_accessions = Counter(), defaultdict(set), defaultdict(set)
    events, excluded = [], Counter()
    for path in sources:
        with path.open(encoding="utf-8", newline="") as handle:
            for line_number, values in enumerate(csv.reader(handle, delimiter="\t"), 1):
                if not values:
                    continue
                row = parse_legacy_hsp(values)
                identifier = row["qseqid"]
                if identifier not in query_rows:
                    excluded[path.name] += 1
                    continue
                if row["qlen"] != len(sequences[identifier]):
                    raise ContractError("Historical HSP query length disagrees with the current FASTA")
                query = query_rows[identifier]
                query["retained_hsp_rows"] += 1
                query["legacy_perfect_rows"] += row["pident"] == 100 and row["qcovs"] == 100 and row["mismatch"] == 0
                strict, tolerant = strict_hsp(row), tolerant_hsp(row)
                query["strict_rows"] += strict
                query["tolerant_rows"] += tolerant
                coverage = 100 * (abs(row["qend"] - row["qstart"]) + 1) / row["qlen"]
                quality[query["group"], min(100, int(row["pident"])), min(100, int(coverage))] += 1
                if not tolerant:
                    continue
                accepted = event(row, metadata[identifier], sequences[identifier], str(path.relative_to(root)), line_number)
                events.append(accepted)
                tolerant_accessions[identifier].add(accepted["accession_version"])
                if strict:
                    strict_accessions[identifier].add(accepted["accession_version"])
    for identifier, row in query_rows.items():
        row["strict_accessions"] = len(strict_accessions[identifier])
        row["tolerant_accessions"] = len(tolerant_accessions[identifier])
    destination.mkdir(parents=True, exist_ok=False)
    write_table(destination / "events.tsv", events, EVENT_SCHEMA)
    write_table(destination / "queries.tsv", list(query_rows.values()), QUERY_SCHEMA)
    write_table(destination / "quality.tsv", [dict(group=key[0], identity_bin=key[1], coverage_bin=key[2], hsp_rows=count)
                                               for key, count in sorted(quality.items())], QUALITY_SCHEMA)
    groups = sorted({row["group"] for row in query_rows.values()})
    summaries = []
    for group in groups:
        selected = [row for row in query_rows.values() if row["group"] == group]
        strict_events = [row for row in events if row["group"] == group and row["strict"] == "yes"]
        summaries.append(dict(group=group, physical_queries=len(selected), condition_canonical_units=len({(row["condition"], row["parent_id"]) for row in selected}),
                              global_canonical_union=len({row["parent_id"] for row in selected}),
                              legacy_positive_ids=sum(row["legacy_perfect_rows"] > 0 for row in selected),
                              strict_positive_ids=sum(row["strict_rows"] > 0 for row in selected), tolerant_positive_ids=sum(row["tolerant_rows"] > 0 for row in selected),
                              strict_rows=sum(row["strict_rows"] for row in selected), tolerant_rows=sum(row["tolerant_rows"] for row in selected),
                              strict_accessions=len({row["accession_version"] for row in strict_events}),
                              strict_orientation_deduplicated_events=len({orientation_key(row) for row in strict_events})))
    write_table(destination / "groups.tsv", summaries, GROUP_SCHEMA)
    derived = derive_event_tables(events, groups, destination)
    if source_hashes != {str(path.relative_to(root)): digest(path) for path in sources}:
        raise ContractError("Historical source changed during reanalysis")
    report = dict(stage="reclassified_retained_output", final_freeze=False, completion_known=False,
                  matched_null_analysed=False, biological_authentication=False, source_hashes=source_hashes,
                  source_fasta_sha256=digest(query_fasta), source_metadata_sha256=digest(query_metadata),
                  parser_sha256=digest(Path(__file__).with_name("matching.py")), analysis_sha256=digest(Path(__file__)),
                  groups=summaries, excluded_noncurrent_rows=dict(excluded), derived=derived,
                  files={path.name: digest(path) for path in destination.glob("*.tsv")})
    save_json(destination / "analysis.json", report)
    return report


def write_completed_analysis(events: list[dict], logical: list[dict], statuses: list[dict], policy: dict,
                             destination: Path, provenance: dict, extra_tables: dict | None = None) -> dict:
    outcomes, contrasts, memberships, eligibility = build_outcomes(events, logical, statuses, policy, provenance["draws_per_parent"],
                                                                   provenance.get("repeated_control_seed_values", 0))
    real_events = project_real_events(events, memberships)
    destination.mkdir(parents=True, exist_ok=False)
    write_table(destination / "events.tsv", events, ANNOTATED_EVENT_SCHEMA)
    write_table(destination / "real_events.tsv", real_events, ANNOTATED_EVENT_SCHEMA)
    write_table(destination / "logical_queries.tsv", logical, LOGICAL_SCHEMA)
    write_table(destination / "search_status.tsv", statuses, STATUS_SCHEMA)
    write_table(destination / "outcomes.tsv", outcomes, OUTCOME_SCHEMA)
    write_table(destination / "contrasts.tsv", contrasts, CONTRAST_SCHEMA)
    from .sensitivity import (COMMON_RULE_SCHEMA, SENSITIVITY_SCHEMA, TESTABILITY_SCHEMA, common_rule_rows, length_sensitivity,
                              length_testability_rows)
    write_table(destination / "length_sensitivity.tsv", length_sensitivity(outcomes, policy, provenance["draws_per_parent"]), SENSITIVITY_SCHEMA)
    write_table(destination / "common_rule_contrasts.tsv", common_rule_rows(outcomes, policy, provenance["draws_per_parent"]), COMMON_RULE_SCHEMA)
    if policy.get("detectability"):
        write_table(destination / "length_testability.tsv", length_testability_rows(policy), TESTABILITY_SCHEMA)
    write_table(destination / "event_memberships.tsv", memberships, MEMBERSHIP_SCHEMA)
    write_table(destination / "eligibility.tsv", eligibility, ELIGIBILITY_SCHEMA)
    save_json(destination / "analysis_policy.json", policy)
    for name, (rows, definition) in (extra_tables or {}).items():
        if Path(name).name != name or not name.endswith(".tsv"):
            raise ContractError("Unsafe extra result table name")
        write_table(destination / name, rows, definition)
    groups = sorted({row["group"] for row in logical if row["role"] == "real"})
    derived = {}
    included = {(row["sequence_id"], row["accession_version"], row["population"]) for row in eligibility if row["state"] == "included"}
    for population in policy["record_populations"]:
        selected = [row for row in real_events if (row["query_id"], row["accession_version"], population) in included]
        derived[population] = derive_event_tables(selected, groups, destination, population + "_", policy)
    report = dict(provenance, stage="completed_reference_screen_analysis", event_contract="annotated_v1", result_contract="reference_screen_v2",
                  final_freeze=False, analysis_policy=policy, logical_queries=len(logical), distinct_search_sequences=len(statuses),
                  canonical_parents=len({row["parent_id"] for row in logical if row["role"] == "real"}),
                  broad_null_adequacy_established=False, derived=derived,
                  event_unit="Canonical search sequence HSP; original parent and group ownership are in event_memberships.tsv and logical_queries.tsv",
                  inference="Fixed-panel equally weighted canonical parents; measure-specific resolved matched subsets with explicit attrition",
                  applicability="Operational policy, not universal sensitivity; original short queries retained",
                  retrieval="Endpoint-specific retained-accession lower bounds, not genomes or exhaustive breadth",
                  mc_caveat="Conditional Monte Carlo errors only; zero empirical variance is not zero uncertainty; seed reuse suppresses independent-draw error claims",
                  analysis_sha256=digest(Path(__file__)), outcomes_sha256=digest(Path(__file__).with_name("outcomes.py")),
                  eligibility_sha256=digest(Path(__file__).with_name("eligibility.py")),
                  files={path.name: digest(path) for path in destination.iterdir() if path.is_file()})
    save_json(destination / "analysis.json", report)
    return report


def analyse_search(panel: Path, search: Path, destination: Path, calibration_migration: Path | None = None) -> dict:
    if (search / "quarantine.json").exists():
        raise ContractError("Search is quarantined; see its quarantine.json and do not use it for scientific results")
    definition = read_json(panel / "panel.json")
    policy = load_policy(panel / "analysis_policy.json")
    if definition.get("analysis_policy") != policy or "analysis_policy.json" not in definition["files"]:
        raise ContractError("Panel lacks a matching, hashed analysis policy; rebuild the panel")
    for name, checksum in definition["files"].items():
        if Path(name).name != name or digest(panel / name) != checksum:
            raise ContractError("Search panel artifact changed")
    queries = read_fasta(panel / "search_queries.fasta")
    logical = read_table(panel / "logical_queries.tsv", LOGICAL_SCHEMA)
    request, complete = read_json(search / "request.json"), read_json(search / "complete.json")
    if request["settings"]["word_size"] != policy["word_size"]:
        raise ContractError("Search word size disagrees with the declared applicability policy")
    if (request["queries_sha256"] != digest(panel / "search_queries.fasta") or complete["search_id"] != object_digest(request)
            or not complete["all_queries_completed"]):
        raise ContractError("Search does not completely cover this panel")
    if request.get("study_binding"):
        declared = request["study_binding"]["analysis"]
        if (policy["position"]["windows"] != declared["position_windows"] or
                policy["position"]["minimum_lengths"] != declared["minimum_record_lengths"]):
            raise ContractError("Analysis policy disagrees with study.json position parameters")
    if request.get("length_calibration"):
        if calibration_migration is not None:
            raise ContractError("Search already carries its runtime calibration; do not supply a migration")
        calibration = request["length_calibration"]
        calibration_source = (dict(kind="bundle_assembly", provenance=request["calibration_provenance"])
                              if request.get("assembly_contract") == "verified_bundle_projection_v1" else dict(kind="runtime_request"))
    elif calibration_migration is not None:
        from .testability_migration import validate_migration
        calibration = validate_migration(calibration_migration, search, request, complete)["length_calibration"]
        calibration_source = dict(kind="post_hoc_fixed_threshold_migration", path=str(calibration_migration.resolve()),
                                  sha256=digest(calibration_migration))
    else:
        raise ContractError("Search lacks length-aware calibration; supply a hash-bound migration from calibrate-search "
                            "for a fixed-threshold historical search, or rerun")
    policy = dict(policy, detectability=calibration["lengths"],
                  sensitivity_baseline_evalue=policy.get("sensitivity_baseline_evalue", request["settings"]["evalue"]))
    events, statuses, seen_queries = [], [], set()
    for batch in complete["batches"]:
        directory = search / batch["directory"]
        if directory.parent.resolve() != search.resolve() or digest(directory / "complete.json") != batch["completion_sha256"]:
            raise ContractError("Search batch completion changed")
        batch_complete = validate_completed_batch(directory, request)
        if batch_complete["search_id"] != complete["search_id"] or set(batch_complete["query_ids"]) & seen_queries:
            raise ContractError("Search batch ownership is duplicated or inconsistent")
        for name, checksum in batch_complete["files"].items():
            if Path(name).name != name or digest(directory / name) != checksum:
                raise ContractError("Search batch artifact changed")
        batch_queries = {name: queries[name] for name in batch_complete["query_ids"]}
        seen_queries.update(batch_queries)
        batch_status = read_table(directory / "status.tsv", STATUS_SCHEMA)
        if len(batch_status) != len(batch_queries) or {row["query_id"] for row in batch_status} != set(batch_queries):
            raise ContractError("Batch status does not cover its declared queries")
        statuses.extend(batch_status)
        for ordinal, row in enumerate(iter_hsps(directory / "blast.tsv", batch_queries), 1):
            if not tolerant_hsp(row):
                continue
            accepted = event(row, {}, queries[row["qseqid"]], directory.name + "/blast.tsv", ordinal)
            accepted.update(annotation_defaults(accepted))
            accepted["parent_id"] = ""
            events.append(accepted)
    if seen_queries != set(queries):
        raise ContractError("Completed search does not cover every submitted sequence")
    provenance = dict(panel_sha256=digest(panel / "panel.json"), search_completion_sha256=digest(search / "complete.json"),
                      search_mode=request["mode"],
                      search_id=complete["search_id"], source_search=str(search.resolve()), draws_per_parent=definition["draws_per_parent"],
                      repeated_control_seed_values=definition["repeated_control_seed_values"], control_algorithm=definition["control_algorithm"],
                      search_settings=request["settings"], length_calibration=calibration,
                      length_calibration_source=calibration_source, study_binding=request.get("study_binding"),
                      search_request_sha256=digest(search / "request.json"), database_acceptance=request.get("database_acceptance"))
    return write_completed_analysis(events, logical, statuses, policy, destination, provenance)
