"""Record cases: one evidence-complete, sequence-verified illustration per declared category.

Candidates are the strict reference events of a completed annotated analysis, one per canonical parent, accession,
footprint and strand, with exact-version archived metadata and record eligibility. Categories are evaluated in the
declared order; within a category the candidates are sorted by the declared keys and attempted in that order (an
accession selected for an earlier category is skipped) until one retrieval verifies, within a bounded number of
attempts. Retrieval is local: up to `flank_bases` on each side with blastdbcmd, checked for accession, range, length
and the aligned subject of the accepted HSP. Failures are recorded and nothing is hand-picked. The cases illustrate
deposited sequence support; they are not biological adjudication.
"""
import csv
import hashlib
import io
import json
import os
import re
import subprocess
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO

from .contracts import ContractError, read_json, read_table, validate_record
from .context import traceback_segment_coverage, validate_flank_record
from .design import digest
from .example import save_json
from .sequences import reverse_complement

CATEGORIES = ("terminal", "synthetic_or_vector", "internal_or_ambiguous")
SORT_KEYS = ("accession_version", "parent_id", "lower", "upper", "canonical_strand")
SPECIFICATION_SCHEMA = {
    "type": "object", "additionalProperties": False,
    "required": ["version", "matching_rule", "category_order", "terminal_window", "terminal_minimum_length", "flank_bases",
                 "maximum_attempts_per_category", "maximum_selected_per_category", "sort_order", "no_repeated_selected_accession",
                 "internal_preferred_to_ambiguous", "synthetic_title_declarations", "synthetic_rule", "incomplete_or_unidentified",
                 "selection_basis"],
    "properties": {
        "version": {"const": "record_triage_v1"}, "matching_rule": {"const": "strict_single_ungapped_full_query"},
        "category_order": {"type": "array", "minItems": 1, "uniqueItems": True, "items": {"enum": list(CATEGORIES)}},
        "terminal_window": {"type": "integer", "minimum": 0}, "terminal_minimum_length": {"type": "integer", "minimum": 1},
        "flank_bases": {"type": "integer", "minimum": 0}, "maximum_attempts_per_category": {"type": "integer", "minimum": 1},
        "maximum_selected_per_category": {"const": 1},
        "sort_order": {"type": "array", "minItems": 1, "uniqueItems": True, "items": {"enum": list(SORT_KEYS)}},
        "no_repeated_selected_accession": {"const": True}, "internal_preferred_to_ambiguous": {"const": True},
        "synthetic_title_declarations": {"type": "array", "minItems": 1, "uniqueItems": True, "items": {"type": "string", "minLength": 1}},
        "synthetic_rule": {"type": "string", "minLength": 1}, "incomplete_or_unidentified": {"type": "string", "minLength": 1},
        "selection_basis": {"type": "string", "minLength": 1}}}
DEFAULT_SPECIFICATION = dict(
    version="record_triage_v1", matching_rule="strict_single_ungapped_full_query", category_order=list(CATEGORIES),
    terminal_window=500, terminal_minimum_length=10000, flank_bases=50, maximum_attempts_per_category=3,
    maximum_selected_per_category=1, sort_order=list(SORT_KEYS), no_repeated_selected_accession=True,
    internal_preferred_to_ambiguous=True,
    synthetic_title_declarations=["cloning vector", "expression vector", "synthetic construct", "synthetic DNA"],
    synthetic_rule="Archived title explicitly declares a cloning/expression vector or synthetic construct; taxon/context alone is insufficient",
    incomplete_or_unidentified="retained_in_candidate_register_not_eligible_for_verified_panel",
    selection_basis="Purposive deterministic triage illustrations; not representative biological validation")
INTERPRETATIONS = {
    "terminal": "Verified terminal sequence support in an operationally eligible linear-genomic record. Review source assembly/library "
                "information; unknown origin is not resolved and contamination is not confirmed.",
    "synthetic_or_vector": "Archived metadata explicitly declares designed vector or synthetic context. The match is compatible with "
                           "intended construct content, not automatically unwanted contamination.",
    "internal_or_ambiguous": "Verified sequence support is internal or contextually ambiguous. Endogenous similarity, shared motifs and "
                             "deposited-record processing remain alternative explanations; no biological attribution is made."}


def selection_specification(settings: dict | None = None) -> dict:
    """The deterministic case-selection rule: the defaults updated by `settings` (the study configuration), validated."""
    specification = dict(DEFAULT_SPECIFICATION, **(settings or {}))
    validate_record(specification, SPECIFICATION_SCHEMA)
    return dict(specification, new_blast_searches=0, external_retrieval_authorized=False)


def declaration_pattern(specification: dict) -> re.Pattern:
    return re.compile("^(?:" + "|".join(re.escape(text) + r"\b" for text in specification["synthetic_title_declarations"]) + ")",
                      flags=re.IGNORECASE)


def category_flags(event: dict, metadata: dict | None, linear_decision: tuple[str, str], completed: bool = True,
                   specification: dict | None = None) -> dict:
    specification = specification or selection_specification()
    exact = (not event["accession_version"].startswith("unresolved:") and event.get("identity_status") == "exact_version"
             and metadata is not None and metadata.get("accessionversion") == event["accession_version"]
             and int(metadata.get("slen", -1)) == event["slen"] and completed)
    title = metadata.get("title", "") if metadata else ""
    declaration = declaration_pattern(specification).match(title)
    genomic = exact and linear_decision[0] == "included"
    window, minimum = specification["terminal_window"], specification["terminal_minimum_length"]
    return dict(evidence_complete=bool(exact),
                terminal=bool(genomic and event["slen"] >= minimum and event["distance_to_end"] <= window),
                synthetic_or_vector=bool(exact and declaration), internal=bool(genomic and event["distance_to_end"] > window),
                ambiguous=bool(exact and linear_decision[0] == "unresolved"),
                synthetic_support="archived_title_explicit_context_declaration" if declaration else "none",
                exclusion_reason="" if exact else "incomplete_or_missing_exact_version_length_metadata")


def choose_cases(candidates: list[dict], retrieve, specification: dict) -> tuple[list, list, list]:
    selected, attempts, categories, used = [], [], [], set()
    order = lambda row: tuple(row[name] for name in specification["sort_order"])
    for category in specification["category_order"]:
        if category == "internal_or_ambiguous":
            field = "internal" if any(row["internal"] and row["accession_version"] not in used for row in candidates) else "ambiguous"
        else:
            field = category
        eligible = sorted((row for row in candidates if row[field]), key=order)
        attempted_accessions, chosen = set(), None
        for rank, row in enumerate(eligible, 1):
            if row["accession_version"] in used or row["accession_version"] in attempted_accessions:
                attempts.append(dict(category=category, candidate_id=row["candidate_id"], rank=rank, attempted=False,
                                     status="skipped_duplicate_accession", reason="selected_earlier_or_attempted_within_category"))
                continue
            if len(attempted_accessions) >= specification["maximum_attempts_per_category"]:
                break
            attempted_accessions.add(row["accession_version"])
            result = retrieve(row, category, len(attempted_accessions))
            attempts.append(dict(category=category, candidate_id=row["candidate_id"], rank=rank, attempted=True,
                                 status=result["status"], reason=result.get("error", "")))
            if result["status"] == "sequence_verified_not_biologically_adjudicated":
                chosen = dict(row, category=category, selection_subcategory=field, **result)
                selected.append(chosen)
                used.add(row["accession_version"])
                break
        categories.append(dict(category=category, selection_subcategory=field, eligible_candidates=len(eligible),
                               attempted_accessions=len(attempted_accessions), selected_candidate=chosen["candidate_id"] if chosen else "",
                               status="selected" if chosen else "unavailable_no_verified_candidate"))
    return selected, attempts, categories


def verified_flanks(record, accession: str, hsp: dict, flank_bases: int = 50) -> dict:
    lower, upper = min(hsp["sstart"], hsp["send"]), max(hsp["sstart"], hsp["send"])
    start, end = max(1, lower - flank_bases), min(hsp["slen"], upper + flank_bases)
    sequence = validate_flank_record(record, accession, start, end)
    aligned = sequence[lower - start:upper - start + 1]
    oriented = aligned if hsp["sstart"] <= hsp["send"] else reverse_complement(aligned)
    if oriented != hsp["sseq"]:
        raise ContractError("Retrieved HSP interval disagrees with the accepted aligned subject")
    left, right = sequence[:lower - start], sequence[upper - start + 1:]
    minus = hsp["sstart"] > hsp["send"]
    return dict(flank_start=start, flank_end=end, forward_sequence=sequence,
                display_left_flank=reverse_complement(right) if minus else left,
                display_right_flank=reverse_complement(left) if minus else right,
                left_flank_clipped=start == 1 and lower - 1 < flank_bases,
                right_flank_clipped=end == hsp["slen"] and hsp["slen"] - upper < flank_bases,
                display_orientation="-" if minus else "+", hsp_subject_sequence=oriented,
                lower=lower, upper=upper, sequence_verified=True)


def retrieve_case(candidate: dict, hsp: dict, metadata: dict, prefix: Path, binary: Path,
                  binary_sha256: str, destination: Path, flank_bases: int = 50) -> dict:
    if digest(binary) != binary_sha256:
        raise ContractError("Pinned blastdbcmd executable changed")
    destination.mkdir(parents=True, exist_ok=False)
    lower, upper = min(hsp["sstart"], hsp["send"]), max(hsp["sstart"], hsp["send"])
    start, end = max(1, lower - flank_bases), min(hsp["slen"], upper + flank_bases)
    command = [str(binary), "-db", str(prefix), "-entry", candidate["accession_version"], "-target_only",
               "-strand", "plus", "-range", f"{start}-{end}", "-outfmt", "%f"]
    result = dict(status="not_checked", error="", command=command, candidate_id=candidate["candidate_id"],
                  hsp=hsp, metadata=metadata, database_prefix=str(prefix), executable_sha256=binary_sha256)
    try:
        process = subprocess.run(command, capture_output=True, text=True, timeout=60,
                                 env=dict(os.environ, LC_ALL="C", BLAST_USAGE_REPORT="false"))
        (destination / "retrieved.fasta").write_text(process.stdout, encoding="utf-8")
        (destination / "stderr.txt").write_text(process.stderr, encoding="utf-8")
        result.update(exit_code=process.returncode)
        if process.returncode:
            raise ContractError("Local exact-version retrieval failed: " + process.stderr.strip())
        record = SeqIO.read(io.StringIO(process.stdout), "fasta")
        result.update(verified_flanks(record, candidate["accession_version"], hsp, flank_bases))
        result["status"] = "sequence_verified_not_biologically_adjudicated"
    except (OSError, ValueError, subprocess.SubprocessError) as error:
        result.update(status="retrieval_or_verification_failed", error=str(error))
    result["files"] = {path.name: digest(path) for path in destination.iterdir() if path.is_file()}
    save_json(destination / "case_evidence.json", result)
    return result

def write_rows(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    fields = fields or list(rows[0])
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def archived_metadata(directory: Path) -> dict[str, dict]:
    """Exact-version ESummary payloads of an archived metadata snapshot, keyed by accession.version, with their source."""
    from .eligibility import read_records

    manifest = read_json(directory / "metadata.json")
    if digest(directory / "subjects.tsv") != manifest["subject_table_sha256"]:
        raise ContractError("Metadata subject table changed")
    subjects = {row["accession_version"]: row for row in read_records(directory / "subjects.tsv")}
    records = {}
    for batch in manifest["batches"]:
        path = directory / batch["file"]
        if Path(batch["file"]).name != batch["file"] or digest(path) != batch["sha256"]:
            raise ContractError("Archived ESummary payload changed")
        result = read_json(path).get("result", {})
        for identifier in result.get("uids", []):
            payload = result[identifier]
            accession = payload.get("accessionversion")
            if accession not in subjects or int(payload.get("slen", -1)) != subjects[accession]["length"]:
                continue
            record = dict(payload=payload, source_file=batch["file"], source_sha256=batch["sha256"],
                          retrieved_utc=batch.get("retrieved_utc", ""), request_url=batch.get("url", ""))
            if accession in records and records[accession]["payload"] != payload:
                raise ContractError("Conflicting archived metadata for one accession.version")
            records[accession] = record
    if set(records) != set(subjects):
        raise ContractError("One or more normalized metadata records lack exact archived source payloads")
    return records


def trace_hsps(search: Path, events: dict[tuple, dict]) -> dict[tuple, dict]:
    """The raw HSP of every event (keyed by its source file and HSP ordinal), checked against the event fields."""
    from .search import iter_hsps, object_digest, validate_completed_batch
    from .sequences import read_fasta

    request, complete = read_json(search / "request.json"), read_json(search / "complete.json")
    wanted = defaultdict(dict)
    for (source_file, ordinal), event in events.items():
        wanted[source_file][ordinal] = event
    traced = {}
    for batch in complete["batches"]:
        folder = search / batch["directory"]
        source_file = folder.name + "/blast.tsv"
        if source_file not in wanted:
            continue
        if digest(folder / "complete.json") != batch["completion_sha256"]:
            raise ContractError("Search batch completion changed: " + batch["directory"])
        completion = validate_completed_batch(folder, request)
        queries = read_fasta(folder / "queries.fasta")
        for ordinal, hsp in enumerate(iter_hsps(folder / "blast.tsv", queries), 1):
            event = wanted[source_file].get(ordinal)
            if event is None:
                continue
            if (any(event[field] != hsp[field] for field in ("qstart", "qend", "qlen", "sstart", "send", "slen", "length", "pident"))
                    or event["raw_subject"] != hsp["sseqid"] or event["query_id"] != hsp["qseqid"]):
                raise ContractError("A raw HSP does not reproduce its analysed event")
            hsp_id = "hsp_" + object_digest(dict(source_file=source_file, source_row=ordinal, source_sha256=completion["files"]["blast.tsv"]))
            traced[source_file, ordinal] = dict(hsp_id=hsp_id, source_file=source_file, source_row=ordinal,
                                                source_sha256=completion["files"]["blast.tsv"], accession_version=event["accession_version"], **hsp)
    if set(traced) != set(events):
        raise ContractError("Not every candidate event was traced to its raw HSP")
    return traced


def candidate_dispositions(candidates: list[dict], categories: list[dict], attempts: list[dict], specification: dict) -> list[dict]:
    """Per category and candidate: eligibility, rank, selection and the reason for the outcome."""
    attempted = {(row["category"], row["candidate_id"]): row for row in attempts}
    by_id = {row["candidate_id"]: row for row in candidates}
    rows, used = [], set()
    for category in categories:
        name, field, chosen = category["category"], category["selection_subcategory"], category["selected_candidate"]
        eligible = sorted((row for row in candidates if row[field]), key=lambda row: tuple(row[key] for key in specification["sort_order"]))
        ranks = {row["candidate_id"]: index for index, row in enumerate(eligible, 1)}
        for candidate in candidates:
            identifier, record = candidate["candidate_id"], attempted.get((name, candidate["candidate_id"]))
            if not candidate[field]:
                state = "criteria_not_met"
                reason = candidate["exclusion_reason"] or candidate["eligibility_reason"] or "category_geometry_or_explicit_context_not_met"
            elif candidate["accession_version"] in used:
                state, reason = "skipped", "accession_selected_in_previous_category"
            elif record:
                state, reason = record["status"], record["reason"] or "first_verified_candidate_in_declared_order"
            elif chosen:
                state, reason = "not_attempted", "selection_stopped_after_first_verified_candidate"
            else:
                state, reason = "not_attempted", f"predeclared_{specification['maximum_attempts_per_category']}_attempt_budget_exhausted"
            rows.append(dict(category=name, candidate_id=identifier, accession_version=candidate["accession_version"],
                             eligible=bool(candidate[field]), rank_in_category=ranks.get(identifier, 0), selected=identifier == chosen,
                             status=state, reason=reason))
        if chosen:
            used.add(by_id[chosen]["accession_version"])
    return rows


def select_record_cases(analysis: Path, metadata: Path, search: Path, physical: dict[str, str], segments: list[dict],
                        acceptance: Path, blastdbcmd: Path, destination: Path, specification: dict) -> dict:
    """Candidates, deterministic selection, local flank retrieval and verification of the record cases of one analysis.

    `physical` maps physical query IDs to sequences; `segments` holds the definition segments of every physical query
    (query_id, segment_index, kind, start, end, element, verification)."""
    from .acceptance import validate_acceptance
    from .analysis import orientation_key, read_events
    from .eligibility import event_eligibility
    from .panels import LOGICAL_SCHEMA
    from .search import STATUS_SCHEMA, object_digest

    report = read_json(analysis / "analysis.json")
    for name, checksum in report["files"].items():
        if Path(name).name != name or digest(analysis / name) != checksum:
            raise ContractError("Analysis artifact changed: " + name)
    if Path(report["source_search"]).resolve() != search.resolve() or digest(search / "complete.json") != report["search_completion_sha256"]:
        raise ContractError("The search is not the source search of the analysis")
    if digest(metadata / "metadata.json") != report.get("metadata_evidence_sha256"):
        raise ContractError("The metadata snapshot is not the one bound to the analysis")
    destination.mkdir(parents=True, exist_ok=False)
    save_json(destination / "case_selection_specification.json",
              dict(specification, source_analysis_sha256=digest(analysis / "analysis.json"),
                   source_metadata_sha256=digest(metadata / "metadata.json")))
    records = archived_metadata(metadata)
    lineage = defaultdict(list)
    if "accession_taxonomy.tsv" in report["files"]:
        with (analysis / "accession_taxonomy.tsv").open(encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                lineage[row["accession_version"]].append(f"{row['rank']}:{row['status']}")
    logical = read_table(analysis / "logical_queries.tsv", LOGICAL_SCHEMA)
    aliases = defaultdict(list)
    for row in logical:
        if row["role"] == "real":
            aliases[row["parent_id"]].append(row)
    events = read_events(analysis, report)
    selected = {(row["source_file"], row["source_row"]): row for row in events if row["query_id"] in aliases and row["strict"] == "yes"}
    traced = trace_hsps(search, selected)
    statuses = {row["query_id"]: row for row in read_table(analysis / "search_status.tsv", STATUS_SCHEMA)}
    decisions = event_eligibility(events, ["linear_genomic"])
    grouped = defaultdict(list)
    for key, event in selected.items():
        grouped[orientation_key(dict(event, parent_id=event["query_id"]))].append(key)
    candidates, sources = [], {}
    for event_key, keys in sorted(grouped.items()):
        source = min(keys)
        event, hsp = selected[source], traced[source]
        record = records.get(event["accession_version"])
        decision = decisions[event["accession_version"], "linear_genomic"]
        flags = category_flags(event, record["payload"] if record else None, decision,
                               statuses[event["query_id"]]["status"] in {"success_hits", "success_zero"}, specification)
        identifier = "candidate_" + object_digest(dict(event_key=list(event_key)))[:20]
        members = sorted(aliases[event["query_id"]], key=lambda row: row["logical_id"])
        candidates.append(dict(candidate_id=identifier, parent_id=event["query_id"], accession_version=event["accession_version"],
                               lower=min(event["sstart"], event["send"]), upper=max(event["sstart"], event["send"]),
                               canonical_strand=event["canonical_strand"], record_length=event["slen"], distance_to_end=event["distance_to_end"],
                               representative_hsp_id=hsp["hsp_id"], supporting_hsp_ids=";".join(traced[key]["hsp_id"] for key in sorted(keys)),
                               physical_aliases=";".join(row["logical_id"] for row in members),
                               groups=";".join(sorted({row["group"] for row in members})), record_type=event["record_type"],
                               molecule=event["molecule"], topology=event["topology"], origin=event["origin"], context=event["context"],
                               lineage_status=";".join(lineage.get(event["accession_version"], [])) or "unknown",
                               taxids=event["record_taxids"], organism=hsp["sscinames"], organism_source="raw_BLAST_sscinames",
                               title=record["payload"].get("title", "") if record else event["title"],
                               metadata_source=record["source_file"] if record else "unresolved", metadata_sha256=record["source_sha256"] if record else "",
                               metadata_retrieved_utc=record["retrieved_utc"] if record else "unknown",
                               linear_eligibility=decision[0], eligibility_reason=decision[1], **flags))
        sources[identifier] = dict(hsp=hsp, metadata=record, aliases=members)
    write_rows(destination / "case_candidates.tsv", candidates, None if candidates else ["candidate_id", "parent_id", "accession_version"])
    accepted = read_json(acceptance)
    evidence = search / "database_evidence.json"
    prefix = Path(accepted["database_prefix"])
    if digest(blastdbcmd) != accepted["binaries"]["blastdbcmd"]["sha256"]:
        raise ContractError("blastdbcmd differs from the executable recorded by the database acceptance")
    validate_acceptance(acceptance, prefix, evidence, accepted["binaries"]["blastn"]["path"])
    binary_sha256 = digest(blastdbcmd)

    def retrieve(candidate, category, attempt):
        source = sources[candidate["candidate_id"]]
        folder = destination / "cases" / f"{category}_attempt_{attempt}"
        result = retrieve_case(candidate, source["hsp"], source["metadata"], prefix, blastdbcmd, binary_sha256, folder,
                               specification["flank_bases"])
        result["evidence_file"] = str(folder.relative_to(destination) / "case_evidence.json")
        return result

    chosen, attempts, categories = choose_cases(candidates, retrieve, specification)
    validate_acceptance(acceptance, prefix, evidence, accepted["binaries"]["blastn"]["path"])
    by_query = defaultdict(list)
    for row in segments:
        by_query[row["query_id"]].append(dict(row, segment_id=f"{row['query_id']}:{row['segment_index']}:{row['kind']}", role=row["kind"]))
    case_rows, case_segments = [], []
    for index, case in enumerate(chosen, 1):
        query = sources[case["candidate_id"]]["aliases"][0]
        hsp = case["hsp"]
        query_segments = by_query[query["logical_id"]]
        if not query_segments:
            raise ContractError("No definition segments for physical query " + query["logical_id"])
        step = 1 if hsp["qend"] >= hsp["qstart"] else -1
        positions = [hsp["qstart"] + offset * step for offset in range(len(hsp["qseq"]))]
        positions = positions if query["orientation"] == "+" else [hsp["qlen"] - position + 1 for position in positions]
        roles = [next(segment["role"] for segment in query_segments if segment["start"] <= position <= segment["end"]) for position in positions]
        sequence = physical[query["logical_id"]]
        case_id = f"case_{index:02d}"
        case_rows.append(dict(case_id=case_id, category=case["category"], selection_subcategory=case["selection_subcategory"],
                              candidate_id=case["candidate_id"], accession_version=case["accession_version"], query_id=query["logical_id"],
                              parent_id=query["parent_id"], physical_aliases=case["physical_aliases"], group=query["group"],
                              condition=query["condition"], query_sequence=sequence,
                              query_sequence_sha256=hashlib.sha256(sequence.encode("ascii")).hexdigest(),
                              alias_orientation=query["orientation"], hsp_id=hsp["hsp_id"], qstart=hsp["qstart"], qend=hsp["qend"],
                              sstart=hsp["sstart"], send=hsp["send"], record_length=case["record_length"],
                              distance_to_end=case["distance_to_end"], qseq=hsp["qseq"], sseq=hsp["sseq"], btop=hsp["btop"],
                              pident=hsp["pident"], span=hsp["length"], evalue=hsp["evalue"], bitscore=hsp["bitscore"],
                              display_segment_roles=";".join(roles), display_physical_positions=";".join(map(str, positions)),
                              **{key: case[key] for key in ("lower", "upper", "flank_start", "flank_end", "forward_sequence",
                                                            "display_left_flank", "display_right_flank", "display_orientation",
                                                            "left_flank_clipped", "right_flank_clipped", "molecule", "topology",
                                                            "record_type", "origin", "context", "lineage_status", "taxids", "organism",
                                                            "organism_source", "title", "metadata_retrieved_utc", "linear_eligibility",
                                                            "eligibility_reason", "evidence_file")},
                              interpretation=INTERPRETATIONS[case["category"]], biological_adjudication="not_performed",
                              source_verification=query_segments[0]["verification"]))
        for segment in query_segments:
            case_segments.append(dict(case_id=case_id, query_id=query["logical_id"], segment_id=segment["segment_id"],
                                      start=segment["start"], end=segment["end"], role=segment["role"], element=segment["element"],
                                      **traceback_segment_coverage(hsp, segment, query["orientation"])))
    write_rows(destination / "record_cases.tsv", case_rows, None if case_rows else ["case_id", "category", "accession_version", "status"])
    write_rows(destination / "record_case_segments.tsv", case_segments, None if case_segments else ["case_id", "segment_id", "role"])
    write_rows(destination / "case_selection_attempts.tsv", attempts, ["category", "candidate_id", "rank", "attempted", "status", "reason"])
    write_rows(destination / "case_category_status.tsv", categories)
    write_rows(destination / "case_candidate_disposition.tsv", candidate_dispositions(candidates, categories, attempts, specification),
               ["category", "candidate_id", "accession_version", "eligible", "rank_in_category", "selected", "status", "reason"])
    save_json(destination / "case_selection.json", dict(specification_sha256=digest(destination / "case_selection_specification.json"),
                                                       candidate_universe_sha256=digest(destination / "case_candidates.tsv"),
                                                       categories=categories, attempts=attempts,
                                                       selected_candidate_ids=[row["candidate_id"] for row in chosen]))
    result = dict(kind="record_cases_v1", analysis_sha256=digest(analysis / "analysis.json"), metadata_sha256=digest(metadata / "metadata.json"),
                  search_completion_sha256=digest(search / "complete.json"), database_acceptance_sha256=digest(acceptance),
                  counts=dict(strict_reference_hsps=len(traced), candidate_events=len(candidates),
                              evidence_complete=sum(row["evidence_complete"] for row in candidates), selected_cases=len(chosen)),
                  category_status=categories, local_retrieval_attempts=sum(row["attempted"] for row in attempts),
                  new_blast_searches=0, external_requests=0,
                  files={str(path.relative_to(destination)): digest(path) for path in sorted(destination.rglob("*")) if path.is_file()})
    (destination / "record_cases.json").write_text(json.dumps(result, indent=1) + "\n", encoding="utf-8")
    return result
