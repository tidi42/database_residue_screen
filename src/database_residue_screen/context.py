import hashlib
import io
import subprocess
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO

from .contracts import ContractError, read_json, read_table, schema, write_table
from .design import digest
from .example import save_json
from .reanalysis import QUERY_SCHEMA, read_events, table_schema
from .search import verify_database
from .sequences import accession_version, read_fasta, reverse_complement


SEGMENT_SCHEMA = table_schema(("query_id", "component_id", "component_kind", "component_orientation", "endpoint", "mismatch_allocation"),
                              ("source_row", "component_start", "component_end", "covered_bases", "component_length", "complete_segment"))
CASE_SCHEMA = table_schema(("case_id", "group", "query_id", "accession_version", "returned_accession", "record_type", "topology", "position_stratum", "status", "error"),
                           ("sstart", "send", "slen", "flank_start", "flank_end"))


def validate_flank_record(record, expected_accession: str, start: int, end: int) -> str:
    returned = accession_version(record.id.split(":", 1)[0])
    if returned != expected_accession:
        raise ContractError(f"Accession lookup mismatch: requested {expected_accession}, returned {returned}")
    if len(record.seq) != end - start + 1:
        raise ContractError("Retrieved flank has an unexpected length")
    if record.id != f"{expected_accession}:{start}-{end}":
        raise ContractError("Retrieved FASTA header does not confirm the requested range")
    return str(record.seq).upper()


def segment_coverage(event: dict, component: dict) -> dict:
    start = max(min(event["qstart"], event["qend"]), component["start"])
    end = min(max(event["qstart"], event["qend"]), component["end"])
    covered = max(0, end - start + 1)
    length = component["end"] - component["start"] + 1
    return dict(query_id=event["query_id"], component_id=component["component_id"], component_kind=component["component_kind"],
                component_orientation=component["strand"], endpoint="strict" if event["strict"] == "yes" else "tolerant_only",
                mismatch_allocation="all_bases_exact" if event["strict"] == "yes" else "unavailable_in_legacy_traceback",
                source_row=event["source_row"], component_start=component["start"], component_end=component["end"],
                covered_bases=covered, component_length=length, complete_segment=int(covered == length))


def traceback_segment_coverage(hsp: dict, component: dict, alias_orientation: str = "+") -> dict:
    length = hsp["qlen"]
    start, end = component["start"], component["end"]
    if alias_orientation not in {"+", "-"} or not 1 <= start <= end <= length:
        raise ContractError("Invalid physical-query segment or alias orientation")
    aligned_query, aligned_subject = hsp["qseq"], hsp["sseq"]
    if (len(aligned_query) != len(aligned_subject) or "-" in aligned_query + aligned_subject
            or len(aligned_query) != abs(hsp["qend"] - hsp["qstart"]) + 1
            or not 1 <= min(hsp["qstart"], hsp["qend"]) <= max(hsp["qstart"], hsp["qend"]) <= length):
        raise ContractError("Component annotation requires a valid single ungapped accepted HSP")
    step = 1 if hsp["qend"] >= hsp["qstart"] else -1
    covered, differences = [], []
    for offset, (query_base, subject_base) in enumerate(zip(aligned_query, aligned_subject)):
        canonical_position = hsp["qstart"] + step * offset
        position = canonical_position if alias_orientation == "+" else length - canonical_position + 1
        if start <= position <= end:
            covered.append(position)
            if query_base != subject_base:
                differences.append(position)
    size, matching = len(covered), len(covered) - len(differences)
    complete = size == end - start + 1
    state = "not_covered" if not size else ("complete" if complete else "partial") + ("_exact" if not differences else "_mismatched")
    return dict(covered_bases=size, matching_bases=matching, mismatches=len(differences),
                component_length=end - start + 1, coverage_fraction=size / (end - start + 1),
                complete_segment=int(complete), support_state=state,
                covered_start=min(covered, default=0), covered_end=max(covered, default=0),
                mismatch_positions=";".join(map(str, sorted(differences))))


def context_evidence(analysis: Path, study_inputs: Path, prefix: Path, evidence: Path, destination: Path, per_stratum: int = 2) -> dict:
    if per_stratum < 1:
        raise ContractError("Case count per stratum must be positive")
    manifest = read_json(analysis / "analysis.json")
    for name in ("events.tsv", "queries.tsv"):
        if digest(analysis / name) != manifest["files"][name]:
            raise ContractError("Source analysis changed")
    verify_database(evidence, prefix)
    events = read_events(analysis, manifest)
    query_rows = {row["query_id"]: row for row in read_table(analysis / "queries.tsv", QUERY_SCHEMA)}
    components = read_table(study_inputs / "components.tsv", schema("components"))
    sequences = read_fasta(study_inputs / "queries.fasta")
    by_query = defaultdict(list)
    for component in components:
        by_query[component["composite_id"]].append(component)
    segment_rows = [segment_coverage(row, component) for row in events for component in by_query[row["query_id"]]]
    strata = defaultdict(dict)
    for row in events:
        if row["strict"] != "yes" or row["accession_version"].startswith("unresolved:"):
            continue
        position = "terminal_le500" if row["distance_to_end"] <= 500 else "internal_gt500"
        key = row["group"], query_rows[row["query_id"]]["construct_class"], row["record_type"], position
        case_key = row["query_id"], row["accession_version"], row["sstart"], row["send"]
        strata[key][case_key] = row
    selected = []
    for key, candidates in sorted(strata.items()):
        ordered = sorted(candidates, key=lambda candidate: hashlib.sha256(repr(candidate).encode("ascii")).hexdigest())
        used_accessions = set()
        for candidate in ordered:
            if candidate[1] in used_accessions:
                continue
            selected.append((key, candidates[candidate]))
            used_accessions.add(candidate[1])
            if len(used_accessions) >= per_stratum:
                break
    destination.mkdir(parents=True, exist_ok=False)
    write_table(destination / "segments.tsv", segment_rows, SEGMENT_SCHEMA)
    save_json(destination / "selection.json", dict(selection="Smallest SHA256 of query/accession/coordinates per group, construct class, record type and position stratum; distinct accessions within stratum",
                                                    per_stratum=per_stratum, strata=len(strata), selected=len(selected),
                                                    limitations="Illustrative sequence checks only; no representative biological validation rate",
                                                    analysis_sha256=digest(analysis / "analysis.json"), database_evidence_sha256=digest(evidence)))
    cases = []
    for index, (key, row) in enumerate(selected, 1):
        case_id = f"case_{index:03d}"
        lower, upper = min(row["sstart"], row["send"]), max(row["sstart"], row["send"])
        start, end = max(1, lower - 100), min(row["slen"], upper + 100)
        command = ["blastdbcmd", "-db", str(prefix.resolve()), "-entry", row["accession_version"], "-range", f"{start}-{end}", "-outfmt", "%f"]
        result = dict(case_id=case_id, group=key[0], query_id=row["query_id"], accession_version=row["accession_version"],
                      returned_accession="",
                      record_type=row["record_type"], topology=row["topology"], position_stratum=key[-1], status="not_checked", error="",
                      sstart=row["sstart"], send=row["send"], slen=row["slen"], flank_start=start, flank_end=end)
        try:
            process = subprocess.run(command, capture_output=True, text=True, check=True, timeout=60)
            with (destination / (case_id + ".fasta")).open("x", encoding="utf-8") as handle:
                handle.write(process.stdout)
            record = SeqIO.read(io.StringIO(process.stdout), "fasta")
            result["returned_accession"] = record.id.split(":", 1)[0]
            sequence = validate_flank_record(record, row["accession_version"], start, end)
            observed = sequence[lower - start:upper - start + 1]
            expected = sequences[row["query_id"]]
            relative_direction = (row["qend"] >= row["qstart"]) == (row["send"] >= row["sstart"])
            if not relative_direction:
                expected = reverse_complement(expected)
            result["status"] = "exact_sequence_corroborated_not_biologically_adjudicated" if observed == expected else "sequence_conflict"
            save_json(destination / (case_id + ".json"), dict(command=command, exit_code=process.returncode, stderr=process.stderr,
                                                            flank_sha256=digest(destination / (case_id + ".fasta")), expected_match=expected,
                                                            observed_match=observed, original_event=row))
        except (OSError, ValueError, subprocess.SubprocessError) as error:
            status = "wrong_accession_returned" if "Accession lookup mismatch" in str(error) else "retrieval_or_validation_failed"
            result.update(status=status, error=str(error))
            save_json(destination / (case_id + "_failure.json"), dict(command=command, error=str(error),
                                                                      stderr=getattr(error, "stderr", None)))
        cases.append(result)
    verify_database(evidence, prefix)
    write_table(destination / "cases.tsv", cases, CASE_SCHEMA)
    report = dict(segment_rows=len(segment_rows), partially_covered_barcode_segments=sum(row["component_kind"] == "barcode" and not row["complete_segment"] for row in segment_rows),
                  tolerant_segment_mismatch_allocation="Unknown; no barcode-specific enrichment claim", selected_cases=len(cases),
                  sequence_corroborated=sum(row["status"].startswith("exact_sequence_corroborated") for row in cases),
                  conflicts=sum(row["status"] == "sequence_conflict" for row in cases), failed=sum(row["status"] == "retrieval_or_validation_failed" for row in cases),
                  wrong_accessions=sum(row["status"] == "wrong_accession_returned" for row in cases),
                  biological_adjudication_completed=False, final_freeze=False,
                  files={path.name: digest(path) for path in destination.iterdir() if path.is_file()})
    save_json(destination / "context.json", report)
    return report