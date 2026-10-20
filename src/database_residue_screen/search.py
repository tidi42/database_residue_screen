import csv
import hashlib
import json
import os
import re
import resource
import shutil
import subprocess
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

from .acceptance import validate_acceptance
from .contracts import ContractError, read_json, read_table, write_table
from .database import verify_inventory, parse_database_info
from .design import digest
from .database_integrity import assert_single_build, volume_evidence
from .example import save_fasta, save_json
from .length_aware import calibrate_lengths
from .matching import BLAST_FIELDS, parse_blast_hsp, strict_hsp, tolerant_hsp
from .sequences import read_fasta


STATUS_SCHEMA = dict(type="object", additionalProperties=False, properties={
    "query_id": {"type": "string"}, "batch_id": {"type": "string"},
    "status": {"type": "string", "enum": ["success_hits", "success_zero"]},
    "hsp_rows": {"type": "integer"}, "unique_targets": {"type": "integer"},
    "strict_targets": {"type": "integer"}, "tolerant_targets": {"type": "integer"},
    "duplicate_hsp_rows": {"type": "integer"}, "target_limit_reached": {"type": "string"},
})


def object_digest(value: dict) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")).hexdigest()


def search_settings(target_limit: int = 500, evalue: float = 10.0, search_space: int | None = None) -> dict:
    if target_limit < 1 or not 0 < evalue < float("inf") or search_space is not None and search_space < 1:
        raise ContractError("Invalid search limit, E-value or effective search space")
    return dict(task="blastn-short", word_size=7, evalue=evalue, dust="no", soft_masking="false", strand="both",
                reward=1, penalty=-3, gapopen=5, gapextend=2, max_target_seqs=target_limit, searchsp=search_space)


def iter_hsps(path: Path, queries: dict[str, str]):
    seen, expected, observed = set(), {}, Counter()
    current = None
    processed = None
    with path.open(encoding="utf-8", newline="") as handle:
        for line_number, line in enumerate(handle, 1):
            if line.startswith("#"):
                if line.startswith("# Query: "):
                    current = line[len("# Query: "):].split()[0]
                    if current not in queries or current in seen or processed is not None:
                        raise ContractError(f"Unexpected or repeated query block at {path.name}:{line_number}")
                    seen.add(current)
                elif match := re.fullmatch(r"# ([0-9]+) hits found\s*", line):
                    if current is None or current in expected:
                        raise ContractError("Repeated or unowned BLAST hit-count marker")
                    expected[current] = int(match.group(1))
                elif match := re.fullmatch(r"# BLAST processed ([0-9]+) queries\s*", line):
                    if processed is not None:
                        raise ContractError("Repeated BLAST completion footer")
                    processed = int(match.group(1))
                continue
            if not line.strip():
                continue
            if processed is not None:
                raise ContractError("HSP data after the BLAST completion footer")
            values = next(csv.reader([line], delimiter="\t", quoting=csv.QUOTE_NONE))
            row = parse_blast_hsp(values, queries)
            if row["qseqid"] != current or current not in expected:
                raise ContractError("HSP row does not belong to its declared query block")
            observed[current] += 1
            yield row
    if (seen != set(queries) or set(expected) != set(queries) or processed != len(queries)
            or any(observed[name] != expected[name] for name in queries)):
        raise ContractError("BLAST output is incomplete or its per-query hit counts disagree with the rows")


def summarise_batch(path: Path, queries: dict[str, str], batch_id: str, target_limit: int) -> list[dict]:
    counts, duplicates = Counter(), Counter()
    targets, strict, tolerant = defaultdict(set), defaultdict(set), defaultdict(set)
    fingerprints, previous_query = set(), None
    for row in iter_hsps(path, queries):
        query_id = row["qseqid"]
        if previous_query != query_id:
            fingerprints.clear()
            previous_query = query_id
        counts[query_id] += 1
        fingerprint = tuple(row[field] for field in BLAST_FIELDS)
        if fingerprint in fingerprints:
            duplicates[query_id] += 1
            continue
        fingerprints.add(fingerprint)
        targets[query_id].add(row["sseqid"])
        if strict_hsp(row):
            strict[query_id].add(row["sseqid"])
        if tolerant_hsp(row):
            tolerant[query_id].add(row["sseqid"])
    return [dict(query_id=query_id, batch_id=batch_id, status="success_hits" if counts[query_id] else "success_zero",
                 hsp_rows=counts[query_id], unique_targets=len(targets[query_id]), strict_targets=len(strict[query_id]),
                 tolerant_targets=len(tolerant[query_id]), duplicate_hsp_rows=duplicates[query_id],
                 target_limit_reached="yes" if len(targets[query_id]) >= target_limit else "no") for query_id in sorted(queries)]


def verify_database(evidence: Path, prefix: Path, hash_content: bool = True) -> dict:
    return verify_inventory(read_json(evidence), prefix, hash_content)


def validate_reuse_request(request: dict) -> None:
    source = request.get("reuse_source")
    if source is None:
        return
    original = source["request"]
    execution_fields = {"threads", "runner_sha256", "reuse_source"}
    if (original.get("reuse_source") is not None or
            {key: value for key, value in original.items() if key not in execution_fields} !=
            {key: value for key, value in request.items() if key not in execution_fields}):
        raise ContractError("Batch reuse requires identical queries, database, parser, executable and scientific settings")


def validate_completed_batch(directory: Path, request: dict) -> dict:
    record = read_json(directory / "complete.json")
    if record.get("projection") is not None:
        from .subset_search import validate_projection_batch

        return validate_projection_batch(directory, request, record)
    required = {"queries.fasta", "blast.tsv", "stderr.txt", "stdout.txt", "status.tsv", "command.json"}
    if record.get("reused_from") is not None:
        validate_reuse_request(request)
        source = request.get("reuse_source")
        original = read_json(directory / "source_complete.json")
        if (source is None or original.get("reused_from") is not None or set(original["files"]) != required or
                original["search_id"] != object_digest(source["request"])):
            raise ContractError("Reused batch has inconsistent source ownership")
        source_checksum = digest(directory / "source_complete.json")
        expected = dict(original, search_id=object_digest(request),
                        files=dict(original["files"], **{"source_complete.json": source_checksum}),
                        reused_from=dict(search_directory=source["directory"], threads=source["request"]["threads"],
                                         completion_sha256=source_checksum))
        if record != expected:
            raise ContractError("Reused batch provenance changed")
        required.add("source_complete.json")
    if record["search_id"] != object_digest(request) or set(record["files"]) != required:
        raise ContractError("Completed batch has inconsistent ownership or an incomplete artifact inventory")
    for name, checksum in record["files"].items():
        if digest(directory / name) != checksum:
            raise ContractError("Completed search artifact changed: " + name)
    return record


def reuse_completed_batch(source: Path, destination: Path, request: dict, batch_id: str, query_ids: list[str]) -> list[Path]:
    attempts = [path for path in source.glob(batch_id + ".attempt_*") if (path / "complete.json").is_file()]
    if not attempts:
        return []
    if len(attempts) != 1 or (source / "quarantine.json").exists():
        raise ContractError("Reused search is quarantined or has duplicate completed batches")
    original = validate_completed_batch(attempts[0], request["reuse_source"]["request"])
    if original["batch_id"] != batch_id or original["query_ids"] != query_ids:
        raise ContractError("Reused batch query membership differs")
    destination.mkdir(exist_ok=False)
    for name in original["files"]:
        shutil.copyfile(attempts[0] / name, destination / name)
    shutil.copyfile(attempts[0] / "complete.json", destination / "source_complete.json")
    source_checksum = digest(destination / "source_complete.json")
    record = dict(original, search_id=object_digest(request),
                  files=dict(original["files"], **{"source_complete.json": source_checksum}),
                  reused_from=dict(search_directory=str(source.resolve()), threads=request["reuse_source"]["request"]["threads"],
                                   completion_sha256=source_checksum))
    save_json(destination / "complete.pending.json", record)
    (destination / "complete.pending.json").rename(destination / "complete.json")
    return [destination]


def run_search(query_fasta: Path, database_prefix: Path, database_evidence: Path, destination: Path,
               settings: dict | None = None, threads: int = 8, batch_size: int = 256, executable: str = "blastn",
               timeout_seconds: float | None = None, mode: str = "production", acceptance: Path | None = None,
               reuse_search: Path | None = None, threshold_mode: str = "fixed",
               maximum_evalue: float | None = None, study: Path | None = None,
               informative_limit: float = .05, testable_limit: float = .01) -> dict:
    if mode not in {"production", "working"} or mode == "production" and acceptance is None or mode == "working" and acceptance is not None:
        raise ContractError("Production search requires --acceptance; working mode must be explicitly selected without production acceptance")
    settings = search_settings() if settings is None else settings
    expected = search_settings(settings["max_target_seqs"], settings["evalue"], settings.get("searchsp"))
    if (settings != expected or threads < 1 or batch_size < 1
            or timeout_seconds is not None and not 0 < timeout_seconds < float("inf")):
        raise ContractError("Unsupported search settings, thread count or batch size")
    queries = read_fasta(query_fasta)
    prefix = database_prefix.resolve()
    database = verify_database(database_evidence, prefix, hash_content=mode == "working")
    accepted = validate_acceptance(acceptance, prefix, database_evidence, executable) if acceptance is not None else None
    assert_single_build(volume_evidence(prefix, database["runtime_metadata"]["volume_names"]))
    binary = shutil.which(executable)
    if binary is None:
        raise ContractError("BLAST executable is unavailable: " + executable)
    binary = str(Path(binary).resolve())
    version = subprocess.run([binary, "-version"], capture_output=True, text=True, check=True, timeout=30).stdout.strip()
    info = subprocess.run([str(Path(binary).with_name("blastdbcmd")), "-db", str(prefix), "-info"],
                          capture_output=True, text=True, check=True, env=dict(os.environ, LC_ALL="C"), timeout=120)
    runtime_metadata = parse_database_info(info.stdout)
    if runtime_metadata != database["runtime_metadata"]:
        raise ContractError("Runtime database dimensions or identity changed since inspection")
    environment_overrides = dict(LC_ALL="C", BLAST_USAGE_REPORT="false")
    if accepted:
        environment_overrides["BLASTDB"] = accepted["taxonomy_support"]["directory"]
    request = dict(kind="accepted_production_search_not_final_freeze" if accepted else "working_search_not_final_freeze",
                   mode=mode, database_acceptance=accepted, environment=environment_overrides,
                   queries_sha256=digest(query_fasta), query_count=len(queries), database_runtime_metadata=runtime_metadata,
                   database_prefix=str(prefix), database_evidence_sha256=digest(database_evidence),
                   database_metadata_conflicts=database["metadata_conflicts"], database_identity_verified=bool(accepted),
                   executable=binary, executable_version=version, executable_sha256=digest(Path(binary)),
                   settings=settings, threads=threads, batch_size=batch_size,
                   runner_sha256=digest(Path(__file__)), parser_sha256=digest(Path(__file__).with_name("matching.py")))
    calibration = calibrate_lengths(binary, [len(sequence) for sequence in queries.values()], database["runtime_metadata"],
                                    settings, threshold_mode, maximum_evalue, informative_limit, testable_limit)
    request["length_calibration"] = calibration
    request["calibration_code_sha256"] = digest(Path(__file__).with_name("length_aware.py"))
    request["failure_policy"] = "abort_run_resume_retries_failed_batches"
    if study is not None:
        from .study_binding import bind_study
        request["study_binding"] = bind_study(study, settings, threshold_mode, maximum_evalue, informative_limit, testable_limit)
    by_length = {row["query_length"]: row for row in calibration["lengths"]}
    if reuse_search is not None:
        reuse_search = reuse_search.resolve()
        if reuse_search == destination.resolve() or (reuse_search / "quarantine.json").exists():
            raise ContractError("Batch reuse needs a distinct, non-quarantined source search")
        request["reuse_source"] = dict(directory=str(reuse_search), request_sha256=digest(reuse_search / "request.json"),
                                      request=read_json(reuse_search / "request.json"))
        validate_reuse_request(request)
    search_id = object_digest(request)
    if destination.exists():
        if read_json(destination / "request.json") != request:
            raise ContractError("Existing search has different inputs, database, code or settings; use a new output directory")
    else:
        destination.mkdir(parents=True, exist_ok=False)
        save_json(destination / "request.json", request)
        shutil.copyfile(database_evidence, destination / "database_evidence.json")
    all_rows, batch_records = [], []
    # Homogeneous length batches ensure each query receives its own statistical
    # threshold and search space. Stable order is part of the request contract.
    batches = []
    for length in sorted(by_length):
        identifiers = sorted(name for name in queries if len(queries[name]) == length)
        batches.extend(identifiers[offset:offset + batch_size] for offset in range(0, len(identifiers), batch_size))
    for batch_index, identifiers in enumerate(batches):
        if accepted and (digest(acceptance) != accepted["sha256"] or (acceptance.parent / "quarantine.json").exists()):
            raise ContractError("Database acceptance changed or was quarantined during search")
        verify_database(database_evidence, prefix, hash_content=False)
        batch_id = f"batch_{batch_index + 1:05d}"
        batch_queries = {query_id: queries[query_id] for query_id in identifiers}
        calibrated = by_length[len(next(iter(batch_queries.values())))]
        batch_settings = dict(settings, evalue=calibrated["evalue"], word_size=calibrated["word_size"],
                              searchsp=calibrated["effective_search_space"])
        attempts = sorted(destination.glob(batch_id + ".attempt_*"))
        completed = [attempt for attempt in attempts if (attempt / "complete.json").is_file()]
        if len(completed) > 1:
            raise ContractError("Multiple completed attempts own the same batch")
        if not completed and reuse_search is not None:
            completed = reuse_completed_batch(reuse_search, destination / f"{batch_id}.attempt_{len(attempts) + 1:05d}",
                                              request, batch_id, list(batch_queries))
        if completed:
            attempt = completed[0]
            record = validate_completed_batch(attempt, request)
            if record["search_id"] != search_id or record["query_ids"] != list(batch_queries):
                raise ContractError("Completed batch ownership disagrees with the search request")
            rows = summarise_batch(attempt / "blast.tsv", batch_queries, batch_id, settings["max_target_seqs"])
        else:
            attempt = destination / f"{batch_id}.attempt_{len(attempts) + 1:05d}"
            attempt.mkdir(exist_ok=False)
            save_fasta(attempt / "queries.fasta", batch_queries)
            command = [binary, "-query", str((attempt / "queries.fasta").resolve()), "-db", str(prefix),
                       "-out", str((attempt / "blast.tsv").resolve()), "-outfmt", "7 " + " ".join(BLAST_FIELDS),
                       "-num_threads", str(threads)]
            for key, value in batch_settings.items():
                if value is not None:
                    command.extend(["-" + key, str(value)])
            save_json(attempt / "command.json", dict(command=command, environment=environment_overrides, search_id=search_id,
                                                        query_ids=list(batch_queries), timeout_seconds=timeout_seconds,
                                                        started_utc=datetime.now(timezone.utc).isoformat()))
            started = time.monotonic()
            try:
                with (attempt / "stdout.txt").open("x", encoding="utf-8") as stdout, (attempt / "stderr.txt").open("x", encoding="utf-8") as stderr:
                    process = subprocess.run(command, stdout=stdout, stderr=stderr, env=dict(os.environ, **environment_overrides),
                                             check=False, timeout=timeout_seconds)
            except subprocess.TimeoutExpired as error:
                save_json(attempt / "incomplete.json", dict(search_id=search_id, query_ids=list(batch_queries), exit_code=None,
                                                             elapsed_seconds=time.monotonic() - started, query_status="incomplete_not_zero",
                                                             error="Runtime limit exceeded; unfinished output is excluded"))
                raise ContractError(f"BLAST batch exceeded its runtime limit: {attempt}") from error
            elapsed = time.monotonic() - started
            if process.returncode:
                save_json(attempt / "failed.json", dict(search_id=search_id, query_ids=list(batch_queries),
                                                         exit_code=process.returncode, elapsed_seconds=elapsed,
                                                         query_status="failed_not_zero", error="See stderr.txt; output is not eligible"))
                raise ContractError(f"BLAST batch failed with exit {process.returncode}: {attempt}")
            try:
                rows = summarise_batch(attempt / "blast.tsv", batch_queries, batch_id, settings["max_target_seqs"])
                verify_database(database_evidence, prefix, hash_content=False)
                if accepted and (digest(acceptance) != accepted["sha256"] or (acceptance.parent / "quarantine.json").exists()):
                    raise ContractError("Database acceptance changed or was quarantined during a batch")
            except (ContractError, OSError, ValueError) as error:
                save_json(attempt / "incomplete.json", dict(search_id=search_id, query_ids=list(batch_queries), exit_code=0,
                                                             elapsed_seconds=elapsed, query_status="incomplete_not_zero", error=str(error)))
                raise
            write_table(attempt / "status.tsv", rows, STATUS_SCHEMA)
            record = dict(search_id=search_id, batch_id=batch_id, query_ids=list(batch_queries), exit_code=0,
                          output_bytes=(attempt / "blast.tsv").stat().st_size,
                          child_process_maxrss=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,
                          memory_scope="Process-lifetime child RSS high-water mark; KiB on Linux, not a per-query measurement",
                          elapsed_seconds=elapsed, completed_utc=datetime.now(timezone.utc).isoformat(),
                          files={name: digest(attempt / name) for name in ("queries.fasta", "blast.tsv", "stderr.txt", "stdout.txt", "status.tsv", "command.json")})
            save_json(attempt / "complete.pending.json", record)
            (attempt / "complete.pending.json").rename(attempt / "complete.json")
            print(f"{batch_id}: {len(rows)} completed queries, {sum(row['hsp_rows'] for row in rows)} HSPs, {elapsed:.1f} seconds", flush=True)
        all_rows.extend(rows)
        batch_records.append(dict(directory=attempt.name, completion_sha256=digest(attempt / "complete.json"),
                      elapsed_seconds=record["elapsed_seconds"], output_bytes=record["output_bytes"],
                      child_process_maxrss=record["child_process_maxrss"]))
    all_rows.sort(key=lambda row: row["query_id"])
    report = dict(search_id=search_id, all_queries_completed=True, final_research_freeze=False,
                  mode=mode, database_acceptance=accepted, output_bytes=sum(row["output_bytes"] for row in batch_records),
                  submitted_queries=len(queries), completed_queries=len(all_rows),
                  censored_queries=sum(not by_length[len(sequence)]["strict_testable"] for sequence in queries.values()),
                  partially_censored_queries=sum(by_length[len(sequence)]["strict_testable"] and not by_length[len(sequence)]["tolerant_fully_testable"] for sequence in queries.values()),
                  testability_counts={endpoint: {state: sum(by_length[len(seq)][endpoint + "_testability"] == state for seq in queries.values())
                                                for state in ("testable", "chance_limited", "censored_untestable")}
                                      for endpoint in ("strict", "tolerant")},
                  success_zero=sum(row["status"] == "success_zero" for row in all_rows),
                  hsp_rows=sum(row["hsp_rows"] for row in all_rows),
                  target_limit_reached=sum(row["target_limit_reached"] == "yes" for row in all_rows),
                  batches=batch_records, interpretation="Completed retrieval at the declared limit; not proof of exhaustive or biologically eligible detection")
    if (destination / "complete.json").is_file():
        if read_json(destination / "complete.json") != report:
            raise ContractError("Search summary disagrees with completed batches")
        if read_table(destination / "status.tsv", STATUS_SCHEMA) != all_rows:
            raise ContractError("Search status table changed since completion")
    else:
        if (destination / "status.tsv").exists():
            if read_table(destination / "status.tsv", STATUS_SCHEMA) != all_rows:
                raise ContractError("Unfinished search summary disagrees with completed batches")
        else:
            write_table(destination / "status.pending.tsv", all_rows, STATUS_SCHEMA)
            (destination / "status.pending.tsv").rename(destination / "status.tsv")
        pending = destination / "complete.pending.json"
        if pending.exists():
            if read_json(pending) != report:
                raise ContractError("Unfinished completion report disagrees with completed batches")
        else:
            save_json(pending, report)
        pending.rename(destination / "complete.json")
    return report