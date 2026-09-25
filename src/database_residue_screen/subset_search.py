import hashlib
import shutil
from collections import Counter
from pathlib import Path

from .contracts import ContractError, read_json, read_table, write_table
from .design import digest
from .example import save_fasta, save_json
from .matching import BLAST_FIELDS
from .panels import LOGICAL_SCHEMA
from .search import STATUS_SCHEMA, iter_hsps, object_digest, summarise_batch, validate_completed_batch
from .sequences import read_fasta


SCIENCE_FIELDS = ("mode", "database_acceptance", "environment", "database_prefix", "database_evidence_sha256",
                  "database_metadata_conflicts", "database_identity_verified", "executable", "executable_version",
                  "executable_sha256", "settings", "parser_sha256")
PROJECTION_FILES = {"queries.fasta", "blast.tsv", "stderr.txt", "stdout.txt", "status.tsv", "command.json",
                    "source_complete.json", "source_batch.json"}


def require_same_science(request, original):
    if any(request.get(name) != original.get(name) for name in SCIENCE_FIELDS):
        raise ContractError("Query reuse requires the same database, executable, parser and scientific settings")


def panel_queries(panel):
    definition = read_json(panel / "panel.json")
    for name, checksum in definition["files"].items():
        if Path(name).name != name or digest(panel / name) != checksum:
            raise ContractError("Reduced panel artifact changed")
    return read_fasta(panel / "search_queries.fasta")


def subset_text(source, query_ids):
    queries = read_fasta(source / "queries.fasta")
    if not set(query_ids) <= queries.keys() or len(query_ids) != len(set(query_ids)):
        raise ContractError("Projected query membership is invalid")
    for _ in iter_hsps(source / "blast.tsv", queries):
        pass
    blocks, current = {}, None
    with (source / "blast.tsv").open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("# BLAST processed "):
                current = None
            elif line.startswith("# Query: "):
                current = line[len("# Query: "):].split()[0]
                blocks[current] = [line]
            elif current is not None:
                blocks[current].append(line)
    return "# Verified projection of completed source query blocks; not a new BLAST execution\n" + "".join(
        "".join(blocks[query_id]) for query_id in query_ids) + f"# BLAST processed {len(query_ids)} queries\n"


def source_records(sources, queries):
    entries, requests, covered, baseline, fingerprints = [], [], set(), None, {}
    for source in sources:
        source = source.resolve()
        if (source / "quarantine.json").exists():
            raise ContractError("Cannot reuse a quarantined search")
        request = read_json(source / "request.json")
        if request.get("assembly_contract"):
            raise ContractError("Use original executions, not a prior query projection, as subset sources")
        if baseline is not None:
            require_same_science(baseline, request)
        baseline = request
        requests.append(dict(directory=str(source), request_sha256=digest(source / "request.json"), request=request))
        source_seen = set()
        for path in sorted(source.glob("batch_*.attempt_*/complete.json")):
            summary = read_json(path)
            selected = sorted(set(summary["query_ids"]) & queries.keys())
            if not selected:
                continue
            record = validate_completed_batch(path.parent, request)
            if record["exit_code"] != 0 or source_seen.intersection(record["query_ids"]):
                raise ContractError("Source search has failed or duplicate completed query ownership")
            source_seen.update(record["query_ids"])
            original_queries = read_fasta(path.parent / "queries.fasta")
            if set(original_queries) != set(record["query_ids"]) or any(original_queries[name] != queries[name] for name in selected):
                raise ContractError("Source query sequence or membership differs from the selected panel")
            observed = {name: Counter() for name in selected}
            for row in iter_hsps(path.parent / "blast.tsv", original_queries):
                if row["qseqid"] in observed:
                    observed[row["qseqid"]][tuple(row[field] for field in BLAST_FIELDS)] += 1
            for name, values in observed.items():
                fingerprint = object_digest({"alignments": sorted((list(key), count) for key, count in values.items())})
                if name in fingerprints and fingerprints[name] != fingerprint:
                    raise ContractError("Overlapping source searches disagree for query: " + name)
                fingerprints[name] = fingerprint
            selected = sorted(set(selected) - covered)
            if selected:
                entries.append(dict(source_index=len(requests) - 1, batch_directory=path.parent.name,
                                    completion_sha256=digest(path), query_ids=selected))
                covered.update(selected)
    return requests, entries, covered


def prepare_subset(panel, sources, destination):
    queries = panel_queries(panel)
    requests, entries, covered = source_records(sources, queries)
    real_ids = {row["sequence_id"] for row in read_table(panel / "logical_queries.tsv", LOGICAL_SCHEMA) if row["role"] == "real"}
    missing = queries.keys() - covered
    destination.mkdir(parents=True, exist_ok=False)
    files, counts = {}, {}
    for name, selected in (("remaining_references", missing & real_ids), ("remaining_controls", missing - real_ids)):
        counts[name] = len(selected)
        if selected:
            path = destination / (name + ".fasta")
            save_fasta(path, {query_id: queries[query_id] for query_id in sorted(selected)})
            files[path.name] = digest(path)
    report = dict(contract="query_subset_plan_v1", panel_sha256=digest(panel / "panel.json"),
                  query_count=len(queries), reusable_queries=len(covered), counts=counts, sources=requests, batches=entries,
                  source_precedence="Earlier listed source wins only if overlapping completed alignments agree",
                  files=files)
    save_json(destination / "subset.json", report)
    return report


def projection_command(source, query_ids):
    return dict(operation="verified_completed_query_projection", source=source, query_ids=query_ids,
                execution="No BLAST executed for this derived batch; source command and completion are preserved")


def projection_record(directory, request, source, original, batch_id, query_ids):
    return dict(search_id=object_digest(request), batch_id=batch_id, query_ids=query_ids, exit_code=0,
                output_bytes=(directory / "blast.tsv").stat().st_size, elapsed_seconds=0.0,
                source_elapsed_seconds=original["elapsed_seconds"], child_process_maxrss=original["child_process_maxrss"],
                memory_scope="Inherited source-process high-water mark; projection did not execute BLAST",
                completed_utc=original["completed_utc"], projection=source,
                files={name: digest(directory / name) for name in sorted(PROJECTION_FILES)})


def validate_projection_batch(directory, request, record):
    if request.get("assembly_contract") != "verified_query_projection_v1" or set(record["files"]) != PROJECTION_FILES:
        raise ContractError("Invalid query-projection contract or file inventory")
    for name, checksum in record["files"].items():
        if digest(directory / name) != checksum:
            raise ContractError("Projected batch artifact changed: " + name)
    source = read_json(directory / "source_batch.json")
    if source != record["projection"] or not 0 <= source["source_index"] < len(request["source_searches"]):
        raise ContractError("Projected batch source mapping changed")
    evidence = request["source_searches"][source["source_index"]]
    root = Path(evidence["directory"])
    original_request = read_json(root / "request.json")
    if ((root / "quarantine.json").exists() or original_request != evidence["request"] or
            digest(root / "request.json") != evidence["request_sha256"] or original_request.get("assembly_contract")):
        raise ContractError("Projected search source changed or was quarantined")
    require_same_science(request, original_request)
    if Path(source["batch_directory"]).name != source["batch_directory"]:
        raise ContractError("Unsafe projected source directory")
    original_directory = root / source["batch_directory"]
    original = validate_completed_batch(original_directory, original_request)
    if (digest(original_directory / "complete.json") != source["completion_sha256"] or
            digest(directory / "source_complete.json") != source["completion_sha256"]):
        raise ContractError("Projected source completion changed")
    query_ids = record["query_ids"]
    source_queries = read_fasta(original_directory / "queries.fasta")
    if not set(query_ids) <= source_queries.keys() or read_fasta(directory / "queries.fasta") != {name: source_queries[name] for name in query_ids}:
        raise ContractError("Projected query sequences differ from their source")
    expected_text = subset_text(original_directory, query_ids)
    if hashlib.sha256(expected_text.encode("utf-8")).hexdigest() != digest(directory / "blast.tsv"):
        raise ContractError("Projected alignments differ from completed source blocks")
    if read_json(directory / "command.json") != projection_command(source, query_ids):
        raise ContractError("Projected execution provenance changed")
    if record != projection_record(directory, request, source, original, record["batch_id"], query_ids):
        raise ContractError("Projected completion record changed")
    statuses = summarise_batch(directory / "blast.tsv", read_fasta(directory / "queries.fasta"), record["batch_id"], request["settings"]["max_target_seqs"])
    if statuses != read_table(directory / "status.tsv", STATUS_SCHEMA):
        raise ContractError("Projected query status differs from recomputed source results")
    return record


def assemble_subset(panel, sources, destination):
    queries = panel_queries(panel)
    requests, entries, covered = source_records(sources, queries)
    if covered != queries.keys():
        raise ContractError(f"Selected panel has {len(queries.keys() - covered)} queries without completed evidence")
    baseline = requests[0]["request"]
    request = {name: value for name, value in baseline.items() if name != "reuse_source"}
    request.update(queries_sha256=digest(panel / "search_queries.fasta"), query_count=len(queries),
                   assembly_contract="verified_query_projection_v1", source_searches=requests,
                   assembler_sha256=digest(Path(__file__)),
                   execution_scope="Source requests preserve actual batch sizes, thread counts and execution times; this is a query-level assembly")
    destination.mkdir(parents=True, exist_ok=False)
    save_json(destination / "request.json", request)
    source_database = Path(requests[0]["directory"]) / "database_evidence.json"
    if digest(source_database) != request["database_evidence_sha256"]:
        raise ContractError("Source database evidence changed")
    shutil.copyfile(source_database, destination / "database_evidence.json")
    statuses, batches = [], []
    for index, entry in enumerate(entries, 1):
        batch_id = f"batch_{index:05d}"
        directory = destination / (batch_id + ".attempt_00001")
        directory.mkdir()
        query_ids = entry["query_ids"]
        source = {name: entry[name] for name in ("source_index", "batch_directory", "completion_sha256")}
        original_directory = Path(requests[source["source_index"]]["directory"]) / source["batch_directory"]
        original = read_json(original_directory / "complete.json")
        save_fasta(directory / "queries.fasta", {name: queries[name] for name in query_ids})
        (directory / "blast.tsv").write_text(subset_text(original_directory, query_ids), encoding="utf-8")
        shutil.copyfile(original_directory / "complete.json", directory / "source_complete.json")
        for name in ("stderr.txt", "stdout.txt"):
            (directory / name).write_text("", encoding="utf-8")
        save_json(directory / "source_batch.json", source)
        save_json(directory / "command.json", projection_command(source, query_ids))
        rows = summarise_batch(directory / "blast.tsv", {name: queries[name] for name in query_ids}, batch_id, request["settings"]["max_target_seqs"])
        write_table(directory / "status.tsv", rows, STATUS_SCHEMA)
        record = projection_record(directory, request, source, original, batch_id, query_ids)
        save_json(directory / "complete.json", record)
        validate_projection_batch(directory, request, record)
        statuses.extend(rows)
        batches.append(dict(directory=directory.name, completion_sha256=digest(directory / "complete.json"),
                            elapsed_seconds=0.0, output_bytes=record["output_bytes"], child_process_maxrss=record["child_process_maxrss"]))
    write_table(destination / "status.tsv", statuses, STATUS_SCHEMA)
    report = dict(search_id=object_digest(request), all_queries_completed=True, final_research_freeze=False,
                  mode=request["mode"], database_acceptance=request.get("database_acceptance"), submitted_queries=len(queries),
                  completed_queries=len(statuses), success_zero=sum(row["status"] == "success_zero" for row in statuses),
                  hsp_rows=sum(row["hsp_rows"] for row in statuses), target_limit_reached=sum(row["target_limit_reached"] == "yes" for row in statuses),
                  output_bytes=sum(row["output_bytes"] for row in batches), batches=batches,
                  interpretation="Verified per-query subset of matching completed searches; retained results, not exhaustive retrieval")
    save_json(destination / "complete.json", report)
    return report


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Assemble a declared query subset from matching, verified completed BLAST results")
    parser.add_argument("panel", type=Path)
    parser.add_argument("--sources", type=Path, nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    assembled = assemble_subset(args.panel, args.sources, args.output)
    print(f"Assembled {assembled['completed_queries']} completed queries from {len(assembled['batches'])} source projections", flush=True)