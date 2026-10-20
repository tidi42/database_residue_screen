import argparse
import csv
import json
import shutil
import time
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from .benchmark_calibration import executable_identities, execute, xml_statistics
from .benchmark_engine import direct_search
from .benchmark_fixtures import make_panels
from .benchmark_truth import independent_hsp_truth, precompute_truth
from .contracts import ContractError, read_json, read_table
from .design import digest
from .example import save_fasta, save_json
from .matching import strict_hsp, tolerant_hsp
from .panels import LOGICAL_SCHEMA
from .search import iter_hsps, search_settings
from .sequences import accession_version, read_fasta


ARMS = ("local_fixture", "production_calibrated", "permissive_diagnostic")


def write_rows(path: Path, rows: list[dict], fields=None):
    fields = fields or list(rows[0])
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def csv_rows(path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def check_hashes(root, files):
    for name, checksum in files.items():
        path = root / name
        if not path.resolve().is_relative_to(root.resolve()) or digest(path) != checksum:
            raise ContractError("Frozen benchmark evidence changed: " + name)


def prepare_benchmark(panel: Path, calibration: Path, destination: Path, extension_panel: Path | None = None,
                      extra_lengths: tuple[int, ...] = ()) -> dict:
    definition, calibrated = read_json(panel / "panel.json"), read_json(calibration / "calibration.json")
    check_hashes(panel, definition["files"])
    check_hashes(calibration, calibrated["files"])
    if calibrated["status"] != "engine_validated_version_matched" or calibrated["route"] != "direct_per_query_XML_and_format7":
        raise ContractError("Current benchmark requires validated direct engine calibration")
    logical = read_table(panel / "logical_queries.tsv", LOGICAL_SCHEMA)
    searched = read_fasta(panel / "search_queries.fasta")
    originals = {row["parent_id"]: searched[row["parent_id"]] for row in logical if row["role"] == "real"}
    lengths = sorted({len(sequence) for sequence in originals.values()})
    extension_logical, extension_originals = [], {}
    if extension_panel is not None:
        extension_definition = read_json(extension_panel / "panel.json")
        check_hashes(extension_panel, extension_definition["files"])
        extension_logical = read_table(extension_panel / "logical_queries.tsv", LOGICAL_SCHEMA)
        extension_search = read_fasta(extension_panel / "search_queries.fasta")
        extension_originals = {row["parent_id"]: extension_search[row["parent_id"]] for row in extension_logical if row["role"] == "real"}
    elif extra_lengths:
        raise ContractError("Extra benchmark lengths belong to a declared extension")
    required = set(lengths) | {len(sequence) for sequence in extension_originals.values()} | set(extra_lengths)
    if not originals or not required <= {row["query_length"] for row in calibrated["search_spaces"]}:
        raise ContractError("Calibration does not cover this panel's query lengths")
    dimensions = {(row["database_bases"], row["database_sequences"]) for row in calibrated["search_spaces"]}
    if len(dimensions) != 1:
        raise ContractError("Calibration mixes target database dimensions")
    bases, sequences = next(iter(dimensions))
    shards, exceptions = make_panels(originals)
    if extension_panel is not None:
        from .benchmark_fixtures import EXTENSION_SEED, REAL_VARIANTS, make_extension_panels
        extension_shards, extension_exceptions = make_extension_panels(extension_originals, tuple(extra_lengths))
        shards, exceptions = shards + extension_shards, exceptions + extension_exceptions
    destination.mkdir(parents=True, exist_ok=False)
    shutil.copytree(calibration, destination / "calibration")
    source_directory = Path(__file__).parent
    source_names = ("current_benchmark.py", "benchmark_fixtures.py", "benchmark_truth.py", "benchmark_calibration.py",
                    "benchmark_engine.py", "matching.py", "search.py", "length_aware.py", "sequences.py", "contracts.py")
    (destination / "code").mkdir()
    for name in source_names:
        shutil.copyfile(source_directory / name, destination / "code" / name)
    query_rows, case_rows, truth_rows, witness_rows, shard_rows = [], [], [], [], []
    aliases = defaultdict(list)
    for row in logical + extension_logical:
        if row["role"] == "real" and row not in aliases[row["parent_id"]]:
            aliases[row["parent_id"]].append(row)
    for index, shard in enumerate(shards, 1):
        folder = destination / "fixtures" / shard["shard_id"]
        folder.mkdir(parents=True)
        save_fasta(folder / "queries.fasta", shard["queries"])
        save_fasta(folder / "subjects.fasta", shard["subjects"])
        truth, witnesses = precompute_truth(shard["queries"], shard["subjects"])
        pairs = [dict(query_id=query_id, subject_id=subject_id, strict=int((query_id, subject_id) in truth["strict"]),
                      tolerant=int((query_id, subject_id) in truth["tolerant"])) for query_id in sorted(shard["queries"]) for subject_id in sorted(shard["subjects"])]
        write_rows(folder / "truth_pairs.tsv", pairs)
        write_rows(folder / "truth_witnesses.tsv", witnesses)
        write_rows(folder / "cases.tsv", shard["cases"])
        for query_id, sequence in shard["queries"].items():
            complexity = next(row["complexity"] for row in shard["cases"] if row["query_id"] == query_id)
            query_rows.append(dict(shard_id=shard["shard_id"], panel=shard["panel"], query_id=query_id, length=len(sequence), complexity=complexity,
                                   aliases=";".join(sorted(row["logical_id"] for row in aliases[query_id])),
                                   groups=";".join(sorted({row["group"] for row in aliases[query_id]})),
                                   conditions=";".join(sorted({row["condition"] for row in aliases[query_id]}))))
        case_rows.extend(dict(shard_id=shard["shard_id"], panel=shard["panel"], **row) for row in shard["cases"])
        truth_rows.extend(dict(shard_id=shard["shard_id"], **row) for row in pairs)
        witness_rows.extend(dict(shard_id=shard["shard_id"], **row) for row in witnesses)
        shard_rows.append(dict(shard_id=shard["shard_id"], panel=shard["panel"], length=shard["length"],
                               queries=len(shard["queries"]), subjects=len(shard["subjects"]), pairs=len(pairs),
                               files={path.name: digest(path) for path in folder.iterdir()}))
        if index % 10 == 0 or index == len(shards):
            print(f"Frozen truth: {index}/{len(shards)} shards, {len(truth_rows)} pairs", flush=True)
    for name, rows in (("query_manifest.tsv", query_rows), ("cases.tsv", case_rows), ("truth_pairs.tsv", truth_rows), ("truth_witnesses.tsv", witness_rows)):
        write_rows(destination / name, rows)
    write_rows(destination / "search_spaces.tsv", calibrated["search_spaces"])
    save_json(destination / "construction_exceptions.json", exceptions)
    extension = None
    if extension_panel is not None:
        extension = dict(seed=EXTENSION_SEED, panel_sha256=digest(extension_panel / "panel.json"), panel_path=str(extension_panel.resolve()),
                         parents=len(extension_originals), extra_lengths=sorted(set(extra_lengths)), real_variants=list(REAL_VARIANTS),
                         combined_variants=["clip_both_at_substitution_at", "clip_both_at_substitution_over",
                                            "reverse_clip_both_at_substitution_at"],
                         panels=dict(C="tolerant and combined variants of every real parent", D="synthetic templates at extra lengths"))
    specification = dict(version="current_benchmark_v3_extension" if extension else "current_benchmark_v2", seed=20260922,
                         extension=extension, panel_sha256=digest(panel / "panel.json"),
                         panel_path=str(panel.resolve()), calibration_sha256=digest(destination / "calibration/calibration.json"),
                         binaries=calibrated["binaries"], source_code={name: digest(source_directory / name) for name in source_names},
                         actual_lengths=lengths, original_parents=len(originals),
                         boundary_templates=sum(len(shard["queries"]) for shard in shards if shard["panel"] == "B"),
                         negative_construction_attempts=100, internal_flank=100, terminal_flank=200,
                         maximum_workers=4, wall_seconds=14400, output_bytes=10 * 1024 ** 3, fixture_timeout_seconds=3600,
                         arms={arm: search_settings(2000, 1e12 if arm == "permissive_diagnostic" else 10) for arm in ARMS},
                         production_database=dict(bases=bases, sequences=sequences),
                         truth_scope="Every query-subject pair within each fixed shard; independent of all BLAST output",
                         shards=shard_rows, query_subject_pairs=len(truth_rows), truth_frozen_before_search=True,
                         files={str(path.relative_to(destination)): digest(path) for path in destination.rglob("*") if path.is_file()})
    save_json(destination / "benchmark_specification.json", specification)
    return specification


def analyse_arm(folder, queries, subjects, pairs, cases, query_info, shard, arm):
    observed = defaultdict(lambda: dict(raw_hsps=0, strict_independent=False, tolerant_independent=False,
                                        strict_accepted=False, tolerant_accepted=False))
    for hsp in iter_hsps(folder / "blast.tsv", queries):
        subject_id = accession_version(hsp["sseqid"])
        if subject_id not in subjects:
            raise ContractError("Returned subject is outside the frozen truth universe")
        checked = independent_hsp_truth(hsp, queries[hsp["qseqid"]], subjects[subject_id])
        entry = observed[hsp["qseqid"], subject_id]
        entry["raw_hsps"] += 1
        for endpoint, predicate in (("strict", strict_hsp), ("tolerant", tolerant_hsp)):
            accepted = predicate(hsp)
            if accepted != checked[endpoint]:
                raise ContractError("Production predicate disagrees with independent returned-HSP truth")
            entry[endpoint + "_independent"] |= checked[endpoint]
            entry[endpoint + "_accepted"] |= accepted
    result = []
    for pair in pairs:
        key = pair["query_id"], pair["subject_id"]
        entry = observed[key]
        case = cases.get(key, {})
        output = dict(shard_id=shard["shard_id"], panel=shard["panel"], arm=arm, query_id=key[0], subject_id=key[1], length=len(queries[key[0]]),
                      complexity=query_info[key[0]]["complexity"], error_class=case.get("error_class", "cross_pair"),
                      context=case.get("context", "cross_pair"), raw_hsps=entry["raw_hsps"])
        for endpoint in ("strict", "tolerant"):
            expected, accepted = int(pair[endpoint]), int(entry[endpoint + "_accepted"])
            if accepted and not expected:
                raise ContractError("Accepted HSP contradicts frozen sequence truth")
            category = "retrieval_or_representation_failure" if expected and not accepted else "recovered" if accepted else "correct_nonpositive"
            output.update({endpoint + "_truth": expected, endpoint + "_independent": int(entry[endpoint + "_independent"]),
                           endpoint + "_accepted": accepted, endpoint + "_classification": category})
        result.append(output)
    return result


def run_shard(destination, specification, shard, spaces, started):
    folder = destination / "fixtures" / shard["shard_id"]
    check_hashes(folder, shard["files"])
    queries = read_fasta(folder / "queries.fasta")
    subjects = read_fasta(folder / "subjects.fasta", allow_ambiguity=True)
    pairs = csv_rows(folder / "truth_pairs.tsv")
    cases = {(row["query_id"], row["subject_id"]): row for row in csv_rows(folder / "cases.tsv")}
    query_info = {row["query_id"]: row for row in csv_rows(destination / "query_manifest.tsv") if row["shard_id"] == shard["shard_id"]}
    work = destination / "searches" / shard["shard_id"]
    work.mkdir(parents=True, exist_ok=False)
    binaries = specification["binaries"]
    save_json(work / "database_build.json", execute([binaries["makeblastdb"]["path"], "-in", str(folder / "subjects.fasta"),
                                                    "-dbtype", "nucl", "-parse_seqids", "-out", str(work / "database")], 60))
    results, completions = [], []
    for arm in ARMS:
        remaining = specification["wall_seconds"] - (time.monotonic() - started)
        if remaining < 2:
            raise ContractError("Predeclared benchmark wall budget exhausted")
        settings = dict(specification["arms"][arm])
        if arm != "local_fixture":
            settings["searchsp"] = spaces[shard["length"]]["effective_search_space"]
        complete = direct_search(folder / "queries.fasta", queries, work / "database", work / arm, binaries, settings,
                                 min(specification["fixture_timeout_seconds"], max(1, int(remaining / 2))))
        for row in xml_statistics(work / arm / "blast.xml"):
            if arm != "local_fixture" and row["effective_search_space"] != settings["searchsp"]:
                raise ContractError("Engine did not use the frozen per-length effective search space")
        results.extend(analyse_arm(work / arm, queries, subjects, pairs, cases, query_info, shard, arm))
        check_hashes(folder, shard["files"])
        completions.append(dict(shard_id=shard["shard_id"], arm=arm, status="complete", queries=len(queries), subjects=len(subjects),
                                pair_denominator=len(pairs), hsps=sum(row["hsp_rows"] for row in complete["statuses"]),
                                completion_sha256=digest(work / arm / "complete.json")))
    write_rows(work / "pair_results.tsv", results)
    return completions


def summarise_benchmark(destination, specification, completions):
    cases = csv_rows(destination / "cases.tsv")
    originals = {row["query_id"]: row for row in csv_rows(destination / "query_manifest.tsv") if row["panel"] == "A"}
    spaces = {int(row["query_length"]): row for row in csv_rows(destination / "search_spaces.tsv")}
    exact_cases = {(row["shard_id"], row["query_id"], row["subject_id"]): row for row in cases if row["panel"] == "A" and row["error_class"] in {"exact", "reverse"}}
    recovery, metrics = defaultdict(dict), defaultdict(Counter)
    combined = destination / "benchmark_pair_results.tsv"
    with combined.open("x", encoding="utf-8", newline="") as handle:
        writer = None
        for shard in specification["shards"]:
            rows = csv_rows(destination / "searches" / shard["shard_id"] / "pair_results.tsv")
            if writer is None:
                writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
                writer.writeheader()
            writer.writerows(rows)
            for row in rows:
                identity = row["shard_id"], row["query_id"], row["subject_id"]
                if identity in exact_cases:
                    recovery[identity][row["arm"]] = row
                scopes = (("all_lengths", "all", 0), ("length", "all", int(row["length"])),
                          ("complexity", row["complexity"], 0), ("variant", row["error_class"], 0), ("panel", row["panel"], 0))
                for endpoint in ("strict", "tolerant"):
                    truth, accepted = int(row[endpoint + "_truth"]), int(row[endpoint + "_accepted"])
                    for population, category, length in scopes:
                        counts = metrics[row["arm"], endpoint, population, category, length]
                        counts["positive_pairs" if truth else "negative_pairs"] += 1
                        counts["true_positive" if accepted and truth else "false_negative" if truth else "false_positive" if accepted else "true_negative"] += 1
        handle.flush()
    metric_rows = []
    for (arm, endpoint, population, category, length), counts in sorted(metrics.items()):
        metric_rows.append(dict(search_space=arm, endpoint=endpoint, population=population, error_class=category, length=length,
                                **{name: counts[name] for name in ("positive_pairs", "negative_pairs", "true_positive", "false_negative", "false_positive", "true_negative")},
                                recall=counts["true_positive"] / counts["positive_pairs"] if counts["positive_pairs"] else "not_estimated",
                                precision=counts["true_positive"] / (counts["true_positive"] + counts["false_positive"]) if counts["true_positive"] + counts["false_positive"] else "not_estimated"))
    detectability = []
    for identity, arms in sorted(recovery.items()):
        case, query = exact_cases[identity], originals[identity[1]]
        space = spaces[int(query["length"])]
        production = int(arms["production_calibrated"]["strict_accepted"])
        diagnostic = int(arms["permissive_diagnostic"]["strict_accepted"])
        attribution = ("evalue_censoring_confirmed" if not production and diagnostic and float(space["exact_evalue"]) > 10 else
                       "retrieval_or_representation_unattributed" if not production else "recovered")
        detectability.append(dict(parent_id=identity[1], aliases=query["aliases"], groups=query["groups"], conditions=query["conditions"],
                                   length=int(query["length"]), context=case["context"], orientation=case["designated_orientation"],
                                   subject_id=identity[2], **{arm + "_recovered": int(arms[arm]["strict_accepted"]) for arm in ARMS},
                                   **{arm + "_raw_hsps": int(arms[arm]["raw_hsps"]) for arm in ARMS},
                                   exact_score=int(query["length"]), calibrated_exact_evalue=space["exact_evalue"],
                                   effective_search_space=space["effective_search_space"], classification=attribution,
                                   zero_interpretation="no_retained_support_despite_recovery_in_tested_context" if production else "no_retained_support_exact_spike_recovery_failed_at_production_scale"))
    write_rows(destination / "metrics.tsv", metric_rows)
    write_rows(destination / "short_query_detectability.tsv", detectability)
    write_rows(destination / "completion.tsv", completions)
    return metric_rows, detectability


def run_benchmark(destination: Path, workers: int = 4) -> dict:
    specification = read_json(destination / "benchmark_specification.json")
    if not 1 <= workers <= specification["maximum_workers"]:
        raise ContractError("Worker count exceeds predeclared resource ceiling")
    check_hashes(destination, specification["files"])
    check_hashes(Path(__file__).parent, specification["source_code"])
    binaries = executable_identities(Path(specification["binaries"]["blastn"]["path"]).parent)
    if binaries != specification["binaries"]:
        raise ContractError("Pinned BLAST binaries changed after fixture freeze")
    spaces = {row["query_length"]: row for row in read_json(destination / "calibration/calibration.json")["search_spaces"]}
    started, completions = time.monotonic(), []
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {executor.submit(run_shard, destination, specification, shard, spaces, started): shard for shard in specification["shards"]}
        try:
            for future in as_completed(futures):
                shard = futures[future]
                completions.extend(future.result())
                total_bytes = sum(path.stat().st_size for path in destination.rglob("*") if path.is_file())
                if total_bytes > specification["output_bytes"] or time.monotonic() - started > specification["wall_seconds"]:
                    raise ContractError("Predeclared benchmark resource ceiling reached")
                print(f"Benchmark: {len(completions) // 3}/{len(futures)} shards complete ({time.monotonic() - started:.1f}s)", flush=True)
        except Exception as error:
            for future in futures:
                future.cancel()
            save_json(destination / "failed.json", dict(status="incomplete_not_zero", error=str(error), completed_arms=completions))
            raise
    metric_rows, detectability = summarise_benchmark(destination, specification, sorted(completions, key=lambda row: (row["shard_id"], row["arm"])))
    check_hashes(destination, specification["files"])
    check_hashes(Path(__file__).parent, specification["source_code"])
    report = dict(contract=specification["version"], status="complete", specification_sha256=digest(destination / "benchmark_specification.json"),
                  current_version="2.17.0+", original_queries=specification["original_parents"], boundary_queries=specification["boundary_templates"], independent_pair_universe=specification["query_subject_pairs"],
                  extension=(dict(specification["extension"], shards=sum(1 for shard in specification["shards"] if shard["panel"] in {"C", "D"}))
                             if specification.get("extension") else None),
                  arms=list(ARMS), completed_shards=len(specification["shards"]), complete_arms=len(completions),
                  truth_independent_of_output=True, predicate_disagreements=0, production_database_searches=0,
                  wall_seconds=time.monotonic() - started, workers=workers, diagnostic_evalue=1e12,
                  short_summary={str(length): dict(parents=len({row["parent_id"] for row in detectability if row["length"] == length}),
                                                       witnesses=sum(row["length"] == length for row in detectability),
                                                       calibrated_recovered=sum(row["production_calibrated_recovered"] for row in detectability if row["length"] == length),
                                                       diagnostic_recovered=sum(row["permissive_diagnostic_recovered"] for row in detectability if row["length"] == length)) for length in specification["actual_lengths"]},
                  scope="Independent software QA and statistical detectability, not biological validation or exhaustive core_nt recall",
                  final_freeze=False, files={str(path.relative_to(destination)): digest(path) for path in destination.rglob("*") if path.is_file()})
    save_json(destination / "benchmark.json", report)
    return report


def main():
    parser = argparse.ArgumentParser(description="Prepare independent truth or run the frozen bounded current-version benchmark")
    parser.add_argument("action", choices=("prepare", "run"))
    parser.add_argument("destination", type=Path)
    parser.add_argument("--panel", type=Path)
    parser.add_argument("--calibration", type=Path)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    if args.action == "prepare":
        if args.panel is None or args.calibration is None:
            parser.error("prepare requires --panel and --calibration")
        report = prepare_benchmark(args.panel, args.calibration, args.destination)
        print(json.dumps({key: report[key] for key in ("original_parents", "boundary_templates", "query_subject_pairs")}, indent=2))
    else:
        report = run_benchmark(args.destination, args.workers)
        print(json.dumps({key: report[key] for key in ("status", "predicate_disagreements", "short_summary", "wall_seconds")}, indent=2))


if __name__ == "__main__":
    main()