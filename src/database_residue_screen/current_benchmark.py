"""Panel-derived planted-sequence benchmark with the production engine, settings and calibrated search spaces.

``prepare_benchmark`` generates one fixture set from the real parents of one panel (plus optional extra template
lengths and declared templates), precomputes the independent sequence truth of every query-subject pair and hashes
fixtures, truth and code before any search. ``run_benchmark`` verifies those hashes, the pinned binaries and the
calibration of one database size, searches each shard in three arms and classifies every pair against the
precomputed truth; the same fixture set is run once per database size.
"""
import argparse
import csv
import hashlib
import json
import re
import shutil
import time
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from .benchmark_calibration import executable_identities, execute, xml_statistics
from .benchmark_engine import direct_search
from .benchmark_fixtures import FIXTURE_SEED, REAL_VARIANTS, VARIANT_SEED, make_panels
from .benchmark_truth import independent_hsp_truth, precompute_truth
from .contracts import ContractError, read_json, read_table
from .design import digest
from .example import save_fasta, save_json
from .matching import strict_hsp, tolerant_hsp
from .panels import LOGICAL_SCHEMA
from .search import iter_hsps, search_settings
from .sequences import accession_version, read_fasta


ARMS = ("local_fixture", "production_calibrated", "permissive_diagnostic")
BENCHMARK_CONTRACT = "panel_benchmark_v1"
PANELS = dict(A="exact and reverse-complement plants of every real parent, internal and terminal, plus a no-seed negative",
              B="boundary templates at the panel lengths with all boundary variants",
              C="tolerant-boundary and combined clip/substitution variants of every real parent",
              D="boundary and declared templates with combined variants at the declared extra lengths")
SOURCE_NAMES = ("current_benchmark.py", "benchmark_fixtures.py", "benchmark_truth.py", "benchmark_calibration.py",
                "benchmark_engine.py", "matching.py", "search.py", "length_aware.py", "sequences.py", "contracts.py")


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
            raise ContractError("Benchmark evidence changed since preparation: " + name)


def load_templates(path: Path) -> list[dict]:
    """Declared templates written by the panel stage: {"templates": [{"name", "sequence", "repeat_from", ...}]}."""
    return read_json(path).get("templates", [])


def checked_templates(templates) -> list[dict]:
    result = []
    for template in templates:
        name, sequence, repeat_from = template.get("name", ""), str(template.get("sequence", "")).upper(), template.get("repeat_from")
        if (not re.fullmatch(r"[A-Za-z0-9_]+", name) or not sequence or set(sequence) - set("ACGT")
                or not isinstance(repeat_from, int) or not 1 <= repeat_from <= len(sequence)):
            raise ContractError(f"Invalid benchmark template {name!r}")
        result.append(dict(name=name, sequence=sequence, repeat_from=repeat_from))
    if len({template["name"] for template in result}) != len(result):
        raise ContractError("Benchmark template names must be unique")
    return result


def prepare_benchmark(panel: Path, destination: Path, extra_lengths: tuple[int, ...] = (), templates=(), shard_limit: int = 0,
                      maximum_workers: int = 4, wall_seconds: int = 14400, output_bytes: int = 10 * 1024 ** 3,
                      fixture_timeout_seconds: int = 3600) -> dict:
    """One hashed fixture set and its truth for the real parents of `panel`. `shard_limit` > 0 keeps the first shards of
    each fixture panel (A-D) for a bounded dry run; the limit and the number of available shards are recorded."""
    definition = read_json(panel / "panel.json")
    check_hashes(panel, definition["files"])
    logical = read_table(panel / "logical_queries.tsv", LOGICAL_SCHEMA)
    searched = read_fasta(panel / "search_queries.fasta")
    originals = {row["parent_id"]: searched[row["parent_id"]] for row in logical if row["role"] == "real"}
    if not originals:
        raise ContractError("The benchmark needs at least one real parent")
    lengths = sorted({len(sequence) for sequence in originals.values()})
    extra_lengths = tuple(sorted(set(extra_lengths)))
    if any(length < 1 for length in extra_lengths) or shard_limit < 0:
        raise ContractError("Extra benchmark lengths must be positive and the shard limit nonnegative")
    if min(maximum_workers, wall_seconds, output_bytes, fixture_timeout_seconds) < 1:
        raise ContractError("Benchmark resource ceilings must be positive")
    templates = checked_templates(templates)
    shards, exceptions = make_panels(originals, extra_lengths, templates=templates)
    available = len(shards)
    if shard_limit:
        kept, per_panel = [], Counter()
        for shard in shards:
            if per_panel[shard["panel"]] < shard_limit:
                kept.append(shard)
                per_panel[shard["panel"]] += 1
        shards = kept
        exceptions = [row for row in exceptions if row["shard_id"] in {shard["shard_id"] for shard in shards}]
    destination.mkdir(parents=True, exist_ok=False)
    source_directory = Path(__file__).parent
    (destination / "code").mkdir()
    for name in SOURCE_NAMES:
        shutil.copyfile(source_directory / name, destination / "code" / name)
    query_rows, case_rows, truth_rows, witness_rows, shard_rows = [], [], [], [], []
    aliases = defaultdict(list)
    for row in logical:
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
            print(f"Precomputed truth: {index}/{len(shards)} shards, {len(truth_rows)} pairs", flush=True)
    for name, rows in (("query_manifest.tsv", query_rows), ("cases.tsv", case_rows), ("truth_pairs.tsv", truth_rows), ("truth_witnesses.tsv", witness_rows)):
        write_rows(destination / name, rows)
    save_json(destination / "construction_exceptions.json", exceptions)
    specification = dict(version=BENCHMARK_CONTRACT, seed=FIXTURE_SEED, variant_seed=VARIANT_SEED, panels=PANELS,
                         extra_lengths=list(extra_lengths), real_variants=list(REAL_VARIANTS),
                         templates=[dict(name=template["name"], repeat_from=template["repeat_from"], length=len(template["sequence"]),
                                         sha256=digest_text(template["sequence"])) for template in templates],
                         combined_variants=["clip_both_at_substitution_at", "clip_both_at_substitution_over",
                                            "reverse_clip_both_at_substitution_at"],
                         panel_sha256=digest(panel / "panel.json"), panel_path=str(panel.resolve()),
                         source_code={name: digest(source_directory / name) for name in SOURCE_NAMES},
                         actual_lengths=lengths, required_lengths=sorted({shard["length"] for shard in shards}),
                         original_parents=len(originals),
                         boundary_templates=sum(len(shard["queries"]) for shard in shards if shard["panel"] == "B"),
                         extra_length_templates=sum(len(shard["queries"]) for shard in shards if shard["panel"] == "D"),
                         shard_limit=shard_limit, shards_available=available,
                         negative_construction_attempts=100, internal_flank=100, terminal_flank=200,
                         maximum_workers=maximum_workers, wall_seconds=wall_seconds, output_bytes=output_bytes,
                         fixture_timeout_seconds=fixture_timeout_seconds,
                         arms={arm: search_settings(2000, 1e12 if arm == "permissive_diagnostic" else 10) for arm in ARMS},
                         truth_scope="Every query-subject pair within each fixed shard; independent of all BLAST output",
                         shards=shard_rows, query_subject_pairs=len(truth_rows), truth_hashed_before_search=True,
                         files={str(path.relative_to(destination)): digest(path) for path in destination.rglob("*") if path.is_file()})
    save_json(destination / "benchmark_specification.json", specification)
    return specification


def digest_text(text: str) -> str:
    return hashlib.sha256(text.encode("ascii")).hexdigest()


def analyse_arm(folder, queries, subjects, pairs, cases, query_info, shard, arm):
    observed = defaultdict(lambda: dict(raw_hsps=0, strict_independent=False, tolerant_independent=False,
                                        strict_accepted=False, tolerant_accepted=False))
    for hsp in iter_hsps(folder / "blast.tsv", queries):
        subject_id = accession_version(hsp["sseqid"])
        if subject_id not in subjects:
            raise ContractError("Returned subject is outside the precomputed truth universe")
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
                raise ContractError("Accepted HSP contradicts the precomputed sequence truth")
            category = "retrieval_or_representation_failure" if expected and not accepted else "recovered" if accepted else "correct_nonpositive"
            output.update({endpoint + "_truth": expected, endpoint + "_independent": int(entry[endpoint + "_independent"]),
                           endpoint + "_accepted": accepted, endpoint + "_classification": category})
        result.append(output)
    return result


def run_shard(fixtures, destination, specification, binaries, shard, spaces, started):
    folder = fixtures / "fixtures" / shard["shard_id"]
    check_hashes(folder, shard["files"])
    queries = read_fasta(folder / "queries.fasta")
    subjects = read_fasta(folder / "subjects.fasta", allow_ambiguity=True)
    pairs = csv_rows(folder / "truth_pairs.tsv")
    cases = {(row["query_id"], row["subject_id"]): row for row in csv_rows(folder / "cases.tsv")}
    query_info = {row["query_id"]: row for row in csv_rows(fixtures / "query_manifest.tsv") if row["shard_id"] == shard["shard_id"]}
    work = destination / "searches" / shard["shard_id"]
    work.mkdir(parents=True, exist_ok=False)
    save_json(work / "database_build.json", execute([binaries["makeblastdb"]["path"], "-in", str(folder / "subjects.fasta"),
                                                    "-dbtype", "nucl", "-parse_seqids", "-out", str(work / "database")], 60))
    results, completions = [], []
    for arm in ARMS:
        remaining = specification["wall_seconds"] - (time.monotonic() - started)
        if remaining < 2:
            raise ContractError("Declared benchmark wall budget exhausted")
        settings = dict(specification["arms"][arm])
        if arm != "local_fixture":
            settings["searchsp"] = spaces[shard["length"]]["effective_search_space"]
        complete = direct_search(folder / "queries.fasta", queries, work / "database", work / arm, binaries, settings,
                                 min(specification["fixture_timeout_seconds"], max(1, int(remaining / 2))))
        for row in xml_statistics(work / arm / "blast.xml"):
            if arm != "local_fixture" and row["effective_search_space"] != settings["searchsp"]:
                raise ContractError("Engine did not use the declared per-length effective search space")
        results.extend(analyse_arm(work / arm, queries, subjects, pairs, cases, query_info, shard, arm))
        check_hashes(folder, shard["files"])
        completions.append(dict(shard_id=shard["shard_id"], arm=arm, status="complete", queries=len(queries), subjects=len(subjects),
                                pair_denominator=len(pairs), hsps=sum(row["hsp_rows"] for row in complete["statuses"]),
                                completion_sha256=digest(work / arm / "complete.json")))
    write_rows(work / "pair_results.tsv", results)
    return completions


def summarise_benchmark(fixtures, destination, specification, completions):
    cases = csv_rows(fixtures / "cases.tsv")
    originals = {row["query_id"]: row for row in csv_rows(fixtures / "query_manifest.tsv") if row["panel"] == "A"}
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
    write_rows(destination / "short_query_detectability.tsv", detectability, None if detectability else ["parent_id", "length", "classification"])
    write_rows(destination / "completion.tsv", completions)
    return metric_rows, detectability


def run_benchmark(fixtures: Path, calibration: Path, destination: Path, workers: int = 4) -> dict:
    """Search the prepared fixture set with the search spaces of one calibrated database size."""
    specification = read_json(fixtures / "benchmark_specification.json")
    if not 1 <= workers <= specification["maximum_workers"]:
        raise ContractError("Worker count exceeds the declared resource ceiling")
    check_hashes(fixtures, specification["files"])
    check_hashes(Path(__file__).parent, specification["source_code"])
    calibrated = read_json(calibration / "calibration.json")
    check_hashes(calibration, calibrated["files"])
    if calibrated["status"] != "engine_validated_version_matched" or calibrated["route"] != "direct_per_query_XML_and_format7":
        raise ContractError("The benchmark requires a validated direct engine calibration")
    spaces = {row["query_length"]: row for row in calibrated["search_spaces"]}
    if not set(specification["required_lengths"]) <= set(spaces):
        raise ContractError("Calibration does not cover the query lengths of the fixture set")
    dimensions = {(row["database_bases"], row["database_sequences"]) for row in calibrated["search_spaces"]}
    if len(dimensions) != 1:
        raise ContractError("Calibration mixes target database dimensions")
    bases, sequences = next(iter(dimensions))
    binaries = executable_identities(Path(calibrated["binaries"]["blastn"]["path"]).parent)
    if binaries != calibrated["binaries"]:
        raise ContractError("Pinned BLAST binaries changed after calibration")
    destination.mkdir(parents=True, exist_ok=False)
    shutil.copytree(calibration, destination / "calibration")
    write_rows(destination / "search_spaces.tsv", calibrated["search_spaces"])
    started, completions = time.monotonic(), []
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {executor.submit(run_shard, fixtures, destination, specification, binaries, shard, spaces, started): shard
                   for shard in specification["shards"]}
        try:
            for future in as_completed(futures):
                completions.extend(future.result())
                total_bytes = sum(path.stat().st_size for path in destination.rglob("*") if path.is_file())
                if total_bytes > specification["output_bytes"] or time.monotonic() - started > specification["wall_seconds"]:
                    raise ContractError("Declared benchmark resource ceiling reached")
                print(f"Benchmark: {len(completions) // 3}/{len(futures)} shards complete ({time.monotonic() - started:.1f}s)", flush=True)
        except Exception as error:
            for future in futures:
                future.cancel()
            save_json(destination / "failed.json", dict(status="incomplete_not_zero", error=str(error), completed_arms=completions))
            raise
    metric_rows, detectability = summarise_benchmark(fixtures, destination, specification,
                                                     sorted(completions, key=lambda row: (row["shard_id"], row["arm"])))
    check_hashes(fixtures, specification["files"])
    check_hashes(Path(__file__).parent, specification["source_code"])
    report = dict(contract=specification["version"], status="complete", specification_sha256=digest(fixtures / "benchmark_specification.json"),
                  fixtures=str(fixtures.resolve()), calibration_sha256=digest(destination / "calibration/calibration.json"),
                  binaries=calibrated["binaries"], production_database=dict(bases=bases, sequences=sequences),
                  current_version="2.17.0+", original_queries=specification["original_parents"], boundary_queries=specification["boundary_templates"],
                  independent_pair_universe=specification["query_subject_pairs"],
                  extra_lengths=specification["extra_lengths"], extra_length_queries=specification["extra_length_templates"],
                  shard_limit=specification["shard_limit"], shards_available=specification["shards_available"],
                  panel_shards={panel: sum(1 for shard in specification["shards"] if shard["panel"] == panel) for panel in sorted(PANELS)},
                  arms=list(ARMS), completed_shards=len(specification["shards"]), complete_arms=len(completions),
                  truth_independent_of_output=True, predicate_disagreements=0, production_database_searches=0,
                  wall_seconds=time.monotonic() - started, workers=workers, diagnostic_evalue=1e12,
                  short_summary={str(length): dict(parents=len({row["parent_id"] for row in detectability if row["length"] == length}),
                                                       witnesses=sum(row["length"] == length for row in detectability),
                                                       calibrated_recovered=sum(row["production_calibrated_recovered"] for row in detectability if row["length"] == length),
                                                       diagnostic_recovered=sum(row["permissive_diagnostic_recovered"] for row in detectability if row["length"] == length)) for length in specification["actual_lengths"]},
                  scope="Independent software QA and statistical detectability, not biological validation or exhaustive database recall",
                  files={str(path.relative_to(destination)): digest(path) for path in destination.rglob("*") if path.is_file()})
    save_json(destination / "benchmark.json", report)
    return report


def main():
    parser = argparse.ArgumentParser(description="Prepare one hashed fixture set with its truth, or run it with one database-size calibration")
    parser.add_argument("action", choices=("prepare", "run"))
    parser.add_argument("fixtures", type=Path)
    parser.add_argument("--panel", type=Path)
    parser.add_argument("--extra-lengths", type=int, nargs="*", default=[])
    parser.add_argument("--templates", type=Path)
    parser.add_argument("--shard-limit", type=int, default=0)
    parser.add_argument("--calibration", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    if args.action == "prepare":
        if args.panel is None:
            parser.error("prepare requires --panel")
        report = prepare_benchmark(args.panel, args.fixtures, tuple(args.extra_lengths), load_templates(args.templates) if args.templates else (),
                                   args.shard_limit)
        print(json.dumps({key: report[key] for key in ("original_parents", "boundary_templates", "query_subject_pairs")}, indent=2))
    else:
        if args.calibration is None or args.output is None:
            parser.error("run requires --calibration and --output")
        report = run_benchmark(args.fixtures, args.calibration, args.output, args.workers)
        print(json.dumps({key: report[key] for key in ("status", "predicate_disagreements", "short_summary", "wall_seconds")}, indent=2))


if __name__ == "__main__":
    main()
