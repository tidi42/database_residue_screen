import math
import random
import subprocess
from pathlib import Path

from .benchmark_truth import precompute_truth, sequence_witnesses
from .contracts import table_schema, write_table
from .database import inspect_database
from .design import digest
from .example import save_fasta, save_json
from .matching import strict_hsp, tolerant_hsp
from .search import iter_hsps, run_search, search_settings
from .sequences import accession_version, reverse_complement


SUPPORTED_MINIMUM = {"strict": 24, "tolerant": 32}
CASE_SCHEMA = dict(type="object", additionalProperties=False, properties={
    "query_id": {"type": "string"}, "subject_id": {"type": "string"}, "error_class": {"type": "string"},
    "length": {"type": "integer"}, "complexity": {"type": "string"},
})
METRIC_SCHEMA = dict(type="object", additionalProperties=False, properties={
    "search_space": {"type": "string"}, "endpoint": {"type": "string"}, "population": {"type": "string"},
    "error_class": {"type": "string"}, "length": {"type": "integer"},
    "positive_pairs": {"type": "integer"}, "true_positive": {"type": "integer"},
    "false_negative": {"type": "integer"}, "negative_pairs": {"type": "integer"}, "false_positive": {"type": "integer"},
    "recall": {"type": "string"}, "precision": {"type": "string"},
})


def exact_reference(query: str, subject: str) -> bool:
    return query in subject or reverse_complement(query) in subject


def gapfree_reference(query: str, subject: str) -> bool:
    return sequence_witnesses(query, subject)["tolerant"] is not None


def make_cases(seed: int = 20260915) -> tuple[dict, dict, list[dict]]:
    generator = random.Random(seed)
    queries, subjects, cases = {}, {}, []
    specifications = [(length, "random", replicate) for length in (6, 8, 12, 16, 20, 24, 32, 48, 80, 128) for replicate in range(3)]
    specifications += [(length, complexity, 0) for length in (24, 32, 80) for complexity in ("homopolymer", "periodic")]
    for length, complexity, replicate in specifications:
        query_id = f"{complexity}_{length:03d}_{replicate}"
        sequence = ("".join(generator.choices("ACGT", k=length)) if complexity == "random" else
                    "A" * length if complexity == "homopolymer" else ("ACG" * length)[:length])
        queries[query_id] = sequence
        middle = length // 2
        replacement = next(base for base in "ACGT" if base != sequence[middle])
        variants = dict(exact=sequence, reverse=reverse_complement(sequence), partial=sequence[:math.floor(0.8 * length)],
                        substitution=sequence[:middle] + replacement + sequence[middle + 1:],
                        indel=sequence[:middle] + replacement + sequence[middle:],
                        no_match="".join(generator.choices("ACGT", k=length)))
        for error_class, variant in variants.items():
            subject_id = f"BENCH{len(subjects) + 1:06d}.1"
            subjects[subject_id] = "N" * 12 + variant + "N" * 12
            cases.append(dict(query_id=query_id, subject_id=subject_id, error_class=error_class,
                              length=length, complexity=complexity))
    return queries, subjects, cases


def metric_row(space: str, endpoint: str, population: str, error_class: str, length: int,
               pairs: set, truth: set, observed: set) -> dict:
    positives = pairs & truth
    true_positive = len(positives & observed)
    false_positive = len((pairs - truth) & observed)
    return dict(search_space=space, endpoint=endpoint, population=population, error_class=error_class, length=length,
                positive_pairs=len(positives), true_positive=true_positive, false_negative=len(positives - observed),
                negative_pairs=len(pairs - truth), false_positive=false_positive,
                recall=str(true_positive / len(positives)) if positives else "not_estimated",
                precision=str(true_positive / (true_positive + false_positive)) if true_positive + false_positive else "not_estimated")


def run_benchmark(destination: Path, production_bases: int, seed: int = 20260915) -> dict:
    if production_bases < 1:
        raise ValueError("Production base count must be positive")
    queries, subjects, cases = make_cases(seed)
    destination.mkdir(parents=True, exist_ok=False)
    save_fasta(destination / "queries.fasta", queries)
    save_fasta(destination / "subjects.fasta", subjects)
    write_table(destination / "cases.tsv", cases, CASE_SCHEMA)
    definition = dict(scope="synthetic_algorithm_evaluation_not_biological_validation", tuning=False, seed=seed,
                      supported_minimum=SUPPORTED_MINIMUM, strict_recall_required=1.0, tolerant_recall_required=0.95,
                      precision_required=0.95, pair_unit="query_subject", production_bases=production_bases,
                      conservative_search_space=production_bases * max(map(len, queries.values())),
                      search_space_caveat="Common unadjusted upper bound using the longest test query; not an exact reproduction of BLAST length adjustments in core_nt",
                      query_sha256=digest(destination / "queries.fasta"), subject_sha256=digest(destination / "subjects.fasta"))
    save_json(destination / "definition.json", definition)
    prefix = destination / "benchmark_db"
    command = ["makeblastdb", "-dbtype", "nucl", "-parse_seqids", "-in", str((destination / "subjects.fasta").resolve()), "-out", str(prefix.resolve())]
    process = subprocess.run(command, capture_output=True, text=True, check=True)
    save_json(destination / "database_build.json", dict(command=command, stdout=process.stdout, stderr=process.stderr, exit_code=process.returncode))
    inspect_database(prefix, destination / "database_evidence")
    all_pairs = {(query_id, subject_id) for query_id in queries for subject_id in subjects}
    truth, witnesses = precompute_truth(queries, subjects)
    truth_rows = [dict(query_id=query_id, subject_id=subject_id, strict=int((query_id, subject_id) in truth["strict"]),
                       tolerant=int((query_id, subject_id) in truth["tolerant"])) for query_id, subject_id in sorted(all_pairs)]
    write_table(destination / "truth_pairs.tsv", truth_rows, table_schema(("query_id", "subject_id"), ("strict", "tolerant")))
    write_table(destination / "truth_witnesses.tsv", witnesses,
                table_schema(("query_id", "subject_id", "endpoint", "strand"), ("qstart", "qend", "sstart", "send", "span", "matches", "mismatches")))
    truth_sha256 = digest(destination / "truth_pairs.tsv")
    metrics, checks = [], []
    for space, search_space in (("small_database", None), ("production_size_bound", definition["conservative_search_space"])):
        directory = destination / space
        completed = run_search(destination / "queries.fasta", prefix, destination / "database_evidence" / "database_identity.json",
                               directory, search_settings(10000, search_space=search_space), threads=1, batch_size=len(queries), mode="working")
        observed = {"strict": set(), "tolerant": set()}
        for batch in completed["batches"]:
            for row in iter_hsps(directory / batch["directory"] / "blast.tsv", queries):
                pair = row["qseqid"], accession_version(row["sseqid"])
                if strict_hsp(row):
                    observed["strict"].add(pair)
                if tolerant_hsp(row):
                    observed["tolerant"].add(pair)
        if digest(destination / "truth_pairs.tsv") != truth_sha256:
            raise ValueError("Independent benchmark truth changed during search")
        for endpoint, expected in truth.items():
            population = all_pairs
            supported = {pair for pair in population if len(queries[pair[0]]) >= SUPPORTED_MINIMUM[endpoint]}
            summary = metric_row(space, endpoint, "prespecified_supported", "all", 0, supported, expected, observed[endpoint])
            metrics.extend([summary, metric_row(space, endpoint, "all_lengths", "all", 0, population, expected, observed[endpoint])])
            if space == "production_size_bound":
                recall = float(summary["recall"]) if summary["recall"] != "not_estimated" else 0
                precision = float(summary["precision"]) if summary["precision"] != "not_estimated" else 0
                checks.append(dict(endpoint=endpoint, recall=recall, precision=precision,
                                   passed=recall >= definition[endpoint + "_recall_required"] and precision >= definition["precision_required"]))
            for length, error_class in sorted({(row["length"], row["error_class"]) for row in cases}):
                subset = {(row["query_id"], row["subject_id"]) for row in cases if row["length"] == length and row["error_class"] == error_class}
                metrics.append(metric_row(space, endpoint, "declared_cases", error_class, length, subset, expected, observed[endpoint]))
    write_table(destination / "metrics.tsv", metrics, METRIC_SCHEMA)
    report = dict(scope=definition["scope"], queries=len(queries), subjects=len(subjects), declared_cases=len(cases),
                  supported_minimum=SUPPORTED_MINIMUM, acceptance=checks, all_supported_checks_passed=all(row["passed"] for row in checks),
                  biological_specificity_validated=False, final_research_freeze=False,
                  truth_independent_of_search=True, truth_sha256=truth_sha256,
                  files={name: digest(destination / name) for name in ("definition.json", "queries.fasta", "subjects.fasta", "cases.tsv", "metrics.tsv", "truth_pairs.tsv", "truth_witnesses.tsv")})
    save_json(destination / "benchmark.json", report)
    return report