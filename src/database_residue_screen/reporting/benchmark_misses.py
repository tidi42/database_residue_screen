"""Benchmark misses: persistent tolerant misses of the panel-derived benchmark at every database size.

A miss is a truth-positive query/subject pair not accepted at the permissive diagnostic arm. Each miss is joined to its
independent tolerant witness (fixed and hashed before any search) and classified:
  hsp_returned_not_qualifying   BLAST returned HSP(s) for the pair but none satisfies the predicate
  no_hsp_seed_present           no HSP although the witness contains an exact run >= the word size
  no_hsp_no_seed                no HSP and no exact run >= the word size in the witness
The witness raw score (reward 1, penalty -3) is reported. Totals are checked against the pair denominators.
"""
from collections import Counter, defaultdict
from pathlib import Path

from ..contracts import ContractError, read_json
from ..matching import BLAST_FIELDS
from ..sequences import accession_version, read_fasta, reverse_complement
from .common import read_tsv, sha256, write_json, write_tsv

WORD = 7
DATABASE_ORDER = ("core_nt_div10", "core_nt", "core_nt_x10")
MISS_FIELDS = ("database", "shard_id", "panel", "query_id", "subject_id", "variant", "context", "length", "complexity", "raw_hsps",
               "witness_span", "witness_mismatches", "witness_raw_score", "witness_longest_exact_run", "reason", "hsp_detail",
               "best_hsp_qspan", "best_hsp_pident", "best_hsp_gapopen")


def database_key(name: str):
    return (DATABASE_ORDER.index(name), name) if name in DATABASE_ORDER else (len(DATABASE_ORDER), name)


def longest_run(query: str, subject: str) -> int:
    best = current = 0
    for a, b in zip(query, subject):
        current = current + 1 if a == b else 0
        best = max(best, current)
    return best


def fixture_set(benchmark: Path, report: dict) -> Path:
    """The fixture set of a benchmark run: the directory named in benchmark.json (checked against the recorded
    specification digest), else the run directory itself."""
    if "fixtures" not in report:
        return benchmark
    root = Path(report["fixtures"])
    if sha256(root / "benchmark_specification.json") != report["specification_sha256"]:
        raise ContractError(f"Benchmark fixture set {root} differs from the one the run searched")
    return root


def benchmark_misses(benchmark: Path, database: str, arm: str = "permissive_diagnostic") -> tuple[list, Counter]:
    report = read_json(benchmark / "benchmark.json")
    if report.get("status") != "complete":
        raise ContractError(f"Benchmark {database} is not complete")
    fixture_root = fixture_set(benchmark, report)
    witnesses = {(row["shard_id"], row["query_id"], row["subject_id"]): row for row in read_tsv(fixture_root / "truth_witnesses.tsv")
                 if row["endpoint"] == "tolerant"}
    pairs = read_tsv(benchmark / "benchmark_pair_results.tsv")
    denominators = Counter()
    misses = []
    for row in pairs:
        if row["arm"] == arm and row["tolerant_truth"] == "1":
            denominators[row["panel"], row["error_class"]] += 1
            if row["tolerant_accepted"] == "0":
                misses.append(row)
    fixtures, hsps, rows = {}, {}, []
    for miss in misses:
        shard = miss["shard_id"]
        if shard not in fixtures:
            folder = fixture_root / "fixtures" / shard
            fixtures[shard] = (read_fasta(folder / "queries.fasta"), read_fasta(folder / "subjects.fasta", allow_ambiguity=True))
        queries, subjects = fixtures[shard]
        witness = witnesses[shard, miss["query_id"], miss["subject_id"]]
        qs, qe, ss, se = (int(witness[key]) for key in ("qstart", "qend", "sstart", "send"))
        query = queries[miss["query_id"]][qs - 1:qe]
        subject = subjects[miss["subject_id"]][min(ss, se) - 1:max(ss, se)]
        if witness["strand"] == "-":
            subject = reverse_complement(subject)
        run = longest_run(query, subject)
        matches, mismatches = int(witness["matches"]), int(witness["mismatches"])
        raw = int(miss["raw_hsps"])
        reason = "hsp_returned_not_qualifying" if raw > 0 else "no_hsp_seed_present" if run >= WORD else "no_hsp_no_seed"
        detail, best = "", None
        if raw > 0:
            key = shard
            if key not in hsps:
                hsps[key] = defaultdict(list)
                for line in (benchmark / "searches" / shard / arm / "blast.tsv").read_text(encoding="utf-8").splitlines():
                    if line and not line.startswith("#"):
                        fields = dict(zip(BLAST_FIELDS, line.split("\t")))
                        hsps[key][fields["qseqid"], accession_version(fields["sseqid"])].append(fields)
            candidates = hsps[key][miss["query_id"], miss["subject_id"]]
            if not candidates:
                raise ContractError(f"Benchmark {database}: recorded raw HSPs of a missed pair are absent from the arm output")
            best = max(candidates, key=lambda fields: float(fields["bitscore"]))
            span = (int(best["qend"]) - int(best["qstart"]) + 1) / int(best["qlen"])
            detail = ("gapped_hsp" if int(best["gapopen"]) > 0 else "hsp_span_below_90pct" if span < 0.9 else
                      "hsp_identity_below_95pct" if float(best["pident"]) < 95 else "hsp_other_predicate")
        rows.append(dict(database=database, shard_id=shard, panel=miss["panel"], query_id=miss["query_id"], subject_id=miss["subject_id"],
                         variant=miss["error_class"], context=miss["context"], length=int(miss["length"]), complexity=miss["complexity"],
                         raw_hsps=raw, witness_span=int(witness["span"]), witness_mismatches=mismatches,
                         witness_raw_score=matches - 3 * mismatches, witness_longest_exact_run=run, reason=reason, hsp_detail=detail,
                         best_hsp_qspan=f"{best['qstart']}-{best['qend']}" if best else "", best_hsp_pident=best["pident"] if best else "",
                         best_hsp_gapopen=best["gapopen"] if best else ""))
    return rows, denominators


def write_benchmark_misses(destination: Path, *, benchmarks: dict[str, Path], arm: str = "permissive_diagnostic") -> dict:
    destination.mkdir(parents=True, exist_ok=False)
    if not benchmarks:
        raise ContractError("No benchmark directories")
    all_rows, denominator_rows, summary_rows, summary = [], [], [], {}
    for database in sorted(benchmarks, key=database_key):
        rows, denominators = benchmark_misses(Path(benchmarks[database]), database, arm)
        if sum(1 for row in rows) > sum(denominators.values()):
            raise ContractError(f"Benchmark {database}: more misses than tolerant truth-positive pairs")
        all_rows.extend(rows)
        denominator_rows.extend(dict(database=database, panel=panel, variant=variant, tolerant_positive_pairs=count,
                                     misses=sum(1 for row in rows if (row["panel"], row["variant"]) == (panel, variant)))
                                for (panel, variant), count in sorted(denominators.items()))
        for keys in (("reason",), ("hsp_detail",), ("variant", "hsp_detail"), ("variant", "reason"), ("panel", "variant"), ("length", "reason"),
                     ("complexity", "reason")):
            counts = Counter(tuple(row[key] for key in keys) for row in rows)
            summary_rows.extend(dict(database=database, dimension="+".join(keys), stratum="|".join(map(str, key)), misses=count)
                                for key, count in sorted(counts.items()))
        summary[database] = dict(arm=arm, misses=len(rows), tolerant_positive_pairs=sum(denominators.values()),
                                 by_reason=dict(sorted(Counter(row["reason"] for row in rows).items())),
                                 by_hsp_detail=dict(sorted(Counter(row["hsp_detail"] for row in rows).items())),
                                 by_variant=dict(sorted(Counter(row["variant"] for row in rows).items())))
    write_tsv(destination / "misses.tsv", all_rows, MISS_FIELDS)
    write_tsv(destination / "miss_summary.tsv", summary_rows, ("database", "dimension", "stratum", "misses"))
    write_tsv(destination / "denominators.tsv", denominator_rows, ("database", "panel", "variant", "tolerant_positive_pairs", "misses"))
    write_json(destination / "summary.json", summary)
    return summary
