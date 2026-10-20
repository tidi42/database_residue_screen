"""Cap-saturation review (protocol v2, R07).

For every query whose search reached the target cap, the weakest retained
subject (lowest best bit score among retained subjects) is compared with the
bit score of the minimum qualifying strict and tolerant raw scores, computed
from the engine's lambda and K for that query length. If the weakest retained
subject still scores at or above a qualifying score, additional qualifying
subjects may lie beyond the cap: the endpoint is cap-censored. Otherwise the
retained list reaches below every qualifying score; that is consistent with
complete retention of qualifying subjects but is not a proof, because BLAST's
target limit is applied during the search. An optional exact scan supplies the
cap-free strict count for comparison.
"""
import math
from collections import defaultdict
from pathlib import Path

from .contracts import ContractError, read_json, table_schema
from .matching import strict_hsp, tolerant_hsp
from .search import iter_hsps, validate_completed_batch
from .sequences import read_fasta

REVIEW_SCHEMA = table_schema(("query_id", "strict_status", "tolerant_status", "interpretation"),
                             ("query_length", "retained_subjects", "strict_subjects", "tolerant_subjects", "target_limit",
                              "minimum_strict_score", "minimum_tolerant_score", "exact_scan_accessions"),
                             ("weakest_retained_bitscore", "strict_threshold_bits", "tolerant_threshold_bits"))


def bits(raw: float, row: dict) -> float:
    return (row["lambda_value"] * raw - math.log(row["kappa"])) / math.log(2)


def review_saturated(search: Path, exact_accessions: dict[str, set] | None = None) -> list[dict]:
    request, complete = read_json(search / "request.json"), read_json(search / "complete.json")
    calibration = {row["query_length"]: row for row in request["length_calibration"]["lengths"]}
    limit = request["settings"]["max_target_seqs"]
    rows = []
    for batch in complete["batches"]:
        directory = search / batch["directory"]
        validate_completed_batch(directory, request)
        queries = read_fasta(directory / "queries.fasta")
        best, strict, tolerant = defaultdict(dict), defaultdict(set), defaultdict(set)
        for hsp in iter_hsps(directory / "blast.tsv", queries):
            query, subject = hsp["qseqid"], hsp["sseqid"]
            best[query][subject] = max(best[query].get(subject, 0.0), float(hsp["bitscore"]))
            if strict_hsp(hsp):
                strict[query].add(subject)
            if tolerant_hsp(hsp):
                tolerant[query].add(subject)
        for query, subjects in best.items():
            if len(subjects) < limit:
                continue
            length = len(queries[query])
            if length not in calibration:
                raise ContractError("No calibration row for a saturated query length")
            row = calibration[length]
            weakest = min(subjects.values())
            strict_bits, tolerant_bits = bits(length, row), bits(row["minimum_tolerant_score"], row)
            status = {name: ("cap_censored_possible" if weakest >= threshold - 0.5 else "qualifying_range_reached")
                      for name, threshold in (("strict", strict_bits), ("tolerant", tolerant_bits))}
            rows.append(dict(query_id=query, query_length=length, retained_subjects=len(subjects), strict_subjects=len(strict[query]),
                             tolerant_subjects=len(tolerant[query]), target_limit=limit, weakest_retained_bitscore=weakest,
                             minimum_strict_score=length, minimum_tolerant_score=row["minimum_tolerant_score"],
                             strict_threshold_bits=strict_bits, tolerant_threshold_bits=tolerant_bits,
                             strict_status=status["strict"], tolerant_status=status["tolerant"],
                             exact_scan_accessions=len(exact_accessions.get(query, set())) if exact_accessions is not None else -1,
                             interpretation="Bit scores are reported rounded by BLAST; a 0.5-bit margin is applied. "
                                            "qualifying_range_reached is consistent with, not proof of, complete retention"))
    return sorted(rows, key=lambda row: row["query_id"])
