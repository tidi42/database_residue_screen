"""Exact-scan validation summary: the exhaustive exact scan against the strict BLAST endpoint.

Reads the exact-scan record (scan.json: per-volume totals against the accepted database dimensions) and the
scan_validation stage tables (per-search concordance and classified discordant occurrences). The validation is accepted
when the scan totals match the database and every strict BLAST HSP of a scanned search is an exact-scan occurrence;
scan-only occurrences are classified (target cap, overlap with a reported HSP, unexplained) and never dropped. The
summary is recomputed from the per-search table and must agree with the stage summary.
"""
from collections import Counter, defaultdict
from pathlib import Path

from ..contracts import ContractError, read_json
from .common import read_tsv, write_json, write_tsv

COUNTS = ("blast_strict_hsps", "scan_occurrences", "matched", "blast_only", "scan_only", "scan_only_target_cap", "scan_only_overlap",
          "scan_only_unexplained")


def write_scan_validation(destination: Path, *, exact_scan: Path, scan_validation: Path, logical: list[dict]) -> dict:
    destination.mkdir(parents=True, exist_ok=False)
    scan = read_json(exact_scan / "scan.json")
    concordance = read_tsv(scan_validation / "scan_concordance_by_query.tsv")
    discordant_path = scan_validation / "scan_discordant_occurrences.tsv"
    discordant = read_tsv(discordant_path) if discordant_path.is_file() else []
    stage_summary = read_json(scan_validation / "scan_summary.json") if (scan_validation / "scan_summary.json").is_file() else {}
    roles = defaultdict(set)
    for row in logical:
        roles[row["sequence_id"]].add(row["role"])
    totals = Counter()
    by_length = defaultdict(Counter)
    for row in concordance:
        role = roles.get(row["query_id"], set())
        key = int(row["length"]), ("reference_and_control" if role == {"real", "control"} else "reference" if "real" in role
                                   else "control" if role else "not_in_panel"), row["target_limit_reached"]
        by_length[key]["searches"] += 1
        by_length[key]["accessions_unequal"] += row["accessions_equal"] != "yes"
        for name in COUNTS:
            by_length[key][name] += int(row[name])
            totals[name] += int(row[name])
    capped = sum(row["target_limit_reached"] == "yes" for row in concordance)
    uncapped_unequal = sum(row["target_limit_reached"] == "no" and row["accessions_equal"] != "yes" for row in concordance)
    recomputed = dict(scanned_search_sequences=len(concordance), capped_searches=capped, uncapped_searches_with_unequal_accessions=uncapped_unequal,
                      **{name: totals[name] for name in COUNTS})
    disagreements = {name: (stage_summary[name], value) for name, value in recomputed.items() if name in stage_summary and stage_summary[name] != value}
    if disagreements:
        raise ContractError(f"Scan validation summary disagrees with its per-search table: {disagreements}")
    if discordant and Counter(row["side"] for row in discordant) != Counter({side: count for side, count in
                                                                            (("blast_only", totals["blast_only"]), ("scan_only", totals["scan_only"])) if count}):
        raise ContractError("Discordant occurrence list disagrees with the per-search counts")
    write_tsv(destination / "concordance_by_length.tsv",
              [dict(query_length=key[0], search_role=key[1], target_limit_reached=key[2], **{name: counts[name] for name in ("searches", *COUNTS,
                                                                                                                         "accessions_unequal")})
               for key, counts in sorted(by_length.items())],
              ("query_length", "search_role", "target_limit_reached", "searches", *COUNTS, "accessions_unequal"))
    classes = Counter((row["side"], row["classification"]) for row in discordant)
    write_tsv(destination / "discordant_classes.tsv", [dict(side=key[0], classification=key[1], occurrences=count)
                                                       for key, count in sorted(classes.items())], ("side", "classification", "occurrences"))
    accepted = totals["blast_only"] == 0 and bool(scan.get("totals_match_database"))
    summary = dict(records=scan.get("records"), bases=scan.get("bases"), expected_records=scan.get("expected_records"),
                   expected_bases=scan.get("expected_bases"), totals_match_database=scan.get("totals_match_database"),
                   minimum_length=scan.get("minimum_length"), scanned_queries=scan.get("scanned_queries"), deflines=scan.get("deflines"),
                   **recomputed, accepted=accepted, all_scan_only_explained=totals["scan_only_unexplained"] == 0,
                   stage_summary_checked=bool(stage_summary))
    rows = [dict(metric=name, value=value) for name, value in summary.items()]
    write_tsv(destination / "scan_validation_summary.tsv", rows, ("metric", "value"))
    write_json(destination / "summary.json", summary)
    return summary
