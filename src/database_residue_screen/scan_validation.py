"""Exact-scan acceptance against a completed production BLAST search (supplementary validation).

For every search sequence covered by the exhaustive exact scan (both strands, length >= the scan minimum):
  1. every strict BLAST HSP (full-length, 100% identity, ungapped) must be an exact scan occurrence at the same record
     coordinates and strand (palindromes: strand-free); otherwise the scan is not accepted;
  2. every scan-only occurrence is classified: 'target_cap' (the search reached the target limit and BLAST did not return
     that record), 'overlap_with_reported_hsp' (BLAST returned a strict HSP of the same query on that record overlapping
     the occurrence) or 'unexplained' (reported as a discrepancy, never dropped).
Accession-level strict counts are compared per query. A scan restricted to some volumes is compared on the records of
those volumes only, and the scope is reported. The BLAST endpoints remain primary; the scan is a second implementation.
"""
import json
from collections import Counter, defaultdict
from pathlib import Path

from .contracts import ContractError, read_json, read_table, table_schema, write_table
from .design import digest
from .exact_scan import read_occurrences, reverse_complement, scanned_accessions
from .matching import strict_hsp
from .search import STATUS_SCHEMA, iter_hsps, object_digest, validate_completed_batch
from .sequences import accession_version, read_fasta

CONCORDANCE_SCHEMA = table_schema(("query_id", "target_limit_reached", "accessions_equal"),
                                  ("length", "blast_strict_hsps", "scan_occurrences", "matched", "blast_only", "scan_only",
                                   "scan_only_target_cap", "scan_only_overlap", "scan_only_unexplained", "blast_strict_accessions",
                                   "scan_accessions", "blast_strict_hsps_outside_scan"))
DISCORDANT_SCHEMA = table_schema(("query_id", "accession_version", "strand", "side", "classification"), ("start", "end"))


def pdb_accession(raw: str) -> str | None:
    """BLAST reports PDB chains as 'gi|N|pdb|4FCY|C'; blastdbcmd %a prints the same record as '4FCY_C'."""
    parts = raw.split("|")
    if "pdb" in parts:
        index = parts.index("pdb")
        if index + 2 < len(parts) and parts[index + 1] and parts[index + 2]:
            return f"{parts[index + 1]}_{parts[index + 2]}"
    return None


def hsp_key(row: dict, palindrome: bool) -> tuple:
    try:
        accession = accession_version(row["sseqid"])
    except ContractError:
        accession = pdb_accession(row["sseqid"]) or "unresolved:" + row["sseqid"]
    strand = "+" if palindrome or row["sstart"] < row["send"] else "-"
    return accession, min(row["sstart"], row["send"]), max(row["sstart"], row["send"]), strand


def compare(query: str, blast: set, scan: set, returned_accessions: set, capped: bool) -> tuple[dict, list[dict]]:
    """Compare one query's strict BLAST HSP keys with its scan occurrence keys (accession, start, end, strand)."""
    details = []
    for key in sorted(blast - scan):
        details.append(dict(query_id=query, accession_version=key[0], start=key[1], end=key[2], strand=key[3], side="blast_only",
                            classification="blast_strict_hsp_not_found_by_scan"))
    by_accession = defaultdict(list)
    for key in blast:
        by_accession[key[0]].append(key)
    counts = Counter()
    for key in sorted(scan - blast):
        if any(other[1] <= key[2] and key[1] <= other[2] for other in by_accession.get(key[0], ())):
            reason = "overlap_with_reported_hsp"
        elif capped and key[0] not in returned_accessions:
            reason = "target_cap"
        else:
            reason = "unexplained"
        counts[reason] += 1
        details.append(dict(query_id=query, accession_version=key[0], start=key[1], end=key[2], strand=key[3], side="scan_only",
                            classification=reason))
    blast_accessions = {key[0] for key in blast}
    scan_accessions = {key[0] for key in scan}
    row = dict(query_id=query, target_limit_reached="yes" if capped else "no", blast_strict_hsps=len(blast), scan_occurrences=len(scan),
               matched=len(blast & scan), blast_only=len(blast - scan), scan_only=len(scan - blast),
               scan_only_target_cap=counts["target_cap"], scan_only_overlap=counts["overlap_with_reported_hsp"],
               scan_only_unexplained=counts["unexplained"], blast_strict_accessions=len(blast_accessions),
               scan_accessions=len(scan_accessions), accessions_equal="yes" if blast_accessions == scan_accessions else "no")
    return row, details


def validate_scan(search: Path, scan: Path, destination: Path) -> dict:
    """Compare a completed search with a completed exact scan of its search sequences; writes the concordance tables."""
    scan_report = read_json(scan / "scan.json")
    if not scan_report.get("totals_match_scanned_volumes", scan_report["totals_match_database"]):
        raise ContractError("Exact scan totals do not match the database evidence of its volumes")
    request, complete = read_json(search / "request.json"), read_json(search / "complete.json")
    if complete["search_id"] != object_digest(request) or not complete["all_queries_completed"]:
        raise ContractError("The search is incomplete or its completion disagrees with its request")
    if scan_report["queries_sha256"] != request["queries_sha256"]:
        raise ContractError("The exact scan covers different search sequences than the search")
    statuses = {row["query_id"]: row for row in read_table(search / "status.tsv", STATUS_SCHEMA)}
    minimum = scan_report["minimum_length"]
    within = scanned_accessions(scan)
    scan_keys, primary = defaultdict(set), {}
    for row in read_occurrences(scan):
        scan_keys[row["query_id"]].add((row["accession_version"], row["start"], row["end"], row["strand"]))
        for alias in filter(None, row["alias_accessions"].split(";")):
            primary[alias] = row["accession_version"]
    remapped, outside = Counter(), Counter()
    blast, returned, queries_seen = defaultdict(set), defaultdict(set), {}
    for batch in complete["batches"]:
        directory = search / batch["directory"]
        if digest(directory / "complete.json") != batch["completion_sha256"]:
            raise ContractError("Changed search batch completion: " + batch["directory"])
        validate_completed_batch(directory, request)
        queries = read_fasta(directory / "queries.fasta")
        if not any(len(sequence) >= minimum for sequence in queries.values()):
            continue
        queries_seen.update({name: sequence for name, sequence in queries.items() if len(sequence) >= minimum})
        palindromes = {name for name, sequence in queries.items() if sequence == reverse_complement(sequence)}
        for row in iter_hsps(directory / "blast.tsv", queries):
            if len(queries[row["qseqid"]]) < minimum:
                continue
            key = hsp_key(row, row["qseqid"] in palindromes)
            if pdb_accession(row["sseqid"]) == key[0]:
                remapped["pdb_chain"] += 1
            if key[0] in primary:
                remapped["alias_to_first_accession"] += 1
                key = (primary[key[0]],) + key[1:]
            if within is not None and key[0] not in within:
                if strict_hsp(row):
                    outside[row["qseqid"]] += 1
                continue
            returned[row["qseqid"]].add(key[0])
            if strict_hsp(row):
                blast[row["qseqid"]].add(key)
    if set(scan_keys) - set(queries_seen):
        raise ContractError("The scan reports occurrences for sequences outside the search")
    rows, details = [], []
    for query in sorted(queries_seen):
        row, detail = compare(query, blast[query], scan_keys[query], returned[query], statuses[query]["target_limit_reached"] == "yes")
        rows.append(dict(row, length=len(queries_seen[query]), blast_strict_hsps_outside_scan=outside[query]))
        details.extend(detail)
    destination.mkdir(parents=True, exist_ok=False)
    write_table(destination / "scan_concordance_by_query.tsv", rows, CONCORDANCE_SCHEMA)
    write_table(destination / "scan_discordant_occurrences.tsv", details, DISCORDANT_SCHEMA)
    uncapped = [row for row in rows if row["target_limit_reached"] == "no"]
    summary = dict(kind="exact_scan_validation_v1", volume_scope=scan_report.get("volume_scope", "all_volumes"),
                   scanned_volumes=scan_report.get("scanned_volumes", len(scan_report["volumes"])),
                   scanned_search_sequences=len(rows), minimum_length=minimum, scan_records=scan_report["records"],
                   scan_bases=scan_report["bases"], scan_deflines=scan_report.get("deflines"),
                   blast_keys_pdb_chain=remapped["pdb_chain"], blast_keys_alias_to_first_accession=remapped["alias_to_first_accession"],
                   blast_strict_hsps_outside_scanned_volumes=sum(outside.values()),
                   blast_strict_hsps=sum(r["blast_strict_hsps"] for r in rows), scan_occurrences=sum(r["scan_occurrences"] for r in rows),
                   matched=sum(r["matched"] for r in rows), blast_only=sum(r["blast_only"] for r in rows),
                   scan_only=sum(r["scan_only"] for r in rows), scan_only_target_cap=sum(r["scan_only_target_cap"] for r in rows),
                   scan_only_overlap=sum(r["scan_only_overlap"] for r in rows), scan_only_unexplained=sum(r["scan_only_unexplained"] for r in rows),
                   capped_searches=len(rows) - len(uncapped),
                   uncapped_searches_with_unequal_accessions=sum(r["accessions_equal"] == "no" for r in uncapped),
                   search_completion_sha256=digest(search / "complete.json"), scan_sha256=digest(scan / "scan.json"),
                   rule="Every strict BLAST HSP must be found by the scan; scan-only occurrences are classified, never dropped")
    summary["accepted"] = summary["blast_only"] == 0
    summary["all_scan_only_explained"] = summary["scan_only_unexplained"] == 0
    summary["files"] = {name: digest(destination / name) for name in ("scan_concordance_by_query.tsv", "scan_discordant_occurrences.tsv")}
    (destination / "scan_validation.json").write_text(json.dumps(summary, indent=1) + "\n", encoding="utf-8")
    return summary
