"""Stage `cap_review`: weakest retained score versus the minimum qualifying scores of every target-limited query."""
import json
from collections import Counter, defaultdict

from ...cap_review import REVIEW_SCHEMA, review_saturated
from ...contracts import read_json, write_table
from ...design import digest
from ...exact_scan import read_occurrences


def cap_review(context):
    """Cap-saturation review of the production search with cap-free strict accession counts from the exact scan."""
    search, scan = context.upstream["search"].output / "search", context.upstream["exact_scan"].output / "scan"
    exact = defaultdict(set)
    for row in read_occurrences(scan):
        exact[row["query_id"]].add(row["accession_version"])
    rows = review_saturated(search, exact)
    write_table(context.output / "cap_review.tsv", rows, REVIEW_SCHEMA)
    scan_report = read_json(scan / "scan.json")
    statuses = Counter(f"{row['strict_status']}|{row['tolerant_status']}" for row in rows)
    report = dict(kind="cap_review_v1", search_completion_sha256=digest(search / "complete.json"), scan_sha256=digest(scan / "scan.json"),
                  target_limit=read_json(search / "request.json")["settings"]["max_target_seqs"], saturated_queries=len(rows),
                  status_counts=dict(statuses), exact_scan_scope=scan_report.get("volume_scope", "all_volumes"),
                  exact_scan_minimum_length=scan_report["minimum_length"],
                  files={"cap_review.tsv": digest(context.output / "cap_review.tsv")})
    (context.output / "cap_review.json").write_text(json.dumps(report, indent=1) + "\n", encoding="utf-8")
    return dict(saturated_queries=len(rows), status_counts=dict(statuses), exact_scan_scope=report["exact_scan_scope"])
