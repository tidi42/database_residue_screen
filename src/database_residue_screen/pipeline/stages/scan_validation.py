"""Stage `scan_validation`: every strict BLAST HSP must be found by the exact scan; scan-only occurrences are classified."""
from ...contracts import ContractError
from ...scan_validation import validate_scan


def scan_validation(context):
    """Concordance of the production search with the exact scan on the scanned search sequences and records."""
    summary = validate_scan(context.upstream["search"].output / "search", context.upstream["exact_scan"].output / "scan",
                            context.output / "scan_validation")
    if not summary["accepted"] and context.config["scan_validation"]["require_acceptance"]:
        raise ContractError(f"Exact scan not accepted: {summary['blast_only']} strict BLAST HSPs were not found by the scan")
    return {key: summary[key] for key in ("volume_scope", "scanned_search_sequences", "blast_strict_hsps", "matched", "blast_only",
                                          "scan_only", "scan_only_target_cap", "scan_only_unexplained", "accepted")}
