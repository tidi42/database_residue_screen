"""Stage `census`: flank census of the retained strict reference windows with the run-time elements and catalogue."""
from pathlib import Path

from ...census import load_elements, run_census
from ...contracts import read_json


def census(context):
    """Retrieve, verify and classify every strict window of the declared census families."""
    panel = context.upstream["panel"].output
    evidence = context.upstream["database"].output / "inventory" / "database_identity.json"
    taxonomy = context.upstream["taxonomy"].output
    report = run_census(taxonomy / "analysis", context.upstream["design"].output / "search_panel", Path(read_json(evidence)["database_prefix"]),
                        context.tool("blastdbcmd"), context.config["census"], load_elements(panel / "census_elements.json"),
                        panel / "census_catalogue.tsv", context.output / "census", taxonomy / "snapshot")
    return {key: report[key] for key in ("windows", "status_counts", "integrity_complete")}
