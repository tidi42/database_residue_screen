"""Stage `taxonomy`: pin the configured taxdump (source URL, retrieval time, hashes) and reconcile taxids."""
from pathlib import Path

from ...taxonomy import annotate_taxonomy, pin_taxonomy


def taxonomy(context):
    """Archive the taxdump snapshot and write the final annotated analysis with exact-rank assignments."""
    section = context.config["taxonomy"]
    pinned = pin_taxonomy(Path(section["directory"]), context.output / "snapshot", section["source_url"], section["retrieved_utc"],
                          section["release_date"])
    report = annotate_taxonomy(context.upstream["events"].output / "analysis", context.output / "snapshot", context.output / "analysis",
                               tuple(section["ranks"]))
    return dict(nodes=pinned["nodes"], retrieved_utc=pinned["retrieved_utc"], ranks=report["taxonomy_ranks"],
                taxonomy_status=report["taxonomy_status"], canonical_parents=report["canonical_parents"])
