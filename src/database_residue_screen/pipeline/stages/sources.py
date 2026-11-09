"""Stage `sources`: verify (and fetch) the pinned panel sources and snapshot them into the attempt."""
from pathlib import Path

from ...constructs import collect_sources


def sources(context):
    """Verify every pinned source by SHA-256 (fetching missing ones) and snapshot it for the panel build."""
    section = context.config["sources"]
    report = collect_sources(Path(section["manifest"]), Path(section["directory"]), context.output / "sources",
                             fetch=section["fetch"], allow_drift=section["allow_source_drift"])
    for record in report["drift"]:
        context.log.warning("source %s drifted: SHA-256 %s instead of the pinned %s (accepted only if the panel reproduces)",
                            record["name"], record["sha256"], record["pinned_sha256"])
    return dict(sources=len(report["sources"]), status_counts=report["status_counts"], origin_counts=report["origin_counts"],
                drift=[record["name"] for record in report["drift"]])
