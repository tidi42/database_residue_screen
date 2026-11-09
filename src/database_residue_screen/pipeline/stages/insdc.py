"""Stage `insdc`: dated, archived GenBank-header snapshot of the records of the positional events."""
from functools import partial

from ...contracts import ContractError
from ...insdc import fetch_batch, fetch_insdc_annotation, technologies_from_metadata
from ...positional import primary_positional_events


def insdc(context):
    """Fetch and archive the structured INSDC annotation (KEYWORDS, /PCR_primers) per exact accession.version."""
    ncbi, section = context.config["ncbi"], context.config["insdc"]
    events, _ = primary_positional_events(context.upstream["taxonomy"].output / "analysis")
    fetch = partial(fetch_batch, email=ncbi["email"], api_key=ncbi.get("api_key") or None)
    report = fetch_insdc_annotation(sorted({row["accession_version"] for row in events}),
                                    technologies_from_metadata(context.upstream["metadata"].output / "metadata"), context.output / "insdc",
                                    section["batch_size"], fetch=fetch, delay=section["delay_seconds"])
    failed = report["status_counts"].get("retrieval_failed", 0)
    if failed and section["require_complete_retrieval"]:
        raise ContractError(f"{failed} INSDC headers could not be retrieved; rerun the stage for a complete snapshot")
    return {key: report[key] for key in ("requested", "status_counts", "targeted_counts", "retrieval_started_utc", "created_utc")}
