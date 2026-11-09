"""Stage `metadata`: fresh exact-version ESummary retrieval, archived once; later stages parse only this snapshot."""
from ...contracts import ContractError
from ...metadata import fetch_metadata


def metadata(context):
    """Retrieve and archive the record summaries of every resolved accession of the analysis."""
    ncbi, section = context.config["ncbi"], context.config["metadata"]
    report = fetch_metadata(context.upstream["analysis"].output / "analysis", context.output / "metadata", section["batch_size"],
                            email=ncbi["email"], api_key=ncbi.get("api_key") or None)
    failed = [row for row in report["problems"] if row.get("reason") == "retrieval_failed"]
    if failed and section["require_complete_retrieval"]:
        raise ContractError(f"{len(failed)} metadata batches could not be retrieved; rerun the stage for a complete snapshot")
    return dict(requested_accessions=report["requested_accessions"], matched_exact_version_records=report["matched_exact_version_records"],
                missing_accessions=len(report["missing_accessions"]), problems=len(report["problems"]), failed_batches=len(failed),
                retrieval_started_utc=report["retrieval_started_utc"], retrieval_completed_utc=report["retrieval_completed_utc"])
