"""Stage `search`: the production BLAST search of every distinct search sequence (references and controls)."""
from .searching import resumable_search


def search(context):
    """Search all distinct search sequences of the search panel in checksummed, resumable batches."""
    section = context.config["search"]
    queries = context.upstream["design"].output / "search_panel" / "search_queries.fasta"
    report, resumed = resumable_search(context, queries, context.output / "search", section, section["target_limit"])
    return dict({key: report[key] for key in ("search_id", "submitted_queries", "completed_queries", "success_zero", "hsp_rows",
                                              "target_limit_reached", "testability_counts", "mode")},
                resumed_batches=len(resumed["copied_batches"]) if resumed else 0)
