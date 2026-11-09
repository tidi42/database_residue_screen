"""Stage `analysis`: reconcile the completed search with the search panel (events, outcomes, contrasts)."""
from ...analysis import analyse_search


def analysis(context):
    """Events, logical outcomes, matched contrasts and derived tables of the completed search."""
    report = analyse_search(context.upstream["design"].output / "search_panel", context.upstream["search"].output / "search",
                            context.output / "analysis")
    return {key: report[key] for key in ("logical_queries", "distinct_search_sequences", "canonical_parents", "search_mode")}
