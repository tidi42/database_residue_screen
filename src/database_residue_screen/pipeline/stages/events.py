"""Stage `events`: join the archived exact-version metadata to the events and recompute the outcome tables."""
from ...metadata import annotate_events


def events(context):
    """Metadata-annotated analysis (record type, topology, molecule, origin) from the archived snapshot."""
    report = annotate_events(context.upstream["analysis"].output / "analysis", context.upstream["metadata"].output / "metadata",
                             context.output / "analysis")
    return {key: report[key] for key in ("logical_queries", "canonical_parents", "taxonomy_status")}
