"""Stage ``metrics``: central publication cohorts and denominators of the final analysis.

Upstream mapping: ``taxonomy`` -> the final (metadata- and taxonomy-annotated) analysis directory; ``design`` -> the prepared
panel (physical IDs); ``search`` -> the production search (archived failed or incomplete batch attempts).
"""
from ...reporting.metrics import compute_metrics
from .report_inputs import directory


def metrics(context):
    """Central publication cohorts and denominators of the final analysis"""
    upstream = context.upstream
    analysis = directory("taxonomy", upstream["taxonomy"].output)
    design = directory("design", upstream["design"].output)
    search = directory("search", upstream["search"].output)
    classes = context.config["metrics"].get("display_classes") or None
    variants = tuple((row["label"], row["testable_chance_limit"], row["informativeness_limit"])
                     for row in context.config["metrics"]["chance_limit_variants"])
    context.log.info("metrics of the final analysis in %s", upstream["taxonomy"].relative)
    summary = compute_metrics(analysis, context.output, design=design, search=search, display_classes=classes,
                              chance_limit_variants=variants)
    return dict(contract=summary["contract"], counts=summary["counts"], zero_searches=summary["zero_searches"]["distinct_searches"],
                display_classes=len(summary["display_classes"]))
