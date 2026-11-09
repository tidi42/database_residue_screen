"""Stage ``figure_sources``: every table read by the R figure scripts (reporting.figure_sources).

Upstream mapping: ``metrics`` and ``reports`` (this analysis), ``taxonomy`` (final analysis), ``design`` (physical queries,
logical queries), ``controls`` (control generator output), ``search`` (raw HSPs for the alignment-quality panels),
``inference`` (precision tables), ``benchmark`` (benchmark runs per database size), ``record_cases`` (selected cases).
"""
from ...reporting.figure_sources import build_figure_sources
from .report_inputs import benchmark_directories, directory


def figure_sources(context):
    """Tables read by the R figure scripts"""
    upstream = {name: attempt.output for name, attempt in context.upstream.items()}
    result = build_figure_sources(
        context.output, metrics=directory("metrics", upstream["metrics"]), reports=directory("reports", upstream["reports"]),
        analysis=directory("taxonomy", upstream["taxonomy"]), design=directory("design", upstream["design"]),
        controls=directory("controls", upstream["controls"]), search=directory("search", upstream["search"]),
        inference=directory("inference", upstream["inference"]), benchmarks=benchmark_directories(upstream["benchmark"]),
        record_cases=directory("record_cases", upstream["record_cases"]))
    context.log.info("figure sources: %d tables", result["tables"])
    return result
