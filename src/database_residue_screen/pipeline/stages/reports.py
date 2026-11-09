"""Stage ``reports``: the publication reports of the final analysis (reporting.reports).

Upstream mapping (stage output layouts in report_inputs.LAYOUT):
  metrics          central metrics directory (must belong to the same analysis)
  taxonomy         final analysis directory (analysis.json)
  design           prepared panel (panel.json, physical_queries.fasta, search_queries.fasta)
  panel            the panel definition (segments column) named in panel.json, checked against its recorded SHA-256
  search           production search (complete.json)
  metadata         archived ESummary snapshot (BioProject IDs)
  insdc            INSDC header snapshot (insdc.json)
  positional       positional-sensitivity table (checked against the recomputation)
  cap_review       cap-review table
  target_limit     target-limit comparison table (target_limit.tsv)
  exact_scan       exact-scan directory (scan.json)
  scan_validation  exact-scan validation tables (scan_concordance_by_query.tsv)
  benchmark        benchmark runs (benchmark.json per database size)
  inference        precision assessment (precision.json) and its precision specification
  record_cases     record-case selection (record_cases.tsv)
  census           flank census (census.json)
  validation       analysis validation record (JSON with "valid")
"""
from ...reporting.reports import build_reports
from .report_inputs import benchmark_directories, directory, file, panel_definition, precision_specification


def reports(context):
    """Publication reports of the final analysis"""
    upstream = {name: attempt.output for name, attempt in context.upstream.items()}
    config = context.config["reports"]
    result = build_reports(
        context.output, metrics=directory("metrics", upstream["metrics"]), analysis=directory("taxonomy", upstream["taxonomy"]),
        design=directory("design", upstream["design"]), search=directory("search", upstream["search"]),
        metadata=directory("metadata", upstream["metadata"]), panel_definition=panel_definition(upstream["panel"]),
        insdc=directory("insdc", upstream["insdc"]), positional_table=file("positional", upstream["positional"]),
        cap_review=file("cap_review", upstream["cap_review"]), target_limit=file("target_limit", upstream["target_limit"]),
        exact_scan=directory("exact_scan", upstream["exact_scan"]), scan_validation=directory("scan_validation", upstream["scan_validation"]),
        benchmarks=benchmark_directories(upstream["benchmark"]), inference=directory("inference", upstream["inference"]),
        specification=precision_specification(upstream["inference"]), record_cases=directory("record_cases", upstream["record_cases"]),
        census=directory("census", upstream["census"]), validation=file("validation", upstream["validation"]),
        benchmark_arm=config["benchmark_arm"], draw_prefixes=tuple(config["draw_prefixes"]))
    context.log.info("reports: %s", ", ".join(result["reports"]))
    return dict(contract=result["contract"], reports=result["reports"], skipped=result["skipped"])
