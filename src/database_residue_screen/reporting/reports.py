"""Publication reports: one directory per report, built from explicit input directories (see docs/reports.md).

  components/             component inventory, within-set distances, barcode-segment readout
  cap_review/             cap-saturated production searches with their review states
  target_limit/           per-parent detection and breadth at 500 and 2,000 targets and in production
  eligibility/            accession eligibility flow, taxonomy context, field population
  positional/             positional sensitivity and record type x taxonomy
  nesting/                per-parent events versus distinct loci
  benchmark_misses/       persistent tolerant misses of the benchmark at every database size
  census/                 flank-census class summary
  exact_scan_validation/  exact scan against the strict BLAST endpoint
  inference/              conditional bounds and control precision
  record_cases/           sequence-level record-case table

Inputs that are not given are listed as skipped in reports.json; the study stage passes every input.
"""
from pathlib import Path

from ..contracts import ContractError, read_json, read_table
from ..outcomes import OUTCOME_SCHEMA
from .benchmark_misses import write_benchmark_misses
from .census_summary import write_census
from .common import load_analysis, new_directory, sha256, verify_manifest, write_json, write_manifest
from .components import write_components
from .eligibility_flow import write_eligibility
from .exact_scan_validation import write_scan_validation
from .inference import write_inference
from .metrics import CONTRACT as METRICS_CONTRACT, read_logical, read_statuses
from .nesting import write_nesting
from .positions import write_positions
from .record_cases import write_record_cases
from .target_limit import write_cap_review, write_target_limit

CONTRACT = "publication_reports_v1"


def build_reports(destination: Path, *, metrics: Path, analysis: Path, design: Path, search: Path, metadata: Path,
                  panel_definition: Path | None = None, insdc: Path | None = None, positional_table: Path | None = None,
                  cap_review: Path | None = None, target_limit: Path | None = None, exact_scan: Path | None = None,
                  scan_validation: Path | None = None, benchmarks: dict | None = None, inference: Path | None = None,
                  specification: Path | None = None, record_cases: Path | None = None, census: Path | None = None,
                  validation: Path | None = None,
                  benchmark_arm: str = "permissive_diagnostic", draw_prefixes=(5, 10, 15)) -> dict:
    analysis, design, search, metadata, metrics = map(Path, (analysis, design, search, metadata, metrics))
    report = load_analysis(analysis)
    if verify_manifest(metrics, METRICS_CONTRACT).get("analysis_sha256") != sha256(analysis / "analysis.json"):
        raise ContractError("The metrics belong to a different analysis")
    policy = report["analysis_policy"]
    logical = read_logical(analysis / "logical_queries.tsv")
    statuses = read_statuses(analysis / "search_status.tsv")
    destination = new_directory(destination)
    outcomes = None

    def analysis_outcomes():
        nonlocal outcomes
        if outcomes is None:
            outcomes = read_table(analysis / "outcomes.tsv", OUTCOME_SCHEMA)
        return outcomes

    summaries, skipped = {}, []
    if panel_definition is not None:
        summaries["components"] = write_components(destination / "components", definition=Path(panel_definition),
                                                   physical_fasta=design / "physical_queries.fasta", logical=logical, search=search, policy=policy)
    else:
        skipped.append("components")
    if cap_review is not None:
        summaries["cap_review"] = write_cap_review(destination / "cap_review", cap_review=Path(cap_review), logical=logical,
                                                   outcomes=analysis_outcomes(), statuses=statuses)
    else:
        skipped.append("cap_review")
    if target_limit is not None:
        summaries["target_limit"] = write_target_limit(destination / "target_limit", table=Path(target_limit), outcomes=analysis_outcomes(),
                                                       logical=logical, production_limit=report["search_settings"]["max_target_seqs"])
    else:
        skipped.append("target_limit")
    summaries["eligibility"] = write_eligibility(destination / "eligibility", analysis=analysis, policy=policy)
    summaries["positional"] = write_positions(destination / "positional", analysis=analysis, metadata=metadata,
                                              insdc=Path(insdc) if insdc is not None else None,
                                              positional_table=Path(positional_table) if positional_table is not None else None)
    summaries["nesting"] = write_nesting(destination / "nesting", analysis=analysis, policy=policy)
    if benchmarks:
        summaries["benchmark_misses"] = write_benchmark_misses(destination / "benchmark_misses",
                                                               benchmarks={name: Path(path) for name, path in benchmarks.items()}, arm=benchmark_arm)
    else:
        skipped.append("benchmark_misses")
    if exact_scan is not None and scan_validation is not None:
        summaries["exact_scan_validation"] = write_scan_validation(destination / "exact_scan_validation", exact_scan=Path(exact_scan),
                                                                   scan_validation=Path(scan_validation), logical=logical)
    else:
        skipped.append("exact_scan_validation")
    if inference is not None:
        summaries["inference"] = write_inference(destination / "inference", inference=Path(inference),
                                                 specification=Path(specification) if specification is not None else None,
                                                 draw_prefixes=draw_prefixes)
    else:
        skipped.append("inference")
    if record_cases is not None:
        summaries["record_cases"] = write_record_cases(destination / "record_cases", record_cases=Path(record_cases),
                                                       physical_fasta=design / "physical_queries.fasta",
                                                       search_fasta=design / "search_queries.fasta", policy=policy)
    else:
        skipped.append("record_cases")
    if census is not None:
        summaries["census"] = write_census(destination / "census", census=Path(census))
    else:
        skipped.append("census")
    validation_record = read_json(Path(validation)) if validation is not None else None
    result = dict(contract=CONTRACT, reports=sorted(summaries), skipped=skipped, summaries=summaries,
                  validation=None if validation_record is None else {key: validation_record[key] for key in ("valid", "result_contract")
                                                                     if key in validation_record})
    write_json(destination / "reports.json", result)
    inputs = {"analysis.json": analysis / "analysis.json", "metrics/manifest.json": metrics / "manifest.json",
              "search/complete.json": search / "complete.json", "design/panel.json": design / "panel.json"}
    for name, path in (("panel_definition", panel_definition), ("cap_review", cap_review), ("target_limit", target_limit),
                       ("positional_table", positional_table), ("specification", specification), ("validation", validation)):
        if path is not None:
            inputs[name] = Path(path)
    for name, directory, marker in (("insdc", insdc, "insdc.json"), ("exact_scan", exact_scan, "scan.json"),
                                    ("inference", inference, "precision.json"), ("record_cases", record_cases, "record_cases.tsv"),
                                    ("census", census, "census.json")):
        if directory is not None and (Path(directory) / marker).is_file():
            inputs[f"{name}/{marker}"] = Path(directory) / marker
    for name, directory in sorted((benchmarks or {}).items()):
        inputs[f"benchmark_{name}/benchmark.json"] = Path(directory) / "benchmark.json"
    write_manifest(destination, CONTRACT, inputs, "Publication reports of one completed analysis; descriptive, no new searches",
                   analysis_sha256=sha256(analysis / "analysis.json"))
    return result
