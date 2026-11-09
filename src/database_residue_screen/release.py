import shutil
import tempfile
from pathlib import Path

from .analysis import ANNOTATED_EVENT_SCHEMA, EVENT_CONTRACT, RESULT_CONTRACT, derive_event_tables, read_events
from .contracts import ContractError, read_json, read_table
from .design import digest
from .eligibility import event_eligibility
from .example import save_json
from .outcomes import CONTRAST_SCHEMA, ELIGIBILITY_SCHEMA, MEMBERSHIP_SCHEMA, OUTCOME_SCHEMA, build_outcomes, project_real_events
from .panels import LOGICAL_SCHEMA
from .policy import load_policy
from .search import STATUS_SCHEMA, object_digest
from .taxonomy import COMPOSITION_SCHEMA, TAXONOMY_SCHEMA, load_taxonomy, taxonomy_composition, taxonomy_rows


def validate_analysis(directory: Path, require_production: bool = False) -> dict:
    report = read_json(directory / "analysis.json")
    if report.get("result_contract") != RESULT_CONTRACT or report.get("event_contract") != EVENT_CONTRACT:
        raise ContractError(f"Expected a completed analysis with result contract {RESULT_CONTRACT} and event contract {EVENT_CONTRACT}")
    acceptance = report.get("database_acceptance")
    if require_production and (report.get("search_mode") != "production" or not acceptance or acceptance.get("full_integrity_passed") is not True):
        raise ContractError("Production result input requires an accepted database search; working fixtures cannot be promoted")
    names = set(report["files"])
    actual = {str(path.relative_to(directory)) for path in directory.rglob("*") if path.is_file()}
    if actual != names | {"analysis.json"}:
        raise ContractError("Result directory contains missing or undeclared files")
    for name, checksum in report["files"].items():
        if Path(name).name != name or digest(directory / name) != checksum:
            raise ContractError("Result artifact checksum mismatch: " + name)
    if report.get("source_search"):
        search = Path(report["source_search"])
        if (search / "quarantine.json").exists():
            raise ContractError("Source search is quarantined")
        if (digest(search / "complete.json") != report["search_completion_sha256"] or
                digest(search / "request.json") != report["search_request_sha256"]):
            raise ContractError("Search completion or request changed since analysis")
        request, complete = read_json(search / "request.json"), read_json(search / "complete.json")
        if (request.get("mode") != report.get("search_mode") or request.get("database_acceptance") != acceptance
                or complete["search_id"] != object_digest(request) or report["search_id"] != complete["search_id"]
                or not complete["all_queries_completed"]):
            raise ContractError("Result provenance disagrees with the original search request")
        if report.get("length_calibration"):
            if (report.get("length_calibration_source") or dict(kind="runtime_request"))["kind"] != "runtime_request":
                raise ContractError("Result testability calibration must come from its search request")
            if request.get("length_calibration") != report["length_calibration"]:
                raise ContractError("Result testability calibration disagrees with its search")
    elif require_production:
        raise ContractError("Production result has no original search evidence")
    if require_production:
        path = Path(acceptance["file"])
        if digest(path) != acceptance["sha256"] or (path.parent / "quarantine.json").exists() or (path.parent / "failed.json").exists():
            raise ContractError("Production database acceptance changed or was quarantined")
        accepted = read_json(path)
        if (accepted.get("accepted") is not True or accepted["inventory_sha256"] != request["database_evidence_sha256"]
                or accepted["database_prefix"] != request["database_prefix"]
                or accepted["binaries"]["blastn"]["sha256"] != request["executable_sha256"]):
            raise ContractError("Production acceptance belongs to a different database or executable")
        for name, checksum in accepted["files"].items():
            if Path(name).name != name or digest(path.parent / name) != checksum:
                raise ContractError("Production database acceptance evidence changed")
    policy = load_policy(directory / "analysis_policy.json")
    if policy != report["analysis_policy"]:
        raise ContractError("Result policy disagrees with its manifest")
    if report.get("length_calibration") and policy.get("detectability") != report["length_calibration"]["lengths"]:
        raise ContractError("Result policy testability rows disagree with the recorded calibration")
    events = read_events(directory, report)
    logical = read_table(directory / "logical_queries.tsv", LOGICAL_SCHEMA)
    statuses = read_table(directory / "search_status.tsv", STATUS_SCHEMA)
    if (report["logical_queries"] != len(logical) or report["distinct_search_sequences"] != len(statuses)
            or report["canonical_parents"] != len({row["parent_id"] for row in logical if row["role"] == "real"})):
        raise ContractError("Result inventory counts disagree with the logical panel")
    outcomes, contrasts, memberships, eligibility = build_outcomes(events, logical, statuses, policy, report["draws_per_parent"],
                                                                   report["repeated_control_seed_values"])
    from .sensitivity import (COMMON_RULE_SCHEMA, SENSITIVITY_SCHEMA, TESTABILITY_SCHEMA, common_rule_rows, length_sensitivity,
                              length_testability_rows)
    if "length_sensitivity.tsv" in names and read_table(directory / "length_sensitivity.tsv", SENSITIVITY_SCHEMA) != length_sensitivity(outcomes, policy, report["draws_per_parent"]):
        raise ContractError("Length sensitivity disagrees with recomputed matched cohorts")
    if "common_rule_contrasts.tsv" in names and read_table(directory / "common_rule_contrasts.tsv", COMMON_RULE_SCHEMA) != common_rule_rows(outcomes, policy, report["draws_per_parent"]):
        raise ContractError("Common strict/tolerant contrasts disagree with recomputed matched cohorts")
    if "length_testability.tsv" in names and read_table(directory / "length_testability.tsv", TESTABILITY_SCHEMA) != length_testability_rows(policy):
        raise ContractError("Per-length testability table disagrees with the analysis policy")
    real_events = project_real_events(events, memberships)
    for name, rows, definition in (("outcomes.tsv", outcomes, OUTCOME_SCHEMA), ("contrasts.tsv", contrasts, CONTRAST_SCHEMA),
                                   ("event_memberships.tsv", memberships, MEMBERSHIP_SCHEMA), ("eligibility.tsv", eligibility, ELIGIBILITY_SCHEMA),
                                   ("real_events.tsv", real_events, ANNOTATED_EVENT_SCHEMA)):
        if read_table(directory / name, definition) != rows:
            raise ContractError(name + " disagrees with recomputed logical outcomes and eligibility")
    groups = sorted({row["group"] for row in logical if row["role"] == "real"})
    decisions = event_eligibility(events, policy["record_populations"])
    with tempfile.TemporaryDirectory() as temporary:
        comparison = Path(temporary)
        for population in policy["record_populations"]:
            selected = [row for row in real_events if decisions[row["accession_version"], population][0] == "included"]
            derived = derive_event_tables(selected, groups, comparison, population + "_", policy)
            if report["derived"][population] != derived:
                raise ContractError("Derived summary disagrees with the event release")
            for suffix in ("density.tsv", "positions.tsv", "overlap.tsv"):
                name = population + "_" + suffix
                if digest(comparison / name) != report["files"][name]:
                    raise ContractError("Derived figure table disagrees with source events: " + name)
    if "taxonomy_snapshot" in report:
        snapshot = report["taxonomy_snapshot"]
        taxonomy, evidence = load_taxonomy(Path(snapshot["directory"]))
        if digest(Path(snapshot["directory"]) / "taxonomy.json") != snapshot["manifest_sha256"] or evidence["files"] != snapshot["files"]:
            raise ContractError("Taxonomy snapshot differs from the annotated release")
        ranks = tuple(report["taxonomy_ranks"])
        assignments = taxonomy_rows(events, taxonomy, ranks, snapshot["manifest_sha256"])
        composition = taxonomy_composition(real_events, assignments, groups, policy, ranks, decisions)
        if read_table(directory / "accession_taxonomy.tsv", TAXONOMY_SCHEMA) != assignments or read_table(directory / "taxonomic_composition.tsv", COMPOSITION_SCHEMA) != composition:
            raise ContractError("Taxonomy assignments or accession denominators disagree with the pinned snapshot")
    return dict(valid=True, result_contract=report["result_contract"], analysis_sha256=digest(directory / "analysis.json"),
                production_evidence_checked=require_production, files=len(names), original_parents=report["canonical_parents"],
                policy=policy, unresolved_outcomes=sum(row["unresolved_subjects"] > 0 for row in outcomes),
                scope="Artifact identity and endpoint arithmetic; not final scientific acceptance, exhaustive retrieval or null adequacy")


def export_figure_inputs(analysis: Path, destination: Path) -> dict:
    validation = validate_analysis(analysis, require_production=True)
    source = read_json(analysis / "analysis.json")
    tables = ["outcomes.tsv", "contrasts.tsv", "eligibility.tsv"]
    tables += [name for name in ("length_sensitivity.tsv", "common_rule_contrasts.tsv", "length_testability.tsv") if name in source["files"]]
    tables += [population + "_" + suffix for population in source["analysis_policy"]["record_populations"]
               for suffix in ("density.tsv", "positions.tsv", "overlap.tsv")]
    if "taxonomy_snapshot" in source:
        tables.append("taxonomic_composition.tsv")
    destination.mkdir(parents=True, exist_ok=False)
    for name in tables:
        shutil.copyfile(analysis / name, destination / name)
    report = dict(contract="reference_figure_inputs_v2", source_analysis_sha256=validation["analysis_sha256"],
                  source_directory=str(analysis.resolve()), policy=source["analysis_policy"],
                  database_acceptance=source["database_acceptance"],
                  interpretation="Production-origin source tables, not automatic acceptance of every claim or panel",
                  uncertainty="Use explicit fixed-panel and Monte Carlo limits; no biological sampling intervals or significance stars supplied",
                  files={name: digest(destination / name) for name in tables})
    save_json(destination / "figure_inputs.json", report)
    return report