"""Figure sources: every table read by the R figure scripts (studies/paper3/figures), contract publication_figure_sources_v1.

Built from the central metrics, the reports and the run's stage outputs only; no searches and no sequences other than
run outputs. Union rows cover every panel condition. Excluded testability states count neither as 0 nor as 1.
Result figures show testable queries only: the alignment-quality tables cover the HSPs of strict-testable reference
sequences, and the control and query-property tables flag strict-testable parents (``strict_testable``).
Tables are written at full precision; booleans as True/False. manifest.json lists the SHA-256 of every table, the
inputs, the counts quoted in the figures and the positional policy.
"""
import math
import shutil
import statistics
from collections import Counter, defaultdict
from pathlib import Path

from ..analysis import orientation_key
from ..contracts import ContractError, read_json, read_table
from ..controls import check_draw, dinucleotides, draw_seed
from ..matching import strict_hsp, tolerant_hsp
from ..outcomes import OUTCOME_SCHEMA, matched_population
from ..positional import primary_positional_events
from ..precision import mean_interval
from ..sequences import canonical, canonical_id, read_fasta
from .common import (GLOBAL, GROUPS, group_key, load_analysis, new_directory, ordered_groups, read_tsv, sha256, verified_files,
                     verify_manifest, write_json, write_tsv)
from .metrics import CONTRACT as METRICS_CONTRACT, read_logical, read_statuses
from .raw_search import search_hsps
from .reports import CONTRACT as REPORTS_CONTRACT

CONTRACT = "publication_figure_sources_v1"
BOOLEAN = ("True", "False")
DETECTION_FIELDS = ("scope", "order", "block", "group", "label", "min_length", "max_length", "conditions", "members", "population", "endpoint",
                    "cohort", "planned_parents", "eligible_parents", "reference_positive", "reference_rate", "control_draws", "control_positive",
                    "control_rate", "censored", "chance_limited", "partially_censored", "unresolved", "incomplete", "other_inapplicable",
                    "target_limited")


def write(path: Path, rows: list[dict], fields=None) -> None:
    write_tsv(path, rows, fields, boolean=BOOLEAN)


def rule_testable_parents(outcomes: list[dict], endpoint: str = "strict") -> set[str]:
    """Parents whose reference search is operational under the rule (a priori testable, all records)."""
    return {row["parent_id"] for row in outcomes if row["role"] == "real" and row["endpoint"] == endpoint
            and row["population"] == "all_records" and row["applicability"] == "operational_retained_search"}


def detection_rows(metrics: Path) -> list[dict]:
    """Table 1 / Figure 3A source: display classes and platform/global unions, both cohorts and populations."""
    rows = [row for row in read_tsv(metrics / "table1_source.tsv") if row["scope"] in {"class", "group", "global"}]
    result = [dict({name: row[name] for name in DETECTION_FIELDS if name in row},
                   block="Reference classes" if row["scope"] == "class" else "Platform unions") for row in rows]
    return sorted(result, key=lambda row: (int(row["order"]), row["population"], row["endpoint"], row["cohort"]))


def paired_breadth(metrics: Path, outcomes: list[dict]) -> list[dict]:
    """Strict all-record breadth of every parent in the native breadth cohort, per platform group."""
    reference, controls = {}, defaultdict(list)
    for row in outcomes:
        if row["endpoint"] == "strict" and row["population"] == "all_records":
            if row["role"] == "real":
                reference[row["parent_id"]] = row
            else:
                controls[row["parent_id"]].append(row)
    rows = []
    for row in read_tsv(metrics / "cohort_membership.tsv"):
        if (row["endpoint"], row["population"], row["measure"], row["state"]) != ("strict", "all_records", "breadth", "eligible"):
            continue
        parent = row["parent_id"]
        draws = controls[parent]
        for group in row["groups"].split(";"):
            rows.append(dict(group=group, parent_id=parent, length=reference[parent]["length"], reference=reference[parent]["retained_accessions"],
                             control_mean=statistics.mean(draw["retained_accessions"] for draw in draws),
                             reference_positive=reference[parent]["positive"], control_positive=sum(draw["positive"] == 1 for draw in draws),
                             reference_cap=reference[parent]["retrieval"] == "target_limit_reached",
                             any_control_cap=any(draw["retrieval"] == "target_limit_reached" for draw in draws)))
    return sorted(rows, key=lambda row: (group_key(row["group"]), row["parent_id"]))


def position_tables(analysis: Path, reports: Path) -> tuple[list, list, dict]:
    events, policy = primary_positional_events(analysis)
    points = {}
    for row in events:
        points[(row["group"], *orientation_key(row))] = dict(group=row["group"], condition=row["condition"], parent_id=row["parent_id"],
                                                             accession_version=row["accession_version"], length=row["length"],
                                                             record_length=row["slen"], distance=row["distance_to_end"])
    points = sorted(points.values(), key=lambda row: (row["group"], row["accession_version"], row["distance"], row["parent_id"]))
    sensitivity = []
    for row in read_tsv(reports / "positional" / "positional_sensitivity.tsv"):
        sensitivity.append(dict(variant=row["variant"], group=row["group"], window=int(row["window"]), minimum_length=int(row["minimum_length"]),
                                events=int(row["events"]), terminal=int(row["terminal_events"]), expected=float(row["expected_n2"]),
                                ratio=None if row["fold_n2"] == "not_estimated" else float(row["fold_n2"]), unit_rule=row["unit_rule"],
                                units_without_project=int(row["units_without_project"])))
    for row in sensitivity:
        if row["variant"] != "all_events" or row["group"] == GLOBAL:
            continue
        subset = [point for point in points if point["group"] == row["group"] and point["record_length"] >= row["minimum_length"]]
        if len(subset) != row["events"] or sum(point["distance"] <= row["window"] for point in subset) != row["terminal"]:
            raise ContractError("Positional events disagree with the positional sensitivity table")
    return points, sensitivity, policy


def control_quality(design: Path, controls: Path, logical: list[dict], testable: set[str]) -> tuple[list, list, list, list]:
    """Control draw checks for every parent. Parent rows, summary counts and query properties flag the strict-testable parents."""
    report = verified_files(controls, read_json(controls / "control_qc.json"))
    originals = read_fasta(design / "physical_queries.fasta")
    parents = {canonical_id(sequence): canonical(sequence) for sequence in originals.values()}
    generated = read_fasta(controls / "controls.fasta")
    draws_table = {row["draw_id"]: row for row in read_tsv(controls / "draws.tsv")}
    memberships = defaultdict(set)
    for row in logical:
        if row["role"] == "real":
            memberships[row["parent_id"]].add(row["group"])
    by_parent, draw_rows = defaultdict(list), []
    for row in logical:
        if row["role"] != "control":
            continue
        identifier, source = row["logical_id"], draws_table[row["logical_id"]]
        parent, sequence = parents[row["parent_id"]], generated[identifier]
        check_draw(parent, sequence)
        if (canonical_id(sequence) != row["sequence_id"] or int(source["seed"]) != row["seed"] or source["parent_id"] != row["parent_id"]
                or int(source["draw_index"]) != row["draw_index"] or draw_seed(report["master_seed"], row["parent_id"], row["draw_index"]) != row["seed"]):
            raise ContractError("Selected control draw sequence or seed mapping changed")
        by_parent[row["parent_id"]].append(sequence)
        draw_rows.append(dict(parent_id=row["parent_id"], draw_index=row["draw_index"], logical_id=identifier, sequence_id=row["sequence_id"],
                              sequence=sequence, seed=row["seed"], invariants_passed=True, parent_identical_exact=sequence == parent,
                              parent_identical_canonical=canonical(sequence) == parent))
    parent_rows = []
    for parent_id, sequences in sorted(by_parent.items()):
        for group in ordered_groups(memberships[parent_id]):
            parent_rows.append(dict(group=group, parent_id=parent_id, strict_testable=parent_id in testable, length=len(parents[parent_id]),
                                    draws=len(sequences),
                                    distinct_exact=len(set(sequences)), distinct_canonical=len({canonical(sequence) for sequence in sequences}),
                                    repeated_exact=len(sequences) - len(set(sequences)),
                                    repeated_canonical=len(sequences) - len({canonical(sequence) for sequence in sequences}),
                                    parent_identical_exact=sum(sequence == parents[parent_id] for sequence in sequences),
                                    parent_identical_canonical=sum(canonical(sequence) == parents[parent_id] for sequence in sequences)))
    repeated = {parent: len(values) - len({canonical(value) for value in values}) for parent, values in by_parent.items()}
    summary = dict(parents=len(by_parent), selected_draws=len(draw_rows), invariants_failed=0, unique_seeds=len({row["seed"] for row in draw_rows}),
                   repeated_canonical=sum(repeated.values()),
                   repeated_exact=sum(len(values) - len(set(values)) for values in by_parent.values()),
                   parent_identical_canonical=sum(row["parent_identical_canonical"] for row in draw_rows),
                   parent_identical_exact=sum(row["parent_identical_exact"] for row in draw_rows),
                   testable_parents=sum(parent in testable for parent in by_parent),
                   testable_draws=sum(row["parent_id"] in testable for row in draw_rows),
                   testable_repeated_canonical=sum(count for parent, count in repeated.items() if parent in testable),
                   testable_parent_identical_canonical=sum(row["parent_identical_canonical"] for row in draw_rows if row["parent_id"] in testable),
                   source_qc_sha256=sha256(controls / "control_qc.json"))
    properties, seen = [], set()
    for row in logical:
        if row["role"] != "real" or (row["parent_id"], row["group"], row["condition"]) in seen:
            continue
        seen.add((row["parent_id"], row["group"], row["condition"]))
        sequence = parents[row["parent_id"]]
        properties.append(dict(group=row["group"], condition=row["condition"], parent_id=row["parent_id"],
                               strict_testable=row["parent_id"] in testable, length=len(sequence),
                               gc_fraction=(sequence.count("G") + sequence.count("C")) / len(sequence),
                               dominant_base=max(Counter(sequence).values()) / len(sequence), distinct_dinucleotides=len(dinucleotides(sequence))))
    return draw_rows, parent_rows, [summary], properties


def cap_tables(statuses: dict, logical: list[dict], outcomes: list[dict]) -> tuple[list, list]:
    by_sequence, parent_groups = defaultdict(list), defaultdict(set)
    for row in logical:
        by_sequence[row["sequence_id"]].append(row)
        if row["role"] == "real":
            parent_groups[row["parent_id"]].add(row["group"])
    lookup = {(row["sequence_id"], row["endpoint"]): row for row in outcomes if row["population"] == "all_records"}
    rows, membership_rows = [], []
    for sequence, status in sorted(statuses.items()):
        if status["target_limit_reached"] != "yes":
            continue
        members = sorted(by_sequence[sequence], key=lambda row: row["logical_id"])
        strict, tolerant = lookup[sequence, "strict"], lookup[sequence, "tolerant"]
        rows.append(dict(sequence_id=sequence, length=members[0]["length"], roles=";".join(sorted({row["role"] for row in members})),
                         groups=";".join(ordered_groups(group for row in members for group in parent_groups[row["parent_id"]])),
                         aliases=";".join(row["logical_id"] for row in members if row["role"] == "real"), logical_memberships=len(members),
                         hsp_rows=status["hsp_rows"], targets=status["unique_targets"], strict_retained=strict["retained_accessions"],
                         tolerant_retained=tolerant["retained_accessions"], strict_resolution=strict["resolution"],
                         tolerant_resolution=tolerant["resolution"]))
        membership_rows.extend(dict(sequence_id=sequence, **{key: row[key] for key in ("logical_id", "role", "parent_id", "draw_index", "seed")})
                               for row in members)
    return rows, membership_rows


def alignment_quality(search: Path, logical: list[dict], statuses: dict, testable: set[str]) -> tuple[list, list, list]:
    """Identity and query-span distributions of every returned HSP of the strict-testable reference sequences, rejected HSPs included."""
    groups = defaultdict(set)
    for row in logical:
        if row["role"] == "real" and row["parent_id"] in testable:
            groups[row["sequence_id"]].add(row["group"])
    marginal, joint, classes, by_query = Counter(), Counter(), Counter(), Counter()
    for _, row in search_hsps(search, set(groups)):
        by_query[row["qseqid"]] += 1
        span = 100 * (abs(row["qend"] - row["qstart"]) + 1) / row["qlen"]
        state = "Strict" if strict_hsp(row) else "Tolerant only" if tolerant_hsp(row) else "Rejected"
        for group in groups[row["qseqid"]]:
            marginal[group, "Identity", row["pident"]] += 1
            marginal[group, "Query span", span] += 1
            joint[group, min(100, math.floor(row["pident"])), min(100, math.floor(span)), row["gapopen"] > 0, state] += 1
            classes[group, state] += 1
    for sequence in groups:
        if sequence in statuses and by_query[sequence] != statuses[sequence]["hsp_rows"]:
            raise ContractError("Raw quality rows do not match the completed search status")
    present = ordered_groups(group for values in groups.values() for group in values)
    marginal_rows = []
    for group in present:
        for measure in ("Identity", "Query span"):
            counts = sorted((value, count) for (name, label, value), count in marginal.items() if name == group and label == measure)
            total, cumulative = sum(count for _, count in counts), 0
            for value, count in counts:
                cumulative += count
                marginal_rows.append(dict(group=group, measure=measure, value=value, hsp_rows=count, total_hsps=total, cumulative=cumulative / total))
    joint_rows = [dict(group=group, identity_bin=identity, span_bin=span, gapped=gapped, acceptance=state, hsp_rows=count)
                  for (group, identity, span, gapped, state), count in sorted(joint.items(), key=lambda item: (group_key(item[0][0]), item[0]))]
    class_rows = [dict(group=group, acceptance=state, hsp_rows=classes[group, state]) for group in present
                  for state in ("Strict", "Tolerant only", "Rejected")]
    return marginal_rows, joint_rows, class_rows


def fixed_prefix_precision(outcomes: list[dict], prefixes: list[int], target_limit: int, draws: int) -> list[dict]:
    """The final matched cohort followed back through the declared control prefixes (diagnostic, no coverage claim)."""
    rows = []
    groups = ordered_groups(row["group"] for row in outcomes if row["role"] == "real")
    alpha = 0.05 / len(prefixes)
    for population in sorted({row["population"] for row in outcomes}):
        for group in groups + [GLOBAL]:
            for endpoint in ("strict", "tolerant"):
                for measure, field, bound in (("detection", "positive", 1), ("breadth", "retained_accessions", target_limit)):
                    selected_rows = [row for row in outcomes if row["endpoint"] == endpoint and row["population"] == population]
                    original, controls, eligible, _ = matched_population(selected_rows, group, "all_conditions", endpoint, population, draws,
                                                                         measure=measure)
                    for prefix in prefixes:
                        selected = [row for parent in sorted(eligible) for row in controls[parent] if row["draw_index"] <= prefix]
                        estimate = statistics.mean(row[field] for row in selected) if selected else None
                        lower, upper, _ = (mean_interval([row[field] for row in selected], [row["seed"] for row in selected], alpha, bound)
                                           if selected else (None, None, None))
                        rows.append(dict(group=group, endpoint=endpoint, population=population, measure=measure, draws_per_parent=prefix,
                                         planned=len(original), matched=len(eligible), estimate=estimate, lower=lower, upper=upper,
                                         status="fixed_final_resolved_subset_diagnostic_no_new_coverage_claim"))
    return rows


def conceptual_examples() -> list[dict]:
    reference = "ACGTACGTACGTACGTACGTACGT"
    definitions = (("Exact full query", reference, 24, 0, 0, "Strict", True, True),
                   ("One substitution", reference[:11] + "A" + reference[12:], 24, 1, 0, "Tolerant", False, True),
                   ("22 of 24 bases", reference[:22] + "..", 22, 0, 0, "Tolerant", False, True),
                   ("21 of 24 bases", reference[:21] + "...", 21, 0, 0, "Rejected", False, False),
                   ("One gapped HSP", reference[:12] + "-" + reference[13:], 24, 0, 1, "Rejected", False, False),
                   ("Two partial HSPs", reference[:12] + " / " + reference[12:], 12, 0, 0, "Rejected", False, False))
    rows = []
    for label, subject, span, mismatches, gaps, decision, expected_strict, expected_tolerant in definitions:
        row = dict(qlen=24, gapopen=gaps, qstart=1, qend=span, sstart=1, send=span - gaps, slen=24, length=span, mismatch=mismatches,
                   pident=100 * (span - mismatches - gaps) / span)
        if strict_hsp(row) != expected_strict or tolerant_hsp(row) != expected_tolerant:
            raise ContractError("Conceptual illustration does not satisfy the implemented predicates")
        rows.append(dict(example=label, reference=reference, subject=subject, decision=decision, strict=expected_strict,
                         tolerant=expected_tolerant, conceptual=True))
    return rows


def benchmark_tables(benchmarks: dict[str, Path]) -> tuple[list, list, dict]:
    metrics, spikes, reports = [], [], {}
    for database in sorted(benchmarks):
        directory = Path(benchmarks[database])
        report = read_json(directory / "benchmark.json")
        if report.get("status") != "complete" or report.get("predicate_disagreements") != 0:
            raise ContractError(f"Benchmark {database} is not complete without predicate disagreements")
        reports[database] = dict(sha256=sha256(directory / "benchmark.json"), contract=report.get("contract"),
                                 pair_universe=report.get("independent_pair_universe"))
        metrics += [dict(database=database, **row) for row in read_tsv(directory / "metrics.tsv")]
        spikes += [dict(database=database, **row) for row in read_tsv(directory / "short_query_detectability.tsv")]
    return metrics, spikes, reports


def copy_table(source: Path | None, destination: Path, fields_if_missing=None) -> None:
    if source is not None and source.is_file():
        shutil.copyfile(source, destination)
    elif fields_if_missing is not None:
        write(destination, [], fields_if_missing)
    else:
        raise ContractError(f"Missing figure-source input: {source}")


def build_figure_sources(destination: Path, *, metrics: Path, reports: Path, analysis: Path, design: Path, controls: Path, search: Path,
                         inference: Path, benchmarks: dict, record_cases: Path | None = None) -> dict:
    metrics, reports, analysis, design, controls, search, inference = map(Path, (metrics, reports, analysis, design, controls, search, inference))
    report = load_analysis(analysis)
    analysis_sha = sha256(analysis / "analysis.json")
    metrics_manifest, reports_manifest = verify_manifest(metrics, METRICS_CONTRACT), verify_manifest(reports, REPORTS_CONTRACT)
    if metrics_manifest.get("analysis_sha256") != analysis_sha or reports_manifest.get("analysis_sha256") != analysis_sha:
        raise ContractError("Metrics or reports belong to a different analysis")
    precision_record = verified_files(inference, read_json(inference / "precision.json"))
    if precision_record.get("source_analysis_sha256") not in (None, analysis_sha):
        raise ContractError("The inference output belongs to a different analysis")
    draws, target_limit = report["draws_per_parent"], report["search_settings"]["max_target_seqs"]
    outcomes = read_table(analysis / "outcomes.tsv", OUTCOME_SCHEMA)
    logical = read_logical(analysis / "logical_queries.tsv")
    statuses = read_statuses(analysis / "search_status.tsv")
    destination = new_directory(destination)
    precision = read_tsv(inference / "precision.tsv")
    prefixes = sorted({int(row["draws_per_parent"]) for row in precision})
    strict_testable = rule_testable_parents(outcomes, "strict")
    draw_rows, parent_rows, qc_summary, properties = control_quality(design, controls, logical, strict_testable)
    caps, cap_memberships = cap_tables(statuses, logical, outcomes)
    quality, joint, accepted = alignment_quality(search, logical, statuses, strict_testable)
    benchmark_metrics, spikes, benchmark_reports = benchmark_tables(benchmarks)
    points, sensitivity, policy = position_tables(analysis, reports)
    flow_scopes = read_tsv(reports / "eligibility" / "eligibility_flow.tsv")
    flow = [dict(state=row["state"], reason=row["reason"], accessions=int(row["accessions"])) for row in flow_scopes
            if row["scope"] == "testable_strict" and row["state"] != "all"]
    totals = {row["reason"]: int(row["accessions"]) for row in flow_scopes if row["scope"] == "testable_strict" and row["state"] == "all"}
    if sum(row["accessions"] for row in flow if row["state"] != "unresolved" or row["reason"] != "unresolved_identifier") != totals["resolved_identifiers"]:
        raise ContractError("Eligibility flow does not partition the resolved strict accessions")
    tables = {"detection": (detection_rows(metrics), DETECTION_FIELDS),
              "paired_breadth": (paired_breadth(metrics, outcomes), ("group", "parent_id", "length", "reference", "control_mean",
                                                                      "reference_positive", "control_positive", "reference_cap", "any_control_cap")),
              "eligibility_flow": (flow, ("state", "reason", "accessions")),
              "eligibility_flow_scopes": (flow_scopes, ("scope", "state", "reason", "accessions")),
              "position_events": (points, ("group", "condition", "parent_id", "accession_version", "length", "record_length", "distance")),
              "position_sensitivity": (sensitivity, ("variant", "group", "window", "minimum_length", "events", "terminal", "expected", "ratio",
                                                     "unit_rule", "units_without_project")),
              "control_draw_qc": (draw_rows, None), "control_parent_qc": (parent_rows, None), "control_qc_summary": (qc_summary, None),
              "query_properties": (properties, None),
              "cap_sequences": (caps, ("sequence_id", "length", "roles", "groups", "aliases", "logical_memberships", "hsp_rows", "targets",
                                       "strict_retained", "tolerant_retained", "strict_resolution", "tolerant_resolution")),
              "cap_memberships": (cap_memberships, ("sequence_id", "logical_id", "role", "parent_id", "draw_index", "seed")),
              "quality_marginal": (quality, ("group", "measure", "value", "hsp_rows", "total_hsps", "cumulative")),
              "quality_joint": (joint, ("group", "identity_bin", "span_bin", "gapped", "acceptance", "hsp_rows")),
              "quality_acceptance": (accepted, ("group", "acceptance", "hsp_rows")),
              "precision": (precision, list(precision[0]) if precision else None),
              "precision_fixed_cohort": (fixed_prefix_precision(outcomes, prefixes, target_limit, draws), None),
              "matching_examples": (conceptual_examples(), None),
              "benchmark_metrics": (benchmark_metrics, None), "short_query_detectability": (spikes, None)}
    for name, (rows, fields) in tables.items():
        write(destination / (name + ".tsv"), rows, fields)
    copy_table(analysis / "length_testability.tsv", destination / "length_testability.tsv")
    copy_table(analysis / "all_records_density.tsv", destination / "density.tsv")
    copy_table(analysis / "logical_queries.tsv", destination / "logical_queries.tsv")
    copy_table(metrics / "testability.tsv", destination / "testability_classes.tsv")
    case_fields = ("case_id", "category", "accession_version", "query_id")
    copy_table(Path(record_cases) / "record_cases.tsv" if record_cases else None, destination / "record_cases.tsv", case_fields)
    copy_table(Path(record_cases) / "record_case_segments.tsv" if record_cases else None, destination / "record_case_segments.tsv",
               ("case_id", "segment_id", "role"))
    copy_table(reports / "record_cases" / "record_case_table.tsv", destination / "record_case_verification.tsv", ("case_id", "category"))
    copy_table(reports / "target_limit" / "target_limit_summary.tsv", destination / "target_limit_summary.tsv", ("group", "endpoint", "population"))
    copy_table(reports / "target_limit" / "target_limit_parents.tsv", destination / "target_limit_parents.tsv", ("parent_id", "endpoint", "population"))
    real = [row for row in logical if row["role"] == "real"]
    counts = dict(physical_references=len(real), canonical_parents=len({row["parent_id"] for row in real}),
                  conditions=len({row["condition"] for row in real}), controls=sum(row["role"] == "control" for row in logical),
                  strict_testable_parents=len(strict_testable), tolerant_testable_parents=len(rule_testable_parents(outcomes, "tolerant")),
                  strict_testable_reference_hsps=sum(statuses[sequence]["hsp_rows"] for sequence in
                                                     {row["sequence_id"] for row in real if row["parent_id"] in strict_testable}
                                                     if sequence in statuses),
                  draws_per_parent=draws, distinct_searches=len(statuses), testable_strict_accessions=totals["resolved_identifiers"],
                  positional_events=len(points), cap_sequences=len(caps))
    files = {path.name: sha256(path) for path in sorted(destination.iterdir()) if path.is_file()}
    manifest = dict(contract=CONTRACT, analysis_sha256=analysis_sha, metrics_manifest_sha256=sha256(metrics / "manifest.json"),
                    reports_manifest_sha256=sha256(reports / "manifest.json"), precision_sha256=sha256(inference / "precision.json"),
                    benchmarks=benchmark_reports, counts=counts, draw_prefixes=prefixes, target_limit=target_limit,
                    groups=ordered_groups(row["group"] for row in real) or list(GROUPS),
                    position_policy={key: policy["position"][key] for key in ("windows", "minimum_lengths", "primary_window", "primary_minimum_length")},
                    rounding="Full precision in all tables; rounding only at display", files=files)
    write_json(destination / "manifest.json", manifest)
    return dict(contract=CONTRACT, tables=len(files), counts=counts)
