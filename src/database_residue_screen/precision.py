import math
import statistics
from collections import Counter
from pathlib import Path

from .contracts import ContractError, read_json, read_table, table_schema, validate_record, write_table
from .design import digest
from .example import save_json
from .outcomes import OUTCOME_SCHEMA, matched_population


PILOT_SCHEMA = {
    "type": "object", "additionalProperties": False,
    "required": ["version", "panel_sha256", "draw_ladder", "maximum_draws", "retention_ladder", "confidence", "detection_half_width",
                 "breadth_absolute_half_width", "breadth_relative_half_width", "stopping", "threads", "batch_size", "timeout_seconds"],
    "properties": {"version": {"const": "1.0.0"}, "panel_sha256": {"type": "string", "pattern": "^([0-9a-f]{64})?$"},
                   "draw_ladder": {"type": "array", "minItems": 1, "uniqueItems": True, "items": {"type": "integer", "minimum": 1}},
                   "maximum_draws": {"type": "integer", "minimum": 1},
                   "retention_ladder": {"type": "array", "minItems": 2, "uniqueItems": True, "items": {"type": "integer", "minimum": 1}},
                   "confidence": {"type": "number", "exclusiveMinimum": 0, "exclusiveMaximum": 1},
                   "detection_half_width": {"type": "number", "exclusiveMinimum": 0, "maximum": 1},
                   "breadth_absolute_half_width": {"type": "number", "exclusiveMinimum": 0},
                   "breadth_relative_half_width": {"type": "number", "exclusiveMinimum": 0, "maximum": 1},
                   "stopping": {"const": "fixed_stages_no_optional_stopping"},
                   "threads": {"type": "integer", "minimum": 1}, "batch_size": {"type": "integer", "minimum": 1},
                   "timeout_seconds": {"type": "integer", "minimum": 1}},
}
PRECISION_SCHEMA = table_schema(("group", "condition", "endpoint", "population", "measure", "estimate", "half_width", "target_half_width",
                                 "lower", "upper", "precision_met", "status"),
                                ("draws_per_parent", "planned_parents", "matched_parents", "control_draws", "unique_seeds", "target_limited_parents", "censored_parents", "partially_censored_parents", "chance_limited_parents", "degenerate_null_parents", "required_independent_draws_for_target"))
CONDITIONAL_SCHEMA = table_schema(("group", "endpoint", "population", "measure", "status", "reference_mean", "control_mean", "difference",
                                   "interval_lower", "interval_upper", "tail_bound", "holm_bound", "observed_reference_total",
                                   "expected_reference_total", "interpretation"),
                                  ("planned_parents", "matched_parents", "control_draws", "unique_seeds", "draws_per_parent",
                                   "a_priori_excluded_parents", "holm_family_size"))
PARENT_PRECISION_SCHEMA = table_schema(("parent_id", "endpoint", "population", "measure", "mean_detection", "mean_breadth", "zero_success_upper", "status"),
                                       ("draws", "completed_draws", "unresolved_draws", "positive_draws", "distinct_sequences", "unique_seeds", "parent_identical_draws"))


def default_pilot_spec() -> dict:
    return dict(version="1.0.0", panel_sha256="", draw_ladder=[10, 50, 100], maximum_draws=100,
                retention_ladder=[500, 2000, 10000], confidence=0.95, detection_half_width=0.01,
                breadth_absolute_half_width=1.0, breadth_relative_half_width=0.05, stopping="fixed_stages_no_optional_stopping",
                threads=8, batch_size=16, timeout_seconds=3600)


def validate_pilot(specification: dict) -> dict:
    validate_record(specification, PILOT_SCHEMA)
    for name in ("draw_ladder", "retention_ladder"):
        if specification[name] != sorted(specification[name]):
            raise ContractError("Pilot ladders must be strictly increasing")
    if specification["draw_ladder"][-1] > specification["maximum_draws"]:
        raise ContractError("Pilot draw ladder exceeds its budget")
    return specification


def write_pilot_spec(panel: Path, destination: Path) -> dict:
    definition = read_json(panel / "panel.json")
    for name, checksum in definition["files"].items():
        if Path(name).name != name or digest(panel / name) != checksum:
            raise ContractError("Panel changed before pilot specification")
    specification = default_pilot_spec()
    specification["panel_sha256"] = digest(panel / "panel.json")
    validate_pilot(specification)
    destination.mkdir(parents=True, exist_ok=False)
    save_json(destination / "pilot.json", specification)
    report = dict(stage="prepared_not_executed", final_freeze=False, pilot_sha256=digest(destination / "pilot.json"),
                  panel_sha256=specification["panel_sha256"], review_before_production=True,
                  assumptions="Distinct seed streams treated as independent; shared seeds form bounded dependence groups. Not biological sampling uncertainty.",
                  precision="Two-sided Hoeffding bounds, simultaneous across the declared draw stages within each reported stratum; not across strata",
                  stopping=specification["stopping"], breadth_target="max(1 retained accession, 5% of the control mean); descriptive censored endpoint",
                  search_started=False, database_accepted=False)
    save_json(destination / "preparation.json", report)
    return report


def grouped_half_width(seeds: list[int], alpha: float, bound: float = 1.0) -> float:
    if not seeds or not 0 < alpha < 1 or bound <= 0 or any(seed <= 0 for seed in seeds):
        raise ContractError("Invalid conditional precision inputs")
    counts = Counter(seeds)
    squared_weights = sum((count / len(seeds)) ** 2 for count in counts.values())
    return bound * math.sqrt(math.log(2 / alpha) * squared_weights / 2)


def precision_rows(outcomes: list[dict], specification: dict, target_limit: int, available_draws: int) -> tuple[list, list]:
    validate_pilot(specification)
    if target_limit < 1:
        raise ContractError("Precision requires the declared retention bound")
    contexts = sorted({(row["group"], "all_conditions") for row in outcomes if row["role"] == "real"} |
                      {(row["group"], row["condition"]) for row in outcomes if row["role"] == "real"} | {("Global union", "all_conditions")})
    strata = sorted({(row["endpoint"], row["population"]) for row in outcomes})
    alpha = (1 - specification["confidence"]) / len(specification["draw_ladder"])
    rows, parent_rows = [], []
    for draws in specification["draw_ladder"]:
        subset = [row for row in outcomes if row["role"] == "real" or row["draw_index"] <= draws]
        for group, condition in contexts:
            for endpoint, population in strata:
                for measure, field, bound in (("detection", "positive", 1), ("breadth", "retained_accessions", target_limit)):
                    real, controls, eligible, attrition = (({}, {}, {}, Counter()) if draws > available_draws else
                                                          matched_population(subset, group, condition, endpoint, population, draws, measure=measure))
                    selected = [row for parent in eligible for row in controls[parent]]
                    seeds = [row["seed"] for row in selected]
                    status = "not_searched" if draws > available_draws else "no_resolved_matched_population" if not selected else "conditional_bound"
                    required = 0
                    values = dict(estimate="not_estimated", half_width="not_estimated", target_half_width="not_estimated", lower="not_estimated", upper="not_estimated", precision_met="not_estimated")
                    if selected:
                        estimate = statistics.mean(row[field] for row in selected)
                        if any(not 0 <= row[field] <= bound for row in selected):
                            raise ContractError("Observed retained outcome exceeds the declared precision bound")
                        radius = grouped_half_width(seeds, alpha, bound)
                        target = specification["detection_half_width"] if measure == "detection" else max(specification["breadth_absolute_half_width"], specification["breadth_relative_half_width"] * estimate)
                        required = math.ceil(math.log(2 / alpha) * bound**2 / (2 * target**2))
                        values.update(estimate=str(estimate), half_width=str(radius), target_half_width=str(target),
                                      lower=str(max(0, estimate - radius)), upper=str(min(bound, estimate + radius)), precision_met="yes" if radius <= target else "no")
                        if len(eligible) != len(real) - attrition["short"] - attrition["censored"] - attrition["partial"] - attrition["chance_limited"]:
                            status = "diagnostic_resolved_subset_no_coverage_claim"
                            values["precision_met"] = "not_estimated"
                    rows.append(dict(group=group, condition=condition, endpoint=endpoint, population=population, measure=measure,
                                     draws_per_parent=draws, planned_parents=len(real), matched_parents=len(eligible), control_draws=len(selected),
                                     unique_seeds=len(set(seeds)), censored_parents=attrition["censored"], partially_censored_parents=attrition["partial"], chance_limited_parents=attrition["chance_limited"],
                                     degenerate_null_parents=attrition["degenerate_null"], required_independent_draws_for_target=required,
                                     target_limited_parents=attrition["target_limited"], status=status, **values))
        if draws <= available_draws:
            for endpoint, population in strata:
                for measure in ("detection", "breadth"):
                    real, controls, eligible, _ = matched_population(subset, "Global union", "all_conditions", endpoint, population, draws, measure=measure)
                    for parent in sorted(real):
                        selected = controls[parent]
                        completed = [row for row in selected if row["completion"] == "successful"]
                        unresolved = sum(row[f"{measure}_resolution"] == "unresolved" for row in selected)
                        positive = sum(row["positive"] == 1 for row in selected)
                        unique_seeds = len({row["seed"] for row in selected})
                        resolved = parent in eligible
                        parent_rows.append(dict(parent_id=parent, endpoint=endpoint, population=population, measure=measure, draws=draws,
                                                completed_draws=len(completed), unresolved_draws=unresolved, positive_draws=positive,
                                                distinct_sequences=len({row["sequence_id"] for row in selected}), unique_seeds=unique_seeds,
                                                parent_identical_draws=sum(row["sequence_id"] == parent for row in selected),
                                                mean_detection=str(statistics.mean(row["positive"] for row in selected)) if resolved and measure == "detection" else "not_estimated",
                                                mean_breadth=str(statistics.mean(row["retained_accessions"] for row in selected)) if resolved and measure == "breadth" else "not_estimated",
                                                zero_success_upper=str(1 - alpha ** (1 / unique_seeds)) if resolved and measure == "detection" and not positive else "not_estimated",
                                                status="conditional_distinct_seed_streams" if resolved else "unresolved_or_inapplicable"))
    return rows, parent_rows


def format_log_probability(value: float) -> str:
    """Print exp(value) without underflow; extreme bounds never become zero."""
    if value >= math.log(1e-300):
        return str(math.exp(value))
    exponent = math.floor(value / math.log(10))
    return f"{10 ** (value / math.log(10) - exponent):.6f}e{exponent}"


def conditional_rows(outcomes: list[dict], specification: dict, target_limit: int,
                     endpoint: str = "strict", population: str = "all_records") -> list[dict]:
    """Post hoc reference-versus-generator contrasts with Holm adjustment (S8.6).

    D = reference mean - control generator mean. The stage-protected Hoeffding
    interval for the generator mean (precision radius at alpha/K) is subtracted
    from the fixed reference mean; p = min(1, 2K exp(-2 D^2 / (B^2 sum w^2))).
    Holm covers every group-by-measure slot, retaining non-estimable slots as
    one. Bounds concern Monte Carlo error of the generator expectation for fixed
    sequences, not exchangeability, contamination or generalisation; expected
    versus observed totals are reported alongside.
    """
    validate_pilot(specification)
    if target_limit < 1:
        raise ContractError("Conditional contrasts require the declared retention bound")
    stages, draws = len(specification["draw_ladder"]), specification["draw_ladder"][-1]
    alpha = (1 - specification["confidence"]) / stages
    subset = [row for row in outcomes if row["role"] == "real" or row["draw_index"] <= draws]
    groups = sorted({row["group"] for row in outcomes if row["role"] == "real"})
    rows, logs = [], []
    for group in groups:
        for measure, field, bound in (("detection", "positive", 1), ("breadth", "retained_accessions", target_limit)):
            real, controls, eligible, attrition = matched_population(subset, group, "all_conditions", endpoint, population, draws, measure=measure)
            selected = [row for parent in eligible for row in controls[parent]]
            excluded = attrition["short"] + attrition["censored"] + attrition["partial"] + attrition["chance_limited"]
            status = ("no_resolved_matched_population" if not selected else "conditional_bound" if len(eligible) == len(real) - excluded
                      else "diagnostic_resolved_subset_no_coverage_claim")
            values = {name: "not_estimated" for name in ("reference_mean", "control_mean", "difference", "interval_lower", "interval_upper",
                                                         "tail_bound", "observed_reference_total", "expected_reference_total")}
            log_tail = 0.0
            if selected:
                if any(not 0 <= row[field] <= bound for row in selected) or any(not 0 <= row[field] <= bound for row in eligible.values()):
                    raise ContractError("Observed retained outcome exceeds the declared conditional bound")
                reference = statistics.mean(row[field] for row in eligible.values())
                control = statistics.mean(statistics.mean(row[field] for row in controls[parent]) for parent in eligible)
                values.update(reference_mean=str(reference), control_mean=str(control), difference=str(reference - control),
                              observed_reference_total=str(sum(row[field] for row in eligible.values())),
                              expected_reference_total=str(sum(statistics.mean(row[field] for row in controls[parent]) for parent in eligible)))
                if status == "conditional_bound":
                    seeds = [row["seed"] for row in selected]
                    radius = grouped_half_width(seeds, alpha, bound)
                    squared = sum((count / len(seeds)) ** 2 for count in Counter(seeds).values())
                    log_tail = min(0.0, math.log(2 * stages) - 2 * (reference - control) ** 2 / (bound ** 2 * squared))
                    values.update(interval_lower=str(reference - min(bound, control + radius)),
                                  interval_upper=str(reference - max(0, control - radius)), tail_bound=format_log_probability(log_tail))
            logs.append(log_tail)
            rows.append(dict(group=group, endpoint=endpoint, population=population, measure=measure, status=status,
                             planned_parents=len(real), matched_parents=len(eligible), control_draws=len(selected),
                             unique_seeds=len({row["seed"] for row in selected}), draws_per_parent=draws, a_priori_excluded_parents=excluded,
                             interpretation="Post hoc exploratory; Monte Carlo error of the shuffle-generator mean for fixed sequences. "
                                            "Not exchangeability, contamination, equivalence or generalisation evidence.", **values))
    family = len(rows)
    order = sorted(range(family), key=lambda index: logs[index])
    running = -math.inf
    for rank, index in enumerate(order):
        running = max(running, min(0.0, math.log(family - rank) + logs[index]))
        rows[index]["holm_bound"] = format_log_probability(running) if rows[index]["status"] == "conditional_bound" else "not_estimated"
        rows[index]["holm_family_size"] = family
    return rows


def assess_precision(analysis: Path, specification_file: Path, destination: Path) -> dict:
    source = read_json(analysis / "analysis.json")
    specification = validate_pilot(read_json(specification_file))
    if source.get("result_contract") != "reference_screen_v2" or not specification["panel_sha256"] or specification["panel_sha256"] != source["panel_sha256"]:
        raise ContractError("Precision specification is not bound to this result panel")
    for name, checksum in source["files"].items():
        if Path(name).name != name or digest(analysis / name) != checksum:
            raise ContractError("Source result changed before precision assessment")
    if source.get("source_search") and (Path(source["source_search"]) / "quarantine.json").exists():
        raise ContractError("Source search is quarantined")
    outcomes = read_table(analysis / "outcomes.tsv", OUTCOME_SCHEMA)
    rows, parents = precision_rows(outcomes, specification, source["search_settings"]["max_target_seqs"], source["draws_per_parent"])
    destination.mkdir(parents=True, exist_ok=False)
    write_table(destination / "precision.tsv", rows, PRECISION_SCHEMA)
    write_table(destination / "parent_precision.tsv", parents, PARENT_PRECISION_SCHEMA)
    if specification["draw_ladder"][-1] <= source["draws_per_parent"]:
        write_table(destination / "conditional_contrasts.tsv",
                    conditional_rows(outcomes, specification, source["search_settings"]["max_target_seqs"]), CONDITIONAL_SCHEMA)
    final_stage = [row for row in rows if row["draws_per_parent"] == specification["draw_ladder"][-1]]
    report = dict(stage="conditional_null_precision_not_retrieval_acceptance", result_contract="measure_precision_v2", final_freeze=False,
                  source_analysis_sha256=digest(analysis / "analysis.json"), specification_sha256=digest(specification_file),
                  all_final_stage_targets_met=bool(final_stage) and all(row["precision_met"] == "yes" and row["status"] == "conditional_bound" for row in final_stage),
                  biological_sampling_uncertainty=False, retrieval_adequacy_established=False,
                  bound="Hoeffding with total weight grouped by seed; simultaneous across fixed stages per stratum, conditional on independence between seed streams",
                  stopping=specification["stopping"],
                  conditional_contrasts="Post hoc exploratory reference-generator contrasts with Holm over all group-by-measure slots (S8.6)",
                  files={path.name: digest(path) for path in sorted(destination.glob("*.tsv"))})
    save_json(destination / "precision.json", report)
    return report