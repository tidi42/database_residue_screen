import statistics
from collections import Counter, defaultdict

from .contracts import ContractError, table_schema
from .eligibility import event_eligibility
from .policy import applicability, condition_role, validate_policy


OUTCOME_SCHEMA = table_schema(("logical_id", "role", "parent_id", "sequence_id", "group", "condition", "endpoint", "population",
                               "applicability", "completion", "resolution", "detection_resolution", "breadth_resolution", "retrieval", "null_state",
                               "panel_version", "condition_role"),
                              ("draw_index", "length", "seed", "retained_accessions", "positive", "unresolved_subjects", "excluded_subjects"))
MEASURES = ("real_detection_rate", "parent_weighted_control_rate", "paired_detection_difference", "real_mean_breadth",
            "parent_weighted_control_breadth", "paired_breadth_difference", "mc_detection_se", "mc_breadth_se", "zero_control_probability_upper95")
CONTRAST_SCHEMA = table_schema(("group", "condition", "condition_role", "panel_version", "endpoint", "population", "measure", "mc_assumptions") + MEASURES,
                               ("planned_parents", "eligible_parents", "excluded_short_parents", "unresolved_parents", "incomplete_parents",
                                "censored_parents", "partially_censored_parents", "chance_limited_parents", "degenerate_null_parents", "real_positive", "control_draws", "draws_per_parent", "target_limited_parents"))
MEMBERSHIP_SCHEMA = table_schema(("sequence_id", "role", "parent_id", "group", "condition"))
ELIGIBILITY_SCHEMA = table_schema(("sequence_id", "accession_version", "population", "state", "reason"))
BUNDLE_COLUMNS = {"panel_version", "condition_role"}
PRIMARY_UNION = "primary_union"
OUTCOME_V1_SCHEMA = table_schema(tuple(name for name, field in OUTCOME_SCHEMA["properties"].items()
                                     if field["type"] == "string" and name not in {"detection_resolution", "breadth_resolution"} | BUNDLE_COLUMNS),
                               tuple(name for name, field in OUTCOME_SCHEMA["properties"].items() if field["type"] == "integer"))
CONTRAST_V1_SCHEMA = table_schema(tuple(name for name, field in CONTRAST_SCHEMA["properties"].items()
                                      if field["type"] == "string" and name not in {"measure"} | BUNDLE_COLUMNS),
                                tuple(name for name, field in CONTRAST_SCHEMA["properties"].items() if field["type"] == "integer"))


def historical_v1_results(outcomes: list[dict], draws_per_parent: int, seed_reuse: int = 0) -> tuple[list, list]:
    legacy = [{name: row[name] for name in OUTCOME_V1_SCHEMA["properties"]} for row in outcomes]
    breadth_gated = [dict(row, detection_resolution=row["breadth_resolution"]) for row in outcomes]
    separated = matched_contrasts(breadth_gated, draws_per_parent, seed_reuse)
    merged = []
    for detection, breadth in zip(separated[::2], separated[1::2]):
        row = {name: detection[name] for name in CONTRAST_V1_SCHEMA["properties"]}
        for name in ("real_mean_breadth", "parent_weighted_control_breadth", "paired_breadth_difference", "mc_breadth_se"):
            row[name] = breadth[name]
        merged.append(row)
    return legacy, merged


def logical_memberships(logical: list[dict], draws_per_parent: int) -> list[dict]:
    if draws_per_parent < 1 or len({row["logical_id"] for row in logical}) != len(logical):
        raise ContractError("Invalid draw count or duplicate logical query identity")
    parents, lengths, draws = defaultdict(set), {}, defaultdict(set)
    for row in logical:
        if row["role"] != "real":
            continue
        parent = row["parent_id"]
        if row["sequence_id"] != parent or row["draw_index"] != 0:
            raise ContractError("Original query does not map to its canonical parent")
        if parent in lengths and lengths[parent] != row["length"]:
            raise ContractError("Physical aliases disagree on parent length")
        lengths[parent] = row["length"]
        parents[parent].add((row["group"], row["condition"]))
    if not parents:
        raise ContractError("No real canonical parents in logical panel")
    memberships = set()
    for row in logical:
        parent = row["parent_id"]
        if parent not in parents or row["length"] != lengths[parent]:
            raise ContractError("Logical query has an unknown parent or inconsistent length")
        if row["role"] == "control":
            if not 0 < row["seed"] <= 2147483647:
                raise ContractError("Control seed is outside the recorded generator range")
            if row["draw_index"] in draws[parent] or not 1 <= row["draw_index"] <= draws_per_parent:
                raise ContractError("Duplicated or out-of-range control draw")
            draws[parent].add(row["draw_index"])
            contexts = parents[parent]
        elif row["role"] == "real":
            contexts = {(row["group"], row["condition"])}
        else:
            raise ContractError("Unknown logical query role")
        for group, condition in contexts:
            memberships.add((row["sequence_id"], row["role"], parent, group, condition))
    if any(draws[parent] != set(range(1, draws_per_parent + 1)) for parent in parents):
        raise ContractError("Logical panel lacks the declared matched draws")
    return [dict(zip(MEMBERSHIP_SCHEMA["properties"], values)) for values in sorted(memberships)]


def project_real_events(events: list[dict], memberships: list[dict]) -> list[dict]:
    contexts = defaultdict(lambda: defaultdict(set))
    for row in memberships:
        if row["role"] == "real":
            contexts[row["sequence_id"]][(row["parent_id"], row["group"])].add(row["condition"])
    return [dict(event, parent_id=parent, group=group, condition=";".join(sorted(conditions)))
            for event in events for (parent, group), conditions in sorted(contexts[event["query_id"]].items())]


def build_outcomes(events: list[dict], logical: list[dict], statuses: list[dict], policy: dict,
                   draws_per_parent: int, seed_reuse: int = 0) -> tuple[list, list, list, list]:
    validate_policy(policy)
    memberships = logical_memberships(logical, draws_per_parent)
    by_sequence = {row["query_id"]: row for row in statuses}
    if len(by_sequence) != len(statuses) or set(by_sequence) != {row["sequence_id"] for row in logical}:
        raise ContractError("Search status membership is duplicated or incomplete; missing is not zero")
    states, reasons = {}, event_eligibility(events, policy["record_populations"])
    for row in events:
        if row["query_id"] not in by_sequence:
            raise ContractError("Event belongs to an undeclared search sequence")
        for population in policy["record_populations"]:
            key = row["accession_version"], population
            for endpoint in ("strict", "tolerant"):
                if row[endpoint] == "yes":
                    states[row["query_id"], endpoint, population, row["accession_version"]] = key
    counts = defaultdict(lambda: defaultdict(set))
    for (sequence, endpoint, population, accession), key in states.items():
        counts[sequence, endpoint, population][reasons[key][0]].add(accession)
    parent_roles = {}
    for row in logical:
        if row["role"] == "real":
            declared = condition_role(policy, row["condition"])
            if parent_roles.setdefault(row["parent_id"], declared) != declared:
                raise ContractError("A canonical parent belongs to conditions with different bundle roles")
    outcomes = []
    for row in logical:
        status = by_sequence[row["sequence_id"]]
        completed = status["status"] in {"success_hits", "success_zero"}
        declared = parent_roles[row["parent_id"]]
        for endpoint in ("strict", "tolerant"):
            for population in policy["record_populations"]:
                selected = counts[row["sequence_id"], endpoint, population]
                included, unresolved = len(selected["included"]), len(selected["unresolved"])
                resolution = ("not_completed" if not completed else "positive_with_unresolved_breadth" if included and unresolved else
                              "unresolved" if unresolved else "resolved_positive" if included else "resolved_zero")
                applicable = applicability(row["length"], endpoint, policy)
                if row["parent_id"] in policy.get("degenerate_parent_ids", []):
                    null_state = "proven_singleton_null"
                else:
                    null_state = "not_proven_singleton"
                censored = applicable in {"censored_untestable", "partially_score_censored", "chance_limited"}
                if censored:
                    resolution = applicable
                unavailable = applicable if censored else "not_completed" if not completed else "inapplicable" if applicable != "operational_retained_search" else ""
                outcomes.append(dict(**{name: row[name] for name in ("logical_id", "role", "parent_id", "sequence_id", "group", "condition", "draw_index", "length", "seed")},
                                     endpoint=endpoint, population=population, applicability=applicable,
                                     completion="successful" if completed else status["status"], resolution=resolution,
                                     detection_resolution=unavailable or ("unresolved" if unresolved and not included else "resolved"),
                                     breadth_resolution=unavailable or ("unresolved" if unresolved else "resolved"),
                                     null_state=null_state, retrieval="target_limit_reached" if status["target_limit_reached"] == "yes" else "retained_not_exhaustive",
                                     panel_version=declared["panel_version"], condition_role=declared["role"],
                                     retained_accessions=included, positive=-1 if censored else 1 if included else -1 if unresolved or not completed or censored else 0,
                                     unresolved_subjects=unresolved, excluded_subjects=len(selected["excluded"])))
    contrasts = matched_contrasts(outcomes, draws_per_parent, seed_reuse)
    membership_keys = {(row["query_id"], row["accession_version"], population) for row in events for population in policy["record_populations"]}
    eligibility = [dict(sequence_id=sequence, accession_version=accession, population=population,
                        state=reasons[accession, population][0], reason=reasons[accession, population][1])
                   for sequence, accession, population in sorted(membership_keys)]
    return outcomes, contrasts, memberships, eligibility


def in_union(row: dict, roles: tuple | None, panel_versions: tuple | None) -> bool:
    # Rows written before bundle roles existed are historical primary v1 rows.
    return ((roles is None or row.get("condition_role", "primary") in roles) and
            (panel_versions is None or row.get("panel_version", "v1") in panel_versions))


def matched_population(outcomes: list[dict], group: str, condition: str, endpoint: str, population: str, draws_per_parent: int, *, measure: str,
                       roles: tuple | None = ("primary",), panel_versions: tuple | None = None):
    """Matched parents and draws for one context.

    Union contexts (condition 'all_conditions', including the Global union) contain only conditions whose bundle role is
    in ``roles`` and whose panel version is in ``panel_versions`` (None = all). Condition-specific contexts are never
    role-filtered, so superseded and alternative families keep their own flagged estimates.
    """
    if measure not in {"detection", "breadth"}:
        raise ContractError("Matched population requires detection or breadth measure")
    selected = [row for row in outcomes if row["endpoint"] == endpoint and row["population"] == population]
    if any(f"{measure}_resolution" not in row for row in selected):
        raise ContractError("Measure-specific outcome contract required")
    real = {row["parent_id"]: row for row in selected if row["role"] == "real"
            and (group == "Global union" or row["group"] == group)
            and (row["condition"] == condition if condition != "all_conditions" else in_union(row, roles, panel_versions))}
    controls = defaultdict(list)
    for row in selected:
        if row["role"] == "control" and row["parent_id"] in real:
            controls[row["parent_id"]].append(row)
    eligible, attrition = {}, Counter()
    for parent, original in real.items():
        draws = controls[parent]
        if len(draws) != draws_per_parent or {row["draw_index"] for row in draws} != set(range(1, draws_per_parent + 1)):
            raise ContractError("Matched outcomes do not contain the declared draw identities")
        family = [original] + draws
        if original.get("null_state") == "proven_singleton_null":
            attrition["degenerate_null"] += 1
        if any(row["applicability"] == "censored_untestable" for row in family):
            attrition["censored"] += 1
        elif any(row["applicability"] == "chance_limited" for row in family):
            attrition["chance_limited"] += 1
        elif any(row["applicability"] == "partially_score_censored" for row in family):
            attrition["partial"] += 1
        elif any(row["applicability"] != "operational_retained_search" for row in family):
            attrition["short"] += 1
        elif any(row["completion"] != "successful" for row in family):
            attrition["incomplete"] += 1
        elif any(row[f"{measure}_resolution"] != "resolved" for row in family):
            attrition["unresolved"] += 1
        else:
            eligible[parent] = original
        if any(row["retrieval"] == "target_limit_reached" for row in family):
            attrition["target_limited"] += 1
    return real, controls, eligible, attrition


def matched_contrasts(outcomes: list[dict], draws_per_parent: int, seed_reuse: int = 0) -> list[dict]:
    real_rows = [row for row in outcomes if row["role"] == "real"]
    declared = {row["condition"]: (row.get("condition_role", "primary"), row.get("panel_version", "v1")) for row in real_rows}
    groups = sorted({row["group"] for row in real_rows if in_union(row, ("primary",), None)})
    contexts = [(group, "all_conditions") for group in groups] + [("Global union", "all_conditions")]
    contexts += sorted({(row["group"], row["condition"]) for row in real_rows})
    union_versions = ";".join(sorted({version for role, version in declared.values() if role == "primary"}))
    strata = sorted({(row["endpoint"], row["population"]) for row in real_rows})
    contrasts = []
    for group, condition in contexts:
        role, version = (PRIMARY_UNION, union_versions) if condition == "all_conditions" else declared[condition]
        for endpoint, population in strata:
            for measure, field in (("detection", "positive"), ("breadth", "retained_accessions")):
                real, controls, eligible, attrition = matched_population(outcomes, group, condition, endpoint, population,
                                                                        draws_per_parent, measure=measure)
                size = len(eligible)
                measures = {name: "not_estimated" for name in MEASURES}
                assumptions = "Conditional independent parent/draw sampling; not biological sampling uncertainty"
                if seed_reuse:
                    assumptions = "Seed reuse recorded; aggregate independent-draw Monte Carlo errors not certified"
                if size:
                    original_mean = statistics.mean(row[field] for row in eligible.values())
                    control_mean = statistics.mean(statistics.mean(row[field] for row in controls[parent]) for parent in eligible)
                    names = (("real_detection_rate", "parent_weighted_control_rate", "paired_detection_difference") if measure == "detection"
                             else ("real_mean_breadth", "parent_weighted_control_breadth", "paired_breadth_difference"))
                    measures.update(zip(names, map(str, (original_mean, control_mean, original_mean - control_mean))))
                    if not seed_reuse:
                        if measure == "detection" and not any(row[field] for parent in eligible for row in controls[parent]):
                            # Independent, possibly heterogeneous Bernoulli draws:
                            # product(1-p_i) <= (1-mean(p_i))**N (AM-GM).
                            measures["zero_control_probability_upper95"] = str(1 - 0.05 ** (1 / (size * draws_per_parent)))
                        if draws_per_parent > 1:
                            measures["mc_" + measure + "_se"] = str(sum(statistics.variance(row[field] for row in controls[parent]) / draws_per_parent
                                                                         for parent in eligible) ** 0.5 / size)
                contrasts.append(dict(group=group, condition=condition, condition_role=role, panel_version=version,
                                      endpoint=endpoint, population=population, measure=measure,
                                      mc_assumptions=assumptions, planned_parents=len(real), eligible_parents=size,
                                      excluded_short_parents=attrition["short"], unresolved_parents=attrition["unresolved"],
                                      incomplete_parents=attrition["incomplete"], censored_parents=attrition["censored"],
                                      partially_censored_parents=attrition["partial"], chance_limited_parents=attrition["chance_limited"], degenerate_null_parents=attrition["degenerate_null"], real_positive=sum(row["positive"] for row in eligible.values()),
                                      control_draws=size * draws_per_parent, draws_per_parent=draws_per_parent,
                                      target_limited_parents=attrition["target_limited"], **measures))
    return contrasts