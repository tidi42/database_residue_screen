from pathlib import Path

from .contracts import ContractError, read_json, validate_record


# Superseded and alternative conditions keep their own condition-level results
# but never enter group or global unions.
CONDITION_ROLE_VALUES = ("primary", "superseded", "alternative")
DEFAULT_CONDITION_ROLE = dict(panel_version="v1", role="primary")

POLICY_SCHEMA = {
    "type": "object", "additionalProperties": False,
    "required": ["version", "profile", "word_size", "minimum_query_length", "record_populations", "unknown_policy", "inference", "position"],
    "properties": {
        "version": {"const": "1.0.0"}, "profile": {"const": "reference_screen"}, "word_size": {"const": 7},
        "minimum_query_length": {"type": "object", "additionalProperties": False, "required": ["strict", "tolerant"],
                                 "properties": {endpoint: {"type": "integer", "minimum": 7} for endpoint in ("strict", "tolerant")}},
        "record_populations": {"type": "array", "minItems": 1, "uniqueItems": True,
                               "items": {"enum": ["all_records", "linear_genomic"]}, "contains": {"const": "all_records"}},
        "unknown_policy": {"const": "retain_unknown_exclude_unresolved_matched_parent"},
        "inference": {"const": "descriptive_retained_lower_bounds_not_biological_prevalence"},
        "sensitivity_baseline_evalue": {"type": "number", "exclusiveMinimum": 0},
        "detectability": {"type": "array", "items": {"type": "object"}},
        "degenerate_parent_ids": {"type": "array", "items": {"type": "string"}, "uniqueItems": True},
        "protocol_sha256": {"type": "string", "pattern": "^[0-9a-f]{64}$"},
        "condition_roles": {"type": "object", "minProperties": 1, "additionalProperties": {
            "type": "object", "additionalProperties": False, "required": ["panel_version", "role"],
            "properties": {"panel_version": {"type": "string", "pattern": "^v[0-9]+$"},
                           "role": {"enum": list(CONDITION_ROLE_VALUES)},
                           "superseded_by": {"type": "array", "items": {"type": "string"}, "uniqueItems": True}}}},
        "position": {"type": "object", "additionalProperties": False,
                     "required": ["primary_window", "primary_minimum_length", "windows", "minimum_lengths"],
                     "properties": {"primary_window": {"type": "integer", "minimum": 0},
                                    "primary_minimum_length": {"type": "integer", "minimum": 1},
                                    "windows": {"type": "array", "minItems": 1, "uniqueItems": True, "items": {"type": "integer", "minimum": 0}},
                                    "minimum_lengths": {"type": "array", "minItems": 1, "uniqueItems": True, "items": {"type": "integer", "minimum": 1}}}},
    },
}


def default_policy() -> dict:
    return dict(version="1.0.0", profile="reference_screen", word_size=7,
                minimum_query_length=dict(strict=7, tolerant=7), record_populations=["all_records", "linear_genomic"],
                unknown_policy="retain_unknown_exclude_unresolved_matched_parent",
                inference="descriptive_retained_lower_bounds_not_biological_prevalence",
                position=dict(primary_window=500, primary_minimum_length=10000,
                              windows=[100, 500, 1000], minimum_lengths=[1000, 10000, 50000]))


def validate_policy(policy: dict) -> dict:
    validate_record(policy, POLICY_SCHEMA)
    rows = policy.get("detectability", [])
    if len({row.get("query_length") for row in rows}) != len(rows):
        raise ContractError("Duplicate calibrated query lengths")
    for row in rows:
        if not {"query_length", "word_size", "strict_testable", "tolerant_fully_testable", "exact_evalue", "tolerant_evalue", "evalue", "exact_witness_recovered", "strict_testability", "tolerant_testability", "strict_expected_chance_bound", "tolerant_expected_chance_bound", "informativeness_limit", "testable_chance_limit"} <= row.keys():
            raise ContractError("Incomplete detectability calibration")
        from .length_aware import verify_calibration_row
        verify_calibration_row(row)
    position = policy["position"]
    if position["primary_window"] not in position["windows"] or position["primary_minimum_length"] not in position["minimum_lengths"]:
        raise ContractError("Primary positional settings must appear in the declared sensitivities")
    roles = policy.get("condition_roles", {})
    for condition, entry in roles.items():
        targets = entry.get("superseded_by", [])
        if (entry["role"] == "superseded") != bool(targets) or any(roles.get(target, {}).get("role") != "primary" for target in targets):
            raise ContractError("Superseded conditions must name primary replacement conditions, and only they may: " + condition)
    return policy


def condition_role(policy: dict, condition: str) -> dict:
    """Return the declared bundle role of a condition; undeclared policies keep the historical all-primary behaviour."""
    roles = policy.get("condition_roles")
    if roles is None:
        return dict(DEFAULT_CONDITION_ROLE)
    if condition not in roles:
        raise ContractError("Condition has no declared bundle role: " + condition)
    return roles[condition]


def load_policy(path: Path | None = None) -> dict:
    return validate_policy(read_json(path) if path is not None else default_policy())


def applicability(length: int, endpoint: str, policy: dict) -> str:
    calibrated = next((row for row in policy.get("detectability", []) if row["query_length"] == length), None)
    if "detectability" in policy and calibrated is None:
        raise ContractError("Calibration does not cover this query length")
    if calibrated is not None:
        state = calibrated.get(endpoint + "_testability")
        if state in {"censored_untestable", "chance_limited"}:
            return state
        if not calibrated["strict_testable"]:
            return "censored_untestable"
        if endpoint == "tolerant" and not calibrated["tolerant_fully_testable"]:
            return "partially_score_censored"
    if length < (calibrated["word_size"] if calibrated else policy["word_size"]):
        return "below_blast_word_size"
    if length < policy["minimum_query_length"][endpoint]:
        return "excluded_by_declared_length_policy"
    return "operational_retained_search"


def primary_event(row: dict, policy: dict) -> bool:
    """True when a projected real event belongs to a primary condition (all conditions of one parent share a role)."""
    if "condition_roles" not in policy:
        return True
    return any(condition_role(policy, condition)["role"] == "primary" for condition in row["condition"].split(";") if condition)
