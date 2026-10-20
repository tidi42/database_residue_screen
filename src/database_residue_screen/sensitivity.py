"""Explicit full inventory and score-calibrated length cohorts (no missing=zero).

Protocol-v2 rates use only a-priori testable, resolved matched parents. The
historical zero-filled fractions (protocol v1: untestable searches counted as
completed zeros) are recomputed alongside, explicitly labelled, so the archived
full-panel and calibrated-length fractions can be reproduced rather than retyped.
"""
import statistics
from collections import defaultdict

from .contracts import table_schema
from .length_aware import score_retained
from .outcomes import matched_population
from .policy import applicability

EXCLUDED_STATES = {"censored_untestable", "chance_limited", "partially_score_censored"}
SENSITIVITY_SCHEMA = table_schema(("cohort", "group", "endpoint", "population", "rate", "interpretation", "historical_zero_filled_rate"),
                                ("panel_parents", "cohort_parents", "eligible_matched_parents", "censored_parents",
                                 "partially_censored_parents", "chance_limited_parents", "observed_positive", "denominator",
                                 "historical_zero_filled_positive", "historical_zero_filled_parents"))
TESTABILITY_SCHEMA = table_schema(("endpoint", "testability", "testability_reason", "score_retainable", "exact_witness_recovered",
                                   "partially_score_censored"),
                                  ("query_length", "minimum_qualifying_score", "word_size", "database_bases", "database_sequences",
                                   "effective_search_space"),
                                  ("evalue_at_minimum_score", "production_evalue", "expected_chance_bound", "informativeness_limit",
                                   "testable_chance_limit"))
COMMON_RULE_SCHEMA = table_schema(("cohort", "group", "condition", "condition_role", "population", "strict_reference_rate", "tolerant_reference_rate",
                                   "strict_control_rate", "tolerant_control_rate", "reference_conversion_rate",
                                   "pooled_control_conversion_rate", "own_control_conversion_rate",
                                   "expected_reference_conversions_under_own_controls", "interpretation"),
                                  ("planned_parents", "strict_eligible_parents", "tolerant_eligible_parents", "common_parents",
                                   "strict_reference_positive", "tolerant_reference_positive", "strict_negative_references",
                                   "converted_references", "lost_strict_references", "control_draws", "strict_control_positive",
                                   "tolerant_control_positive", "strict_negative_control_draws", "converted_control_draws",
                                   "lost_strict_control_draws", "own_strict_negative_control_draws", "own_converted_control_draws",
                                   "references_without_strict_negative_controls"))


def yes(value):
    return "yes" if value else "no"


def historical_zero_filled(outcomes, policy):
    """Re-label a-priori excluded rows with the historical protocol-v1 semantics.

    Retained evidence is unchanged; only the testability gate is removed, so
    completed untestable searches count as zeros exactly as in protocol v1.
    """
    legacy = {key: value for key, value in policy.items() if key != "detectability"}
    rows = []
    for row in outcomes:
        if row["applicability"] in EXCLUDED_STATES:
            completed = row["completion"] == "successful"
            included, unresolved = row["retained_accessions"], row["unresolved_subjects"]
            applicable = applicability(row["length"], row["endpoint"], legacy)
            unavailable = "not_completed" if not completed else "inapplicable" if applicable != "operational_retained_search" else ""
            resolution = ("not_completed" if not completed else "positive_with_unresolved_breadth" if included and unresolved else
                          "unresolved" if unresolved else "resolved_positive" if included else "resolved_zero")
            row = dict(row, applicability=applicable, resolution=resolution,
                       detection_resolution=unavailable or ("unresolved" if unresolved and not included else "resolved"),
                       breadth_resolution=unavailable or ("unresolved" if unresolved else "resolved"),
                       positive=1 if included else -1 if unresolved or not completed else 0)
        rows.append(row)
    return rows


def detection_summaries(outcomes, draws, roles=("primary",), panel_versions=None):
    """All-conditions detection summaries per group and stratum (the matched_contrasts subset used here)."""
    real_rows = [row for row in outcomes if row["role"] == "real"]
    groups = sorted({row["group"] for row in real_rows}) + ["Global union"]
    strata = sorted({(row["endpoint"], row["population"]) for row in real_rows})
    by_stratum = {stratum: [row for row in outcomes if (row["endpoint"], row["population"]) == stratum] for stratum in strata}
    result = {}
    for group in groups:
        for endpoint, population in strata:
            real, _, eligible, attrition = matched_population(by_stratum[endpoint, population], group, "all_conditions", endpoint,
                                                              population, draws, measure="detection", roles=roles,
                                                              panel_versions=panel_versions)
            result[group, endpoint, population] = dict(
                planned_parents=len(real), eligible_parents=len(eligible), censored_parents=attrition["censored"],
                partially_censored_parents=attrition["partial"], chance_limited_parents=attrition["chance_limited"],
                real_positive=sum(row["positive"] for row in eligible.values()),
                real_detection_rate=str(statistics.mean(row["positive"] for row in eligible.values())) if eligible else "not_estimated")
    return result


def length_sensitivity(outcomes, policy, draws):
    """Cohort fractions. Historical inventories (1,288 full / 836 exact-calibrated in the study) are the protocol-v1
    panel, all of its conditions including superseded ones; the primary testable cohort uses primary bundle roles."""
    calibration = {row['query_length']: row for row in policy.get('detectability', [])}
    baseline = policy.get('sensitivity_baseline_evalue', 10.)
    historical_scope = (None, ("v1",))
    cohorts = dict(full_panel_inventory=(outcomes, historical_scope))
    if calibration:
        lengths = {length for length, row in calibration.items() if row['exact_evalue'] <= baseline and row['exact_witness_recovered']}
        cohorts['baseline_exact_calibrated_length_panel'] = ([row for row in outcomes if row['length'] in lengths], historical_scope)
        cohorts['primary_testable_parents'] = ([row for row in outcomes
                                                 if calibration[row['length']][row['endpoint'] + '_testability'] == 'testable'],
                                                (("primary",), None))
    rows = []
    for cohort, (selected, (roles, versions)) in cohorts.items():
        full = detection_summaries(outcomes, draws, roles, versions)
        full_counts = {group: summary['planned_parents'] for (group, _, _), summary in full.items()}
        current = detection_summaries(selected, draws, roles, versions)
        historical = detection_summaries(historical_zero_filled(selected, policy), draws, roles, versions)
        for key, contrast in current.items():
            group, endpoint, population = key
            rows.append(dict(cohort=cohort, group=group, endpoint=endpoint, population=population,
                             panel_parents=full_counts[group], cohort_parents=contrast['planned_parents'],
                             eligible_matched_parents=contrast['eligible_parents'], censored_parents=contrast['censored_parents'],
                             partially_censored_parents=contrast['partially_censored_parents'], chance_limited_parents=contrast['chance_limited_parents'], observed_positive=contrast['real_positive'],
                             denominator=contrast['eligible_parents'], rate=contrast['real_detection_rate'],
                             historical_zero_filled_positive=historical[key]['real_positive'],
                             historical_zero_filled_parents=historical[key]['eligible_parents'],
                             historical_zero_filled_rate=historical[key]['real_detection_rate'],
                             interpretation='Rate among resolved, a-priori testable matched parents; full inventory is not a zero-filled '
                                            'denominator. historical_zero_filled_* reproduce protocol v1 (untestable counted as zero) '
                                            'for sensitivity only'))
        for group, endpoint, population in sorted(set(full) - set(current)):
            rows.append(dict(cohort=cohort, group=group, endpoint=endpoint, population=population,
                             panel_parents=full_counts[group], cohort_parents=0, eligible_matched_parents=0,
                             censored_parents=0, partially_censored_parents=0, chance_limited_parents=0, observed_positive=0,
                             denominator=0, rate='not_estimated', historical_zero_filled_positive=0, historical_zero_filled_parents=0,
                             historical_zero_filled_rate='not_estimated', interpretation='Empty calibrated-length cohort; no rate estimated'))
    return rows


def length_testability_rows(policy):
    """Per-length/rule a-priori classes with the inputs that determine them."""
    rows = []
    for row in sorted(policy.get('detectability', []), key=lambda value: value['query_length']):
        for endpoint, score, evalue, retained in (("strict", row['query_length'], row['exact_evalue'], row['strict_testable']),
                                                  ("tolerant", row['minimum_tolerant_score'], row['tolerant_evalue'],
                                                   row['tolerant_fully_testable'])):
            rows.append(dict(query_length=row['query_length'], endpoint=endpoint, testability=row[endpoint + '_testability'],
                             testability_reason=row.get(endpoint + '_testability_reason', 'not_recorded'),
                             score_retainable=yes(score_retained(row, endpoint)), exact_witness_recovered=yes(row['exact_witness_recovered']),
                             partially_score_censored=yes(row['partially_score_censored']) if endpoint == "tolerant" else "not_applicable",
                             minimum_qualifying_score=score, word_size=row['word_size'], database_bases=row['database_bases'],
                             database_sequences=row['database_sequences'], effective_search_space=row['effective_search_space'],
                             evalue_at_minimum_score=evalue, production_evalue=row['evalue'],
                             expected_chance_bound=row[endpoint + '_expected_chance_bound'],
                             informativeness_limit=row['informativeness_limit'], testable_chance_limit=row['testable_chance_limit']))
    return rows


def rate(numerator, denominator):
    return str(numerator / denominator) if denominator else "not_estimated"


def common_rule_rows(outcomes, policy, draws):
    """Strict/tolerant comparison on parents resolved and eligible for BOTH rules.

    Conversion counts strict-negative units that gain tolerant support. Each
    strict-negative reference is compared with its own controls' conversion
    fraction, retaining draw identities, so parent clustering is preserved.
    """
    rows = []
    for cohort, selected in (("primary_testable", outcomes),
                             ("historical_zero_filled_diagnostic", historical_zero_filled(outcomes, policy))):
        real_rows = [row for row in selected if row["role"] == "real"]
        declared = {row["condition"]: row.get("condition_role", "primary") for row in real_rows}
        contexts = sorted({(row["group"], "all_conditions") for row in real_rows if declared[row["condition"]] == "primary"})
        contexts += [("Global union", "all_conditions")]
        contexts += sorted({(row["group"], row["condition"]) for row in real_rows})
        for group, condition in contexts:
            for population in sorted({row["population"] for row in real_rows}):
                real, controls_s, eligible_s, _ = matched_population(selected, group, condition, "strict", population, draws, measure="detection")
                _, controls_t, eligible_t, _ = matched_population(selected, group, condition, "tolerant", population, draws, measure="detection")
                common = sorted(set(eligible_s) & set(eligible_t))
                counts = defaultdict(int)
                expected, without = 0.0, 0
                for parent in common:
                    strict, tolerant = eligible_s[parent]["positive"], eligible_t[parent]["positive"]
                    counts["strict_reference_positive"] += strict
                    counts["tolerant_reference_positive"] += tolerant
                    counts["strict_negative_references"] += strict == 0
                    counts["converted_references"] += strict == 0 and tolerant == 1
                    counts["lost_strict_references"] += strict == 1 and tolerant == 0
                    by_draw = {row["draw_index"]: row["positive"] for row in controls_t[parent]}
                    pairs = [(row["positive"], by_draw[row["draw_index"]]) for row in controls_s[parent]]
                    negative = sum(s == 0 for s, _ in pairs)
                    converted = sum(s == 0 and t == 1 for s, t in pairs)
                    counts["strict_control_positive"] += sum(s for s, _ in pairs)
                    counts["tolerant_control_positive"] += sum(t for _, t in pairs)
                    counts["strict_negative_control_draws"] += negative
                    counts["converted_control_draws"] += converted
                    counts["lost_strict_control_draws"] += sum(s == 1 and t == 0 for s, t in pairs)
                    if strict == 0:
                        counts["own_strict_negative_control_draws"] += negative
                        counts["own_converted_control_draws"] += converted
                        if negative:
                            expected += converted / negative
                        else:
                            without += 1
                size = len(common)
                rows.append(dict(cohort=cohort, group=group, condition=condition,
                                 condition_role="primary_union" if condition == "all_conditions" else declared[condition],
                                 population=population, planned_parents=len(real),
                                 strict_eligible_parents=len(eligible_s), tolerant_eligible_parents=len(eligible_t), common_parents=size,
                                 control_draws=size * draws, references_without_strict_negative_controls=without,
                                 **{name: counts[name] for name in ("strict_reference_positive", "tolerant_reference_positive",
                                                                     "strict_negative_references", "converted_references",
                                                                     "lost_strict_references", "strict_control_positive",
                                                                     "tolerant_control_positive", "strict_negative_control_draws",
                                                                     "converted_control_draws", "lost_strict_control_draws",
                                                                     "own_strict_negative_control_draws", "own_converted_control_draws")},
                                 strict_reference_rate=rate(counts["strict_reference_positive"], size),
                                 tolerant_reference_rate=rate(counts["tolerant_reference_positive"], size),
                                 strict_control_rate=rate(counts["strict_control_positive"], size * draws),
                                 tolerant_control_rate=rate(counts["tolerant_control_positive"], size * draws),
                                 reference_conversion_rate=rate(counts["converted_references"], counts["strict_negative_references"]),
                                 pooled_control_conversion_rate=rate(counts["converted_control_draws"], counts["strict_negative_control_draws"]),
                                 own_control_conversion_rate=rate(counts["own_converted_control_draws"], counts["own_strict_negative_control_draws"]),
                                 expected_reference_conversions_under_own_controls=(str(expected) if counts["strict_negative_references"] and not without
                                                                                  else "not_estimated"),
                                 interpretation=("Common strict/tolerant cohort; descriptive matched-background comparison, no significance test. "
                                                 + ("Protocol-v2 testable parents only." if cohort == "primary_testable" else
                                                    "Historical zero-filled semantics: a diagnostic, never a primary estimate."))))
    return rows
