"""Central publication metrics: one computation of every cohort and denominator used by reports, figures and text.

Inputs are a completed, taxonomy-annotated analysis (outcomes, logical queries, search status, policy with the
runtime length calibration), optionally the prepared design (physical IDs) and the search directory (archived
failed or incomplete batch attempts). Nothing is searched or re-annotated.

Cohorts (per endpoint, record population and measure; matched families of a reference parent and all its control
draws, exactly as :func:`database_residue_screen.outcomes.matched_population`):

* endpoint-native cohort: parents whose reference and every control draw are a-priori testable, completed and
  resolved for that endpoint and measure;
* common strict/tolerant cohort: parents in the native cohorts of both rules (same population and measure).

Excluded states (censored_untestable, chance_limited, partially_score_censored, incomplete, unresolved) are counted
neither as 0 nor as 1; they appear only as attrition. Sensitivity panels re-gate the same evidence: the three
declared chance-limit pairs, the ungated panel (no testability gate) and the exact-score-retainable length panel;
in every panel incomplete or unresolved evidence stays excluded and is never counted as zero.
"""
from collections import Counter, defaultdict
from pathlib import Path

from ..contracts import ContractError, read_json, read_table
from ..length_aware import score_retained, testability, testability_reason
from ..outcomes import OUTCOME_SCHEMA, matched_population
from ..policy import applicability
from ..sensitivity import ungated_outcomes
from .common import (EXCLUDED_STATES, GLOBAL, NE, load_analysis, mean, new_directory, ordered_groups, rate, ratio, read_tsv,
                     sha256, verify_manifest, write_json, write_manifest, write_tsv)

CONTRACT = "publication_metrics_v1"
OPERATIONAL = "operational_retained_search"
# (group, label, members): a member is a condition, optionally qualified by a construct class ("condition:class").
# Labels are unique within a group; figures show them as "group / label".
DISPLAY_CLASSES = (
    ("ONT", "Native barcodes", ("ont_barcode",)),
    ("ONT", "K40 contexts", ("ont_native_barcode_context_k40",)),
    ("ONT", "K46 contexts", ("ont_native_barcode_context_k46",)),
    ("ONT", "Rapid barcode contexts", ("ont_rapid_barcode_context_90",)),
    ("ONT", "Adapters and primers", ("ont_ligation_adapter", "ont_rapid_adapter", "ont_pcr_adapter", "ont_cdna_adapter")),
    ("Illumina", "Bare i7 indexes", ("illumina_i7_index",)),
    ("Illumina", "Bare i5 indexes", ("illumina_i5_index",)),
    ("Illumina", "Embedded i7 indexes", ("illumina_i7_index_adapter",)),
    ("Illumina", "Embedded i5 indexes", ("illumina_i5_index_adapter",)),
    ("Illumina", "Adapters and primers", ("illumina_truseq", "illumina_nextera", "illumina_smallrna")),
    ("BGI/DNBSEQ", "Adapters", ("illumina_bgi",)),
    ("PacBio", "Barcodes", ("pacbio_barcode",)),
    ("PacBio", "Barcode-hairpin queries", ("pacbio_barcode_adapter",)),
    ("PacBio", "SMRTbell adapter", ("pacbio_smrtbell_adapter:adapter",)),
    ("PacBio", "C2 primer", ("pacbio_smrtbell_adapter:primer", "pacbio_c2_primer")),
)
DEFAULT_CHANCE_LIMIT_VARIANTS = (("stricter", 0.001, 0.01), ("looser", 0.05, 0.10))


def chance_limit_pairs(policy: dict, variants=DEFAULT_CHANCE_LIMIT_VARIANTS) -> tuple:
    """(label, testable limit, informativeness limit): the baseline pair of the analysis policy, then the variants."""
    baseline = {(row["testable_chance_limit"], row["informativeness_limit"]) for row in policy["detectability"]}
    if len(baseline) != 1:
        raise ContractError("The analysis policy must use one chance-limit pair for every calibrated length")
    pairs = (("baseline", *baseline.pop()),) + tuple((label, float(t), float(i)) for label, t, i in variants)
    labels = [label for label, _, _ in pairs]
    if len(set(labels)) != len(labels) or any(not 0 < t <= i < 1 for _, t, i in pairs):
        raise ContractError("Chance-limit variants need unique labels other than baseline and 0 < testable <= informativeness < 1")
    return pairs
STATES = ("eligible", "censored", "chance_limited", "partial", "short", "incomplete", "unresolved")
ATTRITION = (("censored", "censored"), ("chance_limited", "chance_limited"), ("partially_censored", "partial"),
             ("other_inapplicable", "short"), ("incomplete", "incomplete"), ("unresolved", "unresolved"))
VALUE = {"detection": "positive", "breadth": "retained_accessions"}
TABLE1_FIELDS = ("scope", "order", "group", "label", "conditions", "members", "min_length", "max_length", "population", "endpoint",
                 "cohort", "planned_parents", "eligible_parents", "reference_positive", "reference_rate", "control_draws",
                 "control_positive", "control_rate", "censored", "chance_limited", "partially_censored", "unresolved", "incomplete",
                 "other_inapplicable", "target_limited", "degenerate_null", "native_eligible_parents")


def parse_display_classes(classes=None) -> list[dict]:
    """Normalise display classes given as (group, label, members) tuples or {group, label, members} mappings."""
    result = []
    for number, item in enumerate(DISPLAY_CLASSES if classes is None else classes, 1):
        group, label, members = (item["group"], item["label"], item["members"]) if isinstance(item, dict) else item
        selectors = []
        for member in members:
            condition, _, construct = str(member).partition(":")
            if not condition:
                raise ContractError(f"Display class {label!r}: empty member")
            selectors.append((condition, construct or None))
        result.append(dict(order=number, group=group, label=label, selectors=tuple(selectors), declared=True))
    if len({(row["group"], row["label"]) for row in result}) != len(result):
        raise ContractError("Display class labels must be unique within a group")
    return result


def read_logical(path: Path) -> list[dict]:
    rows = read_tsv(path)
    for row in rows:
        for name in ("draw_index", "length", "seed"):
            row[name] = int(row[name])
    return rows


def read_statuses(path: Path) -> dict[str, dict]:
    rows = read_tsv(path)
    for row in rows:
        row["hsp_rows"] = int(row["hsp_rows"])
        row["unique_targets"] = int(row["unique_targets"])
    statuses = {row["query_id"]: row for row in rows}
    if len(statuses) != len(rows):
        raise ContractError("Duplicate search status rows")
    return statuses


def family_state(original: dict, draws: list[dict], measure: str) -> str:
    """Attrition state of one matched family, in the order of outcomes.matched_population."""
    family = [original] + draws
    for state, test in (("censored", lambda row: row["applicability"] == "censored_untestable"),
                        ("chance_limited", lambda row: row["applicability"] == "chance_limited"),
                        ("partial", lambda row: row["applicability"] == "partially_score_censored"),
                        ("short", lambda row: row["applicability"] != OPERATIONAL),
                        ("incomplete", lambda row: row["completion"] != "successful"),
                        ("unresolved", lambda row: row[measure + "_resolution"] != "resolved")):
        if any(test(row) for row in family):
            return state
    return "eligible"


class Cohorts:
    """Per-parent matched-family states of one outcome set, for every endpoint, population and measure."""

    def __init__(self, outcomes: list[dict], draws: int):
        self.draws = draws
        self.real, self.controls = defaultdict(dict), defaultdict(lambda: defaultdict(list))
        for row in outcomes:
            stratum = row["endpoint"], row["population"]
            if row["role"] == "real":
                self.real[stratum][row["parent_id"]] = row
            elif row["role"] == "control":
                self.controls[stratum][row["parent_id"]].append(row)
        self.strata = sorted(self.real)
        self.states = {}
        for stratum in self.strata:
            for parent in self.real[stratum]:
                controls = self.controls[stratum][parent]
                if len(controls) != draws or {row["draw_index"] for row in controls} != set(range(1, draws + 1)):
                    raise ContractError("Matched outcomes do not contain the declared draw identities")
            for measure in VALUE:
                self.states[stratum + (measure,)] = {parent: family_state(original, self.controls[stratum][parent], measure)
                                                     for parent, original in self.real[stratum].items()}

    def parents(self) -> set[str]:
        return set().union(*(set(rows) for rows in self.real.values())) if self.real else set()

    def eligible(self, endpoint, population, measure, parents) -> set[str]:
        states = self.states[endpoint, population, measure]
        return {parent for parent in parents if states.get(parent) == "eligible"}

    def summary(self, endpoint, population, measure, parents, selected=None) -> dict:
        """Counts of one scope; ``selected`` (default: the native cohort) are the parents whose values enter."""
        stratum, states = (endpoint, population), self.states[endpoint, population, measure]
        planned = [parent for parent in parents if parent in self.real[stratum]]
        counts = Counter(states[parent] for parent in planned)
        selected = sorted(self.eligible(endpoint, population, measure, planned) if selected is None else selected)
        field = VALUE[measure]
        references = [self.real[stratum][parent] for parent in selected]
        controls = [row for parent in selected for row in self.controls[stratum][parent]]
        return dict(planned_parents=len(planned), eligible_parents=len(selected), native_eligible_parents=counts["eligible"],
                    reference_positive=sum(row["positive"] == 1 for row in references),
                    reference_total=sum(row[field] for row in references), reference_mean=mean(row[field] for row in references),
                    control_draws=len(controls), control_positive=sum(row["positive"] == 1 for row in controls),
                    control_total=sum(row[field] for row in controls), control_mean=mean(row[field] for row in controls),
                    target_limited=sum(any(row["retrieval"] == "target_limit_reached"
                                           for row in [self.real[stratum][parent]] + self.controls[stratum][parent]) for parent in planned),
                    degenerate_null=sum(self.real[stratum][parent]["null_state"] == "proven_singleton_null" for parent in planned),
                    **{name: counts[state] for name, state in ATTRITION})

    def check(self, outcomes: list[dict]) -> None:
        """The per-parent states reproduce outcomes.matched_population for the global union of every stratum."""
        for endpoint, population in self.strata:
            stratum_rows = [row for row in outcomes if row["endpoint"] == endpoint and row["population"] == population]
            for measure in VALUE:
                real, _, eligible, attrition = matched_population(stratum_rows, GLOBAL, "all_conditions", endpoint, population,
                                                                  self.draws, measure=measure)
                states = Counter(self.states[endpoint, population, measure].values())
                expected = Counter(eligible=len(eligible), **{state: attrition[state] for state in STATES if state != "eligible"})
                if len(real) != sum(states.values()) or +expected != +states:
                    raise ContractError(f"Cohort states disagree with matched_population ({endpoint}, {population}, {measure})")


def regate(row: dict, policy: dict) -> dict:
    """Re-derive one outcome row's applicability, resolution and value under another policy (as build_outcomes)."""
    completed = row["completion"] == "successful"
    included, unresolved = row["retained_accessions"], row["unresolved_subjects"]
    resolution = ("not_completed" if not completed else "positive_with_unresolved_breadth" if included and unresolved else
                  "unresolved" if unresolved else "resolved_positive" if included else "resolved_zero")
    applicable = applicability(row["length"], row["endpoint"], policy)
    censored = applicable in EXCLUDED_STATES
    if censored:
        resolution = applicable
    unavailable = applicable if censored else "not_completed" if not completed else "inapplicable" if applicable != OPERATIONAL else ""
    return dict(row, applicability=applicable, resolution=resolution,
                detection_resolution=unavailable or ("unresolved" if unresolved and not included else "resolved"),
                breadth_resolution=unavailable or ("unresolved" if unresolved else "resolved"),
                positive=-1 if censored else 1 if included else -1 if unresolved or not completed else 0)


def limit_policy(policy: dict, testable_limit: float, informative_limit: float) -> dict:
    """The policy with every length re-classified under another (testable, informativeness) chance-limit pair."""
    rows = []
    for row in policy["detectability"]:
        new = dict(row, testable_chance_limit=testable_limit, informativeness_limit=informative_limit)
        for endpoint, retained in (("strict", "strict_testable"), ("tolerant", "tolerant_fully_testable")):
            chance = row[endpoint + "_expected_chance_bound"]
            new[endpoint + "_testability"] = testability(row[retained], chance, informative_limit, testable_limit)
            new[endpoint + "_testability_reason"] = testability_reason(score_retained(row, endpoint), row["exact_witness_recovered"],
                                                                       chance, informative_limit, testable_limit)
        rows.append(new)
    return dict(policy, detectability=rows)


def exact_score_lengths(policy: dict) -> set[int]:
    """Lengths whose exact full-length score is retained at the production E-value with a recovered exact witness."""
    baseline = policy.get("sensitivity_baseline_evalue", 10.0)
    return {row["query_length"] for row in policy.get("detectability", [])
            if row["exact_evalue"] <= baseline and row["exact_witness_recovered"]}


def scopes(logical: list[dict], classes: list[dict]) -> list[dict]:
    """Display classes, conditions, platform groups and the global union, with their member parents."""
    real = [row for row in logical if row["role"] == "real"]
    if not real:
        raise ContractError("No reference queries in the logical panel")
    covered, result = set(), []

    def scope(kind, order, group, label, selected, members, declared=True):
        parents = sorted({row["parent_id"] for row in selected})
        lengths = sorted({row["length"] for row in selected})
        return dict(scope=kind, order=order, group=group, label=label, parents=parents, declared=declared,
                    conditions=";".join(sorted({row["condition"] for row in selected})), members=members,
                    physical_ids=len(selected), min_length=lengths[0] if lengths else 0, max_length=lengths[-1] if lengths else 0)

    for item in classes:
        selected = [row for row in real if row["group"] == item["group"] and any(
            row["condition"] == condition and (construct is None or row.get("construct_class") == construct)
            for condition, construct in item["selectors"])]
        if not selected:
            continue
        covered.update(row["logical_id"] for row in selected)
        members = ";".join(condition + (":" + construct if construct else "") for condition, construct in item["selectors"])
        result.append(scope("class", item["order"], item["group"], item["label"], selected, members))
    uncovered = sorted({(row["group"], row["condition"], row.get("construct_class", "")) for row in real if row["logical_id"] not in covered})
    order = len(classes)
    for group, condition, construct in uncovered:
        order += 1
        selected = [row for row in real if (row["group"], row["condition"], row.get("construct_class", "")) == (group, condition, construct)]
        result.append(scope("class", order, group, condition + (f" ({construct})" if construct else ""), selected,
                            f"{condition}:{construct}" if construct else condition, declared=False))
    groups = ordered_groups(row["group"] for row in real)
    for number, group in enumerate(groups, 1):
        result.append(scope("group", 100 + number, group, group, [row for row in real if row["group"] == group], "all_conditions"))
    result.append(scope("global", 100 + len(groups) + 1, GLOBAL, GLOBAL, real, "all_conditions"))
    conditions = sorted({(row["group"], row["condition"]) for row in real}, key=lambda item: (groups.index(item[0]), item[1]))
    for number, (group, condition) in enumerate(conditions, 201):
        result.append(scope("condition", number, group, condition,
                            [row for row in real if (row["group"], row["condition"]) == (group, condition)], condition))
    return result


def cohort_rows(cohorts: Cohorts, scope_list: list[dict]) -> list[dict]:
    rows = []
    for scope in scope_list:
        for endpoint, population in cohorts.strata:
            for measure in VALUE:
                common = (cohorts.eligible("strict", population, measure, scope["parents"])
                          & cohorts.eligible("tolerant", population, measure, scope["parents"]))
                for cohort, selected in (("native", None), ("common_strict_tolerant", common)):
                    summary = cohorts.summary(endpoint, population, measure, scope["parents"], selected)
                    rows.append(dict(scope=scope["scope"], order=scope["order"], group=scope["group"], label=scope["label"],
                                     conditions=scope["conditions"], members=scope["members"], min_length=scope["min_length"],
                                     max_length=scope["max_length"], population=population, endpoint=endpoint, measure=measure,
                                     cohort=cohort, **summary,
                                     reference_rate=ratio(summary["reference_positive"], summary["eligible_parents"]),
                                     control_rate=ratio(summary["control_positive"], summary["control_draws"])))
    return rows


def table1_rows(cohort_table: list[dict]) -> list[dict]:
    """Table 1 source: detection rows of every scope, both cohorts, one column per cohort and endpoint."""
    return [{name: row[name] for name in TABLE1_FIELDS} for row in cohort_table if row["measure"] == "detection"]


def membership_rows(cohorts: Cohorts, logical: list[dict]) -> list[dict]:
    real = defaultdict(list)
    for row in logical:
        if row["role"] == "real":
            real[row["parent_id"]].append(row)
    rows = []
    for endpoint, population in cohorts.strata:
        for measure, field in VALUE.items():
            states = cohorts.states[endpoint, population, measure]
            common = (cohorts.eligible("strict", population, measure, states) & cohorts.eligible("tolerant", population, measure, states))
            for parent in sorted(states):
                members = real[parent]
                reference = cohorts.real[endpoint, population][parent]
                controls = cohorts.controls[endpoint, population][parent]
                rows.append(dict(parent_id=parent, groups=";".join(ordered_groups(row["group"] for row in members)),
                                 conditions=";".join(sorted({row["condition"] for row in members})),
                                 physical_ids=";".join(sorted(row["logical_id"] for row in members)), length=reference["length"],
                                 endpoint=endpoint, population=population, measure=measure, state=states[parent],
                                 common_strict_tolerant=parent in common, reference_value=reference[field],
                                 control_mean=mean(row[field] for row in controls) if states[parent] == "eligible" else None,
                                 target_limited=any(row["retrieval"] == "target_limit_reached" for row in [reference] + controls),
                                 degenerate_null=reference["null_state"] == "proven_singleton_null"))
    return rows


def panel_rows(label: str, cohorts: Cohorts, scope_list: list[dict], parents: set[str], gating: str, **extra) -> list[dict]:
    rows = []
    for scope in scope_list:
        if scope["scope"] not in {"group", "global"}:
            continue
        members = sorted(set(scope["parents"]) & parents)
        for endpoint, population in cohorts.strata:
            summary = cohorts.summary(endpoint, population, "detection", members)
            rows.append(dict(panel=label, gating=gating, **extra, group=scope["group"], endpoint=endpoint, population=population,
                             scope_parents=len(scope["parents"]), panel_parents=summary["planned_parents"],
                             eligible_parents=summary["eligible_parents"], reference_positive=summary["reference_positive"],
                             reference_rate=ratio(summary["reference_positive"], summary["eligible_parents"]),
                             control_draws=summary["control_draws"], control_positive=summary["control_positive"],
                             control_rate=ratio(summary["control_positive"], summary["control_draws"]),
                             **{name: summary[name] for name, _ in ATTRITION}))
    return rows


def conversion_rows(cohorts: Cohorts, scope_list: list[dict], label: str) -> list[dict]:
    """Strict -> tolerant conversion in the common detection cohort against each parent's own controls."""
    rows = []
    populations = sorted({population for _, population in cohorts.strata})
    for scope in scope_list:
        for population in populations:
            if ("strict", population) not in cohorts.real or ("tolerant", population) not in cohorts.real:
                continue
            strict = cohorts.eligible("strict", population, "detection", scope["parents"])
            tolerant = cohorts.eligible("tolerant", population, "detection", scope["parents"])
            common = sorted(strict & tolerant)
            counts, expected, without = Counter(), 0.0, 0
            for parent in common:
                s, t = (cohorts.real[endpoint, population][parent]["positive"] for endpoint in ("strict", "tolerant"))
                counts["strict_reference_positive"] += s
                counts["tolerant_reference_positive"] += t
                counts["strict_negative_references"] += s == 0
                counts["converted_references"] += s == 0 and t == 1
                counts["lost_strict_references"] += s == 1 and t == 0
                by_draw = {row["draw_index"]: row["positive"] for row in cohorts.controls["tolerant", population][parent]}
                pairs = [(row["positive"], by_draw[row["draw_index"]]) for row in cohorts.controls["strict", population][parent]]
                negative, converted = sum(a == 0 for a, _ in pairs), sum(a == 0 and b == 1 for a, b in pairs)
                counts["strict_control_positive"] += sum(a for a, _ in pairs)
                counts["tolerant_control_positive"] += sum(b for _, b in pairs)
                counts["strict_negative_control_draws"] += negative
                counts["converted_control_draws"] += converted
                counts["lost_strict_control_draws"] += sum(a == 1 and b == 0 for a, b in pairs)
                if s == 0:
                    counts["own_strict_negative_control_draws"] += negative
                    counts["own_converted_control_draws"] += converted
                    if negative:
                        expected += converted / negative
                    else:
                        without += 1
            size, draws = len(common), cohorts.draws
            rows.append(dict(cohort=label, scope=scope["scope"], group=scope["group"], label=scope["label"], population=population,
                             planned_parents=len([p for p in scope["parents"] if p in cohorts.real["strict", population]]),
                             strict_eligible_parents=len(strict), tolerant_eligible_parents=len(tolerant), common_parents=size,
                             control_draws=size * draws, references_without_strict_negative_controls=without,
                             **{name: counts[name] for name in (
                                 "strict_reference_positive", "tolerant_reference_positive", "strict_negative_references",
                                 "converted_references", "lost_strict_references", "strict_control_positive", "tolerant_control_positive",
                                 "strict_negative_control_draws", "converted_control_draws", "lost_strict_control_draws",
                                 "own_strict_negative_control_draws", "own_converted_control_draws")},
                             strict_reference_rate=ratio(counts["strict_reference_positive"], size),
                             tolerant_reference_rate=ratio(counts["tolerant_reference_positive"], size),
                             strict_control_rate=ratio(counts["strict_control_positive"], size * draws),
                             tolerant_control_rate=ratio(counts["tolerant_control_positive"], size * draws),
                             reference_conversion_rate=ratio(counts["converted_references"], counts["strict_negative_references"]),
                             pooled_control_conversion_rate=ratio(counts["converted_control_draws"], counts["strict_negative_control_draws"]),
                             own_control_conversion_rate=ratio(counts["own_converted_control_draws"], counts["own_strict_negative_control_draws"]),
                             expected_reference_conversions_under_own_controls=(expected if counts["strict_negative_references"] and not without
                                                                                else None)))
    return rows


def breadth_lower_bound_rows(outcomes: list[dict], scope_list: list[dict]) -> list[dict]:
    """Testable, completed references: resolved breadth and a lower bound adding detection-resolved parents with
    unresolved breadth at their retained resolved accessions."""
    rows = []
    references = defaultdict(dict)
    for row in outcomes:
        if row["role"] == "real" and row["applicability"] == OPERATIONAL and row["completion"] == "successful":
            references[row["endpoint"], row["population"]][row["parent_id"]] = row
    for endpoint, population in sorted(references):
        for scope in scope_list:
            if scope["scope"] not in {"group", "global"}:
                continue
            real = [references[endpoint, population][parent] for parent in scope["parents"] if parent in references[endpoint, population]]
            if not real:
                continue
            resolved = [row for row in real if row["breadth_resolution"] == "resolved"]
            partial = sorted((row for row in real if row["detection_resolution"] == "resolved" and row["breadth_resolution"] != "resolved"),
                             key=lambda row: (-row["retained_accessions"], row["parent_id"]))
            rows.append(dict(group=scope["group"], endpoint=endpoint, population=population, testable_completed_parents=len(real),
                             breadth_resolved_parents=len(resolved), unresolved_breadth_parents=len(partial),
                             resolved_mean_breadth=mean(row["retained_accessions"] for row in resolved),
                             lower_bound_mean_breadth_including_unresolved=mean(row["retained_accessions"] for row in resolved + partial),
                             unresolved_parents=";".join(f"{row['parent_id']}:{row['retained_accessions']}+{row['unresolved_subjects']}unresolved"
                                                         for row in partial)))
    return rows


def testability_rows(policy: dict, logical: list[dict], statuses: dict, limit_pairs=None) -> list[dict]:
    """A-priori classes per query length and rule, under every chance-limit pair, with the panel at that length."""
    physical, parents, controls, searches = Counter(), defaultdict(set), Counter(), defaultdict(set)
    for row in logical:
        if row["role"] == "real":
            physical[row["length"]] += 1
            parents[row["length"]].add(row["parent_id"])
        else:
            controls[row["length"]] += 1
        searches[row["length"]].add(row["sequence_id"])
    pairs = {label: limit_policy(policy, testable, informative)
             for label, testable, informative in (limit_pairs or chance_limit_pairs(policy))}
    retainable = exact_score_lengths(policy)
    rows = []
    for row in sorted(policy["detectability"], key=lambda item: item["query_length"]):
        length = row["query_length"]
        for endpoint, score, evalue, retained in (("strict", length, row["exact_evalue"], "strict_testable"),
                                                  ("tolerant", row["minimum_tolerant_score"], row["tolerant_evalue"], "tolerant_fully_testable")):
            classes = {label: next(item for item in pairs[label]["detectability"] if item["query_length"] == length)[endpoint + "_testability"]
                       for label in pairs}
            if classes["baseline"] != row[endpoint + "_testability"]:
                raise ContractError("The baseline chance-limit pair does not reproduce the calibrated classes")
            rows.append(dict(query_length=length, endpoint=endpoint, testability=row[endpoint + "_testability"],
                             testability_reason=row.get(endpoint + "_testability_reason", "not_recorded"),
                             score_retainable=score_retained(row, endpoint), exact_witness_recovered=row["exact_witness_recovered"],
                             fully_retained=row[retained],
                             partially_score_censored=row["partially_score_censored"] if endpoint == "tolerant" else "not_applicable",
                             exact_score_panel=length in retainable, minimum_qualifying_score=score, evalue_at_minimum_score=evalue,
                             production_evalue=row["evalue"], effective_search_space=row["effective_search_space"],
                             expected_chance_bound=row[endpoint + "_expected_chance_bound"],
                             **{f"class_{label}": classes[label] for label in classes},
                             physical_ids=physical[length], canonical_parents=len(parents[length]), control_draws=controls[length],
                             distinct_searches=len(searches[length] & set(statuses))))
    return rows


def search_roles(logical: list[dict]) -> dict[str, str]:
    roles = defaultdict(set)
    for row in logical:
        roles[row["sequence_id"]].add(row["role"])
    return {sequence: "reference_and_control" if value == {"real", "control"} else "reference" if "real" in value else "control"
            for sequence, value in roles.items()}


def search_state(status: dict | None, forced: bool) -> str:
    if status is None:
        return "incomplete"
    if status["status"] == "success_zero":
        return "score_forced_zero" if forced else "ordinary_zero"
    return "capped" if status["target_limit_reached"] == "yes" else "with_hits_uncapped"


def archived_attempts(search: Path | None) -> dict:
    """Failed or incomplete batch attempts archived in a search directory (each later rerun or still open)."""
    if search is None:
        return dict(searched=False)
    counts = Counter()
    for attempt in sorted(Path(search).glob("batch_*.attempt_*")):
        for name, state in (("failed.json", "failed"), ("incomplete.json", "incomplete")):
            if (attempt / name).is_file():
                queries = read_json(attempt / name).get("query_ids", [])
                counts[state + "_attempts"] += 1
                counts[state + "_attempt_queries"] += len(queries) if isinstance(queries, list) else 0
    return dict(searched=True, **{key: counts[key] for key in ("failed_attempts", "failed_attempt_queries", "incomplete_attempts",
                                                               "incomplete_attempt_queries")})


def zero_search_rows(policy: dict, logical: list[dict], statuses: dict) -> tuple[list, list, dict]:
    calibration = {row["query_length"]: row for row in policy["detectability"]}
    forced = {length for length, row in calibration.items() if row["exact_evalue"] > row["evalue"]}
    roles, lengths = search_roles(logical), {row["sequence_id"]: row["length"] for row in logical}
    by_length = defaultdict(Counter)
    for sequence, role in roles.items():
        length = lengths[sequence]
        by_length[length, role][search_state(statuses.get(sequence), length in forced)] += 1
    states = ("with_hits_uncapped", "capped", "score_forced_zero", "ordinary_zero", "incomplete")
    table = []
    for (length, role), counts in sorted(by_length.items()):
        table.append(dict(query_length=length, search_role=role, distinct_searches=sum(counts.values()),
                          with_hits=counts["with_hits_uncapped"] + counts["capped"], **{state: counts[state] for state in states},
                          zero=counts["score_forced_zero"] + counts["ordinary_zero"], score_forced=length in forced,
                          exact_evalue=calibration[length]["exact_evalue"], production_evalue=calibration[length]["evalue"]))
    units = Counter()
    for row in logical:
        unit = "reference_physical_id" if row["role"] == "real" else "control_draw"
        units[unit, "score_forced" if row["length"] in forced else "retainable",
              search_state(statuses.get(row["sequence_id"]), row["length"] in forced)] += 1
    unit_rows = [dict(unit=unit, length_class=length_class, state=state, count=count) for (unit, length_class, state), count in sorted(units.items())]
    totals = Counter()
    for row in table:
        for state in states:
            totals[state] += row[state]
    summary = dict(distinct_searches=sum(row["distinct_searches"] for row in table), score_forced_lengths=sorted(forced),
                   with_hits=totals["with_hits_uncapped"] + totals["capped"], zero=totals["score_forced_zero"] + totals["ordinary_zero"],
                   **{state: totals[state] for state in states},
                   score_forced_share_of_zero=ratio(totals["score_forced_zero"], totals["score_forced_zero"] + totals["ordinary_zero"]),
                   score_forced_with_hits=sum(row["with_hits"] for row in table if row["score_forced"]))
    return table, unit_rows, summary


def attrition_rows(logical: list[dict], statuses: dict, cohorts: Cohorts, scope_list: list[dict], physical_ids: int | None) -> list[dict]:
    """physical IDs -> canonical parents -> logical queries -> distinct searches -> completion -> testability -> resolution."""
    rows = []
    roles = search_roles(logical)
    by_parent = defaultdict(list)
    for row in logical:
        by_parent[row["parent_id"]].append(row)
    for scope in scope_list:
        if scope["scope"] not in {"group", "global"}:
            continue
        parents = set(scope["parents"])
        members = [row for parent in sorted(parents) for row in by_parent[parent]]
        real = [row for row in members if row["role"] == "real" and (scope["scope"] == "global" or row["group"] == scope["group"])]
        controls = [row for row in members if row["role"] == "control"]
        searches = {row["sequence_id"] for row in members}

        def add(step, count, unit, endpoint="", population="", measure="", detail=""):
            rows.append(dict(group=scope["group"], step=step, endpoint=endpoint, population=population, measure=measure, unit=unit,
                             count=count, detail=detail))
        add("physical_ids", len(real), "reference physical query IDs")
        add("canonical_parents", len(parents), "identical or reverse-complement physical sequences counted once")
        add("parent_condition_memberships", len({(row["parent_id"], row["group"], row["condition"]) for row in real}), "parent x condition")
        add("logical_queries", len(real) + len(controls), "reference physical IDs plus control draws",
            detail=f"references={len(real)};control_draws={len(controls)}")
        add("distinct_searches", len(searches), "distinct searched sequences (identical and reverse-complement sequences once)",
            detail=";".join(f"{role}={count}" for role, count in sorted(Counter(roles[sequence] for sequence in searches).items())))
        completion = Counter("incomplete" if sequence not in statuses else "completed_with_hits" if statuses[sequence]["status"] == "success_hits"
                             else "completed_zero" for sequence in searches)
        for state in ("completed_with_hits", "completed_zero", "incomplete"):
            add(state, completion[state], "distinct searches")
        add("capped", sum(sequence in statuses and statuses[sequence]["target_limit_reached"] == "yes" for sequence in searches),
            "distinct searches at the target limit (subset of completed_with_hits)")
        populations = sorted({population for _, population in cohorts.strata})
        first = "all_records" if "all_records" in populations else populations[0]
        for endpoint in sorted({endpoint for endpoint, _ in cohorts.strata}):
            # Testability and completion are properties of the searches, not of the record population.
            states = cohorts.states[endpoint, first, "detection"]
            planned = [parent for parent in parents if parent in states]
            counts = Counter(states[parent] for parent in planned)
            testable = len(planned) - counts["censored"] - counts["chance_limited"] - counts["partial"] - counts["short"]
            for state, name in (("testable", None), ("censored", "censored"), ("chance_limited", "chance_limited"),
                                ("partially_censored", "partial"), ("other_inapplicable", "short")):
                add("testability_" + state, testable if name is None else counts[name], "parents (matched families)", endpoint)
            add("testable_incomplete", counts["incomplete"], "testable parents with an incomplete family member", endpoint)
            for population in populations:
                for measure in VALUE:
                    measure_counts = Counter(cohorts.states[endpoint, population, measure][parent] for parent in planned)
                    add("resolved", measure_counts["eligible"], "testable, complete and resolved parents (native cohort)",
                        endpoint, population, measure)
                    add("unresolved", measure_counts["unresolved"], "testable, complete parents with unresolved evidence",
                        endpoint, population, measure)
    if physical_ids is not None:
        global_physical = next(row["count"] for row in rows if row["group"] == GLOBAL and row["step"] == "physical_ids")
        if global_physical != physical_ids:
            raise ContractError(f"Design physical IDs ({physical_ids}) differ from the analysed reference queries ({global_physical})")
    return rows


def compute_metrics(analysis: Path, destination: Path, *, design: Path | None = None, search: Path | None = None,
                    display_classes=None, chance_limit_variants=DEFAULT_CHANCE_LIMIT_VARIANTS) -> dict:
    """Write the central metrics tables of one completed analysis into a new directory."""
    analysis = Path(analysis)
    report = load_analysis(analysis)
    policy, draws = report["analysis_policy"], report["draws_per_parent"]
    if not policy.get("detectability"):
        raise ContractError("Publication metrics need the runtime length calibration (analysis_policy.detectability)")
    outcomes = read_table(analysis / "outcomes.tsv", OUTCOME_SCHEMA)
    logical = read_logical(analysis / "logical_queries.tsv")
    statuses = read_statuses(analysis / "search_status.tsv")
    physical_ids = None
    if design is not None:
        physical_ids = len(read_tsv(Path(design) / "physical_queries.tsv"))
    classes = parse_display_classes(display_classes)
    scope_list = scopes(logical, classes)
    base = Cohorts(outcomes, draws)
    base.check(outcomes)
    destination = new_directory(destination)

    cohort_table = cohort_rows(base, scope_list)
    write_tsv(destination / "cohorts.tsv", cohort_table)
    write_tsv(destination / "table1_source.tsv", table1_rows(cohort_table))
    write_tsv(destination / "cohort_membership.tsv", membership_rows(base, logical))
    write_tsv(destination / "display_classes.tsv", [
        dict(order=scope["order"], group=scope["group"], label=scope["label"], members=scope["members"], conditions=scope["conditions"],
             declared=scope["declared"], physical_ids=scope["physical_ids"], canonical_parents=len(scope["parents"]),
             min_length=scope["min_length"], max_length=scope["max_length"]) for scope in scope_list if scope["scope"] == "class"])

    # Chance-limit pairs: the same evidence re-gated; the baseline must reproduce the analysed outcomes exactly.
    limit_rows, testable_lengths = [], {}
    limit_pairs = chance_limit_pairs(policy, chance_limit_variants)
    for label, testable_limit, informative_limit in limit_pairs:
        modified = limit_policy(policy, testable_limit, informative_limit)
        relabelled = [regate(row, modified) for row in outcomes]
        if label == "baseline" and relabelled != outcomes:
            raise ContractError("Re-gating with the baseline chance limits does not reproduce the analysed outcomes")
        testable_lengths[label] = {endpoint: sorted(row["query_length"] for row in modified["detectability"]
                                                    if row[endpoint + "_testability"] == "testable") for endpoint in ("strict", "tolerant")}
        cohorts = base if label == "baseline" else Cohorts(relabelled, draws)
        for row in panel_rows(label, cohorts, scope_list, cohorts.parents(), "testable",
                              testable_chance_limit=testable_limit, informativeness_limit=informative_limit):
            row["testable_lengths"] = ";".join(map(str, testable_lengths[label][row["endpoint"]]))
            limit_rows.append(row)
    write_tsv(destination / "chance_limits.tsv", limit_rows)

    # Sensitivity panels: ungated (no testability gate) and the exact-score-retainable length panel.
    ungated = Cohorts(ungated_outcomes(outcomes, policy), draws)
    retainable_lengths = exact_score_lengths(policy)
    exact_parents = {parent for parent, row in base.real[base.strata[0]].items() if row["length"] in retainable_lengths}
    panel_table = (panel_rows("testable_cohort", base, scope_list, base.parents(), "testable")
                   + panel_rows("ungated_panel", ungated, scope_list, ungated.parents(), "ungated")
                   + panel_rows("exact_score_panel", base, scope_list, exact_parents, "testable")
                   + panel_rows("exact_score_panel", ungated, scope_list, exact_parents, "ungated"))
    write_tsv(destination / "sensitivity_panels.tsv", panel_table)
    write_tsv(destination / "conversion.tsv", conversion_rows(base, scope_list, "testable") + conversion_rows(ungated, scope_list, "ungated"))
    write_tsv(destination / "breadth_lower_bound.tsv", breadth_lower_bound_rows(outcomes, scope_list))
    write_tsv(destination / "testability.tsv", testability_rows(policy, logical, statuses, limit_pairs))
    attrition = attrition_rows(logical, statuses, base, scope_list, physical_ids)
    write_tsv(destination / "attrition.tsv", attrition)
    zero_table, zero_units, zero_summary = zero_search_rows(policy, logical, statuses)
    write_tsv(destination / "zero_searches.tsv", zero_table)
    write_tsv(destination / "zero_search_units.tsv", zero_units)

    headline = {f"{row['endpoint']}|{row['population']}|{row['cohort']}": dict(
        planned=row["planned_parents"], eligible=row["eligible_parents"], positive=row["reference_positive"],
        rate=rate(row["reference_positive"], row["eligible_parents"]), control=f"{row['control_positive']}/{row['control_draws']}")
        for row in cohort_table if row["scope"] == "global" and row["measure"] == "detection"}
    real = [row for row in logical if row["role"] == "real"]
    summary = dict(
        contract=CONTRACT, draws_per_parent=draws, target_limit=report["search_settings"]["max_target_seqs"],
        production_evalue=report["search_settings"]["evalue"],
        counts=dict(physical_ids=len(real), canonical_parents=len({row["parent_id"] for row in real}),
                    conditions=len({row["condition"] for row in real}), groups=ordered_groups(row["group"] for row in real),
                    logical_queries=len(logical), control_draws=sum(row["role"] == "control" for row in logical),
                    distinct_searches=len(statuses)),
        display_classes=[dict(order=scope["order"], group=scope["group"], label=scope["label"], members=scope["members"],
                              declared=scope["declared"]) for scope in scope_list if scope["scope"] == "class"],
        chance_limit_pairs=[dict(label=label, testable_chance_limit=t, informativeness_limit=i, testable_lengths=testable_lengths[label])
                            for label, t, i in limit_pairs],
        exact_score_panel=dict(lengths=sorted(retainable_lengths), parents=len(exact_parents)),
        zero_searches=zero_summary, archived_search_attempts=archived_attempts(search), global_detection=headline,
        cohort_rule="Native: reference and every control draw testable, complete and resolved for the endpoint and measure. "
                    "Common: native for both rules. Excluded states are neither 0 nor 1.")
    write_json(destination / "metrics.json", summary)
    inputs = {"analysis.json": analysis / "analysis.json", "outcomes.tsv": analysis / "outcomes.tsv",
              "logical_queries.tsv": analysis / "logical_queries.tsv", "search_status.tsv": analysis / "search_status.tsv"}
    if design is not None:
        inputs["design/physical_queries.tsv"] = Path(design) / "physical_queries.tsv"
    if search is not None and (Path(search) / "complete.json").is_file():
        inputs["search/complete.json"] = Path(search) / "complete.json"
    write_manifest(destination, CONTRACT, inputs, "Central cohorts and denominators of one completed analysis; no new searches",
                   analysis_sha256=sha256(analysis / "analysis.json"))
    return summary


def read_metrics(directory: Path) -> dict:
    """metrics.json of a verified metrics directory."""
    verify_manifest(directory, CONTRACT)
    return read_json(Path(directory) / "metrics.json")


__all__ = ["DEFAULT_CHANCE_LIMIT_VARIANTS", "chance_limit_pairs", "CONTRACT", "Cohorts", "DISPLAY_CLASSES", "NE", "compute_metrics", "exact_score_lengths",
           "limit_policy", "parse_display_classes", "read_metrics", "regate", "scopes"]
