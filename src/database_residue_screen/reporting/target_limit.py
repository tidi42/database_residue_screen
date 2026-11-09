"""Target-cap reports: cap-review summary and the 500 versus 2,000 target-limit comparison.

cap_review_searches.tsv joins every cap-saturated search of the production search (cap-review stage) to its logical
memberships and production outcomes; cap_review_counts.tsv counts strict/tolerant review states by search role.

The target-limit comparison reads the target_limit stage table (long format, one row per canonical reference parent,
target limit, endpoint and record population; see docs/reports.md) and adds the production values (2,000 targets,
mixed reference/control batches) of the final analysis. Per parent it reports detection and breadth at every limit and
whether they agree; the summary counts agreeing parents per group. Incomplete or unresolved values are never zero.
"""
from collections import Counter, defaultdict
from pathlib import Path

from ..contracts import ContractError
from .common import GLOBAL, mean, ordered_groups, ratio, read_tsv, write_json, write_tsv

LIMIT_COLUMNS = ("target_limit", "max_target_seqs", "limit")
COMPLETED = {"successful", "success_hits", "success_zero", "completed"}


def read_cap_review(path: Path) -> list[dict]:
    rows = read_tsv(path)
    required = {"query_id", "query_length", "retained_subjects", "strict_subjects", "tolerant_subjects", "strict_status", "tolerant_status"}
    if rows and not required <= set(rows[0]):
        raise ContractError(f"{path.name}: not a cap-review table")
    return rows


def write_cap_review(destination: Path, *, cap_review: Path, logical: list[dict], outcomes: list[dict], statuses: dict) -> dict:
    destination.mkdir(parents=True, exist_ok=False)
    reviewed = {row["query_id"]: row for row in read_cap_review(cap_review)}
    capped = sorted(sequence for sequence, status in statuses.items() if status["target_limit_reached"] == "yes")
    if set(reviewed) - set(statuses):
        raise ContractError("Cap review lists searches outside the analysed search")
    members, parent_groups = defaultdict(list), defaultdict(set)
    for row in logical:
        members[row["sequence_id"]].append(row)
        if row["role"] == "real":
            parent_groups[row["parent_id"]].add(row["group"])
    production = {(row["sequence_id"], row["endpoint"]): row for row in outcomes if row["population"] == "all_records"}
    rows = []
    for sequence in sorted(set(capped) | set(reviewed)):
        review, owners = reviewed.get(sequence, {}), members[sequence]
        roles = sorted({row["role"] for row in owners})
        strict, tolerant = production.get((sequence, "strict"), {}), production.get((sequence, "tolerant"), {})
        rows.append(dict(sequence_id=sequence, length=owners[0]["length"] if owners else review.get("query_length", ""),
                         search_role="reference_and_control" if roles == ["control", "real"] else "reference" if roles == ["real"] else "control",
                         groups=";".join(ordered_groups(group for row in owners for group in parent_groups[row["parent_id"]])),
                         physical_ids=";".join(sorted(row["logical_id"] for row in owners if row["role"] == "real")),
                         control_draws=";".join(f"{row['parent_id']}:{row['draw_index']}" for row in sorted(owners, key=lambda row: row["logical_id"])
                                                if row["role"] == "control"),
                         logical_memberships=len(owners), hsp_rows=statuses[sequence]["hsp_rows"], targets=statuses[sequence]["unique_targets"],
                         reviewed=bool(review), retained_subjects=review.get("retained_subjects", ""),
                         strict_subjects=review.get("strict_subjects", ""), tolerant_subjects=review.get("tolerant_subjects", ""),
                         weakest_retained_bitscore=review.get("weakest_retained_bitscore", ""),
                         strict_threshold_bits=review.get("strict_threshold_bits", ""), tolerant_threshold_bits=review.get("tolerant_threshold_bits", ""),
                         strict_status=review.get("strict_status", "not_reviewed"), tolerant_status=review.get("tolerant_status", "not_reviewed"),
                         exact_scan_accessions=review.get("exact_scan_accessions", ""),
                         strict_retained=strict.get("retained_accessions", ""), tolerant_retained=tolerant.get("retained_accessions", ""),
                         strict_resolution=strict.get("resolution", ""), tolerant_resolution=tolerant.get("resolution", "")))
    fields = ("sequence_id", "length", "search_role", "groups", "physical_ids", "control_draws", "logical_memberships", "hsp_rows", "targets",
              "reviewed", "retained_subjects", "strict_subjects", "tolerant_subjects", "weakest_retained_bitscore", "strict_threshold_bits",
              "tolerant_threshold_bits", "strict_status", "tolerant_status", "exact_scan_accessions", "strict_retained", "tolerant_retained",
              "strict_resolution", "tolerant_resolution")
    write_tsv(destination / "cap_review_searches.tsv", rows, fields)
    counts = Counter((row["search_role"], row["strict_status"], row["tolerant_status"]) for row in rows)
    write_tsv(destination / "cap_review_counts.tsv", [dict(search_role=key[0], strict_status=key[1], tolerant_status=key[2], searches=count)
                                                      for key, count in sorted(counts.items())],
              ("search_role", "strict_status", "tolerant_status", "searches"))
    summary = dict(capped_searches=len(capped), reviewed_searches=len(reviewed), unreviewed_capped_searches=len(set(capped) - set(reviewed)),
                   reviewed_not_capped=len(set(reviewed) - set(capped)),
                   cap_censored_possible=sum("cap_censored_possible" in (row["strict_status"], row["tolerant_status"]) for row in rows),
                   accepted_support_in_capped_searches=sum(bool(row["strict_retained"]) or bool(row["tolerant_retained"]) for row in rows
                                                           if row["strict_retained"] != "" and row["tolerant_retained"] != ""))
    write_json(destination / "summary.json", summary)
    return summary


def _limit(row: dict) -> int:
    for name in LIMIT_COLUMNS:
        if row.get(name, "") != "":
            return int(str(row[name]).removeprefix("limit_"))
    raise ContractError("Target-limit table row without a target limit column (" + ", ".join(LIMIT_COLUMNS) + ")")


def normalise_limit_row(row: dict) -> dict:
    """One target-limit observation in comparison form; incomplete or unresolved evidence is never zero."""
    parent = row.get("parent_id") or row.get("sequence_id")
    if not parent:
        raise ContractError("Target-limit table row without parent_id/sequence_id")
    completed = row.get("completion", row.get("status", "successful")) in COMPLETED
    retained, unresolved = int(row.get("retained_accessions", 0) or 0), int(row.get("unresolved_subjects", 0) or 0)
    capped = row.get("target_limit_reached", "") == "yes" or row.get("retrieval", "") == "target_limit_reached"
    detection = None if not completed or (unresolved and not retained) else int(retained > 0)
    breadth = None if not completed or unresolved else retained
    return dict(parent_id=parent, endpoint=row["endpoint"], population=row["population"], limit=_limit(row), completed=completed,
                retained=retained, unresolved=unresolved, capped=capped, detection=detection, breadth=breadth)


def production_values(outcomes: list[dict]) -> dict:
    values = {}
    for row in outcomes:
        if row["role"] != "real":
            continue
        completed = row["completion"] == "successful"
        retained, unresolved = row["retained_accessions"], row["unresolved_subjects"]
        values[row["parent_id"], row["endpoint"], row["population"]] = dict(
            completed=completed, retained=retained, unresolved=unresolved, capped=row["retrieval"] == "target_limit_reached",
            detection=None if not completed or (unresolved and not retained) else int(retained > 0),
            breadth=None if not completed or unresolved else retained, applicability=row["applicability"])
    return values


def write_target_limit(destination: Path, *, table: Path, outcomes: list[dict], logical: list[dict], production_limit: int) -> dict:
    destination.mkdir(parents=True, exist_ok=False)
    observations = [normalise_limit_row(row) for row in read_tsv(table)]
    limits = sorted({row["limit"] for row in observations})
    if len(limits) < 2:
        raise ContractError("The target-limit comparison needs at least two target limits")
    arms = defaultdict(dict)
    for row in observations:
        key = row["parent_id"], row["endpoint"], row["population"]
        if row["limit"] in arms[key]:
            raise ContractError(f"Duplicate target-limit observation: {key} at {row['limit']}")
        arms[key][row["limit"]] = row
    production = production_values(outcomes)
    groups = defaultdict(set)
    for row in logical:
        if row["role"] == "real":
            groups[row["parent_id"]].add(row["group"])
    parents = sorted({key[0] for key in arms})
    unknown = [parent for parent in parents if parent not in groups]
    if unknown:
        raise ContractError(f"Target-limit parents outside the analysed panel: {', '.join(unknown[:5])}")
    missing = sorted({key[0] for key in production} - set(parents))
    rows = []
    for key in sorted(arms):
        parent, endpoint, population = key
        observed, reference = arms[key], production.get(key)
        if reference is None:
            raise ContractError(f"Target-limit row without a production outcome: {key}")
        if set(observed) != set(limits):
            raise ContractError(f"Parent {parent} lacks a target-limit arm ({endpoint}, {population})")
        detections = [observed[limit]["detection"] for limit in limits]
        breadths = [observed[limit]["breadth"] for limit in limits]
        resolved = None not in detections
        row = dict(parent_id=parent, groups=";".join(ordered_groups(groups[parent])), endpoint=endpoint, population=population,
                   applicability=reference["applicability"])
        for limit in limits:
            item = observed[limit]
            row.update({f"completed_{limit}": item["completed"], f"capped_{limit}": item["capped"], f"retained_{limit}": item["retained"],
                        f"unresolved_{limit}": item["unresolved"], f"detection_{limit}": item["detection"], f"breadth_{limit}": item["breadth"]})
        row.update(completed_production=reference["completed"], capped_production=reference["capped"], retained_production=reference["retained"],
                   unresolved_production=reference["unresolved"], detection_production=reference["detection"],
                   breadth_production=reference["breadth"],
                   detection_resolved_all=resolved and reference["detection"] is not None,
                   detection_agrees=None if not resolved else len(set(detections)) == 1,
                   breadth_agrees=None if None in breadths else len(set(breadths)) == 1,
                   highest_limit_matches_production=(None if observed[limits[-1]]["breadth"] is None or reference["breadth"] is None
                                                     else observed[limits[-1]]["breadth"] == reference["breadth"]))
        rows.append(row)
    write_tsv(destination / "target_limit_parents.tsv", rows)
    summary_rows = []
    for endpoint, population in sorted({(key[1], key[2]) for key in arms}):
        selected = [row for row in rows if (row["endpoint"], row["population"]) == (endpoint, population)]
        for group in ordered_groups(group for row in selected for group in row["groups"].split(";")) + [GLOBAL]:
            members = [row for row in selected if group == GLOBAL or group in row["groups"].split(";")]
            for gate in ("all_parents", "testable"):
                chosen = [row for row in members if gate == "all_parents" or row["applicability"] == "operational_retained_search"]
                summary = dict(group=group, endpoint=endpoint, population=population, gate=gate, parents=len(chosen))
                for limit in limits:
                    breadth = [row[f"breadth_{limit}"] for row in chosen if row[f"breadth_{limit}"] is not None]
                    summary.update({f"capped_{limit}": sum(row[f"capped_{limit}"] for row in chosen),
                                    f"detected_{limit}": sum(row[f"detection_{limit}"] == 1 for row in chosen),
                                    f"breadth_resolved_{limit}": len(breadth), f"mean_breadth_{limit}": mean(breadth)})
                agreeing = [row for row in chosen if row["detection_agrees"] is not None]
                breadth_pairs = [row for row in chosen if row["breadth_agrees"] is not None]
                summary.update(detection_compared=len(agreeing), detection_changed=sum(not row["detection_agrees"] for row in agreeing),
                               breadth_compared=len(breadth_pairs), breadth_changed=sum(not row["breadth_agrees"] for row in breadth_pairs),
                               breadth_unchanged_fraction=ratio(sum(row["breadth_agrees"] for row in breadth_pairs), len(breadth_pairs)),
                               production_compared=sum(row["highest_limit_matches_production"] is not None for row in chosen),
                               production_differs=sum(row["highest_limit_matches_production"] is False for row in chosen))
                summary_rows.append(summary)
    write_tsv(destination / "target_limit_summary.tsv", summary_rows)
    result = dict(limits=limits, production_limit=production_limit, parents=len(parents), production_parents_without_arms=len(missing),
                  rows=len(rows), detection_changed=sum(row["detection_agrees"] is False for row in rows),
                  breadth_changed=sum(row["breadth_agrees"] is False for row in rows),
                  production_differs=sum(row["highest_limit_matches_production"] is False for row in rows),
                  statement="Dedicated reference-only arms with identical settings except the target limit; production values from the "
                            "final analysis. Incomplete or unresolved values are never zero; every declared comparison is reported.")
    write_json(destination / "summary.json", result)
    return result
