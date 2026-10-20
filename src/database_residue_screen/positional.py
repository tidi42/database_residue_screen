"""Positional robustness for strict terminal events (protocol v2, R09 and D-D).

All variants start from the primary event population: strict, a-priori testable,
primary-condition events in eligible linear genomic records. Variants:
  all_events                 the primary positional definition
  one_per_accession          one event per accession.version
  one_per_bioproject         one event per BioProject; records without a project are their own unit (counted, never dropped)
  excluding_targeted_insdc   records flagged targeted from structured INSDC annotation removed (secondary)
The kept event of a unit is the one with the smallest sha256(accession|parent|start|end|strand), which is independent of
position. Windows are the declared policy windows plus the protocol's 25-nt diagnostic window.
"""
import glob
import hashlib
import json
from pathlib import Path

from .contracts import ContractError, read_json, read_table, table_schema
from .design import digest
from .eligibility import event_eligibility
from .policy import applicability, primary_event
from .reanalysis import ANNOTATED_EVENT_SCHEMA, position_summary

DIAGNOSTIC_WINDOW = 25
SENSITIVITY_SCHEMA = table_schema(("variant", "group", "population", "fold_n2", "unit_rule"),
                                  ("window", "minimum_length", "events", "terminal_events", "units_without_project"),
                                  ("expected_n2",))


def selection_key(row: dict) -> str:
    text = "|".join(str(row[key]) for key in ("accession_version", "parent_id", "sstart", "send", "canonical_strand"))
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def one_event_per(events: list[dict], unit) -> list[dict]:
    chosen = {}
    for row in events:
        key = unit(row), row["group"]
        if key not in chosen or selection_key(row) < selection_key(chosen[key]):
            chosen[key] = row
    return list(chosen.values())


def primary_positional_events(analysis: Path) -> tuple[list[dict], dict]:
    report = read_json(analysis / "analysis.json")
    for name, checksum in report["files"].items():
        if digest(analysis / name) != checksum:
            raise ContractError("Analysis artifact changed: " + name)
    policy = report["analysis_policy"]
    events = read_table(analysis / "real_events.tsv", ANNOTATED_EVENT_SCHEMA)
    decisions = event_eligibility(events, ["linear_genomic"])
    selected = [row for row in events if row["strict"] == "yes" and primary_event(row, policy)
                and ("detectability" not in policy or applicability(row["qlen"], "strict", policy) == "operational_retained_search")
                and decisions[row["accession_version"], "linear_genomic"][0] == "included"]
    return selected, policy


def projects_from_metadata(metadata: Path) -> dict[str, str]:
    """Exact-version BioProject IDs from archived ESummary payloads (projectid; empty when absent)."""
    projects = {}
    for path in sorted(glob.glob(str(metadata / "*esummary_*.json"))):
        payload = json.loads(Path(path).read_text())
        for key, value in payload.get("result", {}).items():
            if isinstance(value, dict) and value.get("accessionversion"):
                project = str(value.get("projectid") or "").strip()
                projects[value["accessionversion"]] = "" if project in {"", "0"} else project
    return projects


def positional_sensitivity(analysis: Path, projects: dict[str, str], targeted: dict[str, dict] | None = None) -> list[dict]:
    events, policy = primary_positional_events(analysis)
    groups = sorted({row["group"] for row in events})
    variants = {"all_events": (events, "every primary strict event"),
                "one_per_accession": (one_event_per(events, lambda row: row["accession_version"]), "one event per accession.version"),
                "one_per_bioproject": (one_event_per(events, lambda row: projects.get(row["accession_version"]) or "record:" + row["accession_version"]),
                                       "one event per BioProject; records without project are their own unit")}
    if targeted is not None:
        variants["excluding_targeted_insdc"] = ([row for row in events if targeted.get(row["accession_version"], {}).get("targeted_insdc") != "yes"],
                                                "records flagged targeted from structured INSDC annotation removed")
    without_project = len({row["accession_version"] for row in events if not projects.get(row["accession_version"])})
    windows = sorted(set(policy["position"]["windows"]) | {DIAGNOSTIC_WINDOW})
    rows = []
    for variant, (selected, rule) in variants.items():
        pooled = [dict(row, group="Global union") for row in selected]
        for group in groups + ["Global union"]:
            source = pooled if group == "Global union" else selected
            for window in windows:
                for minimum in policy["position"]["minimum_lengths"]:
                    summary = position_summary(source, group, window, minimum, True)
                    rows.append(dict(variant=variant, group=group, population=summary["population"], window=window, minimum_length=minimum,
                                     events=summary["events"], terminal_events=summary["terminal_events"], expected_n2=summary["expected_n2"],
                                     fold_n2=summary["fold_n2"], unit_rule=rule, units_without_project=without_project))
    return rows
