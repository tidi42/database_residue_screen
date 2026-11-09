"""Target-limit sensitivity: qualifying accession sets of the same reference sequences searched at different target limits.

Each arm is a complete search of the canonical reference search sequences of one panel with identical FASTA, batching,
threads and per-length search spaces; only the target limit differs. Per search sequence the strict and tolerant
qualifying accession sets are built with the matching rules of the analysis: every HSP passing the tolerant predicate
is an event, the strict predicate marks strict events, and subjects resolve to their accession.version (subjects that
do not resolve are kept apart as unresolved, never counted as accessions). A comparison table states per canonical
parent and endpoint the detection state and breadth in both searches and the accession-set differences; every
differing accession is listed. `target_limit.tsv` lists per target-limit arm, canonical parent and endpoint the completion,
the target-limit flag and the retained and unresolved accession counts (all records). Controls are not part of the
comparison.
"""
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

from .contracts import ContractError, read_json, read_table, table_schema, write_table
from .design import digest
from .eligibility import record_eligibility
from .matching import strict_hsp, tolerant_hsp
from .panels import LOGICAL_SCHEMA
from .search import STATUS_SCHEMA, iter_hsps, object_digest, validate_completed_batch
from .sequences import accession_version, read_fasta

ENDPOINTS = ("strict", "tolerant")
COMPARISON_SCHEMA = table_schema(
    ("comparison", "left", "right", "parent_id", "groups", "conditions", "endpoint", "testability", "left_status", "right_status",
     "left_detection", "right_detection", "detection_agrees", "left_target_limit_reached", "right_target_limit_reached",
     "accession_sets_identical", "left_set_sha256", "right_set_sha256"),
    ("length", "left_breadth", "right_breadth", "shared", "left_only", "right_only", "left_unresolved", "right_unresolved"))
DIFFERENCE_SCHEMA = table_schema(("comparison", "parent_id", "endpoint", "accession_version", "present_in"))
ARM_SCHEMA = table_schema(("search", "parent_id", "endpoint", "population", "completion", "target_limit_reached"),
                          ("target_limit", "retained_accessions", "unresolved_subjects"))


def qualifying_sets(search: Path, sequence_ids: set[str]) -> dict[str, dict]:
    """Per search sequence of `sequence_ids`: completion status, target-limit flag and the strict/tolerant accession sets."""
    request, complete = read_json(search / "request.json"), read_json(search / "complete.json")
    if complete["search_id"] != object_digest(request) or not complete["all_queries_completed"]:
        raise ContractError(f"{search}: search incomplete or its completion disagrees with its request")
    statuses = {row["query_id"]: row for row in read_table(search / "status.tsv", STATUS_SCHEMA)}
    missing = sorted(sequence_ids - set(statuses))
    if missing:
        raise ContractError(f"{search}: {len(missing)} reference sequences were not searched")
    result = {name: dict(status=statuses[name]["status"], target_limit_reached=statuses[name]["target_limit_reached"],
                         **{endpoint: set() for endpoint in ENDPOINTS}, **{endpoint + "_unresolved": set() for endpoint in ENDPOINTS})
              for name in sequence_ids}
    for batch in complete["batches"]:
        directory = search / batch["directory"]
        if digest(directory / "complete.json") != batch["completion_sha256"]:
            raise ContractError("Changed search batch completion: " + batch["directory"])
        record = validate_completed_batch(directory, request)
        if not set(record["query_ids"]) & sequence_ids:
            continue
        queries = read_fasta(directory / "queries.fasta")
        for row in iter_hsps(directory / "blast.tsv", queries):
            if row["qseqid"] not in sequence_ids or not tolerant_hsp(row):
                continue
            try:
                accession = accession_version(row["sseqid"])
            except ContractError:
                accession = "unresolved:" + row["sseqid"]
            state = record_eligibility(dict(accession_version=accession), "all_records")[0]
            entry = result[row["qseqid"]]
            for endpoint in ENDPOINTS:
                if endpoint == "tolerant" or strict_hsp(row):
                    entry[endpoint if state == "included" else endpoint + "_unresolved"].add(accession)
    return result


def detection(entry: dict, endpoint: str) -> str:
    """Detection state with the outcome rules of the analysis (before testability gating)."""
    if entry["status"] not in {"success_hits", "success_zero"}:
        return "not_completed"
    if entry[endpoint]:
        return "positive"
    return "unresolved" if entry[endpoint + "_unresolved"] else "zero"


def set_digest(values: set[str]) -> str:
    return hashlib.sha256("".join(value + "\n" for value in sorted(values)).encode("utf-8")).hexdigest()


def compare_searches(panel: Path, searches: dict[str, Path], comparisons: list[tuple[str, str]], destination: Path,
                     testability: dict[int, dict] | None = None, arms=None) -> dict:
    """Per-parent comparison of reference qualifying sets between named searches (`comparisons`: (left, right) labels).
    `arms` names the searches listed in target_limit.tsv (default: every search)."""
    logical = read_table(panel / "logical_queries.tsv", LOGICAL_SCHEMA)
    references = defaultdict(list)
    for row in logical:
        if row["role"] == "real":
            references[row["parent_id"]].append(row)
    if not references:
        raise ContractError("The panel has no reference sequences")
    sets = {label: qualifying_sets(path, set(references)) for label, path in searches.items()}
    rows, differences, totals = [], [], defaultdict(Counter)
    for left, right in comparisons:
        name = f"{left}_vs_{right}"
        for parent, members in sorted(references.items()):
            length = members[0]["length"]
            for endpoint in ENDPOINTS:
                a, b = sets[left][parent], sets[right][parent]
                shared, only_left, only_right = a[endpoint] & b[endpoint], a[endpoint] - b[endpoint], b[endpoint] - a[endpoint]
                state = (testability or {}).get(length, {}).get(endpoint + "_testability", "not_calibrated")
                left_detection, right_detection = detection(a, endpoint), detection(b, endpoint)
                rows.append(dict(comparison=name, left=left, right=right, parent_id=parent,
                                 groups=";".join(sorted({row["group"] for row in members})),
                                 conditions=";".join(sorted({row["condition"] for row in members})), endpoint=endpoint, testability=state,
                                 length=length, left_status=a["status"], right_status=b["status"], left_detection=left_detection,
                                 right_detection=right_detection, detection_agrees="yes" if left_detection == right_detection else "no",
                                 left_breadth=len(a[endpoint]), right_breadth=len(b[endpoint]), shared=len(shared),
                                 left_only=len(only_left), right_only=len(only_right),
                                 left_unresolved=len(a[endpoint + "_unresolved"]), right_unresolved=len(b[endpoint + "_unresolved"]),
                                 left_target_limit_reached=a["target_limit_reached"], right_target_limit_reached=b["target_limit_reached"],
                                 accession_sets_identical="yes" if a[endpoint] == b[endpoint] and a[endpoint + "_unresolved"] == b[endpoint + "_unresolved"] else "no",
                                 left_set_sha256=set_digest(a[endpoint]), right_set_sha256=set_digest(b[endpoint])))
                differences.extend(dict(comparison=name, parent_id=parent, endpoint=endpoint, accession_version=accession, present_in=label)
                                   for label, values in ((left, only_left), (right, only_right)) for accession in sorted(values))
                counter = totals[name, endpoint]
                counter["parents"] += 1
                counter["identical_sets"] += a[endpoint] == b[endpoint]
                counter["detection_discordant"] += left_detection != right_detection
                counter["breadth_changed"] += len(a[endpoint]) != len(b[endpoint])
                counter["left_only_accessions"] += len(only_left)
                counter["right_only_accessions"] += len(only_right)
                if state == "testable":
                    counter["testable_parents"] += 1
                    counter["testable_detection_discordant"] += left_detection != right_detection
                    counter["testable_breadth_changed"] += len(a[endpoint]) != len(b[endpoint])
    destination.mkdir(parents=True, exist_ok=False)
    write_table(destination / "target_limit_comparison.tsv", rows, COMPARISON_SCHEMA)
    write_table(destination / "target_limit_differences.tsv", differences, DIFFERENCE_SCHEMA)
    limits = {label: read_json(path / "request.json")["settings"]["max_target_seqs"] for label, path in searches.items()}
    arms = [dict(search=label, target_limit=limits[label], parent_id=parent, endpoint=endpoint, population="all_records",
                 completion=sets[label][parent]["status"], target_limit_reached=sets[label][parent]["target_limit_reached"],
                 retained_accessions=len(sets[label][parent][endpoint]), unresolved_subjects=len(sets[label][parent][endpoint + "_unresolved"]))
            for label in (searches if arms is None else arms) for parent in sorted(references) for endpoint in ENDPOINTS]
    write_table(destination / "target_limit.tsv", arms, ARM_SCHEMA)
    report = dict(kind="target_limit_comparison_v1", panel_sha256=digest(panel / "panel.json"), reference_parents=len(references),
                  searches={label: dict(directory=str(path.resolve()), completion_sha256=digest(path / "complete.json"),
                                        target_limit=limits[label],
                                        target_limit_reached=sum(entry["target_limit_reached"] == "yes" for entry in sets[label].values()))
                            for label, path in searches.items()},
                  comparisons={f"{name} {endpoint}": dict(counter) for (name, endpoint), counter in sorted(totals.items())},
                  matching="Analysis rules: tolerant HSP events, strict predicate, accession.version subjects; unresolved subjects apart",
                  scope="Canonical reference search sequences only; detection before testability gating, testability class reported",
                  files={name: digest(destination / name) for name in ("target_limit_comparison.tsv", "target_limit_differences.tsv",
                                                                         "target_limit.tsv")})
    (destination / "target_limit.json").write_text(json.dumps(report, indent=1) + "\n", encoding="utf-8")
    return report
