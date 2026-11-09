"""Stage `target_limit`: all reference search sequences re-searched in dedicated arms that differ only in the target limit."""
from ...contracts import ContractError, read_json
from ...subset_search import write_query_sets
from ...target_limit import compare_searches
from .searching import resumable_search

IDENTICAL = ("queries_sha256", "query_count", "database_prefix", "database_evidence_sha256", "executable_sha256", "threads", "batch_size",
             "runner_sha256", "parser_sha256", "calibration_code_sha256")


def target_limit(context):
    """Search the references once per configured target limit and compare the qualifying accession sets per parent:
    consecutive arms with each other and the arm at the production limit with the production search."""
    section, arms = context.config["search"], sorted(context.config["target_limit"]["arms"])
    if section["target_limit"] not in arms:
        raise ContractError("The target-limit arms must include the production target limit")
    panel = context.upstream["design"].output / "search_panel"
    sets = write_query_sets(panel, context.output / "query_sets")
    references = context.output / "query_sets" / "references.fasta"
    searches, requests, resumed = {}, {}, {}
    for limit in arms:
        label = f"limit_{limit}"
        destination = context.output / f"search_{limit}"
        report, restored = resumable_search(context, references, destination, section, limit)
        searches[label], requests[label] = destination, read_json(destination / "request.json")
        resumed[label] = len(restored["copied_batches"]) if restored else 0
    first = requests[f"limit_{arms[0]}"]
    for label, request in requests.items():
        differing = [key for key in IDENTICAL if request.get(key) != first.get(key)]
        settings = {key: value for key, value in request["settings"].items() if key != "max_target_seqs"}
        if differing or settings != {key: value for key, value in first["settings"].items() if key != "max_target_seqs"}:
            raise ContractError(f"Target-limit arm {label} differs from the first arm beyond its target limit: {', '.join(differing) or 'settings'}")
        spaces = [(row["query_length"], row["effective_search_space"], row["evalue"]) for row in request["length_calibration"]["lengths"]]
        if spaces != [(row["query_length"], row["effective_search_space"], row["evalue"]) for row in first["length_calibration"]["lengths"]]:
            raise ContractError(f"Target-limit arm {label} uses different per-length search spaces")
    production = context.upstream["search"].output / "search"
    searches["production"] = production
    comparisons = [(f"limit_{left}", f"limit_{right}") for left, right in zip(arms, arms[1:])]
    comparisons.append((f"limit_{section['target_limit']}", "production"))
    calibration = read_json(production / "request.json")["length_calibration"]["lengths"]
    report = compare_searches(panel, searches, comparisons, context.output / "comparison", {row["query_length"]: row for row in calibration},
                              arms=[f"limit_{limit}" for limit in arms])
    return dict(reference_parents=report["reference_parents"], references=sets["counts"]["references"], arms=arms,
                comparisons=[f"{left}_vs_{right}" for left, right in comparisons], resumed_batches=resumed,
                identical_sets={name: values["identical_sets"] for name, values in report["comparisons"].items()})
