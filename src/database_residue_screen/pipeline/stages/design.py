"""Stage `design`: bind the study panel, its control draws and the analysis policy into the search panel."""
from ...panels import prepare_panel
from ...policy import default_policy, validate_policy


def analysis_policy(section: dict) -> dict:
    policy = default_policy()
    policy.update(record_populations=list(section["record_populations"]), minimum_query_length=dict(section["minimum_query_length"]),
                  position=dict(primary_window=section["primary_window"], primary_minimum_length=section["primary_minimum_length"],
                                windows=list(section["windows"]), minimum_lengths=list(section["minimum_lengths"])))
    return validate_policy(policy)


def design(context):
    """Physical and logical identities, distinct search sequences and the hashed analysis policy."""
    panel = context.upstream["panel"].output
    report = prepare_panel(panel / "queries.fasta", panel / "queries.tsv", context.upstream["controls"].output / "controls",
                           context.output / "search_panel", context.config["controls"]["draws"],
                           policy=analysis_policy(context.config["analysis_policy"]))
    return {key: report[key] for key in ("source_physical_queries", "physical_real_queries", "canonical_parents", "control_draws",
                                         "draws_per_parent", "logical_queries", "distinct_search_sequences")}
