"""Stage `controls`: pinned dinucleotide-shuffled control draws of every canonical parent of the study panel."""
from ...contracts import ContractError
from ...controls import generate_controls


def controls(context):
    """Generate and audit the matched control draws (uShuffle, k=2, per-parent seeds)."""
    section = context.config["controls"]
    if section["k"] != 2:
        raise ContractError("The control generator implements k=2 (dinucleotide) shuffles only")
    report = generate_controls(context.upstream["panel"].path("queries.fasta"), context.output / "controls", section["draws"], section["seed"])
    return {key: report[key] for key in ("canonical_parents", "physical_queries", "draws_per_parent", "total_draws", "all_invariants_valid",
                                         "parent_identical_draws", "zero_diversity_parents")}
