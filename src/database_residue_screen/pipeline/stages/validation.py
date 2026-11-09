"""Stage `validation`: recompute and validate the final annotated analysis (production evidence required)."""
from ...example import save_json
from ...release import validate_analysis


def validation(context):
    """Hashes, policy, logical weights, eligibility, derived tables and taxonomy of the final analysis, recomputed."""
    result = validate_analysis(context.upstream["taxonomy"].output / "analysis", context.config["validation"]["require_production"])
    save_json(context.output / "validation.json", result)
    return {key: result[key] for key in ("valid", "production_evidence_checked", "files", "original_parents", "unresolved_outcomes")}
