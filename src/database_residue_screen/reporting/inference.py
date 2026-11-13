"""Inference summary: conditional reference-versus-generator bounds and control precision.

Conditional bounds: the group-by-measure slots (strict, all records; detection and breadth per platform group) of the conditional
relative-entropy Hoeffding contrasts with Holm adjustment over the full family (non-estimable slots stay in the family). The bound is
stage-protected over the declared draw prefixes (5, 10 and 15 control draws per parent). These are Monte Carlo bounds
of the shuffle-generator mean for fixed reference sequences, not contamination, exchangeability or generalisation tests.
Control precision: every fixed-prefix precision row (groups, conditions, global union, endpoints, populations) with its status and the
draws needed for the declared precision target; targets are those of the specification.
"""
from collections import Counter
from pathlib import Path

from ..contracts import ContractError, read_json
from .common import GLOBAL, group_key, ordered_groups, read_tsv, verified_files, write_json, write_tsv

INTERPRETATION = ("Monte Carlo error of the shuffle-generator mean for fixed reference sequences; stage-protected over the declared draw "
                  "prefixes, Holm-adjusted over all group-by-measure slots. Not a contamination, exchangeability, equivalence or "
                  "generalisation test.")
SLOT_FIELDS = ("group", "measure", "endpoint", "population", "status", "planned_parents", "matched_parents", "a_priori_excluded_parents",
               "control_draws", "unique_seeds", "draws_per_parent", "reference_mean", "control_mean", "difference", "interval_lower",
               "interval_upper", "tail_bound", "holm_bound", "holm_family_size", "observed_reference_total", "expected_reference_total")


def write_inference(destination: Path, *, inference: Path, specification: Path | None = None, draw_prefixes=(5, 10, 15)) -> dict:
    destination.mkdir(parents=True, exist_ok=False)
    record = verified_files(inference, read_json(inference / "precision.json"))
    precision = read_tsv(inference / "precision.tsv")
    prefixes = sorted({int(row["draws_per_parent"]) for row in precision})
    declared = list(read_json(specification)["draw_ladder"]) if specification is not None else prefixes
    if prefixes != declared or (draw_prefixes is not None and declared != list(draw_prefixes)):
        raise ContractError(f"Precision draw prefixes {prefixes} differ from the declared {declared} / required {list(draw_prefixes or [])}")
    conditional = read_tsv(inference / "conditional_contrasts.tsv")
    groups = ordered_groups(row["group"] for row in precision if row["condition"] == "all_conditions" and row["group"] != GLOBAL)
    families = {int(row["holm_family_size"]) for row in conditional}
    slots = {(row["group"], row["measure"]) for row in conditional}
    if families != {len(conditional)} or len(slots) != len(conditional) or slots != {(group, measure) for group in groups
                                                                                       for measure in ("detection", "breadth")}:
        raise ContractError(f"Unexpected conditional family: sizes {sorted(families)}, {len(conditional)} slots for groups {groups}")
    stage_alpha = None
    if specification is not None:
        stage_alpha = (1 - read_json(specification)["confidence"]) / len(declared)
    rows = [dict({name: row.get(name, "") for name in SLOT_FIELDS}, draw_prefixes=";".join(map(str, declared)), stage_alpha=stage_alpha,
                 interpretation=INTERPRETATION)
            for row in sorted(conditional, key=lambda row: (group_key(row["group"]), row["measure"] != "detection"))]
    write_tsv(destination / "conditional_slots.tsv", rows)
    fields = list(precision[0]) if precision else ["group"]
    write_tsv(destination / "precision.tsv", precision, fields)
    status = Counter((int(row["draws_per_parent"]), "all_conditions" if row["condition"] == "all_conditions" else "condition",
                      row["status"], row["precision_met"]) for row in precision)
    write_tsv(destination / "precision_status.tsv",
              [dict(draws_per_parent=key[0], scope=key[1], status=key[2], precision_met=key[3], rows=count) for key, count in sorted(status.items())],
              ("draws_per_parent", "scope", "status", "precision_met", "rows"))
    parents_path = inference / "parent_precision.tsv"
    parent_status = Counter()
    if parents_path.is_file():
        for row in read_tsv(parents_path):
            parent_status[int(row["draws"]), row["endpoint"], row["population"], row["measure"], row["status"]] += 1
    write_tsv(destination / "parent_precision_status.tsv",
              [dict(draws=key[0], endpoint=key[1], population=key[2], measure=key[3], status=key[4], parents=count)
               for key, count in sorted(parent_status.items())], ("draws", "endpoint", "population", "measure", "status", "parents"))
    final = [row for row in precision if int(row["draws_per_parent"]) == declared[-1]]
    summary = dict(draw_prefixes=declared, holm_family_size=len(conditional), stage_alpha=stage_alpha,
                   conditional_slots=[{key: row[key] for key in ("group", "measure", "status", "holm_bound", "interval_lower", "interval_upper",
                                                                   "observed_reference_total", "expected_reference_total")} for row in rows],
                   precision_rows=len(precision), final_stage_rows=len(final),
                   final_stage_precision_met=dict(sorted(Counter(row["precision_met"] for row in final).items())),
                   final_stage_status=dict(sorted(Counter(row["status"] for row in final).items())),
                   all_final_stage_targets_met=record.get("all_final_stage_targets_met"), interpretation=INTERPRETATION)
    write_json(destination / "summary.json", summary)
    return summary
