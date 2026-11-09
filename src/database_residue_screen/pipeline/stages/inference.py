"""Stage `inference`: control-precision specification with explicit draw and retention ladders and its assessment."""
from ...precision import assess_precision, write_precision_spec

SETTINGS = ("confidence", "detection_half_width", "breadth_absolute_half_width", "breadth_relative_half_width", "threads", "batch_size",
            "timeout_seconds")


def inference(context):
    """Bind the precision specification to the search panel and assess the declared draw prefixes on the final analysis."""
    section = context.config["inference"]
    prepared = write_precision_spec(context.upstream["design"].output / "search_panel", context.output / "specification", section["draw_ladder"],
                                section["retention_ladder"], section["maximum_draws"], **{key: section[key] for key in SETTINGS})
    report = assess_precision(context.upstream["taxonomy"].output / "analysis", context.output / "specification" / "precision_specification.json",
                              context.output / "precision")
    return dict(draw_ladder=prepared["draw_ladder"], retention_ladder=prepared["retention_ladder"],
                all_final_stage_targets_met=report["all_final_stage_targets_met"], files=sorted(report["files"]))
