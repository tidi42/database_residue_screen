"""Registration of the publication reporting stages in the study pipeline (see docs/reports.md).

``register(pipeline)`` declares the configuration sections ``[metrics]``, ``[reports]`` and ``[figures]`` and the stages

  metrics          central cohorts and denominators of the final analysis
  reports          the publication reports
  figure_sources   every table read by the R figure scripts
  figures          the rendered figures (tool ``Rscript``; the R scripts directory is a hashed stage input)

on top of the acquisition and analysis stages named below, which must already be registered. Stage code lives in
:mod:`database_residue_screen.pipeline.stages`; this module only registers it, so its own code is in no stage closure.
"""
from .stages import figure_sources, figures, metrics, reports

METRICS_INPUTS = ("design", "search", "taxonomy")
REPORTS_INPUTS = ("metrics", "panel", "design", "search", "metadata", "taxonomy", "insdc", "positional", "record_cases", "census",
                  "exact_scan", "scan_validation", "cap_review", "target_limit", "benchmark", "inference", "validation")
FIGURE_SOURCE_INPUTS = ("metrics", "reports", "taxonomy", "design", "controls", "search", "inference", "benchmark", "record_cases")
FIGURE_INPUTS = ("figure_sources",)

DISPLAY_CLASS = {"type": "object", "additionalProperties": False, "required": ["group", "label", "members"],
                 "properties": {"group": {"type": "string", "minLength": 1}, "label": {"type": "string", "minLength": 1},
                                "members": {"type": "array", "minItems": 1, "items": {"type": "string", "pattern": "^[^:]+(:[^:]+)?$"},
                                            "description": "condition or condition:construct_class"}}}
SECTIONS = {
    "metrics": {"type": "object", "properties": {
        "display_classes": {"type": "array", "items": DISPLAY_CLASS, "default": [],
                            "description": "Table 1 display classes in order; empty: the built-in Paper 3 classes"},
        "chance_limit_variants": {
            "type": "array", "uniqueItems": True,
            "default": [{"label": "stricter", "testable_chance_limit": 0.001, "informativeness_limit": 0.01},
                        {"label": "looser", "testable_chance_limit": 0.05, "informativeness_limit": 0.10}],
            "items": {"type": "object", "additionalProperties": False,
                      "required": ["label", "testable_chance_limit", "informativeness_limit"],
                      "properties": {"label": {"type": "string", "pattern": "^[a-z][a-z0-9_]*$"},
                                     "testable_chance_limit": {"type": "number", "exclusiveMinimum": 0, "exclusiveMaximum": 1},
                                     "informativeness_limit": {"type": "number", "exclusiveMinimum": 0, "exclusiveMaximum": 1}}},
            "description": "Sensitivity pairs re-gating the same evidence; the baseline pair comes from the analysis policy"}}},
    "reports": {"type": "object", "properties": {
        "benchmark_arm": {"type": "string", "default": "permissive_diagnostic",
                          "description": "Benchmark arm whose persistent tolerant misses are decomposed"},
        "draw_prefixes": {"type": "array", "items": {"type": "integer", "minimum": 1}, "minItems": 1, "default": [5, 10, 15],
                          "description": "Required control draw prefixes of the inference stage"}}},
    "figures": {"type": "object", "properties": {
        "scripts": {"type": "string", "format": "path", "default": "figures",
                    "description": "Directory with build_figures.sh and the R figure scripts (relative to the configuration)"}}},
}


def register(pipeline):
    """Add the reporting configuration sections and the metrics, reports, figure_sources and figures stages."""
    for name, schema in SECTIONS.items():
        pipeline.section(name, schema)
    pipeline.stage(depends_on=METRICS_INPUTS, config_sections=("metrics",))(metrics.metrics)
    pipeline.stage(depends_on=REPORTS_INPUTS, config_sections=("reports",))(reports.reports)
    pipeline.stage(depends_on=FIGURE_SOURCE_INPUTS)(figure_sources.figure_sources)
    pipeline.stage(depends_on=FIGURE_INPUTS, config_sections=("figures",), tools=("Rscript",),
                   input_files=lambda config: [config["figures"]["scripts"]])(figures.figures)
    return pipeline
