"""The `screen` pipeline: a complete reference-panel screen of one BLAST nucleotide database from one study configuration.

This module only declares the configuration sections and registers the stages; every stage function lives in its own
module under `stages` (see docs/study-workflow.md for the stage list, data flow, configuration and outputs).
"""
from pathlib import Path

from ..census import SPECIFICATION_SCHEMA as CENSUS_SCHEMA
from ..materials import CATALOGUE_CHECK_SCHEMA, CATALOGUE_RULE_SCHEMA, ELEMENT_RULE_SCHEMA, TEMPLATE_RULE_SCHEMA
from ..record_context import SPECIFICATION_SCHEMA as CASE_SCHEMA
from ..study import Pipeline
from .stages import (analysis, benchmark, cap_review, census, controls, database, design, events, exact_scan, inference, insdc, metadata,
                     panel, positional, record_cases, scan_validation, search, sources, target_limit, taxonomy, validation)

PATH = {"type": "string", "format": "path"}
POSITIVE = {"type": "integer", "minimum": 1}
COUNT = {"type": "integer", "minimum": 0}
PROBABILITY = {"type": "number", "exclusiveMinimum": 0, "exclusiveMaximum": 1}
LADDER = {"type": "array", "minItems": 1, "uniqueItems": True, "items": POSITIVE}


def closed(properties: dict, required=(), **extra) -> dict:
    return dict(type="object", additionalProperties=False, required=list(required), properties=properties, **extra)


PIPELINE = Pipeline("screen", "Reference-panel screen: panel, controls, BLAST search, analysis, metadata, taxonomy, positional, "
                              "record cases, census, exact scan, target limits, benchmark, inference and validation",
                    secrets=("NCBI_API_KEY",))

PIPELINE.section("sources", closed({
    "manifest": dict(PATH, description="Pinned source manifest (URL, retrieval date, SHA-256; no content)"),
    "directory": dict(PATH, description="Local source directory: matching user-supplied copies are used, missing files are fetched into it"),
    "fetch": {"type": "boolean", "default": True, "description": "Download sources missing from the directory"},
    "allow_source_drift": {"type": "boolean", "default": False,
                           "description": "Accept a changed file only if every derived panel sequence still matches its definition hash"},
}, ("manifest", "directory")))
PIPELINE.section("panel", closed({
    "definition": dict(PATH, description="Panel definition TSV (rules, locators and SHA-256; no sequences)"),
    "census_elements": {"type": "array", "default": [], "items": ELEMENT_RULE_SCHEMA},
    "census_catalogue": {"type": "array", "default": [], "items": CATALOGUE_RULE_SCHEMA},
    "census_catalogue_checks": {"type": "array", "default": [], "items": CATALOGUE_CHECK_SCHEMA},
    "benchmark_templates": {"type": "array", "default": [], "items": TEMPLATE_RULE_SCHEMA},
}, ("definition",)))
PIPELINE.section("subset", closed({
    "parents": dict(COUNT, default=0, description="Dry run: number of canonical parents with the smallest parent IDs (0: whole panel)"),
    "query_ids": {"type": "array", "default": [], "uniqueItems": True, "items": {"type": "string", "pattern": "^[A-Za-z0-9_.:+-]+$"},
                  "description": "Dry run: explicit physical query IDs (alternative to parents)"},
}))
PIPELINE.section("controls", closed({
    "draws": dict(POSITIVE, description="Matched control draws per canonical parent"),
    "seed": dict(COUNT, description="Master seed of the per-parent draw seeds"),
    "k": {"const": 2, "default": 2, "description": "Preserved k-let size (uShuffle dinucleotide shuffle)"},
}, ("draws", "seed")))
PIPELINE.section("analysis_policy", closed({
    "record_populations": {"type": "array", "minItems": 1, "uniqueItems": True, "items": {"enum": ["all_records", "linear_genomic"]}},
    "minimum_query_length": closed({"strict": {"type": "integer", "minimum": 7}, "tolerant": {"type": "integer", "minimum": 7}},
                                   ("strict", "tolerant")),
    "primary_window": COUNT, "primary_minimum_length": POSITIVE,
    "windows": {"type": "array", "minItems": 1, "uniqueItems": True, "items": COUNT},
    "minimum_lengths": {"type": "array", "minItems": 1, "uniqueItems": True, "items": POSITIVE},
}, ("record_populations", "minimum_query_length", "primary_window", "primary_minimum_length", "windows", "minimum_lengths")))
PIPELINE.section("database", closed({
    "prefix": dict(PATH, description="BLAST database prefix"),
    "taxonomy_support": dict(PATH, description="Directory with the BLAST taxonomy files taxdb.bti and taxdb.btd"),
    "download_evidence": dict(PATH, description="Pinned download evidence directory (required for core_nt)"),
    "hash_all": {"type": "boolean", "default": True, "description": "Hash every database file (required for acceptance)"},
    "check_timeout_seconds": {"type": "number", "exclusiveMinimum": 0},
    "expected": closed({"sequences": POSITIVE, "bases": POSITIVE, "volumes": POSITIVE, "build_date": {"type": "string"}}, default={}),
    "probes": {"type": "array", "minItems": 1, "maxItems": 100, "items": closed({
        "accession_version": {"type": "string", "pattern": "^[A-Za-z][A-Za-z0-9_]*[.][0-9]+$"}, "expected_length": POSITIVE,
        "allow_absent": {"type": "boolean"}, "basis": {"type": "string", "minLength": 1}},
        ("accession_version", "expected_length", "allow_absent", "basis"))},
}, ("prefix", "taxonomy_support", "probes")))
PIPELINE.section("search", closed({
    "threads": POSITIVE, "batch_size": POSITIVE, "target_limit": POSITIVE, "evalue": {"type": "number", "exclusiveMinimum": 0},
    "threshold_mode": {"enum": ["fixed", "adaptive"], "default": "fixed"},
    "timeout_seconds": {"type": "number", "exclusiveMinimum": 0},
    "mode": {"enum": ["production", "working"], "default": "production",
             "description": "production requires the accepted database of the database stage"},
    "informativeness_limit": PROBABILITY, "testable_chance_limit": PROBABILITY,
    "resume": {"type": "boolean", "default": True, "description": "Continue from completed batches of an unfinished identical attempt"},
    "scoring": closed({"task": {"type": "string"}, "word_size": POSITIVE, "dust": {"type": "string"}, "soft_masking": {"type": "string"},
                       "strand": {"type": "string"}, "reward": {"type": "integer"}, "penalty": {"type": "integer"},
                       "gapopen": {"type": "integer"}, "gapextend": {"type": "integer"}},
                      ("task", "word_size", "dust", "soft_masking", "strand", "reward", "penalty", "gapopen", "gapextend")),
}, ("threads", "batch_size", "target_limit", "evalue", "informativeness_limit", "testable_chance_limit", "scoring")))
PIPELINE.section("ncbi", closed({
    "email": {"type": "string", "minLength": 3, "description": "Contact address sent with E-utilities requests (not a secret)"},
    "api_key": {"type": "string", "description": "E-utilities API key; use ${NCBI_API_KEY} (recorded as <redacted>)"},
}, ("email",)))
PIPELINE.section("metadata", closed({
    "batch_size": {"type": "integer", "minimum": 1, "maximum": 100, "default": 100},
    "require_complete_retrieval": {"type": "boolean", "default": True},
}))
PIPELINE.section("taxonomy", closed({
    "directory": dict(PATH, description="Local taxdump directory (nodes.dmp, names.dmp, merged.dmp, delnodes.dmp)"),
    "source_url": {"type": "string", "pattern": "^https://"}, "retrieved_utc": {"type": "string", "minLength": 1},
    "release_date": {"type": "string", "default": ""},
    "ranks": {"type": "array", "minItems": 1, "uniqueItems": True, "items": {"type": "string", "minLength": 1}},
}, ("directory", "source_url", "retrieved_utc", "ranks")))
PIPELINE.section("insdc", closed({
    "batch_size": {"type": "integer", "minimum": 1, "maximum": 200, "default": 100},
    "delay_seconds": {"type": "number", "minimum": 0, "default": 0.4},
    "require_complete_retrieval": {"type": "boolean", "default": True},
}))
PIPELINE.section("positional", closed({"diagnostic_window": dict(COUNT, description="Diagnostic window added to the policy windows")},
                                      ("diagnostic_window",)))
PIPELINE.section("record_cases", CASE_SCHEMA)
PIPELINE.section("census", CENSUS_SCHEMA)
PIPELINE.section("exact_scan", closed({
    "minimum_length": POSITIVE, "workers": POSITIVE,
    "volume_limit": dict(COUNT, default=0, description="Dry run: scan only the first N database volumes (0: all)"),
    "resume": {"type": "boolean", "default": True},
}, ("minimum_length", "workers")))
PIPELINE.section("scan_validation", closed({"require_acceptance": {"type": "boolean", "default": True}}))
PIPELINE.section("target_limit", closed({"arms": {"type": "array", "minItems": 2, "uniqueItems": True, "items": POSITIVE}}, ("arms",)))
PIPELINE.section("benchmark", closed({
    "source_archive": dict(PATH, description="NCBI BLAST+ source tarball of the pinned version (length-adjustment routine)"),
    "extra_lengths": {"type": "array", "uniqueItems": True, "items": POSITIVE},
    "sizes": {"type": "array", "minItems": 1, "items": closed({"label": {"type": "string", "pattern": "^[A-Za-z0-9_]+$"},
                                                               "factor": {"type": "number", "exclusiveMinimum": 0}}, ("label", "factor"))},
    "run_sizes": {"type": "array", "default": [], "uniqueItems": True, "items": {"type": "string"},
                  "description": "Dry run: labels of the sizes to run (empty: all)"},
    "shard_limit": dict(COUNT, default=0, description="Dry run: keep the first N shards of each fixture panel (0: all)"),
    "workers": POSITIVE, "maximum_workers": dict(POSITIVE, default=4), "wall_seconds": dict(POSITIVE, default=14400),
    "output_bytes": dict(POSITIVE, default=10 * 1024 ** 3), "fixture_timeout_seconds": dict(POSITIVE, default=3600),
}, ("source_archive", "extra_lengths", "sizes", "workers")))
PIPELINE.section("inference", closed({
    "draw_ladder": LADDER, "retention_ladder": dict(LADDER, minItems=2), "maximum_draws": POSITIVE, "confidence": PROBABILITY,
    "detection_half_width": {"type": "number", "exclusiveMinimum": 0, "maximum": 1},
    "breadth_absolute_half_width": {"type": "number", "exclusiveMinimum": 0},
    "breadth_relative_half_width": {"type": "number", "exclusiveMinimum": 0, "maximum": 1},
    "threads": POSITIVE, "batch_size": POSITIVE, "timeout_seconds": POSITIVE,
}, ("draw_ladder", "retention_ladder", "maximum_draws", "confidence", "detection_half_width", "breadth_absolute_half_width",
    "breadth_relative_half_width", "threads", "batch_size", "timeout_seconds")))
PIPELINE.section("validation", closed({"require_production": {"type": "boolean", "default": True}}))


def database_inputs(config: dict) -> list[str]:
    """Small identity files of the database configuration: BLAST taxonomy support, `.njs` metadata, download evidence."""
    section = config["database"]
    files = [str(Path(section["taxonomy_support"], name)) for name in ("taxdb.bti", "taxdb.btd")]
    if Path(section["prefix"] + ".njs").is_file():
        files.append(section["prefix"] + ".njs")
    if section.get("download_evidence"):
        files.append(section["download_evidence"])
    return files


def taxdump_inputs(config: dict) -> list[str]:
    return [str(Path(config["taxonomy"]["directory"], name)) for name in ("nodes.dmp", "names.dmp", "merged.dmp", "delnodes.dmp")]


S = "subset"
PIPELINE.stage(config_sections=("sources",), input_files=lambda config: [config["sources"]["manifest"]])(sources.sources)
PIPELINE.stage(depends_on=("sources",), config_sections=("panel", S), tools=("pdftotext",),
               input_files=lambda config: [config["panel"]["definition"]])(panel.panel)
PIPELINE.stage(depends_on=("panel",), config_sections=("controls", S))(controls.controls)
PIPELINE.stage(depends_on=("panel", "controls"), config_sections=("controls", "analysis_policy", S))(design.design)
PIPELINE.stage(config_sections=("database",), tools=("blastn", "blastdbcmd", "blastdbcheck"), input_files=database_inputs)(database.database)
PIPELINE.stage(depends_on=("design", "database"), config_sections=("search", S), tools=("blastn",))(search.search)
PIPELINE.stage(depends_on=("design", "search"), config_sections=(S,))(analysis.analysis)
PIPELINE.stage(depends_on=("analysis",), config_sections=("ncbi", "metadata", S))(metadata.metadata)
PIPELINE.stage(depends_on=("analysis", "metadata"), config_sections=(S,))(events.events)
PIPELINE.stage(depends_on=("events",), config_sections=("taxonomy", S), input_files=taxdump_inputs)(taxonomy.taxonomy)
PIPELINE.stage(depends_on=("taxonomy", "metadata"), config_sections=("ncbi", "insdc", S))(insdc.insdc)
PIPELINE.stage(depends_on=("taxonomy", "metadata", "insdc"), config_sections=("positional", S))(positional.positional)
PIPELINE.stage(depends_on=("taxonomy", "metadata", "search", "database", "panel"), config_sections=("record_cases", S),
               tools=("blastdbcmd",))(record_cases.record_cases)
PIPELINE.stage(depends_on=("taxonomy", "design", "panel", "database"), config_sections=("census", S), tools=("blastdbcmd",))(census.census)
PIPELINE.stage(depends_on=("design", "database"), config_sections=("exact_scan", S), tools=("blastdbcmd",))(exact_scan.exact_scan)
PIPELINE.stage(depends_on=("exact_scan", "search"), config_sections=("scan_validation", "exact_scan", S))(scan_validation.scan_validation)
PIPELINE.stage(depends_on=("search", "exact_scan"), config_sections=("exact_scan", S))(cap_review.cap_review)
PIPELINE.stage(depends_on=("design", "database", "search"), config_sections=("target_limit", "search", S),
               tools=("blastn",))(target_limit.target_limit)
PIPELINE.stage(depends_on=("design", "panel", "database"), config_sections=("benchmark", S),
               tools=("blastn", "makeblastdb", "blast_formatter", "blastdbcmd", "cc"),
               input_files=lambda config: [config["benchmark"]["source_archive"]])(benchmark.benchmark)
PIPELINE.stage(depends_on=("design", "taxonomy"), config_sections=("inference", S))(inference.inference)
PIPELINE.stage(depends_on=("taxonomy",), config_sections=("validation", S))(validation.validation)

from .reporting import register as _register_reporting  # noqa: E402

_register_reporting(PIPELINE)
