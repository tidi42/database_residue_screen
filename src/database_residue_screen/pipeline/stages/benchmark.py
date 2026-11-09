"""Stage `benchmark`: one hashed fixture set from the search panel, calibrated and run at each configured database size."""
import json
from pathlib import Path

from ...benchmark_engine import calibrate_direct
from ...contracts import ContractError, read_json
from ...current_benchmark import load_templates, prepare_benchmark, run_benchmark
from ...design import digest
from ...example import save_json

BLAST_TOOLS = ("blastn", "makeblastdb", "blast_formatter", "blastdbcmd")


def benchmark(context):
    """Fixtures and truth hashed before any search; per database size a source-verified calibration and a run."""
    section = context.config["benchmark"]
    fixtures = context.output / "fixtures"
    specification = prepare_benchmark(context.upstream["design"].output / "search_panel", fixtures, tuple(section["extra_lengths"]),
                                      load_templates(context.upstream["panel"].output / "benchmark_templates.json"), section["shard_limit"],
                                      section["maximum_workers"], section["wall_seconds"], section["output_bytes"],
                                      section["fixture_timeout_seconds"])
    binaries = Path(context.tool("blastn")).parent
    if any(Path(context.tool(name)).parent != binaries for name in BLAST_TOOLS):
        raise ContractError("blastn, makeblastdb, blast_formatter and blastdbcmd must be configured from one BLAST+ directory")
    labels = [size["label"] for size in section["sizes"]]
    unknown = sorted(set(section["run_sizes"]) - set(labels))
    if unknown or len(set(labels)) != len(labels):
        raise ContractError("Benchmark sizes need unique labels; run_sizes names unknown sizes: " + ", ".join(unknown))
    selected = [size for size in section["sizes"] if not section["run_sizes"] or size["label"] in section["run_sizes"]]
    evidence = context.upstream["database"].output / "inventory" / "database_identity.json"
    runtime = read_json(evidence)["runtime_metadata"]
    runs = {}
    for size in selected:
        label, factor = size["label"], size["factor"]
        bases, sequences = int(round(runtime["base_count"] * factor)), int(round(runtime["sequence_count"] * factor))
        (context.output / "dimensions").mkdir(exist_ok=True)
        save_json(context.output / "dimensions" / f"{label}.json",
                  dict(runtime_metadata=dict(runtime, base_count=bases, sequence_count=sequences), scale_factor=factor,
                       source_evidence_sha256=digest(evidence),
                       note="Accepted database dimensions" if factor == 1 else "Emulated database dimensions for benchmark search spaces only"))
        calibration = context.output / "calibration" / label
        calibrate_direct(Path(section["source_archive"]), binaries, calibration, bases, sequences, tuple(specification["required_lengths"]),
                         context.tool("cc"))
        report = run_benchmark(fixtures, calibration, context.output / "runs" / label, section["workers"])
        runs[label] = dict(factor=factor, bases=bases, sequences=sequences, completed_shards=report["completed_shards"],
                           benchmark_sha256=digest(context.output / "runs" / label / "benchmark.json"))
        context.log.info("benchmark size %s: %d shards complete", label, report["completed_shards"])
    summary = dict(kind="benchmark_sizes_v1", fixtures_specification_sha256=digest(fixtures / "benchmark_specification.json"),
                   fixture_shards=len(specification["shards"]), shards_available=specification["shards_available"],
                   shard_limit=specification["shard_limit"], query_subject_pairs=specification["query_subject_pairs"],
                   required_lengths=specification["required_lengths"], extra_lengths=specification["extra_lengths"],
                   templates=specification["templates"], sizes=runs, configured_sizes=labels, workers=section["workers"])
    (context.output / "benchmark.json").write_text(json.dumps(summary, indent=1) + "\n", encoding="utf-8")
    return dict(fixture_shards=summary["fixture_shards"], shards_available=summary["shards_available"], sizes=sorted(runs),
                query_subject_pairs=summary["query_subject_pairs"])
