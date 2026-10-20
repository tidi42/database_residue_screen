import argparse
import json
import subprocess
import sys
from collections import Counter
from pathlib import Path

from . import __version__
from .acceptance import accept_database
from .artifacts import freeze_design, validate_outputs, write_preflight
from .benchmark import run_benchmark
from .contracts import ContractError, read_table, schema
from .context import context_evidence
from .controls import generate_controls
from .database import inspect_database
from .database_integrity import audit_integrity
from .design import check_design
from .example import create_example
from .legacy import import_legacy
from .metadata import annotate_events, fetch_metadata
from .panels import prepare_panel, select_original_queries
from .precision import assess_precision, write_pilot_spec
from .preparation import prepare_reference_inputs
from .reanalysis import analyse_search, reanalyse_legacy
from .release import export_figure_inputs, validate_analysis
from .search import run_search, search_settings
from .taxonomy import annotate_taxonomy, pin_taxonomy


def parser() -> argparse.ArgumentParser:
    command = argparse.ArgumentParser(description="Database residue screening: pre-freeze input and output validation")
    command.add_argument("--version", action="version", version=__version__)
    subcommands = command.add_subparsers(dest="command", required=True)
    example = subcommands.add_parser("example", help="Generate synthetic software-test inputs; not a research study")
    example.add_argument("destination", type=Path)
    preflight = subcommands.add_parser("preflight", help="Validate input contracts and emit a checksummed preflight bundle")
    preflight.add_argument("study", type=Path)
    preflight.add_argument("--output", type=Path)
    validate = subcommands.add_parser("validate-output", help="Check preflight structure, hashes, rows and source consistency")
    validate.add_argument("study", type=Path)
    validate.add_argument("directory", type=Path)
    freeze = subcommands.add_parser("freeze", help="Write a G0-only design lock if all required design evidence is complete")
    freeze.add_argument("study", type=Path)
    freeze.add_argument("--output", type=Path, required=True)
    legacy = subcommands.add_parser("import-legacy", help="Read-only import of the historical core_nt analysis archive")
    legacy.add_argument("root", type=Path)
    legacy.add_argument("--output", type=Path, required=True)
    database = subcommands.add_parser("inspect-database", help="Capture BLAST database identity evidence and conflicts, without modifying it")
    database.add_argument("prefix", type=Path)
    database.add_argument("--output", type=Path, required=True)
    database.add_argument("--blastdbcmd", default="blastdbcmd")
    database.add_argument("--study-root", type=Path)
    database.add_argument("--hash-all", action="store_true", help="Hash all database files; may read hundreds of GB")
    current = subcommands.add_parser("current-benchmark", help="Run or prepare frozen current-engine fixtures for any supplied panel")
    current.add_argument("action", choices=("prepare", "run"))
    current.add_argument("destination", type=Path)
    current.add_argument("--panel", type=Path)
    current.add_argument("--calibration", type=Path)
    current.add_argument("--workers", type=int, default=4)
    current.add_argument("--extension-panel", type=Path, help="Panel whose real parents receive the separately seeded B01 extension shards")
    current.add_argument("--extra-lengths", type=int, nargs="*", default=[], help="Extension template lengths (panel D)")
    calibration = subcommands.add_parser("calibrate-benchmark", help="Source-verified NCBI calibration for explicit panel lengths and DB evidence")
    calibration.add_argument("archive", type=Path)
    calibration.add_argument("binaries", type=Path)
    calibration.add_argument("evidence", type=Path)
    calibration.add_argument("queries", type=Path)
    calibration.add_argument("--output", type=Path, required=True)
    table = subcommands.add_parser("validate-table", help="Validate one table's row contract; not a complete search-result certification")
    table.add_argument("schema_name", choices=("queries", "sources", "components", "subjects", "validation_cases", "search_status"))
    table.add_argument("table", type=Path)
    controls = subcommands.add_parser("controls", help="Generate and audit pinned dinucleotide controls; not a research freeze")
    controls.add_argument("queries", type=Path)
    controls.add_argument("--output", type=Path, required=True)
    controls.add_argument("--draws", type=int, required=True)
    controls.add_argument("--seed", type=int, default=42)
    search = subcommands.add_parser("search", help="Run or resume checksummed BLAST batches with explicit query completion")
    search.add_argument("queries", type=Path)
    search.add_argument("database", type=Path)
    search.add_argument("--evidence", type=Path, required=True)
    search.add_argument("--output", type=Path, required=True)
    search.add_argument("--threads", type=int, default=8)
    search.add_argument("--batch-size", type=int, default=256)
    search.add_argument("--targets", type=int, default=500)
    search.add_argument("--evalue", type=float, default=10.0)
    search.add_argument("--search-space", type=int)
    search.add_argument("--threshold-mode", choices=("adaptive", "fixed"), default="fixed",
                        help="fixed: production default; adaptive: diagnostic chance-saturated output only")
    search.add_argument("--informativeness-limit", type=float, default=.05)
    search.add_argument("--testable-chance-limit", type=float, default=.01)
    search.add_argument("--maximum-evalue", type=float, help="Optional random-HSP budget; may censor short lengths")
    search.add_argument("--study", type=Path, help="Validate and hash the declared study settings")
    search.add_argument("--blastn", default="blastn")
    search.add_argument("--timeout-seconds", type=float)
    search.add_argument("--mode", choices=("production", "working"), default="production")
    search.add_argument("--acceptance", type=Path, help="Accepted database artifact; required for production")
    search.add_argument("--reuse-search", type=Path, help="Reuse completed batches with identical science, retaining original execution provenance")
    acceptance = subcommands.add_parser("accept-database", help="Verify full hashes, volumes, integrity and exact-ID probes before production; can take hours")
    acceptance.add_argument("database", type=Path)
    acceptance.add_argument("--evidence", type=Path, required=True)
    acceptance.add_argument("--probes", type=Path, required=True)
    acceptance.add_argument("--taxonomy-support", type=Path, required=True, help="Directory containing verified taxdb.bti and taxdb.btd for blastdbcheck")
    acceptance.add_argument("--download-evidence", type=Path)
    acceptance.add_argument("--blastn", default="blastn")
    acceptance.add_argument("--blastdbcmd", default="blastdbcmd")
    acceptance.add_argument("--blastdbcheck", default="blastdbcheck")
    acceptance.add_argument("--timeout-seconds", type=float)
    acceptance.add_argument("--output", type=Path, required=True)
    benchmark = subcommands.add_parser("benchmark", help="Measure synthetic strict/tolerant detection at production-scale significance")
    benchmark.add_argument("--output", type=Path, required=True)
    benchmark.add_argument("--production-bases", type=int, required=True)
    benchmark.add_argument("--seed", type=int, default=20260915)
    panel = subcommands.add_parser("prepare-panel", help="Record real/control identities and deduplicate only the physical search work")
    panel.add_argument("queries", type=Path)
    panel.add_argument("metadata", type=Path)
    panel.add_argument("controls", type=Path)
    panel.add_argument("--output", type=Path, required=True)
    panel.add_argument("--draws", type=int, required=True)
    panel.add_argument("--pilot-parents", type=int, default=0)
    panel.add_argument("--selection", choices=("original", "all"), default="original")
    panel.add_argument("--policy", type=Path, help="Versioned analysis policy; defaults do not claim exhaustive query recovery")
    originals = subcommands.add_parser("select-originals", help="Export original author queries, excluding added reference fragments, with hashes")
    originals.add_argument("queries", type=Path)
    originals.add_argument("metadata", type=Path)
    originals.add_argument("--output", type=Path, required=True)
    pilot = subcommands.add_parser("pilot-spec", help="Prepare a panel-bound fixed-stage control precision and retention specification without searching")
    pilot.add_argument("panel", type=Path)
    pilot.add_argument("--output", type=Path, required=True)
    precision = subcommands.add_parser("assess-precision", help="Report conditional fixed-stage null uncertainty, preserving repeated draws and shared seeds")
    precision.add_argument("analysis", type=Path)
    precision.add_argument("specification", type=Path)
    precision.add_argument("--output", type=Path, required=True)
    preparation = subcommands.add_parser("prepare-reference-inputs", help="Verify original author panels and prepare matched controls and pilot inputs without database access")
    preparation.add_argument("root", type=Path)
    preparation.add_argument("imported_inputs", type=Path)
    preparation.add_argument("--draws", type=int, default=100)
    preparation.add_argument("--seed", type=int, default=42)
    preparation.add_argument("--output", type=Path, required=True)
    retrospective = subcommands.add_parser("reanalyse-retained", help="Reclassify historical real HSPs; never infer successful zeros or a valid null")
    retrospective.add_argument("root", type=Path)
    retrospective.add_argument("queries", type=Path)
    retrospective.add_argument("metadata", type=Path)
    retrospective.add_argument("--output", type=Path, required=True)
    analyse = subcommands.add_parser("analyse-search", help="Reconcile completed real/control searches on preserved parent and draw identities")
    analyse.add_argument("panel", type=Path)
    analyse.add_argument("search", type=Path)
    analyse.add_argument("--output", type=Path, required=True)
    analyse.add_argument("--calibration-migration", type=Path,
                         help="calibration_migration.json from calibrate-search; only for fixed-threshold searches without runtime calibration")
    insdc = subcommands.add_parser("fetch-insdc", help="Dated GenBank-header snapshot of primary eligible linear-genomic records for the "
                                   "structured targeted/amplicon sensitivity (no titles)")
    insdc.add_argument("analysis", type=Path)
    insdc.add_argument("--metadata", type=Path, required=True, help="Archived ESummary directory (tech field)")
    insdc.add_argument("--output", type=Path, required=True)
    positional = subcommands.add_parser("positional-sensitivity", help="Terminal-event robustness: one event per accession/BioProject, "
                                        "25-nt diagnostic window, optional INSDC targeted exclusion")
    positional.add_argument("analysis", type=Path)
    positional.add_argument("--metadata", type=Path, required=True, help="Archived ESummary directory (projectid field)")
    positional.add_argument("--insdc", type=Path, help="fetch-insdc output directory")
    positional.add_argument("--output", type=Path, required=True, help="TSV output path")
    caps = subcommands.add_parser("cap-review", help="Weakest retained subject score versus minimum qualifying scores for every "
                                  "cap-saturated query of a completed search")
    caps.add_argument("search", type=Path)
    caps.add_argument("--exact-scan", type=Path, help="exact-scan output for cap-free strict counts")
    caps.add_argument("--output", type=Path, required=True, help="TSV output path")
    census = subcommands.add_parser("flank-census", help="Descriptive flank classes of retained strict hits for declared families; "
                                    "integrity-checked local retrieval, no new search")
    census.add_argument("analysis", type=Path)
    census.add_argument("--panel", type=Path, required=True)
    census.add_argument("--database", type=Path, required=True)
    census.add_argument("--catalogue", type=Path, required=True, help="TSV set/name/sequence of catalogued barcodes")
    census.add_argument("--taxonomy", type=Path, help="Pinned taxdump snapshot for the phage Mu lineage class")
    census.add_argument("--blastdbcmd", default="blastdbcmd")
    census.add_argument("--output", type=Path, required=True)
    scan = subcommands.add_parser("exact-scan", help="Exhaustive exact occurrences of queries on both strands in every database volume; "
                                  "independent of BLAST seeds, E-values and caps (needs pyahocorasick)")
    scan.add_argument("queries", type=Path)
    scan.add_argument("database", type=Path)
    scan.add_argument("--evidence", type=Path, required=True, help="Database evidence with runtime volume names and totals")
    scan.add_argument("--output", type=Path, required=True)
    scan.add_argument("--blastdbcmd", default="blastdbcmd")
    scan.add_argument("--workers", type=int, default=8)
    scan.add_argument("--minimum-length", type=int, default=23)
    scan.add_argument("--resume", action="store_true", help="Continue an interrupted scan with an identical frozen request")
    bundle = subcommands.add_parser("assemble-bundle", help="One analysable search from a migrated historical search and new calibrated "
                                    "searches; verified projections, nothing re-run")
    bundle.add_argument("panel", type=Path, help="Bundle panel whose search_queries.fasta all sources together cover")
    bundle.add_argument("--historical", type=Path, required=True, help="Completed historical fixed-threshold search")
    bundle.add_argument("--migration", type=Path, required=True, help="calibration_migration.json of the historical search")
    bundle.add_argument("--new-search", type=Path, nargs="+", required=True, help="Completed runtime-calibrated searches (bridge first)")
    bundle.add_argument("--protocol", type=Path, help="Frozen protocol file to bind by hash")
    bundle.add_argument("--output", type=Path, required=True)
    migration = subcommands.add_parser("calibrate-search", help="A-priori testability for a completed fixed-threshold search that predates "
                                       "runtime calibration; reads no HSPs and repeats no search")
    migration.add_argument("search", type=Path)
    migration.add_argument("--queries", type=Path, required=True, help="The searched FASTA (panel search_queries.fasta)")
    migration.add_argument("--output", type=Path, required=True)
    migration.add_argument("--blastn", default="blastn", help="Same BLAST+ version as the historical search")
    migration.add_argument("--informativeness-limit", type=float, default=.05)
    migration.add_argument("--testable-chance-limit", type=float, default=.01)
    validate_analysis_command = subcommands.add_parser("validate-analysis", help="Recompute and validate result hashes, policy, parent weights, eligibility and source tables")
    validate_analysis_command.add_argument("analysis", type=Path)
    validate_analysis_command.add_argument("--require-production", action="store_true")
    figure_inputs = subcommands.add_parser("figure-inputs", help="Export one validated production result release for figures; no historical fallback or plotting")
    figure_inputs.add_argument("analysis", type=Path)
    figure_inputs.add_argument("--output", type=Path, required=True)
    metadata = subcommands.add_parser("fetch-metadata", help="Archive NCBI exact-version record summaries for an event release")
    metadata.add_argument("analysis", type=Path)
    metadata.add_argument("--output", type=Path, required=True)
    metadata.add_argument("--reuse-metadata", type=Path, help="Verified exact-version metadata archive; fetch only uncached matching records")
    annotation = subcommands.add_parser("annotate-events", help="Join exact-version metadata and regenerate endpoint-specific event tables")
    annotation.add_argument("analysis", type=Path)
    annotation.add_argument("metadata", type=Path)
    annotation.add_argument("--output", type=Path, required=True)
    taxonomy = subcommands.add_parser("pin-taxonomy", help="Archive and hash a local NCBI taxdump with explicit source and retrieval date")
    taxonomy.add_argument("directory", type=Path)
    taxonomy.add_argument("--source-url", required=True)
    taxonomy.add_argument("--retrieved-utc", required=True)
    taxonomy.add_argument("--release-date", default="")
    taxonomy.add_argument("--output", type=Path, required=True)
    lineage = subcommands.add_parser("annotate-taxonomy", help="Reconcile dated taxids and recompute accession-level summaries and eligible outcomes")
    lineage.add_argument("analysis", type=Path)
    lineage.add_argument("snapshot", type=Path)
    lineage.add_argument("--ranks", nargs="+", default=["domain", "superkingdom", "genus"])
    lineage.add_argument("--output", type=Path, required=True)
    context = subcommands.add_parser("context-evidence", help="Export segment coverage and a stratified sequence/flank case series")
    context.add_argument("analysis", type=Path)
    context.add_argument("inputs", type=Path)
    context.add_argument("database", type=Path)
    context.add_argument("--evidence", type=Path, required=True)
    context.add_argument("--output", type=Path, required=True)
    integrity = subcommands.add_parser("audit-database-integrity", help="Check every volume build date and archive exact-accession lookup probes")
    integrity.add_argument("database", type=Path)
    integrity.add_argument("--evidence", type=Path, required=True)
    integrity.add_argument("--output", type=Path, required=True)
    integrity.add_argument("--probe", action="append", default=[], help="Explicit exact accession.version probe; repeat as needed")
    return command


def main(arguments: list[str] | None = None) -> int:
    args = parser().parse_args(arguments)
    try:
        if args.command == "example":
            result = dict(study=str(create_example(args.destination)), scope="synthetic_example")
        elif args.command == "preflight":
            report = write_preflight(args.study, args.output) if args.output else check_design(args.study)
            result = dict(study_id=report["study_id"], scope=report["scope"], ready_for_freeze=report["ready_for_freeze"],
                          counts=report["counts"], blocker_counts=dict(Counter(issue["code"] for issue in report["issues"] if issue["severity"] == "blocker")),
                          warning_counts=dict(Counter(issue["code"] for issue in report["issues"] if issue["severity"] == "warning")))
            if args.output:
                result["report"] = str(args.output / "preflight.json")
            print(json.dumps(result, indent=2))
            return 0 if report["ready_for_freeze"] else 2
        elif args.command == "validate-output":
            result = validate_outputs(args.study, args.directory)
        elif args.command == "freeze":
            result = freeze_design(args.study, args.output)
        elif args.command == "import-legacy":
            study, report = import_legacy(args.root, args.output)
            result = dict(study=str(study), legacy_report=str(args.output / "legacy_inventory.json"),
                          original_queries=report["original_physical_queries"], missing_files=report["missing_files"],
                          derived_component_queries=len(report["derived_standalone_queries"]), research_design_ready=False)
        elif args.command == "inspect-database":
            report = inspect_database(args.prefix, args.output, args.blastdbcmd, args.hash_all, args.study_root)
            result = dict(runtime_metadata=report["runtime_metadata"], metadata_conflicts=report["metadata_conflicts"],
                          missing_volume_files=report["missing_required_volume_files"], total_bytes=report["total_bytes"],
                          all_database_files_hashed=report["all_database_files_hashed"], report=str(args.output / "database_identity.json"))
        elif args.command == "controls":
            result = generate_controls(args.queries, args.output, args.draws, args.seed)
        elif args.command == "search":
            result = run_search(args.queries, args.database, args.evidence, args.output,
                                search_settings(args.targets, args.evalue, args.search_space), args.threads, args.batch_size,
                                args.blastn, args.timeout_seconds, args.mode, args.acceptance, args.reuse_search,
                                args.threshold_mode, args.maximum_evalue, args.study, args.informativeness_limit, args.testable_chance_limit)
        elif args.command == "accept-database":
            result = accept_database(args.database, args.evidence, args.probes, args.output, args.blastn, args.download_evidence,
                                     args.timeout_seconds, args.taxonomy_support, args.blastdbcmd, args.blastdbcheck)
            result = {key: value for key, value in result.items() if key not in {"inventory", "download"}}
        elif args.command == "current-benchmark":
            from .current_benchmark import prepare_benchmark, run_benchmark as run_current
            if args.action == "prepare":
                if args.panel is None or args.calibration is None:
                    raise ContractError("prepare requires --panel and --calibration")
                result = prepare_benchmark(args.panel, args.calibration, args.destination, args.extension_panel, tuple(args.extra_lengths))
            else:
                result = run_current(args.destination, args.workers)
        elif args.command == "calibrate-benchmark":
            from .benchmark_engine import calibrate_direct
            from .contracts import read_json
            from .sequences import read_fasta
            dimensions = read_json(args.evidence)["runtime_metadata"]
            lengths = tuple(sorted({len(s) for s in read_fasta(args.queries).values()}))
            result = calibrate_direct(args.archive, args.binaries, args.output, dimensions["base_count"], dimensions["sequence_count"], lengths)
        elif args.command == "benchmark":
            result = run_benchmark(args.output, args.production_bases, args.seed)
        elif args.command == "prepare-panel":
            result = prepare_panel(args.queries, args.metadata, args.controls, args.output, args.draws, args.pilot_parents, args.selection, args.policy)
        elif args.command == "select-originals":
            result = select_original_queries(args.queries, args.metadata, args.output)
            result = {key: value for key, value in result.items() if key != "queries"}
        elif args.command == "pilot-spec":
            result = write_pilot_spec(args.panel, args.output)
        elif args.command == "assess-precision":
            result = assess_precision(args.analysis, args.specification, args.output)
        elif args.command == "prepare-reference-inputs":
            result = prepare_reference_inputs(args.root, args.imported_inputs, args.output, args.draws, args.seed)
        elif args.command == "reanalyse-retained":
            result = reanalyse_legacy(args.root, args.queries, args.metadata, args.output)
        elif args.command == "analyse-search":
            result = analyse_search(args.panel, args.search, args.output, args.calibration_migration)
        elif args.command == "fetch-insdc":
            from .insdc import fetch_insdc_annotation, technologies_from_metadata
            from .positional import primary_positional_events
            events, _ = primary_positional_events(args.analysis)
            report = fetch_insdc_annotation(sorted({row["accession_version"] for row in events}),
                                            technologies_from_metadata(args.metadata), args.output)
            result = {key: report[key] for key in ("requested", "status_counts", "targeted_counts")}
        elif args.command == "positional-sensitivity":
            from .contracts import write_table as _write_table
            from .insdc import read_annotation
            from .positional import SENSITIVITY_SCHEMA, positional_sensitivity, projects_from_metadata
            rows = positional_sensitivity(args.analysis, projects_from_metadata(args.metadata),
                                          read_annotation(args.insdc) if args.insdc else None)
            _write_table(args.output, rows, SENSITIVITY_SCHEMA)
            result = dict(rows=len(rows), output=str(args.output))
        elif args.command == "cap-review":
            from .cap_review import REVIEW_SCHEMA, review_saturated
            from .contracts import write_table as _write_table
            exact = None
            if args.exact_scan:
                from .exact_scan import read_occurrences
                exact = {}
                for row in read_occurrences(args.exact_scan):
                    exact.setdefault(row["query_id"], set()).add(row["accession_version"])
            rows = review_saturated(args.search, exact)
            _write_table(args.output, rows, REVIEW_SCHEMA)
            result = dict(saturated_queries=len(rows), statuses=dict(Counter((r["strict_status"], r["tolerant_status"]) for r in rows).most_common()))
        elif args.command == "flank-census":
            import shutil as _shutil
            from .census import run_census
            binary = _shutil.which(args.blastdbcmd)
            if binary is None:
                raise ContractError("blastdbcmd is unavailable: " + args.blastdbcmd)
            report = run_census(args.analysis, args.panel, args.database, str(Path(binary).resolve()), args.catalogue, args.output,
                                args.taxonomy)
            result = {key: report[key] for key in ("windows", "status_counts", "integrity_complete")}
        elif args.command == "exact-scan":
            import shutil as _shutil
            from .exact_scan import run_exact_scan
            binary = _shutil.which(args.blastdbcmd)
            if binary is None:
                raise ContractError("blastdbcmd is unavailable: " + args.blastdbcmd)
            report = run_exact_scan(args.queries, args.database, args.evidence, args.output, str(Path(binary).resolve()),
                                    args.workers, args.minimum_length, args.resume)
            result = {key: report[key] for key in ("scanned_queries", "records", "bases", "totals_match_database", "occurrences", "reused_volumes")}
        elif args.command == "assemble-bundle":
            from .bundle import assemble_bundle
            report = assemble_bundle(args.panel, args.historical, args.migration, args.new_search, args.output, args.protocol)
            result = {key: report[key] for key in ("search_id", "submitted_queries", "completed_queries", "success_zero", "hsp_rows",
                                                   "target_limit_reached")}
        elif args.command == "calibrate-search":
            from .testability_migration import calibrate_existing_search
            record = calibrate_existing_search(args.search, args.queries, args.output, args.blastn,
                                               args.informativeness_limit, args.testable_chance_limit)
            result = dict(kind=record["kind"], source_search_id=record["source_search_id"],
                          lengths=[dict(query_length=row["query_length"], strict=row["strict_testability"], tolerant=row["tolerant_testability"])
                                   for row in record["length_calibration"]["lengths"]],
                          migration=str(args.output / "calibration_migration.json"))
        elif args.command == "validate-analysis":
            result = validate_analysis(args.analysis, args.require_production)
        elif args.command == "figure-inputs":
            result = export_figure_inputs(args.analysis, args.output)
        elif args.command == "fetch-metadata":
            result = fetch_metadata(args.analysis, args.output, reuse_metadata=args.reuse_metadata)
            result = {key: value for key, value in result.items() if key != "batches"}
        elif args.command == "annotate-events":
            result = annotate_events(args.analysis, args.metadata, args.output)
        elif args.command == "pin-taxonomy":
            result = pin_taxonomy(args.directory, args.output, args.source_url, args.retrieved_utc, args.release_date)
        elif args.command == "annotate-taxonomy":
            result = annotate_taxonomy(args.analysis, args.snapshot, args.output, tuple(args.ranks))
        elif args.command == "context-evidence":
            result = context_evidence(args.analysis, args.inputs, args.database, args.evidence, args.output)
            result = {key: value for key, value in result.items() if key != "files"}
        elif args.command == "audit-database-integrity":
            result = audit_integrity(args.database, args.evidence, args.output, tuple(args.probe))
            result = {key: value for key, value in result.items() if key not in {"volumes", "accession_probes"}}
        else:
            rows = read_table(args.table, schema(args.schema_name))
            result = dict(valid_row_contract=True, rows=len(rows), scope="Row validation only; not full execution or result validation")
        print(json.dumps(result, indent=2))
        return 0
    except (ContractError, OSError, ValueError, subprocess.SubprocessError) as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
