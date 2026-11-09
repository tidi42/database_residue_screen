import argparse
import json
import subprocess
import sys
from collections import Counter
from pathlib import Path

from . import __version__
from .acceptance import accept_database
from .analysis import analyse_search
from .artifacts import validate_outputs, write_preflight
from .contracts import ContractError, read_table, schema
from .controls import generate_controls
from .database import inspect_database
from .database_integrity import audit_integrity
from .design import check_design
from .example import create_example
from .metadata import annotate_events, fetch_metadata
from .panels import prepare_panel
from .precision import assess_precision, write_precision_spec
from .release import export_figure_inputs, validate_analysis
from .search import run_search, search_settings
from .taxonomy import annotate_taxonomy, pin_taxonomy


def parser() -> argparse.ArgumentParser:
    command = argparse.ArgumentParser(description="Database residue screening: one reference panel with matched controls against "
                                                  "a BLAST nucleotide database, from input contracts to validated result tables")
    command.add_argument("--version", action="version", version=__version__)
    subcommands = command.add_subparsers(dest="command", required=True)
    example = subcommands.add_parser("example", help="Generate synthetic software-test inputs; not a research study")
    example.add_argument("destination", type=Path)
    preflight = subcommands.add_parser("preflight", help="Validate study input contracts and write checksummed preflight outputs")
    preflight.add_argument("study", type=Path)
    preflight.add_argument("--output", type=Path)
    validate = subcommands.add_parser("validate-output", help="Check preflight structure, hashes, rows and source consistency")
    validate.add_argument("study", type=Path)
    validate.add_argument("directory", type=Path)
    database = subcommands.add_parser("inspect-database", help="Capture BLAST database identity evidence and conflicts, without modifying it")
    database.add_argument("prefix", type=Path)
    database.add_argument("--output", type=Path, required=True)
    database.add_argument("--blastdbcmd", default="blastdbcmd")
    database.add_argument("--study-root", type=Path)
    database.add_argument("--hash-all", action="store_true", help="Hash all database files; may read hundreds of GB")
    benchmark = subcommands.add_parser("benchmark", help="Prepare one hashed planted-sequence fixture set and its truth from a panel's "
                                       "real parents, or run it with the production settings and the search spaces of one calibration")
    benchmark.add_argument("action", choices=("prepare", "run"))
    benchmark.add_argument("fixtures", type=Path, help="Fixture-set directory (written by prepare, read by run)")
    benchmark.add_argument("--panel", type=Path)
    benchmark.add_argument("--extra-lengths", type=int, nargs="*", default=[],
                           help="Additional boundary-template lengths (panel D); each must be covered by the calibration of a run")
    benchmark.add_argument("--templates", type=Path, help="Declared templates (JSON written by the study panel stage)")
    benchmark.add_argument("--shard-limit", type=int, default=0, help="Keep the first N shards of each fixture panel (0 = all)")
    benchmark.add_argument("--calibration", type=Path, help="calibrate-benchmark output of one database size (run)")
    benchmark.add_argument("--output", type=Path, help="Run directory (run)")
    benchmark.add_argument("--workers", type=int, default=4)
    calibration = subcommands.add_parser("calibrate-benchmark", help="Source-verified NCBI calibration for explicit panel lengths and DB evidence")
    calibration.add_argument("archive", type=Path)
    calibration.add_argument("binaries", type=Path)
    calibration.add_argument("evidence", type=Path)
    calibration.add_argument("queries", type=Path)
    calibration.add_argument("--output", type=Path, required=True)
    table = subcommands.add_parser("validate-table", help="Validate one table's row contract; not a complete search-result certification")
    table.add_argument("schema_name", choices=("queries", "sources", "components", "subjects", "validation_cases", "search_status"))
    table.add_argument("table", type=Path)
    controls = subcommands.add_parser("controls", help="Generate and audit pinned dinucleotide-shuffled controls")
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
    panel = subcommands.add_parser("prepare-panel", help="Bind every panel query and its control draws; deduplicate only the physical search work")
    panel.add_argument("queries", type=Path)
    panel.add_argument("metadata", type=Path)
    panel.add_argument("controls", type=Path)
    panel.add_argument("--output", type=Path, required=True)
    panel.add_argument("--draws", type=int, required=True)
    panel.add_argument("--sample-parents", type=int, default=0, help="Optional stratified subset of canonical parents (0 = all)")
    panel.add_argument("--policy", type=Path, help="Versioned analysis policy; defaults do not claim exhaustive query recovery")
    precision_spec = subcommands.add_parser("precision-spec", help="Write a panel-bound control precision and retention specification")
    precision_spec.add_argument("panel", type=Path)
    precision_spec.add_argument("--output", type=Path, required=True)
    precision_spec.add_argument("--draw-ladder", type=int, nargs="+", required=True, help="Increasing draw prefixes per parent, e.g. 5 10 15")
    precision_spec.add_argument("--retention-ladder", type=int, nargs="+", required=True, help="Increasing target limits, e.g. 500 2000")
    precision_spec.add_argument("--maximum-draws", type=int, help="Draw budget per parent; defaults to the largest draw prefix")
    precision = subcommands.add_parser("assess-precision", help="Conditional null precision at the declared draw prefixes and Holm-adjusted "
                                       "reference-generator contrasts, preserving repeated draws and shared seeds")
    precision.add_argument("analysis", type=Path)
    precision.add_argument("specification", type=Path)
    precision.add_argument("--output", type=Path, required=True)
    analyse = subcommands.add_parser("analyse-search", help="Reconcile completed real/control searches on preserved parent and draw identities")
    analyse.add_argument("panel", type=Path)
    analyse.add_argument("search", type=Path)
    analyse.add_argument("--output", type=Path, required=True)
    insdc = subcommands.add_parser("fetch-insdc", help="Dated GenBank-header snapshot of the eligible linear-genomic records of the "
                                   "positional analysis, for the structured targeted/amplicon sensitivity (no titles)")
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
    census.add_argument("--specification", type=Path, required=True, help="Census families and class rules (JSON)")
    census.add_argument("--elements", type=Path, required=True, help="Flank elements derived at run time (JSON)")
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
    scan.add_argument("--resume", action="store_true", help="Continue an interrupted scan with an identical request")
    validate_analysis_command = subcommands.add_parser("validate-analysis", help="Recompute and validate result hashes, policy, parent weights, eligibility and source tables")
    validate_analysis_command.add_argument("analysis", type=Path)
    validate_analysis_command.add_argument("--require-production", action="store_true")
    figure_inputs = subcommands.add_parser("figure-inputs", help="Export one validated production result release for figures; no plotting")
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
    integrity = subcommands.add_parser("audit-database-integrity", help="Check every volume build date and archive exact-accession lookup probes")
    integrity.add_argument("database", type=Path)
    integrity.add_argument("--evidence", type=Path, required=True)
    integrity.add_argument("--output", type=Path, required=True)
    integrity.add_argument("--probe", action="append", default=[], help="Explicit exact accession.version probe; repeat as needed")
    fetch = subcommands.add_parser("fetch-sources", help="Download pinned panel sources and verify their SHA-256; a mismatching "
                                   "local file is never replaced")
    fetch.add_argument("manifest", type=Path, help="Source manifest (e.g. studies/<study>/sources.json)")
    fetch.add_argument("--directory", type=Path, required=True)
    build = subcommands.add_parser("build-panel", help="Derive every panel query from verified local sources and check it against "
                                   "the definition hashes (manifest: sources.json beside the definition)")
    build.add_argument("definition", type=Path)
    build.add_argument("--sources", type=Path, required=True, help="Directory with the pinned source files")
    build.add_argument("--output", type=Path, required=True)
    run_study_command = subcommands.add_parser("run-study", help="Run or resume a whole study from one TOML configuration into an "
                                               "immutable, provenance-recording run directory")
    run_study_command.add_argument("config", type=Path)
    run_study_command.add_argument("--output", type=Path, required=True, help="Run directory")
    stage_selection = run_study_command.add_mutually_exclusive_group()
    stage_selection.add_argument("--until", metavar="STAGE", help="Run STAGE and every stage it depends on")
    stage_selection.add_argument("--only", metavar="STAGE", nargs="+", action="extend", default=[],
                                 help="Run only these stages; their other dependencies must already be complete")
    run_study_command.add_argument("--working", action="store_true",
                                   help="Working mode: allow a dirty or untracked code tree, recorded in every stage record")
    run_study_command.add_argument("--rerun-changed", action="store_true",
                                   help="Give stages whose request changed a new attempt (earlier attempts are kept); selected "
                                        "descendants rerun because their upstream changes")
    study_status_command = subcommands.add_parser("study-status", help="Stage status, attempts, durations and request digests of a run")
    study_status_command.add_argument("run", type=Path)
    study_status_command.add_argument("--verify", action="store_true", help="Re-hash the outputs of every current attempt")
    study_status_command.add_argument("--verify-code", action="store_true",
                                      help="Check every current attempt against the stage code, tools and dependencies now "
                                           "installed (run before tagging a release)")
    return command


def main(arguments: list[str] | None = None) -> int:
    args = parser().parse_args(arguments)
    try:
        if args.command == "example":
            result = dict(study=str(create_example(args.destination)), scope="synthetic_example")
        elif args.command == "preflight":
            report = write_preflight(args.study, args.output) if args.output else check_design(args.study)
            result = dict(study_id=report["study_id"], scope=report["scope"], design_ready=report["design_ready"],
                          counts=report["counts"], blocker_counts=dict(Counter(issue["code"] for issue in report["issues"] if issue["severity"] == "blocker")),
                          warning_counts=dict(Counter(issue["code"] for issue in report["issues"] if issue["severity"] == "warning")))
            if args.output:
                result["report"] = str(args.output / "preflight.json")
            print(json.dumps(result, indent=2))
            return 0 if report["design_ready"] else 2
        elif args.command == "validate-output":
            result = validate_outputs(args.study, args.directory)
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
        elif args.command == "benchmark":
            from .current_benchmark import load_templates, prepare_benchmark, run_benchmark
            if args.action == "prepare":
                if args.panel is None:
                    raise ContractError("prepare requires --panel")
                report = prepare_benchmark(args.panel, args.fixtures, tuple(args.extra_lengths),
                                           load_templates(args.templates) if args.templates else (), args.shard_limit)
                result = {key: report[key] for key in ("original_parents", "boundary_templates", "required_lengths", "query_subject_pairs")}
            else:
                if args.calibration is None or args.output is None:
                    raise ContractError("run requires --calibration and --output")
                report = run_benchmark(args.fixtures, args.calibration, args.output, args.workers)
                result = {key: report[key] for key in ("status", "production_database", "completed_shards", "short_summary")}
        elif args.command == "calibrate-benchmark":
            from .benchmark_engine import calibrate_direct
            from .contracts import read_json
            from .sequences import read_fasta
            dimensions = read_json(args.evidence)["runtime_metadata"]
            lengths = tuple(sorted({len(s) for s in read_fasta(args.queries).values()}))
            result = calibrate_direct(args.archive, args.binaries, args.output, dimensions["base_count"], dimensions["sequence_count"], lengths)
        elif args.command == "prepare-panel":
            result = prepare_panel(args.queries, args.metadata, args.controls, args.output, args.draws, args.sample_parents, args.policy)
        elif args.command == "precision-spec":
            result = write_precision_spec(args.panel, args.output, args.draw_ladder, args.retention_ladder, args.maximum_draws)
        elif args.command == "assess-precision":
            result = assess_precision(args.analysis, args.specification, args.output)
        elif args.command == "analyse-search":
            result = analyse_search(args.panel, args.search, args.output)
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
            from .census import load_elements
            from .contracts import read_json as _read_json
            report = run_census(args.analysis, args.panel, args.database, str(Path(binary).resolve()), _read_json(args.specification),
                                load_elements(args.elements), args.catalogue, args.output, args.taxonomy)
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
        elif args.command == "audit-database-integrity":
            result = audit_integrity(args.database, args.evidence, args.output, tuple(args.probe))
            result = {key: value for key, value in result.items() if key not in {"volumes", "accession_probes"}}
        elif args.command == "fetch-sources":
            from .constructs import fetch_sources
            result = fetch_sources(args.manifest, args.directory)
        elif args.command == "build-panel":
            from .constructs import build_panel
            report = build_panel(args.definition, args.sources, args.output)
            result = {key: report[key] for key in ("queries", "canonical_parents", "by_group", "by_family", "by_verification", "files")}
        elif args.command == "run-study":
            from .study import run_study
            result = run_study(args.config, args.output, until=args.until, only=args.only, rerun_changed=args.rerun_changed,
                               mode="working" if args.working else "production", console=sys.stderr)
        elif args.command == "study-status":
            from .study import format_status, study_status
            report = study_status(args.run, verify=args.verify, verify_code=args.verify_code)
            print(format_status(report))
            failed = any(row.get("verified") is False for row in report["stages"]) or not report.get("code", {}).get("ok", True)
            return 1 if failed else 0
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
