"""Hash-bound a-priori testability for completed fixed-threshold searches.

Protocol-v2 classes depend only on query lengths, the accepted database
dimensions, the fixed blastn-short scoring contract and the production E-value
threshold, never on BLAST output. A completed search that predates runtime
calibration but used that contract can therefore be classified afterwards.
The migration is a separate, hashed artefact: it never edits the historical
request, batches or outputs and it reads no HSP of the historical search.
"""
import math
import shutil
import subprocess
from pathlib import Path

from .contracts import ContractError, read_json
from .design import digest
from .example import save_json
from .length_aware import calibrate_lengths, verify_calibration_row
from .search import object_digest
from .sequences import read_fasta

MIGRATION_KIND = "post_hoc_fixed_threshold_testability_migration_v1"
MIGRATION_FILE = "calibration_migration.json"
SCORING_CONTRACT = dict(task="blastn-short", word_size=7, dust="no", soft_masking="false", strand="both",
                        reward=1, penalty=-3, gapopen=5, gapextend=2)


def completed_search(search: Path) -> tuple[dict, dict]:
    if (search / "quarantine.json").exists():
        raise ContractError("Search is quarantined; do not derive testability for scientific results")
    request, complete = read_json(search / "request.json"), read_json(search / "complete.json")
    if complete["search_id"] != object_digest(request) or not complete["all_queries_completed"]:
        raise ContractError("Search request and completion disagree or the search is incomplete")
    return request, complete


def database_dimensions(search: Path, request: dict) -> dict:
    evidence = search / "database_evidence.json"
    if digest(evidence) != request["database_evidence_sha256"]:
        raise ContractError("Search database evidence changed since the search")
    runtime = read_json(evidence)["runtime_metadata"]
    if request.get("database_runtime_metadata") not in (None, runtime):
        raise ContractError("Recorded runtime database dimensions disagree with the search evidence")
    return runtime


def calibrate_existing_search(search: Path, queries: Path, destination: Path, executable: str = "blastn",
                              informative_limit: float = .05, testable_limit: float = .01) -> dict:
    request, complete = completed_search(search)
    if request.get("length_calibration"):
        raise ContractError("Search already carries its runtime calibration; analyse it without a migration")
    if digest(queries) != request["queries_sha256"]:
        raise ContractError("Query FASTA is not the one searched")
    settings = request["settings"]
    if any(settings.get(name) != value for name, value in SCORING_CONTRACT.items()):
        raise ContractError("Historical settings differ from the validated fixed scoring contract; recalibrate and validate")
    if not 0 < settings["evalue"] < math.inf:
        raise ContractError("Historical E-value threshold is not finite and positive")
    runtime = database_dimensions(search, request)
    binary = shutil.which(executable)
    if binary is None:
        raise ContractError("BLAST executable is unavailable: " + executable)
    binary = str(Path(binary).resolve())
    version = subprocess.run([binary, "-version"], capture_output=True, text=True, check=True, timeout=30).stdout.strip()
    historical = str(request.get("executable_version", "")).strip()
    if not historical or version.splitlines()[0] != historical.splitlines()[0]:
        raise ContractError("Installed BLAST version differs from the version that ran the search")
    lengths = sorted({len(sequence) for sequence in read_fasta(queries).values()})
    calibration = calibrate_lengths(binary, lengths, runtime, settings, "fixed", None, informative_limit, testable_limit)
    record = dict(kind=MIGRATION_KIND, source_search=str(search.resolve()), source_request_sha256=digest(search / "request.json"),
                  source_completion_sha256=digest(search / "complete.json"), source_search_id=complete["search_id"],
                  queries_sha256=request["queries_sha256"], database_evidence_sha256=request["database_evidence_sha256"],
                  database_runtime_metadata=runtime, historical_executable_version=historical, executable=binary,
                  executable_version=version, executable_sha256=digest(Path(binary)),
                  calibration_code_sha256=digest(Path(__file__).with_name("length_aware.py")),
                  migration_code_sha256=digest(Path(__file__)), length_calibration=calibration,
                  statement="A-priori classes from query lengths, recorded database dimensions and the historical fixed "
                            "E threshold; no historical BLAST output was read and no search was repeated.")
    destination.mkdir(parents=True, exist_ok=False)
    save_json(destination / MIGRATION_FILE, record)
    return record


def validate_migration(path: Path, search: Path, request: dict, complete: dict) -> dict:
    migration = read_json(path)
    if migration.get("kind") != MIGRATION_KIND:
        raise ContractError("Not a fixed-threshold testability migration")
    if request.get("length_calibration"):
        raise ContractError("Search already carries its runtime calibration; a migration is not allowed")
    if (migration["source_request_sha256"] != digest(search / "request.json") or migration["source_search_id"] != complete["search_id"]
            or migration["source_completion_sha256"] != digest(search / "complete.json")
            or migration["queries_sha256"] != request["queries_sha256"]
            or migration["database_evidence_sha256"] != request["database_evidence_sha256"]):
        raise ContractError("Testability migration belongs to a different search")
    runtime = database_dimensions(search, request)
    calibration = migration["length_calibration"]
    if (migration["database_runtime_metadata"] != runtime or calibration["threshold_mode"] != "fixed"
            or calibration["maximum_evalue"] is not None):
        raise ContractError("Testability migration does not describe the fixed historical search")
    for row in calibration["lengths"]:
        if (row["evalue"] != request["settings"]["evalue"] or row["database_bases"] != runtime["base_count"]
                or row["database_sequences"] != runtime["sequence_count"]
                or request["settings"].get("searchsp") not in (None, row["effective_search_space"])):
            raise ContractError("Migrated calibration disagrees with the historical threshold or database")
        verify_calibration_row(row)
    return migration
