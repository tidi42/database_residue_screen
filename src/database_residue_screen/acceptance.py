import os
import re
import shutil
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path

from .contracts import ContractError, read_json
from .database import metadata_conflicts, parse_database_info, verify_inventory
from .database_integrity import assert_single_build, volume_evidence
from .design import digest
from .example import save_json


PROBE_SCHEMA = {"type": "object", "additionalProperties": False, "required": ["version", "probes"],
                "properties": {"version": {"const": "1.0.0"}, "probes": {"type": "array", "minItems": 1, "maxItems": 100,
                    "items": {"type": "object", "additionalProperties": False,
                              "required": ["accession_version", "expected_length", "allow_absent", "basis"],
                              "properties": {"accession_version": {"type": "string", "pattern": "^[A-Za-z][A-Za-z0-9_]*[.][0-9]+$"},
                                             "expected_length": {"type": "integer", "minimum": 1}, "allow_absent": {"type": "boolean"},
                                             "basis": {"type": "string", "minLength": 1}}}}}}


def executable_record(executable: str) -> dict:
    located = shutil.which(executable)
    if located is None:
        raise ContractError("Required executable is unavailable: " + executable)
    path = Path(located).resolve()
    result = subprocess.run([str(path), "-version"], capture_output=True, text=True, check=True, timeout=30)
    return dict(path=str(path), sha256=digest(path), version=result.stdout.strip())


def probe_result(requested: dict, process: subprocess.CompletedProcess) -> dict:
    lines = process.stdout.strip().splitlines()
    fields = lines[0].split("\t") if len(lines) == 1 else []
    matched = (process.returncode == 0 and len(fields) == 2 and fields[0] == requested["accession_version"] and
               fields[1].isdigit() and int(fields[1]) == requested["expected_length"])
    absent = process.returncode != 0 and not process.stdout.strip() and "entry not found" in process.stderr.lower()
    status = "exact_identity_and_length" if matched else "absent_recorded" if absent and requested["allow_absent"] else "failed"
    return dict(requested=requested, status=status, passed=status != "failed", exit_code=process.returncode,
                stdout=process.stdout, stderr=process.stderr)


def verify_download_evidence(directory: Path, runtime: dict) -> dict:
    marker = read_json(directory / "DOWNLOAD_VERIFIED.json")
    metadata = read_json(directory / "core_nt-nucl-metadata.json")
    conflicts = [conflict for conflict in metadata_conflicts(runtime, metadata)
                 if conflict["field"] != "database_version" or "db-version" in metadata]
    if (marker.get("status") != "all_archives_checksum_verified" or marker.get("remote_release_metadata_unchanged") is not True
            or marker.get("archives") != metadata.get("number-of-volumes") or conflicts):
        raise ContractError("Pinned download completion or release metadata do not match the extracted database")
    if digest(directory / "metadata-after-download.json") != digest(directory / "core_nt-nucl-metadata.json"):
        raise ContractError("Download release metadata changed")
    root = directory.parent.resolve()
    request_hashes = {}
    for line in (directory / "request.sha256").read_text(encoding="utf-8").splitlines():
        fields = line.split(maxsplit=1)
        if len(fields) != 2 or not re.fullmatch(r"[0-9a-f]{64}", fields[0]):
            raise ContractError("Invalid download request checksum manifest")
        name = fields[1].lstrip("*")
        path = (root / name).resolve()
        if name in request_hashes or not path.is_relative_to(directory.resolve()) or digest(path) != fields[0]:
            raise ContractError("Pinned download request changed or contains an unsafe path")
        request_hashes[name] = fields[0]
    required = {"evidence/core_nt-nucl-metadata.json", "evidence/archives.txt", "evidence/archives.md5", "evidence/download_core_nt_snapshot.sh"}
    if set(request_hashes) != required:
        raise ContractError("Download request inventory is incomplete")
    archives = (directory / "archives.txt").read_text(encoding="utf-8").splitlines()
    expected = [Path(name).name for name in metadata["files"]]
    if archives != expected or len(set(archives)) != len(archives):
        raise ContractError("Archive membership differs from the pinned release")
    log_lines = (directory / "archive_checksum_verification.log").read_text(encoding="utf-8").splitlines()
    if sorted(log_lines) != sorted(name + ": OK" for name in archives):
        raise ContractError("Archive checksum report is incomplete")
    exit_record = read_json(directory / "last_exit.json")
    if exit_record.get("exit_code") != 0:
        raise ContractError("Download service did not record successful completion")
    return dict(directory=str(directory.resolve()), marker=marker, release_metadata=metadata, request_hashes=request_hashes,
                files={name: digest(directory / name) for name in ("DOWNLOAD_VERIFIED.json", "request.sha256", "last_exit.json",
                                                                   "metadata-after-download.json", "archive_checksum_verification.log")},
                scope="Pinned launcher's archive checksum evidence; extracted content is independently hashed and checked below")


def accept_database(prefix: Path, inventory_file: Path, probes_file: Path, destination: Path,
                    blastn: str = "blastn", download_evidence: Path | None = None, timeout_seconds: float | None = None,
                    taxonomy_support: Path | None = None, blastdbcmd: str = "blastdbcmd", blastdbcheck: str = "blastdbcheck") -> dict:
    prefix = prefix.resolve()
    inventory = read_json(inventory_file)
    if (inventory.get("metadata_consistent") is not True or inventory.get("metadata_conflicts") or
            inventory.get("all_database_files_hashed") is not True or not all(record.get("sha256") for record in inventory.get("files", []))):
        raise ContractError("Acceptance requires consistent metadata and complete content hashes")
    probes = read_json(probes_file, PROBE_SCHEMA)["probes"]
    if len({row["accession_version"] for row in probes}) != len(probes):
        raise ContractError("Accession probes must be unique")
    if prefix.name == "core_nt" and download_evidence is None:
        raise ContractError("core_nt acceptance requires the pinned download evidence")
    if taxonomy_support is None:
        raise ContractError("Full blastdbcheck requires an explicit taxonomy support directory containing taxdb.bti and taxdb.btd")
    support = dict(directory=str(taxonomy_support.resolve()), files={name: digest(taxonomy_support / name) for name in ("taxdb.bti", "taxdb.btd")})
    environment = dict(os.environ, BLASTDB=support["directory"], LC_ALL="C")
    binaries = {name: executable_record(value) for name, value in (("blastn", blastn), ("blastdbcmd", blastdbcmd), ("blastdbcheck", blastdbcheck))}
    verify_inventory(inventory, prefix)
    observed = subprocess.run([binaries["blastdbcmd"]["path"], "-db", str(prefix), "-info"], capture_output=True, text=True, check=True, timeout=120, env=environment)
    runtime = parse_database_info(observed.stdout)
    if runtime != inventory["runtime_metadata"] or metadata_conflicts(runtime, inventory["supplied_json_metadata"]):
        raise ContractError("Live database metadata differs from the complete inventory")
    volumes = volume_evidence(prefix, runtime["volume_names"], binaries["blastdbcmd"]["path"], environment)
    assert_single_build(volumes)
    if (len(volumes) != runtime["volume_count"] or len({row["volume"] for row in volumes}) != len(volumes)
            or sum(row["metadata"]["sequence_count"] for row in volumes) != runtime["sequence_count"]
            or sum(row["metadata"]["base_count"] for row in volumes) != runtime["base_count"]
            or any(row["metadata"]["build_date"] != runtime["build_date"] or row["metadata"]["database_version"] != runtime["database_version"] for row in volumes)):
        raise ContractError("Volume counts, versions or build dates disagree with aggregate metadata")
    download = verify_download_evidence(download_evidence, runtime) if download_evidence is not None else None
    destination.mkdir(parents=True, exist_ok=False)
    command = [binaries["blastdbcheck"]["path"], "-db", str(prefix), "-dbtype", "nucl", "-full", "-verbosity", "3"]
    save_json(destination / "request.json", dict(database_prefix=str(prefix), inventory_sha256=digest(inventory_file), probes_sha256=digest(probes_file),
                                                binaries=binaries, taxonomy_support=support, environment=dict(BLASTDB=support["directory"], LC_ALL="C"),
                                                command=command, started_utc=datetime.now(timezone.utc).isoformat()))
    started = time.monotonic()
    try:
        with (destination / "integrity_stdout.txt").open("x", encoding="utf-8") as stdout, (destination / "integrity_stderr.txt").open("x", encoding="utf-8") as stderr:
            checked = subprocess.run(command, stdout=stdout, stderr=stderr, check=False, timeout=timeout_seconds, env=environment)
        probe_rows = []
        for requested in probes:
            probe_command = [binaries["blastdbcmd"]["path"], "-db", str(prefix), "-entry", requested["accession_version"], "-target_only", "-outfmt", "%a\t%l"]
            process = subprocess.run(probe_command, capture_output=True, text=True, check=False, timeout=60, env=environment)
            probe_rows.append(dict(probe_result(requested, process), command=probe_command))
        save_json(destination / "checks.json", dict(volumes=volumes, accession_probes=probe_rows, full_integrity_exit=checked.returncode,
                                                   full_integrity_elapsed_seconds=time.monotonic() - started))
        if checked.returncode or not all(row["passed"] for row in probe_rows) or not any(row["status"] == "exact_identity_and_length" for row in probe_rows):
            raise ContractError("Full integrity check or exact-accession probes failed; no acceptance artifact written")
        verify_inventory(inventory, prefix, hash_content=False)
        if any(digest(taxonomy_support / name) != checksum for name, checksum in support["files"].items()):
            raise ContractError("BLAST taxonomy support changed during acceptance")
    except (ContractError, OSError, subprocess.SubprocessError) as error:
        save_json(destination / "failed.json", dict(error=str(error), accepted=False, elapsed_seconds=time.monotonic() - started))
        raise
    report = dict(contract="database_acceptance_v1", accepted=True, database_prefix=str(prefix), inventory_sha256=digest(inventory_file),
                  inventory=inventory, full_integrity_passed=True, exact_accession_probes_passed=True, binaries=binaries,
                  download=download, taxonomy_support=support, accepted_utc=datetime.now(timezone.utc).isoformat(), final_research_freeze=False,
                  identity_scope="Full extracted-file SHA256 inventory, consistent volume metadata, full blastdbcheck and declared exact-ID/length probes",
                  runtime_mutation_guard="Membership, size and nanosecond mtime checked at startup and between batches; keep accepted files immutable",
                  files={path.name: digest(path) for path in destination.iterdir() if path.is_file()})
    save_json(destination / "accepted.pending.json", report)
    (destination / "accepted.pending.json").rename(destination / "accepted.json")
    return report


def validate_acceptance(path: Path, prefix: Path, evidence: Path, blastn: str) -> dict:
    report = read_json(path)
    if (report.get("contract") != "database_acceptance_v1" or report.get("accepted") is not True
            or report.get("full_integrity_passed") is not True or report.get("exact_accession_probes_passed") is not True
            or Path(report["database_prefix"]).resolve() != prefix.resolve() or report["inventory_sha256"] != digest(evidence)):
        raise ContractError("Production search requires matching accepted database evidence")
    if (path.parent / "failed.json").exists() or (path.parent / "quarantine.json").exists():
        raise ContractError("Database acceptance has been failed or quarantined")
    for name, checksum in report["files"].items():
        if Path(name).name != name or digest(path.parent / name) != checksum:
            raise ContractError("Database acceptance evidence changed")
    if set(report["files"]) != {"request.json", "checks.json", "integrity_stdout.txt", "integrity_stderr.txt"}:
        raise ContractError("Database acceptance artifact inventory is incomplete")
    if report["binaries"]["blastn"] != executable_record(blastn):
        raise ContractError("BLAST executable differs from the accepted version or binary")
    support = report["taxonomy_support"]
    if set(support["files"]) != {"taxdb.bti", "taxdb.btd"} or any(digest(Path(support["directory"]) / name) != checksum for name, checksum in support["files"].items()):
        raise ContractError("Accepted BLAST taxonomy support changed")
    if report["inventory"] != read_json(evidence):
        raise ContractError("Accepted inventory differs from supplied database evidence")
    verify_inventory(report["inventory"], prefix, hash_content=False)
    return dict(file=str(path.resolve()), sha256=digest(path), contract=report["contract"],
                accepted_utc=report["accepted_utc"], full_integrity_passed=True, taxonomy_support=support)